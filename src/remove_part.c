
#include "allvars.h"
#include "proto.h"

#define SIZE_FILL 96
#define SKIP {my_fwrite(&blksize, sizeof(int), 1, file);}

int NumPartCut, NumGasNew, NumPartNew, TotNumPart, *PartFlag;
double MTot, CMVel[3];
FILE *file;

void remove_part_block(int block);
void write_part_block(int block);


void remove_part()
{
  char snapdir[MAX_STRING_LEN], fill[SIZE_FILL];
  int i, block, blksize, dummy;
  double mtot, cmvel[3], fac;

  /*box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
    }

  All.CutSizeHalf = All.CutSize / 2;

  All.BoxSize = All.CutSize;*/

  NumPartCut = NumGas;

  /*if(All.CutDM)
    NumPartCut += NumPart[1];*/

  PartFlag = mymalloc_movable(&PartFlag, "PartFlag", NumPartCut * sizeof(int));

  memset(PartFlag, 0, NumPartCut * sizeof(int));

  for(i = NumPartNew = 0; i < NumPartCut; i++)
    if(PV.ID[i] != 0 && PV.Mass[i] != 0)
      {
	PartFlag[i] = 1;

	if(i < NumGas)
	  NumGasNew++;

	NumPartNew++;
      }

  if(All.RemoveCosm)
    {
      for(i = 0; i < 3; i++)
	CMVel[i] = 0;
    }
  else
    {
      for(i = mtot = cmvel[0] = cmvel[1] = cmvel[2] = 0; i < NumPartCut; i++)
	if(PartFlag[i])
	  {
	    mtot += PV.Mass[i];

	    cmvel[0] += PV.Mass[i] * PV.VX[i];
	    cmvel[1] += PV.Mass[i] * PV.VY[i];
	    cmvel[2] += PV.Mass[i] * PV.VZ[i];
	  }

      MPI_Allreduce(&mtot, &MTot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(cmvel, CMVel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < 3; i++)
	CMVel[i] /= MTot;
    }

  for(block = 0; block < NUM_BLOCKS; block++)
    remove_part_block(block);

  memset(NumPart, 0, NUM_PART_TYPES * sizeof(int));

  NumPart[0] = NumGasNew;
  //NumPart[1] = NumPartNew - NumGasNew;

  for(i = TotNumPart = 0; i < NUM_PART_TYPES; i++)
    TotNumPart += All.GlobNumPart[i];

  mpi_printf("Old particle count: %d (gas) / %d (total)\n", All.GlobNumPart[0], TotNumPart);

  MPI_Allreduce(NumPart, All.GlobNumPart, NUM_PART_TYPES, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = TotNumPart = 0; i < NUM_PART_TYPES; i++)
    TotNumPart += All.GlobNumPart[i];

  mpi_printf("New particle count: %d (gas) / %d (total)\n", All.GlobNumPart[0], TotNumPart);

  memset(All.Mass, 0, NUM_PART_TYPES * sizeof(double));

  if(All.LengthUnit == 1)
    fac = 1e-3;
  else if(All.LengthUnit == 2)
    fac = ASTRONOMICAL_UNIT / UNIT_LENGTH;
  else
    fac = 1;

  if(!All.ComovingUnits && All.RemoveCosm)
    fac *= (1 + All.RedShift);

  if(All.ComovingUnits && !All.RemoveCosm)
    fac /= (1 + All.RedShift);

  if(All.RemoveCosm)
    fac *= All.HubbleParam;

  All.BoxSize *= fac;

  mpi_printf("New Box Size: %g\n", All.BoxSize);

  All.NumFiles = NTask;

  if(!All.RemoveCosm)
    {
      All.Time = All.RedShift = All.OmegaM = All.OmegaLambda = 0;

      All.HubbleParam = 1;
    }

  dummy = 0;

  memset(fill, 0, SIZE_FILL * sizeof(char));

  sprintf(snapdir, "%s/%s", All.Path, All.Base);

  mkdir(snapdir, 0755);

  if(All.RemoveCosm)
    sprintf(All.SnapFile, "%s/%s/%s_ic", All.Path, All.Base, All.Base);
  else
    sprintf(All.SnapFile, "%s/%s/%s_900", All.Path, All.Base, All.Base);

  if(NTask > 1)
    sprintf(All.SnapFile, "%s.%d", All.SnapFile, ThisTask);

  if(!(file = fopen(All.SnapFile, "w")))
    {
      printf("Can't open file for writing snapshot!");

      terminate("");
    }

  blksize = 256;

  SKIP;
  my_fwrite(NumPart, sizeof(int), NUM_PART_TYPES, file);
  my_fwrite(All.Mass, sizeof(double), NUM_PART_TYPES, file);
  my_fwrite(&All.Time, sizeof(double), 1, file);
  my_fwrite(&All.RedShift, sizeof(double), 1, file);
  my_fwrite(&dummy, sizeof(int), 1, file);
  my_fwrite(&dummy, sizeof(int), 1, file);
  my_fwrite(All.GlobNumPart, sizeof(int), NUM_PART_TYPES, file);
  my_fwrite(&dummy, sizeof(int), 1, file);
  my_fwrite(&All.NumFiles, sizeof(int), 1, file);
  my_fwrite(&All.BoxSize, sizeof(double), 1, file);
  my_fwrite(&All.OmegaM, sizeof(double), 1, file);
  my_fwrite(&All.OmegaLambda, sizeof(double), 1, file);
  my_fwrite(&All.HubbleParam, sizeof(double), 1, file);
  my_fwrite(fill, sizeof(char), SIZE_FILL, file);
  SKIP;

  for(block = 0; block < NUM_BLOCKS; block++)
    write_part_block(block);

  myfree(PartFlag);
}


void remove_part_block(int block)
{
  int i, j;
  double fac;

  if(BlockFlag[block])
    {
      if(block == IO_POS)
	{
	  if(All.LengthUnit == 1)
	    fac = 1e-3;
	  else if(All.LengthUnit == 2)
	    fac = ASTRONOMICAL_UNIT / UNIT_LENGTH;
	  else
	    fac = 1;

	  if(!All.ComovingUnits && All.RemoveCosm)
	    fac *= (1 + All.RedShift);

	  if(All.ComovingUnits && !All.RemoveCosm)
	    fac /= (1 + All.RedShift);

	  if(All.RemoveCosm)
	    fac *= All.HubbleParam;

	  for(i = j = 0; i < NumPartCut; i++)
	    if(PartFlag[i])
	      {
		PV.X[j] *= fac;
		PV.Y[j] *= fac;
		PV.Z[j] *= fac;

		j++;
	      }
	}
      else if(block == IO_VEL)
	{
	  if(All.RemoveCosm)
	    fac = sqrt(1 + All.RedShift);
	  else
	    fac = 1;

	  for(i = j = 0; i < NumPartCut; i++)
	    if(PartFlag[i])
	      {
		PV.VX[j] = fac * (PV.VX[i] - CMVel[0]);
		PV.VY[j] = fac * (PV.VY[i] - CMVel[1]);
		PV.VZ[j] = fac * (PV.VZ[i] - CMVel[2]);

		j++;
	      }
	}
      else if(block == IO_ID)
	{
	  for(i = j = 0; i < NumPartCut; i++)
	    if(PartFlag[i])
	      {
		PV.ID[j] = PV.ID[i];

		j++;
	      }
	}
      else if(block == IO_MASS)
	{
	  fac = SOLAR_MASS / UNIT_MASS;

	  if(All.RemoveCosm)
	    fac *= All.HubbleParam;

	  for(i = j = 0; i < NumPartCut; i++)
	    if(PartFlag[i])
	      {
		PV.Mass[j] = fac * PV.Mass[i];

		j++;
	      }
	}
      else if(block == IO_U)
	{
	  fac = UNIT_MASS / UNIT_ENERGY;

	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.U[j] = fac * PV.U[i];

		j++;
	      }
	}
      else if(block == IO_RHO)
	{
	  fac = pow(UNIT_LENGTH, 3) / UNIT_MASS;

	  if(All.RemoveCosm)
	    fac /= pow(All.HubbleParam, 2) * pow(1 + All.RedShift, 3);

	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.Rho[j] = fac * PV.Rho[i];

		j++;
	      }
	}
/*      else if(block == IO_VOL)
	{
	  if(All.LengthUnit == 1)
	    fac = 1e-9;
	  else if(All.LengthUnit == 2)
	    fac = pow(ASTRONOMICAL_UNIT / UNIT_LENGTH, 3);
	  else
	    fac = 1;

	  if(!All.ComovingUnits && All.CutCosm)
	    fac *= pow(1 + All.RedShift, 3);

	  if(All.ComovingUnits && !All.CutCosm)
	    fac /= pow(1 + All.RedShift, 3);

	  if(All.CutCosm)
	    fac *= pow(All.HubbleParam, 3);

	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.Vol[j] = fac * PV.Vol[i];

		j++;
	      }
	}
      else if(block == IO_DELAUNAY)
	{
	  if(All.LengthUnit == 1)
	    fac = 1e-3;
	  else if(All.LengthUnit == 2)
	    fac = ASTRONOMICAL_UNIT / UNIT_LENGTH;
	  else
	    fac = 1;

	  if(!All.ComovingUnits && All.CutCosm)
	    fac *= (1 + All.RedShift);

	  if(All.ComovingUnits && !All.CutCosm)
	    fac /= (1 + All.RedShift);

	  if(All.CutCosm)
	    fac *= All.HubbleParam;

	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.Delaunay[j] = fac * PV.Delaunay[i];

		j++;
	      }
	}*/
      else if(block == IO_GRAVACC)
	{
	  fac = pow(1 + All.RedShift, 2);

	  for(i = j = 0; i < NumPartCut; i++)
	    if(PartFlag[i])
	      {
		PV.GravAccX[j] = fac * PV.GravAccX[i];
		PV.GravAccY[j] = fac * PV.GravAccY[i];
		PV.GravAccZ[j] = fac * PV.GravAccZ[i];

		j++;
	      }
	}
      else if(block == IO_GRADP)
	{
	  fac = pow(1 + All.RedShift, -4) * pow(All.HubbleParam, -3);

	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.GradPX[j] = fac * PV.GradPX[i];
		PV.GradPY[j] = fac * PV.GradPY[i];
		PV.GradPZ[j] = fac * PV.GradPZ[i];

		j++;
	      }
	}
      else if(block == IO_CHEM)
	{
	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
#ifdef OLDCHEM
		PV.AbH2[j] = PV.AbH2[i];
		PV.AbHII[j] = PV.AbHII[i];
		PV.AbDII[j] = PV.AbDII[i];
		PV.AbHD[j] = PV.AbHD[i];
		PV.AbHeII[j] = PV.AbHeII[i];
		PV.AbHeIII[j] = PV.AbHeIII[i];
#else
		PV.AbHM[j] = PV.AbHM[i];
		PV.AbH2[j] = PV.AbH2[i];
		PV.AbHII[j] = PV.AbHII[i];
#endif

		j++;
	      }
	}
      else if(block == IO_GAMMA)
	{
	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.Gamma[j] = PV.Gamma[i];

		j++;
	      }
	}
      else if(block == IO_RATES)
	{
	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.PdVRate[j] = PV.PdVRate[i];
		PV.H2Rate[j] = PV.H2Rate[i];
		PV.CIERate[j] = PV.CIERate[i];
		PV.ChemRate[j] = PV.ChemRate[i];

		j++;
	      }
	}
      else if(block == IO_ALLOWREF)
	{
	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.AllowRef[j] = PV.AllowRef[i];

		j++;
	      }
	}
      else if(block == IO_DIVVEL)
	{
	  for(i = j = 0; i < NumGas; i++)
	    if(PartFlag[i])
	      {
		PV.DivVel[j] = PV.DivVel[i];

		j++;
	      }
	}
      else{
	terminate("Unknown block!");
      }
    }
}


void write_part_block(int block)
{
  int i, blksize;

  if(BlockFlag[block])
    {
      if(block == IO_POS)
	{
	  blksize = 3 * NumPartNew;

	  SKIP;
	  for(i = 0; i < NumPartNew; i++)
	    {
	      my_fwrite(&PV.X[i], sizeof(double), 1, file);
	      my_fwrite(&PV.Y[i], sizeof(double), 1, file);
	      my_fwrite(&PV.Z[i], sizeof(double), 1, file);
	    }
	  SKIP;
	}
      else if(block == IO_VEL)
	{
	  blksize = 3 * NumPartNew;

	  SKIP;
	  for(i = 0; i < NumPartNew; i++)
	    {
	      my_fwrite(&PV.VX[i], sizeof(double), 1, file);
	      my_fwrite(&PV.VY[i], sizeof(double), 1, file);
	      my_fwrite(&PV.VZ[i], sizeof(double), 1, file);
	    }
	  SKIP;
	}
      else if(block == IO_ID)
	{
	  blksize = NumPartNew;

	  SKIP;
	  my_fwrite(PV.ID, sizeof(int), NumPartNew, file);
	  SKIP;
	}
      else if(block == IO_MASS)
	{
	  blksize = NumPartNew;

	  SKIP;
	  my_fwrite(PV.Mass, sizeof(double), NumPartNew, file);
	  SKIP;
	}
      else if(block == IO_U)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.U, sizeof(double), NumGasNew, file);
	  SKIP;
	}
      else if(block == IO_RHO)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.Rho, sizeof(double), NumGasNew, file);
	  SKIP;
	}
/*      else if(block == IO_VOL)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.Vol, sizeof(double), NumGasNew, file);
	  SKIP;
	}
      else if(block == IO_DELAUNAY)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.Delaunay, sizeof(double), NumGasNew, file);
	  SKIP;
	}*/
      else if(block == IO_GRAVACC)
	{
	  blksize = 3 * NumPartNew;

	  SKIP;
	  for(i = 0; i < NumPartNew; i++)
	    {
	      my_fwrite(&PV.GravAccX[i], sizeof(double), 1, file);
	      my_fwrite(&PV.GravAccY[i], sizeof(double), 1, file);
	      my_fwrite(&PV.GravAccZ[i], sizeof(double), 1, file);
	    }
	  SKIP;
	}
      else if(block == IO_GRADP)
	{
	  blksize = 3 * NumGasNew;

	  SKIP;
	  for(i = 0; i < NumGasNew; i++)
	    {
	      my_fwrite(&PV.GradPX[i], sizeof(double), 1, file);
	      my_fwrite(&PV.GradPY[i], sizeof(double), 1, file);
	      my_fwrite(&PV.GradPZ[i], sizeof(double), 1, file);
	    }
	  SKIP;
	}
      else if(block == IO_CHEM)
	{
	  blksize = 6 * NumGasNew;

	  SKIP;
	  for(i = 0; i < NumGasNew; i++)
	    {
#ifdef OLDCHEM
	      my_fwrite(&PV.AbH2[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbHII[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbDII[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbHD[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbHeII[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbHeIII[i], sizeof(double), 1, file);
#else
	      my_fwrite(&PV.AbHM[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbH2[i], sizeof(double), 1, file);
	      my_fwrite(&PV.AbHII[i], sizeof(double), 1, file);
#endif
	    }
	  SKIP;
	}
      else if(block == IO_GAMMA)
	{	
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.Gamma, sizeof(double), NumGasNew, file);
	  SKIP;
	}
      else if(block == IO_RATES)
	{
	  blksize = 4 * NumGasNew;

	  SKIP;
	  for(i = 0; i < NumGasNew; i++)
	    {
	      my_fwrite(&PV.PdVRate[i], sizeof(double), 1, file);
	      my_fwrite(&PV.H2Rate[i], sizeof(double), 1, file);
	      my_fwrite(&PV.CIERate[i], sizeof(double), 1, file);
	      my_fwrite(&PV.ChemRate[i], sizeof(double), 1, file);
	    }
	  SKIP;
	}
      else if(block == IO_ALLOWREF)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.AllowRef, sizeof(int), NumGasNew, file);
	  SKIP;
	}
      else if(block == IO_DIVVEL)
	{
	  blksize = NumGasNew;

	  SKIP;
	  my_fwrite(PV.DivVel, sizeof(double), NumGasNew, file);
	  SKIP;
	}
      else
	terminate("Unknown block!");
    }
}
