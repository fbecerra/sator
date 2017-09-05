
#include "allvars.h"
#include "proto.h"

#define CONVERT(x, y) (data_type == 1 ? *((float *) x + y) : *((double *) x + y));

static int BlockDataType, BlockBytesPerElement, BlockNumDims;
static size_t BlockOffset;
static int DataDataType[NUM_BLOCKS], DataBytesPerElement[NUM_BLOCKS];
static int DataNumDims[NUM_BLOCKS], DataNumElements[NUM_BLOCKS];
static size_t DataOffset[NUM_BLOCKS];

void get_snap_path(void);
void get_part_distribution(void);
void get_part_list(void);
void read_blocks(void);
int block_present(int block, int type);
void get_block_info(int type, int snap, int block, int list_start);
void get_block_value(int block);
void calc_temp(void);


void read_snap()
{
  All.SnapNum = All.SnapNumStart + All.Iter;

  get_snap_path();

  get_part_distribution();

  get_part_list();

  read_blocks();

  calc_temp();

  /*
  int i;
  int vali, vglobmini, vmini, vglobmaxi, vmaxi;

  vmini = INT_MAX;
  vmaxi = INT_MIN;

  for(i = 0; i < NumGas; i++)
    {
      vali = PV.SIDM_NumTotalScatter[i];

      vmini = dmin(vali, vmini);
      vmaxi = dmax(vali, vmaxi);
    }

  MPI_Allreduce(&vmini, &vglobmini, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&vmaxi, &vglobmaxi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("min = %d, max = %d\n", vglobmini, vglobmaxi);
  */
  /*
  int j;
  double vald, vmind, vglobmind, vmaxd, vglobmaxd;

  vmind = DBL_MAX;
  vmaxd = -DBL_MAX;

  for(j = 0; j < NumGas; j++)
    {
      vald = PV.DivVel[j];

      vmind = dmin(vald, vmind);
      vmaxd = dmax(vald, vmaxd);
    }

  MPI_Allreduce(&vmind, &vglobmind, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&vmaxd, &vglobmaxd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("min = %g, max = %g\n", vglobmind, vglobmaxd);
  */
}


void get_snap_path()
{
  if(ThisTask == 0)
    {
      char *systype;
      FILE *fd;

      systype = getenv("SYSTYPE");

      if(!strcmp(systype, "odyssey-opteron"))
//	strcpy(All.Path, "/n/hernquistfs3/fbecerra");
	strcpy(All.Path, "/n/hernquistfs2/fbecerra");
      else if(!strcmp(systype, "odin"))
	strcpy(All.Path, "/ptmp/mpa/tgreif");
      else if(!strcmp(systype, "stampede"))
//	strcpy(All.Path, "/scratch/02563/fbecerra");
	strcpy(All.Path, "/scratch/00025/tgreif");
      else if(!strcmp(systype, "lonestar"))
	strcpy(All.Path, "/scratch/00025/tgreif");
      else
	terminate("Unknown system type!");

      if(All.SnapNum == -1)
	sprintf(All.SnapNumString, "ic");
      else
	sprintf(All.SnapNumString, "%03d", All.SnapNum);

      sprintf(All.SnapFile, "%s/%s/%s_%s", All.Path, All.Base, All.Base, All.SnapNumString);

      strcpy(All.SnapBase, All.SnapFile);

      if(!(fd = fopen(All.SnapFile, "r")))
	{
	  sprintf(All.SnapFile, "%s.0", All.SnapBase);

	  if(!(fd = fopen(All.SnapFile, "r")))
	    {
	      if(!strcmp(All.SnapNumString, "ic"))
		{
		  terminate("Could not find snapshot file!");
		}
	      else
		{
		  sprintf(All.SnapBase, "%s/%s/snapdir_%s/%s_%s", All.Path, All.Base, All.SnapNumString, All.Base, All.SnapNumString);

		  sprintf(All.SnapFile, "%s.0", All.SnapBase);

		  if(!(fd = fopen(All.SnapFile, "r")))
		    terminate("Could not find snapshot file!");
		}
	    }
	}

      fclose(fd);
    }

  MPI_Bcast(&All.Path, MAX_STRING_LEN * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.SnapNumString, MAX_STRING_LEN * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.SnapBase, MAX_STRING_LEN * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.SnapFile, MAX_STRING_LEN * sizeof(char), MPI_BYTE, 0, MPI_COMM_WORLD);
}


void get_part_distribution()
{
  int task, type, snap;
  FILE *fd;

  if(ThisTask == 0)
    {
      if(!(fd = fopen(All.SnapFile, "r")))
	terminate("Could not read snapshot file!");

      fseek(fd, 28, SEEK_SET);
      my_fread(&All.Mass, sizeof(double), NUM_PART_TYPES, fd);
      my_fread(&All.Time, sizeof(double), 1, fd);
      my_fread(&All.RedShift, sizeof(double), 1, fd);
      fseek(fd, 8, SEEK_CUR);
      my_fread(&All.GlobNumPart[0], sizeof(int), NUM_PART_TYPES, fd);
      fseek(fd, 4, SEEK_CUR);
      my_fread(&All.NumFiles, sizeof(int), 1, fd);
      my_fread(&All.BoxSize, sizeof(double), 1, fd);
      my_fread(&All.OmegaM, sizeof(double), 1, fd);
      my_fread(&All.OmegaLambda, sizeof(double), 1, fd);
      my_fread(&All.HubbleParam, sizeof(double), 1, fd);

      // Debug
      mpi_printf("Time: %.10f \n", All.Time);

      fclose(fd);
    }

  MPI_Bcast(All.Mass, NUM_PART_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.Time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.RedShift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(All.GlobNumPart, NUM_PART_TYPES, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.NumFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.BoxSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.OmegaM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.OmegaLambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.HubbleParam, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  All.SnapNumPart = mymalloc("SnapNumPart", All.NumFiles * sizeof(struct PartList));
  All.TaskNumPart = mymalloc("TaskNumPart", NTask * sizeof(struct PartList));

  /*
  if(ThisTask == 0)
    for(type = 0; type < NUM_PART_TYPES; type++)
      printf("%d\n", All.GlobNumPart[type]);
  */

  if(ThisTask == 0)
    {
      for(snap = 0; snap < All.NumFiles; snap++)
	{
	  if(All.NumFiles > 1)
	    sprintf(All.SnapFile, "%s.%d", All.SnapBase, snap);

	  if(!(fd = fopen(All.SnapFile, "r")))
	    terminate("Could not read snapshot file!");

	  fseek(fd, 4, SEEK_SET);

	  my_fread(&All.SnapNumPart[snap].type[0], sizeof(int), NUM_PART_TYPES, fd);

	  fclose(fd);
	}

      for(task = 0; task < NTask; task++)
	for(type = 0; type < NUM_PART_TYPES; type++)
	  {
	    All.TaskNumPart[task].type[type] = All.GlobNumPart[type] / NTask;

	    if(task < All.GlobNumPart[type] % NTask)
	      All.TaskNumPart[task].type[type] += 1;
	  }
    }

  MPI_Bcast(All.SnapNumPart, All.NumFiles * sizeof(struct PartList), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(All.TaskNumPart, NTask * sizeof(struct PartList), MPI_BYTE, 0, MPI_COMM_WORLD);

  for(type = TotNumPart = 0; type < NUM_PART_TYPES; type++)
    {
      NumPart[type] = All.TaskNumPart[ThisTask].type[type];

      TotNumPart += NumPart[type];
    }

  NumGas = NumPart[0];
  NumSinks = NumPart[5];

  MPI_Allreduce(&NumSinks, &All.TotNumSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  NumSinksIter[All.Iter] = All.TotNumSinks;
}


void get_part_list()
{
  int i, j;
  int partcount, rank_idx;
  int count[NUM_PART_TYPES];

  for(i = 0; i < NUM_PART_TYPES; i++)
    {
      partcount = rank_idx = 0;

      SnapListStart[i] = 0;
      PartListStart[i] = 0;

      SnapListEnd[i] = 0;
      PartListEnd[i] = -1;

      for(j = 0; j < All.NumFiles; j++)
	{
	  partcount += All.SnapNumPart[j].type[i];

	  while(partcount > All.TaskNumPart[rank_idx].type[i])
	    {
	      partcount -= All.TaskNumPart[rank_idx].type[i];

	      if(ThisTask == rank_idx)
		{
		  SnapListEnd[i] = j;
		  PartListEnd[i] = All.SnapNumPart[j].type[i] - partcount - 1;
		}

	      rank_idx += 1;

	      if(ThisTask == rank_idx)
		{
		  SnapListStart[i] = j;
		  PartListStart[i] = All.SnapNumPart[j].type[i] - partcount;
		}
	    }

	  if(ThisTask == rank_idx)
	    {
	      SnapListEnd[i] = j;
	      PartListEnd[i] = All.SnapNumPart[j].type[i] - 1;
	    }
	}
    }

  for(i = 0; i < NUM_PART_TYPES; i++)
    {
      count[i] = 0;

      for(j = SnapListStart[i]; j <= SnapListEnd[i]; j++)
	{
	  if(j == SnapListEnd[i])
	    {
	      if(j == SnapListStart[i])
		count[i] += PartListEnd[i] - PartListStart[i] + 1;
	      else
		count[i] += PartListEnd[i] + 1;
	    }
	  else if(j == SnapListStart[i])
	    count[i] += All.SnapNumPart[j].type[i] - PartListStart[i];
	  else
	    count[i] += All.SnapNumPart[j].type[i];
	}

      if(count[i] != NumPart[i])
	{
	  printf("task %d: type = %d, count = %d, desnumpart = %d\n", ThisTask, i, count[i], NumPart[i]); 

	  terminate("Particle numbers don't match!");
	}
    }
  /*
  if(ThisTask == 0)
    for(i = 0; i < All.NumFiles; i++)
      for(j = 0; j < NUM_PART_TYPES; j++)
	{
	  printf("%d %d\n", i, All.SnapNumPart[i].type[j]);

	  fflush(stdout);
	}

  MPI_Barrier(MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
	for(j = 0; j < NUM_PART_TYPES; j++)
	  {
	    printf("task %d: %d\n", i, NumPart[j]);
	    printf("task %d: %d\n", i, SnapListStart[j]);
	    printf("task %d: %d\n", i, SnapListEnd[j]);
	    printf("task %d: %d\n", i, PartListStart[j]);
	    printf("task %d: %d\n\n", i, PartListEnd[j]);
	  }

      MPI_Barrier(MPI_COMM_WORLD);

    }

  MPI_Barrier(MPI_COMM_WORLD);
  */
  //dump_memory_table();

}


void read_blocks()
{
  int i, type, snap, block;
  int list_start, list_end;
  int num_elements, ttype;
  char buf[MAXCHARS - 1];
  FILE *fd;

  for(block = 0; block < NUM_BLOCKS; block++)
    DataNumElements[block] = DataOffset[block] = 0;

  for(block = 0; block < NUM_BLOCKS; block++)
    if(BlockFlag[block])
      for(type = 0; type < NUM_PART_TYPES; type++)
	if(block_present(block, type))
	  //if(NumPart[type] > 0)
	    for(snap = SnapListStart[type]; snap <= SnapListEnd[type]; snap++)
	      {
		if(snap == SnapListEnd[type])
		  {
		    if(snap == SnapListStart[type])
		      {
			list_start = PartListStart[type];
			list_end = PartListEnd[type];
		      }
		    else
		      {
			list_start = 0;
			list_end = PartListEnd[type];
		      }
		  }
		else if(snap == SnapListStart[type])
		  {
		    list_start = PartListStart[type];
		    list_end = All.SnapNumPart[snap].type[type] - 1;
		  }
		else
		  {
		    list_start = 0;
		    list_end = All.SnapNumPart[snap].type[type] - 1;
		  }

		get_block_info(type, snap, block, list_start);

		DataDataType[block] = BlockDataType;
		DataBytesPerElement[block] = BlockBytesPerElement;
		DataNumDims[block] = BlockNumDims;

		num_elements = BlockNumDims * (list_end - list_start + 1);

		DataNumElements[block] += num_elements;

		//if((block == IO_SIDM_HSML || block == IO_SIDM_DENSITY) && type == 1)
		//printf("bla = %d %d %d %d %d %d %d %d %d %d %llu\n", block, type, ThisTask, snap, list_start, list_end, BlockDataType, BlockBytesPerElement, BlockNumDims, num_elements, DataOffset[block]);

		sprintf(buf, "PD[%d]", block);

		if(PD[block] == 0)
		  PD[block] = mymalloc_movable(&PD[block], buf, DataNumElements[block] * BlockBytesPerElement);
		else
		  PD[block] = myrealloc_movable(PD[block], DataNumElements[block] * BlockBytesPerElement);

		if(block == IO_MASS && All.Mass[type] != 0)
		  {
		    for(i = 0; i < num_elements / BlockNumDims; i++)
		      memcpy(PD[block] + DataOffset[block] + BlockBytesPerElement * i, &All.Mass[type], sizeof(double));
		  }
		else
		  {
		    if(All.NumFiles > 1)
		      sprintf(All.SnapFile, "%s.%d", All.SnapBase, snap);
		    else
		      strcpy(All.SnapFile, All.SnapBase);

		    if(!(fd = fopen(All.SnapFile, "r")))
		      terminate("Could not read snapshot file!");

		    fseek(fd, BlockOffset, SEEK_SET);

		    my_fread(PD[block] + DataOffset[block], num_elements * BlockBytesPerElement, 1, fd);

		    fclose(fd);
		  }

		//if(block == IO_SIDM_DENSITY && type == 1)
		//printf("task = %d, off = %llu, %d, %g %g %g %g\n", ThisTask, BlockOffset, num_elements, *((double*) (PD[block] + DataOffset[block])), *((double*) (PD[block] + DataOffset[block] + BlockBytesPerElement)), *((double*) (PD[block] + DataOffset[block] + BlockBytesPerElement * (num_elements - 2))), *((double*) (PD[block] + DataOffset[block] + BlockBytesPerElement * (num_elements - 1))));

		DataOffset[block] += num_elements * BlockBytesPerElement;
	      }

  for(block = 0; block < NUM_BLOCKS; block++)
    if(BlockFlag[block])
      get_block_value(block);
}


int block_present(int block, int type)
{
  if(type == 0 && block != IO_SIDM_PSUM && block != IO_SIDM_NUMNGB && block != IO_SIDM_NUMTOTALSCATTER && block != IO_SIDM_HSML && block != IO_SIDM_DENSITY && block != IO_SIDM_VELDISP)
    return 1;

  if(type == 1 && (block == IO_SIDM_PSUM || block == IO_SIDM_NUMNGB || block == IO_SIDM_NUMTOTALSCATTER || block == IO_SIDM_HSML || block == IO_SIDM_DENSITY || block == IO_SIDM_VELDISP))
    return 1;

  if(block == IO_POS || block == IO_VEL || block == IO_ID || block == IO_MASS || block == IO_GRAVACC)
    return 1;

  return 0;
}


void get_block_info(int type, int snap, int block, int list_start)
{
  int i, j, block_type, num_elements;

  BlockOffset = 264;

  for(i = 0; i <= block; i++)
    if(BlockPresent[i])
      {
	BlockOffset += 4;

	// number of dimensions

	if(i == IO_POS || i == IO_VEL || i == IO_GRAVACC || i == IO_GRADP)
	  BlockNumDims = 3;
	else if(i == IO_CHEM)
#ifdef OLDCHEM
	  BlockNumDims = 6;
#else
	  BlockNumDims = 3;
#endif
	else if(i == IO_RATES)
	  BlockNumDims = 4;
	else
	  BlockNumDims = 1;

	// number of elements

	num_elements = 0;

	if(i == block)
	  block_type = type;
	else
	  block_type = NUM_PART_TYPES;

	if(i == IO_POS || i == IO_VEL || i == IO_ID || i == IO_GRAVACC)
	  {
	    for(j = 0; j < block_type; j++)
	      num_elements += All.SnapNumPart[snap].type[j];
	  }
	else if(i == IO_MASS)
	  {
	    for(j = 0; j < block_type; j++)
	      if(All.Mass[j] == 0)
		num_elements += All.SnapNumPart[snap].type[j];
	  }
	else if(i == IO_SIDM_PSUM || i == IO_SIDM_NUMNGB || i == IO_SIDM_NUMTOTALSCATTER || i == IO_SIDM_HSML || i == IO_SIDM_DENSITY || i == IO_SIDM_VELDISP)
	  {
	    if(i == block)
	      num_elements = 0;
	    else
	      num_elements = All.SnapNumPart[snap].type[1];
	  }
	else
	  {
	    if(i == block)
	      num_elements = 0;
	    else
	      num_elements = All.SnapNumPart[snap].type[0];
	  }

	if(i == block)
	  num_elements += list_start;

	// type and bytes

	if(i == IO_ID || i == IO_ALLOWREF || i == IO_SIDM_NUMTOTALSCATTER)
	  {
	    BlockDataType = 0;

	    BlockBytesPerElement = 4;
	  }
	else
	  {
	    if(All.FlagDouble == 0)
	      {
		BlockDataType = 1;

		BlockBytesPerElement = 4;
	      }
	    else
	      {
		BlockDataType = 2;

		BlockBytesPerElement = 8;
	      }
	  }

	// offset in bytes

	BlockOffset += BlockNumDims * num_elements * BlockBytesPerElement;

	if(i != block)
	  BlockOffset += 4;
      }
}


void get_block_value(int block)
{
  int i;
  double fac, mu, gamma, abe;
  int data_type = DataDataType[block];
  int bytes_per_element = DataBytesPerElement[block];
  int num_dims = DataNumDims[block];
  int num_elements = DataNumElements[block];
  int bytes = num_elements * bytes_per_element;

  if(block == IO_POS)
    {
      if(All.LengthUnit == 0)
	fac = 1;
      else if(All.LengthUnit == 1)
	fac = 1e3;
      else if(All.LengthUnit == 2)
	fac = UNIT_LENGTH / ASTRONOMICAL_UNIT;
      else
	terminate("Length unit not implemented!");

      if(!All.ComovingUnits)
	fac /= (1 + All.RedShift);

      fac /= All.HubbleParam;

      PV.X = mymalloc_movable(&PV.X, "PV.X", bytes / num_dims);
      PV.Y = mymalloc_movable(&PV.Y, "PV.Y", bytes / num_dims);
      PV.Z = mymalloc_movable(&PV.Z, "PV.Z", bytes / num_dims);

      for(i = 0; i < num_elements / num_dims; i++)
	{
	  PV.X[i] = fac * CONVERT(PD[block], num_dims * i);
	  PV.Y[i] = fac * CONVERT(PD[block], num_dims * i + 1);
	  PV.Z[i] = fac * CONVERT(PD[block], num_dims * i + 2);
	}
    }
  else if(block == IO_VEL)
    {
      fac = 1 / sqrt(1 + All.RedShift);

      PV.VX = mymalloc_movable(&PV.VX, "PV.VX", bytes / num_dims);
      PV.VY = mymalloc_movable(&PV.VY, "PV.VY", bytes / num_dims);
      PV.VZ = mymalloc_movable(&PV.VZ, "PV.VZ", bytes / num_dims);

      for(i = 0; i < num_elements / num_dims; i++)
	{
	  PV.VX[i] = fac * CONVERT(PD[block], num_dims * i);
	  PV.VY[i] = fac * CONVERT(PD[block], num_dims * i + 1);
	  PV.VZ[i] = fac * CONVERT(PD[block], num_dims * i + 2);
	}
    }
  else if(block == IO_ID)
    {
      PV.ID = mymalloc_movable(&PV.ID, "PV.ID", bytes);

      for(i = 0; i < num_elements; i++){
	PV.ID[i] = *((int *) PD[block] + i);
        if (PV.ID[i] == 1018821992)
           mpi_printf("id = %i, idx = %i, ", PV.ID[i], i);
      }
    }
  else if(block == IO_MASS)
    {
      PV.Mass = mymalloc_movable(&PV.Mass, "PV.Mass", bytes);

      fac = UNIT_MASS / SOLAR_MASS / All.HubbleParam;

      for(i = 0; i < num_elements; i++){
	PV.Mass[i] = fac * CONVERT(PD[block], i);
        if (PV.ID[i] == 1018821992)
           mpi_printf("id = %i, mass = %.10f \n", PV.ID[i], PV.Mass[i]);
      }
    }
  else if(block == IO_U)
    {
      PV.U = mymalloc_movable(&PV.U, "PV.U", bytes);
      PV.Temp = mymalloc_movable(&PV.Temp, "PV.Temp", bytes);
      fac = UNIT_ENERGY / UNIT_MASS;

      for(i = 0; i < num_elements; i++)
	PV.U[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_RHO)
    {
      PV.Rho = mymalloc_movable(&PV.Rho, "PV.Rho", bytes);
      PV.NH = mymalloc_movable(&PV.NH, "PV.NH", bytes);

      fac = UNIT_MASS / pow(UNIT_LENGTH, 3) * pow(1 + All.RedShift, 3) * pow(All.HubbleParam, 2);

      for(i = 0; i < num_elements; i++)
	PV.Rho[i] = fac * CONVERT(PD[block], i);

      fac = HYDROGEN_MASSFRAC / PROTONMASS;

      for(i = 0; i < num_elements; i++)
	PV.NH[i] = fac * PV.Rho[i];
    }
  else if(block == IO_VOL)
    {
      PV.Vol = mymalloc_movable(&PV.Vol, "PV.Vol", bytes);
      PV.Hsml = mymalloc_movable(&PV.Hsml, "PV.Hsml", bytes);

      if(All.LengthUnit == 0)
	fac = 1;
      else if(All.LengthUnit == 1)
	fac = 1e9;
      else if(All.LengthUnit == 2)
	fac = pow(UNIT_LENGTH / ASTRONOMICAL_UNIT, 3);
      else
	terminate("Length unit not implemented!");

      if(!All.ComovingUnits)
	fac /= pow(1 + All.RedShift, 3);

      fac /= pow(All.HubbleParam, 3);

      for(i = 0; i < num_elements; i++)
	PV.Vol[i] = fac * CONVERT(PD[block], i);

      fac = pow(3. / 4 / M_PI, 1. / 3);

      for(i = 0; i < num_elements; i++){
	PV.Hsml[i] = fac * pow(PV.Vol[i], 1. / 3);
      }
    }
/*  else if(block == IO_DELAUNAY)
    {
      PV.Delaunay = mymalloc_movable(&PV.Delaunay, "PV.Delaunay", bytes);

      if(All.LengthUnit == 0)
	fac = 1;
      else if(All.LengthUnit == 1)
	fac = 1e3;
      else if(All.LengthUnit == 2)
	fac = UNIT_LENGTH / ASTRONOMICAL_UNIT;
      else
	terminate("Length unit not implemented!");

      if(!All.ComovingUnits)
	fac /= (1 + All.RedShift);

      fac /= All.HubbleParam;

      for(i = 0; i < num_elements; i++)
	PV.Delaunay[i] = fac * CONVERT(PD[block], i);
    }*/
  else if(block == IO_GRAVACC)
    {
      PV.GravAccX = mymalloc_movable(&PV.GravAccX, "PV.GravAccX", bytes / num_dims);
      PV.GravAccY = mymalloc_movable(&PV.GravAccY, "PV.GravAccY", bytes / num_dims);
      PV.GravAccZ = mymalloc_movable(&PV.GravAccZ, "PV.GravAccZ", bytes / num_dims);

      fac = pow(1 + All.RedShift, -2);

      for(i = 0; i < num_elements / num_dims; i++)
	{
	  PV.GravAccX[i] = fac * CONVERT(PD[block], num_dims * i);
	  PV.GravAccY[i] = fac * CONVERT(PD[block], num_dims * i + 1);
	  PV.GravAccZ[i] = fac * CONVERT(PD[block], num_dims * i + 2);
	}
    }
  else if(block == IO_GRADP)
    {
      PV.GradPX = mymalloc_movable(&PV.GradPX, "PV.GradPX", bytes / num_dims);
      PV.GradPY = mymalloc_movable(&PV.GradPY, "PV.GradPY", bytes / num_dims);
      PV.GradPZ = mymalloc_movable(&PV.GradPZ, "PV.GradPZ", bytes / num_dims);

      fac = pow(1 + All.RedShift, 4) * pow(All.HubbleParam, 3);

      for(i = 0; i < num_elements / num_dims; i++)
	{
	  PV.GradPX[i] = fac * CONVERT(PD[block], num_dims * i);
	  PV.GradPY[i] = fac * CONVERT(PD[block], num_dims * i + 1);
	  PV.GradPZ[i] = fac * CONVERT(PD[block], num_dims * i + 2);
	}
    }
  else if(block == IO_CHEM)
    {
#ifdef OLDCHEM
      PV.AbH2 = mymalloc_movable(&PV.AbH2, "PV.AbH2", bytes / num_dims);
      PV.AbHII = mymalloc_movable(&PV.AbHII, "PV.AbHII", bytes / num_dims);
      PV.AbDII = mymalloc_movable(&PV.AbDII, "PV.AbDII", bytes / num_dims);
      PV.AbHD = mymalloc_movable(&PV.AbHD, "PV.AbHD", bytes / num_dims);
      PV.AbHeII = mymalloc_movable(&PV.AbHeII, "PV.AbHeII", bytes / num_dims);
      PV.AbHeIII = mymalloc_movable(&PV.AbHeIII, "PV.AbHeIII", bytes / num_dims);
#else
      PV.AbHM = mymalloc_movable(&PV.AbHM, "PV.AbHM", bytes / num_dims);
      PV.AbH2 = mymalloc_movable(&PV.AbH2, "PV.AbH2", bytes / num_dims);
      PV.AbHII = mymalloc_movable(&PV.AbHII, "PV.AbHII", bytes / num_dims);
#endif
      PV.Mu = mymalloc_movable(&PV.Mu, "PV.Mu", bytes / num_dims);

      for(i = 0; i < num_elements / num_dims; i++)
	{
#ifdef OLDCHEM
	  PV.AbH2[i] = CONVERT(PD[block], num_dims * i);
	  PV.AbHII[i] = CONVERT(PD[block], num_dims * i + 1);
	  PV.AbDII[i] = CONVERT(PD[block], num_dims * i + 2);
	  PV.AbHD[i] = CONVERT(PD[block], num_dims * i + 3);
	  PV.AbHeII[i] = CONVERT(PD[block], num_dims * i + 4);
	  PV.AbHeIII[i] = CONVERT(PD[block], num_dims * i + 5);

	  abe = PV.AbHII[i] + PV.AbDII[i] + PV.AbHeII[i] + 2. * PV.AbHeIII[i];
#else
	  PV.AbHM[i] = CONVERT(PD[block], num_dims * i);
	  PV.AbH2[i] = CONVERT(PD[block], num_dims * i + 1);
	  PV.AbHII[i] = CONVERT(PD[block], num_dims * i + 2);

	  abe = PV.AbHII[i];
#endif
	  PV.Mu[i] = (1. + 4. * HE_ABUND) / (1 + HE_ABUND - PV.AbH2[i] + abe);
	}
    }
  else if(block == IO_GAMMA)
    {
      PV.Gamma = mymalloc_movable(&PV.Gamma, "PV.Gamma", bytes);

      for(i = 0; i < num_elements; i++)
	PV.Gamma[i] = CONVERT(PD[block], i);
    }
  else if(block == IO_RATES)
    {
      PV.PdVRate = mymalloc_movable(&PV.PdVRate, "PV.PdVRate", bytes / num_dims);
      PV.H2Rate = mymalloc_movable(&PV.H2Rate, "PV.H2Rate", bytes / num_dims);
      PV.CIERate = mymalloc_movable(&PV.CIERate, "PV.CIERate", bytes / num_dims);
      PV.ChemRate = mymalloc_movable(&PV.ChemRate, "PV.ChemRate", bytes / num_dims);

      for(i = 0; i < num_elements / num_dims; i++)
	{
	  PV.PdVRate[i] = CONVERT(PD[block], num_dims * i);
	  PV.H2Rate[i] = CONVERT(PD[block], num_dims * i + 1);
	  PV.CIERate[i] = CONVERT(PD[block], num_dims * i + 2);
	  PV.ChemRate[i] = CONVERT(PD[block], num_dims * i + 3);
	}
    }
  else if(block == IO_ALLOWREF)
    {
      PV.AllowRef = mymalloc_movable(&PV.AllowRef, "PV.AllowRef", bytes);

      for(i = 0; i < num_elements; i++)
	PV.AllowRef[i] = *((int *) PD[block] + i);
    }
  else if(block == IO_DIVVEL)
    {
      PV.DivVel = mymalloc_movable(&PV.DivVel, "PV.DivVel", bytes);

      fac = (1 + All.RedShift) * All.HubbleParam / (UNIT_LENGTH / UNIT_VELOCITY);

      for(i = 0; i < num_elements; i++)
	PV.DivVel[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_SIDM_PSUM)
    {
      PV.SIDM_PSum = mymalloc_movable(&PV.SIDM_PSum, "PV.SIDM_PSum", bytes);

      fac = 1.;

      for(i = 0; i < num_elements; i++)
	PV.SIDM_PSum[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_SIDM_NUMNGB)
    {
      PV.SIDM_NumNgb = mymalloc_movable(&PV.SIDM_NumNgb, "PV.SIDM_NumNgb", bytes);

      fac = 1.;

      for(i = 0; i < num_elements; i++)
	PV.SIDM_NumNgb[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_SIDM_NUMTOTALSCATTER)
    {
      PV.SIDM_NumTotalScatter = mymalloc_movable(&PV.SIDM_NumTotalScatter, "PV.SIDM_NumTotalScatter", bytes);

      for(i = 0; i < num_elements; i++)
	PV.SIDM_NumTotalScatter[i] = *((int *) PD[block] + i);
    }
  else if(block == IO_SIDM_HSML)
    {
      PV.SIDM_Hsml = mymalloc_movable(&PV.SIDM_Hsml, "PV.SIDM_Hsml", bytes);

      if(All.LengthUnit == 0)
	fac = 1;
      else if(All.LengthUnit == 1)
	fac = 1e3;
      else if(All.LengthUnit == 2)
	fac = UNIT_LENGTH / ASTRONOMICAL_UNIT;
      else
	terminate("Length unit not implemented!");

      if(!All.ComovingUnits)
	fac /= (1 + All.RedShift);

      fac /= All.HubbleParam;

      for(i = 0; i < num_elements; i++)
	PV.SIDM_Hsml[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_SIDM_DENSITY)
    {
      PV.SIDM_Density = mymalloc_movable(&PV.SIDM_Density, "PV.SIDM_Density", bytes);

      fac = UNIT_MASS / pow(UNIT_LENGTH, 3) * pow(1 + All.RedShift, 3) * pow(All.HubbleParam, 2);

      fac *= HYDROGEN_MASSFRAC / PROTONMASS;

      for(i = 0; i < num_elements; i++)
	PV.SIDM_Density[i] = fac * CONVERT(PD[block], i);
    }
  else if(block == IO_SIDM_VELDISP)
    {
      PV.SIDM_VelDisp = mymalloc_movable(&PV.SIDM_VelDisp, "PV.SIDM_VelDisp", bytes);

      fac = 1.;

      for(i = 0; i < num_elements; i++)
	PV.SIDM_VelDisp[i] = fac * CONVERT(PD[block], i);
    }
  else
    terminate("Unknown IO block!");

  myfree_movable(PD[block]);

  PD[block] = 0;
}


void calc_temp()
{
  int i;
  double mu, gamma;

  if(BlockFlag[IO_U])
    {
      for(i = 0; i < NumGas; i++)
	{
	  if(BlockFlag[IO_CHEM])
	    mu = PV.Mu[i];
	  else
	    mu = MU_NEUTRAL;

	  if(BlockFlag[IO_GAMMA])
	    gamma = PV.Gamma[i];
	  else
	    gamma = GAMMA_ADB;

	  PV.Temp[i] = mu * (gamma - 1) * PROTONMASS / BOLTZMANN * PV.U[i];
	}
    }
}


void clear_snap()
{
  int block;

  for(block = 0; block < NUM_BLOCKS; block++)
    if(BlockFlag[block])
      {
	if(block == IO_POS)
	  {
	    myfree_movable(PV.X);
	    myfree_movable(PV.Y);
	    myfree_movable(PV.Z);

	    PV.X = PV.Y = PV.Z = 0;
	  }
	else if(block == IO_VEL)
	  {
	    myfree_movable(PV.VX);
	    myfree_movable(PV.VY);
	    myfree_movable(PV.VZ);

	    PV.VX = PV.VY = PV.VZ = 0;
	  }
	else if(block == IO_ID)
	  {
	    myfree_movable(PV.ID);

	    PV.ID = 0;
	  }
	else if(block == IO_MASS)
	  {
	    myfree_movable(PV.Mass);

	    PV.Mass = 0;
	  }
	else if(block == IO_U)
	  {
	    myfree_movable(PV.U);
	    myfree_movable(PV.Temp);

	    PV.U = PV.Temp = 0;
	  }
	else if(block == IO_RHO)
	  {
	    myfree_movable(PV.Rho);
	    myfree_movable(PV.NH);

	    PV.Rho = PV.NH = 0;
	  }
	else if(block == IO_VOL)
	  {
	    myfree_movable(PV.Vol);
	    myfree_movable(PV.Hsml);

	    PV.Vol = PV.Hsml = 0;
	  }
/*	else if(block == IO_DELAUNAY)
	  {
	    myfree_movable(PV.Delaunay);

	    PV.Delaunay = 0;
	  }*/
	else if(block == IO_GRAVACC)
	  {
	    myfree_movable(PV.GravAccX);
	    myfree_movable(PV.GravAccY);
	    myfree_movable(PV.GravAccZ);

	    PV.GravAccX = PV.GravAccY = PV.GravAccZ = 0;
	  }
	else if(block == IO_GRADP)
	  {
	    myfree_movable(PV.GradPX);
	    myfree_movable(PV.GradPY);
	    myfree_movable(PV.GradPZ);

	    PV.GradPX = PV.GradPY = PV.GradPZ = 0;
	  }
	else if(block == IO_CHEM)
	  {
#ifdef OLDCHEM
	    myfree_movable(PV.AbH2);
	    myfree_movable(PV.AbHII);
	    myfree_movable(PV.AbDII);
	    myfree_movable(PV.AbHD);
	    myfree_movable(PV.AbHeII);
	    myfree_movable(PV.AbHeIII);

	    PV.AbH2 = PV.AbHII = PV.AbDII = PV.AbHD = PV.AbHeII = PV.AbHeIII = 0;
#else
	    myfree_movable(PV.AbHM);
	    myfree_movable(PV.AbH2);
	    myfree_movable(PV.AbHII);

	    PV.AbHM = PV.AbH2 = PV.AbHII = 0;
#endif
	    myfree_movable(PV.Mu);

	    PV.Mu = 0;
	  }
	else if(block == IO_GAMMA)
	  {
	    myfree_movable(PV.Gamma);

	    PV.Gamma = 0;
	  }
	else if(block == IO_RATES)
	  {
	    myfree_movable(PV.PdVRate);
	    myfree_movable(PV.H2Rate);
	    myfree_movable(PV.CIERate);
	    myfree_movable(PV.ChemRate);

	    PV.PdVRate = PV.H2Rate = PV.CIERate = PV.ChemRate = 0;
	  }
	else if(block == IO_ALLOWREF)
	  {
	    myfree_movable(PV.AllowRef);

	    PV.AllowRef = 0;
	  }
	else if(block == IO_DIVVEL)
	  {
	    myfree_movable(PV.DivVel);

	    PV.DivVel = 0;
	  }
	else if(block == IO_SIDM_PSUM)
	  {
	    myfree_movable(PV.SIDM_PSum);

	    PV.SIDM_PSum = 0;
	  }
	else if(block == IO_SIDM_NUMNGB)
	  {
	    myfree_movable(PV.SIDM_NumNgb);

	    PV.SIDM_NumNgb = 0;
	  }
	else if(block == IO_SIDM_NUMTOTALSCATTER)
	  {
	    myfree_movable(PV.SIDM_NumTotalScatter);

	    PV.SIDM_NumTotalScatter = 0;
	  }
	else if(block == IO_SIDM_HSML)
	  {
	    myfree_movable(PV.SIDM_Hsml);

	    PV.SIDM_Hsml = 0;
	  }
	else if(block == IO_SIDM_DENSITY)
	  {
	    myfree_movable(PV.SIDM_Density);

	    PV.SIDM_Density = 0;
	  }
	else if(block == IO_SIDM_VELDISP)
	  {
	    myfree_movable(PV.SIDM_VelDisp);

	    PV.SIDM_VelDisp = 0;
	  }
	else
	  terminate("Unknown IO block!");
      }

  myfree(All.TaskNumPart);
  myfree(All.SnapNumPart);
}
