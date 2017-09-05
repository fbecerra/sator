#include "allvars.h"
#include "proto.h"


void disk()
{
  char disk_file[MAX_STRING_LEN], label[MAX_STRING_LEN];
  int block_toomre, block_gammie;
  int i, j, max_string_len, *r_idx, num_blocks;
  double fac, *r, *rad, rfac, loc_r_min, r_min, log_r_min, loc_r_max, r_max, log_r_max;
  double *loc_mass, *mass, *enc, mass_be, *loc_mbe, *mbe, *tot_mbe, *mbe_ratio, dist;
  double tot_mass, loc_vel_cm_x, vel_cm_x, loc_vel_cm_y, vel_cm_y, loc_vel_cm_z, vel_cm_z;
  double *loc_vel_clump_x, *loc_vel_clump_y, *loc_vel_clump_z;
  double *rad_vel_clump_x, *rad_vel_clump_y, *rad_vel_clump_z;
  double *vel_clump_x, *vel_clump_y, *vel_clump_z;
  double vrad, tacc, tff, rho_avg, *acc_ratio;
  FILE *file;

  max_string_len = MAX_STRING_LEN;

  if(!BlockFlag[IO_POS] || !BlockFlag[IO_VEL] || !BlockFlag[IO_MASS] || !BlockFlag[IO_ID] || !BlockFlag[IO_U] || !BlockFlag[IO_RHO] || !BlockFlag[IO_VOL])
    terminate("Required block not read!");

  block_toomre = 1;
  block_gammie = 0;

  if(All.LengthUnit == 0)
    fac = UNIT_LENGTH;
  else if(All.LengthUnit == 1)
    fac = UNIT_LENGTH / 1e3;
  else if(All.LengthUnit == 2)
    fac = ASTRONOMICAL_UNIT;
  else
    terminate("Length unit not implemented!");

  if(All.ComovingUnits)
    fac /= (1 + All.RedShift);

  box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
    }

  if(All.FlagRotate)
    box_rotate();

  num_blocks = 0;

  if(block_toomre)
     num_blocks += 4;

  if(block_gammie)
     num_blocks += 4;

  r = (double *) mymalloc_movable(&r, "r", TotNumPart * sizeof(double));

  for(i = loc_r_max = 0, loc_r_min = DBL_MAX; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
        continue;

      r[i] = sqrt(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i]);

      loc_r_min = dmin(PV.Hsml[i], loc_r_min);
      loc_r_max = dmax(r[i], loc_r_max);
    }

  MPI_Allreduce(&loc_r_min, &r_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_r_max, &r_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(All.DiskMin)
    r_min = All.DiskMin;

  log_r_min = log10(r_min);
  log_r_max = log10(r_max);

  rfac = All.DiskBins / (log_r_max - log_r_min);

  r_idx = (int *) mymalloc_movable(&r_idx, "r_idx", NumGas * sizeof(int));

  for(i = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
        continue;

      if(r[i] == 0)
        r_idx[i] = -1;
      else
        r_idx[i] = rfac * (log10(r[i]) - log_r_min);
    }

  rad = (double *) mymalloc_movable(&rad, "rad", All.DiskBins * sizeof(double));

  for(i = 0; i < All.DiskBins; i++)
    rad[i] = pow(10., log_r_min + i / rfac);


  /* Define variable */

  loc_mass = (double *) mymalloc_movable(&loc_mass, "loc_mass", All.DiskBins * sizeof(double));
  mass = (double *) mymalloc_movable(&mass, "mass", All.DiskBins * sizeof(double));

/*  if(block_toomre){
  }

  if(block_gammie){
  } */


  for(i = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
        continue;

      j = r_idx[i];

      if(j < 0 || j > All.DiskBins - 1)
        continue;

      loc_mass[j] += PV.Mass[i];

      loc_vel_cm_x += PV.Mass[i] * PV.VX[i];
      loc_vel_cm_y += PV.Mass[i] * PV.VY[i];
      loc_vel_cm_z += PV.Mass[i] * PV.VZ[i];

      loc_nh[j] += PV.Mass[i] * log10(PV.NH[i]);

      loc_temp[j] += PV.Mass[i] * log10(PV.Temp[i]);

      if(block_toomre)
        {
        }

      if(block_gammie)
        {
        }
    }

  MPI_Allreduce(loc_mass, mass, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_vel_cm_x, vel_cm_x, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_vel_cm_y, vel_cm_y, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_vel_cm_z, vel_cm_z, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_nh, nh, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_temp, temp, All.DiskBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(BlockFlag[IO_VEL])
    {
      if(All.LengthUnit == 0)
        fac = UNIT_LENGTH;
      else if(All.LengthUnit == 1)
        fac = UNIT_LENGTH / 1e3;
      else if(All.LengthUnit == 2)
        fac = ASTRONOMICAL_UNIT;
      else
        terminate("Length unit not implemented!");

      if(All.ComovingUnits)
        fac /= (1 + All.RedShift);

      for(i = 0; i < NumGas; i++)
        {
          if(PV.Mass[i] == 0 && PV.ID[i] == 0)
            continue;

          j = r_idx[i];

          if(j < 0 || j > All.RadBins - 1)
            continue;

          vradx = PV.X[i] * (PV.VX[i] - vel_cm_x[j]);
          vrady = PV.Y[i] * (PV.VY[i] - vel_cm_y[j]);
          vradz = PV.Z[i] * (PV.VZ[i] - vel_cm_z[j]);

          loc_vel_rad[j] += PV.Mass[i] * (vradx + vrady + vradz) / r[i];

          loc_angmom_x[j] += PV.Mass[i] * (PV.Y[i] * (PV.VZ[i] - vel_cm_z[j]) - PV.Z[i] * (PV.VY[i] - vel_cm_y[j]));
          loc_angmom_y[j] += PV.Mass[i] * (PV.Z[i] * (PV.VX[i] - vel_cm_x[j]) - PV.X[i] * (PV.VZ[i] - vel_cm_z[j]));
          loc_angmom_z[j] += PV.Mass[i] * (PV.X[i] * (PV.VY[i] - vel_cm_y[j]) - PV.Y[i] * (PV.VX[i] - vel_cm_x[j]));

          loc_csnd[j] += PV.Mass[i] * log10(sqrt(PV.Gamma[i] * BOLTZMANN * PV.Temp[i] / PV.Mu[i] / PROTONMASS) / 1e5);
        }

      MPI_Allreduce(loc_vel_rad, vel_rad, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(loc_angmom_x, angmom_x, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_angmom_y, angmom_y, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_angmom_z, angmom_z, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(loc_csnd, csnd, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
        {
          vel_ang[i] = sqrt(angmom_x[i] * angmom_x[i] + angmom_y[i] * angmom_y[i] + angmom_z[i] * angmom_z[i]) / rad[i] / rad[i] / mass[i];

          t_orb[i] = 2 * M_PI / vel_ang[i];

          
        }
