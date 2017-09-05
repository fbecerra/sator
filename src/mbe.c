#include "allvars.h"
#include "proto.h"


void mbe()
{
  char mbe_file[MAX_STRING_LEN], label[MAX_STRING_LEN];
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

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] += All.ShiftX;
      PV.Y[i] += All.ShiftY;
      PV.Z[i] += All.ShiftZ;
    }

  dist = sqrt(pow(All.ShiftX , 2) + pow(All.ShiftY, 2) + pow(All.ShiftZ, 2));

  num_blocks = 2;

  r = (double *) mymalloc_movable(&r, "r", NumGas * sizeof(double));

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

  if(All.MBEMin)
    r_min = All.MBEMin;

  log_r_min = log10(r_min);
  log_r_max = log10(r_max);

  rfac = All.MBEBins / (log_r_max - log_r_min);

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

  rad = (double *) mymalloc_movable(&rad, "rad", All.MBEBins * sizeof(double));

  for(i = 0; i < All.MBEBins; i++)
    rad[i] = pow(10., log_r_min + i / rfac);

  loc_mass = (double *) mymalloc_movable(&loc_mass, "loc_mass", All.MBEBins * sizeof(double));
  mass = (double *) mymalloc_movable(&mass, "mass", All.MBEBins * sizeof(double));
  enc = (double *) mymalloc_movable(&enc, "enc", All.MBEBins * sizeof(double));
  loc_mbe = (double *) mymalloc_movable(&loc_mbe, "loc_mbe", All.MBEBins * sizeof(double));
  mbe = (double *) mymalloc_movable(&mbe, "mbe", All.MBEBins * sizeof(double));
  tot_mbe = (double *) mymalloc_movable(&tot_mbe, "tot_mbe", All.MBEBins * sizeof(double));
  mbe_ratio = (double *) mymalloc_movable(&mbe_ratio, "mbe_ratio", All.MBEBins * sizeof(double));
  loc_vel_clump_x = (double *) mymalloc_movable(&loc_vel_clump_x, "loc_vel_clump_x", All.MBEBins * sizeof(double));
  rad_vel_clump_x = (double *) mymalloc_movable(&rad_vel_clump_x, "rad_vel_clump_x", All.MBEBins * sizeof(double));
  vel_clump_x = (double *) mymalloc_movable(&vel_clump_x, "vel_clump_x", All.MBEBins * sizeof(double));
  loc_vel_clump_y = (double *) mymalloc_movable(&loc_vel_clump_y, "loc_vel_clump_y", All.MBEBins * sizeof(double));
  rad_vel_clump_y = (double *) mymalloc_movable(&rad_vel_clump_y, "rad_vel_clump_y", All.MBEBins * sizeof(double));
  vel_clump_y = (double *) mymalloc_movable(&vel_clump_y, "vel_clump_y", All.MBEBins * sizeof(double));
  loc_vel_clump_z = (double *) mymalloc_movable(&loc_vel_clump_z, "loc_vel_clump_z", All.MBEBins * sizeof(double));
  rad_vel_clump_z = (double *) mymalloc_movable(&rad_vel_clump_z, "rad_vel_clump_z", All.MBEBins * sizeof(double));
  vel_clump_z = (double *) mymalloc_movable(&vel_clump_z, "vel_clump_z", All.MBEBins * sizeof(double));
  acc_ratio = (double *) mymalloc_movable(&acc_ratio, "acc_ratio", All.MBEBins * sizeof(double));

  memset(loc_mass, 0, All.MBEBins * sizeof(double));
  memset(loc_mbe, 0, All.MBEBins * sizeof(double));
  memset(loc_vel_clump_x, 0, All.MBEBins * sizeof(double));
  memset(loc_vel_clump_y, 0, All.MBEBins * sizeof(double));
  memset(loc_vel_clump_z, 0, All.MBEBins * sizeof(double));

  tot_mass = loc_vel_cm_x = loc_vel_cm_y = loc_vel_cm_z = 0;

  for(i = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      j = r_idx[i];

      if(j < 0 || j > All.MBEBins - 1)
	continue;

      loc_mass[j] += PV.Mass[i];

      mass_be = 2.6e-11 * pow(PV.Temp[i], 3. / 2) * pow(PV.Rho[i], -1. / 2) * pow(PV.Mu[i], -3. / 2) * pow(PV.Gamma[i], 2);

      loc_mbe[j] += PV.Mass[i] * mass_be;

      tot_mass += PV.Mass[i];

      loc_vel_cm_x += PV.Mass[i] * PV.VX[i];
      loc_vel_cm_y += PV.Mass[i] * PV.VY[i];
      loc_vel_cm_z += PV.Mass[i] * PV.VZ[i];
    }

  MPI_Allreduce(loc_mass, mass, All.MBEBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_mbe, mbe, All.MBEBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_vel_cm_x, &vel_cm_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_vel_cm_y, &vel_cm_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_vel_cm_z, &vel_cm_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  vel_cm_x /= tot_mass;
  vel_cm_y /= tot_mass;
  vel_cm_z /= tot_mass;

  for(i = 0; i < All.MBEBins; i++)
    if(!mass[i])
      mpi_printf("Warning! Bin %d without mass!\n", i);

  for(i = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      j = r_idx[i];

      if(j < 0 || j > All.MBEBins - 1)
	continue;

      loc_vel_clump_x[j] += PV.Mass[i] * (PV.VX[i] - vel_cm_x);
      loc_vel_clump_y[j] += PV.Mass[i] * (PV.VY[i] - vel_cm_y);
      loc_vel_clump_z[j] += PV.Mass[i] * (PV.VZ[i] - vel_cm_z);
    }

  MPI_Allreduce(loc_vel_clump_x, rad_vel_clump_x, All.MBEBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_vel_clump_y, rad_vel_clump_y, All.MBEBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_vel_clump_z, rad_vel_clump_z, All.MBEBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < All.MBEBins; i++)
    {
      if(i == 0)
	{
	  enc[i] = mass[i];
	  tot_mbe[i] = mbe[i];
	  vel_clump_x[i] = rad_vel_clump_x[i];
	  vel_clump_y[i] = rad_vel_clump_y[i];
	  vel_clump_z[i] = rad_vel_clump_z[i];
	}
      else
	{
	  enc[i] = enc[i - 1] + mass[i];
	  tot_mbe[i] = tot_mbe[i - 1] + mbe[i];
	  vel_clump_x[i] = vel_clump_x[i - 1] + rad_vel_clump_x[i];
	  vel_clump_y[i] = vel_clump_y[i - 1] + rad_vel_clump_y[i];
	  vel_clump_z[i] = vel_clump_z[i - 1] + rad_vel_clump_z[i];
	}

      //mbe_ratio[i] = enc[i] / (mbe[i] / mass[i]);
      mbe_ratio[i] = enc[i] / (tot_mbe[i] / enc[i]);

      rho_avg = 3. / 4. / M_PI * enc[i] * SOLAR_MASS / pow(fac * rad[i], 3);

      tff = sqrt(3 * M_PI / 32 / GRAVITY / rho_avg);

      vrad = dabs(All.ShiftX * vel_clump_x[i] + All.ShiftY * vel_clump_y[i] + All.ShiftZ * vel_clump_z[i]) / enc[i] / dist;

      tacc = fac * dist / (vrad * 1e5);

      acc_ratio[i] = tff / tacc;
    }

  dist = log10(dist);

  for(i = 0; i < All.MBEBins; i++)
    {
      rad[i] = log10(rad[i]);
      enc[i] = log10(enc[i]);
      mbe_ratio[i] = log10(mbe_ratio[i]);
      acc_ratio[i] = log10(acc_ratio[i]);
    }

  if(ThisTask == 0)
    {
      sprintf(mbe_file, "%s/%s/%s_%s.mbe", All.Path, All.Base, All.Base, All.SnapNumString);

      if(!(file = fopen(mbe_file, "w")))
	terminate("Could not write mbe file!");

      fwrite(&num_blocks, sizeof(int), 1, file);
      fwrite(&max_string_len, sizeof(int), 1, file);
      fwrite(&All.MBEBins, sizeof(int), 1, file);
      fwrite(&All.LengthUnit, sizeof(int), 1, file);
      fwrite(&All.ComovingUnits, sizeof(int), 1, file);
      fwrite(&dist, sizeof(double), 1, file);
      fwrite(&log_r_min, sizeof(double), 1, file);
      fwrite(&rad[All.MBEBins - 1], sizeof(double), 1, file);
      fwrite(rad, All.MBEBins * sizeof(double), 1, file);

      memset(label, 0, MAX_STRING_LEN * sizeof(char));
      strcpy(label, "mbe_ratio");
      fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
      fwrite(mbe_ratio, All.MBEBins * sizeof(double), 1, file);

      memset(label, 0, MAX_STRING_LEN * sizeof(char));
      strcpy(label, "acc_ratio");
      fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
      fwrite(acc_ratio, All.MBEBins * sizeof(double), 1, file);

      fclose(file);
    }

  myfree_movable(acc_ratio);
  myfree_movable(vel_clump_z);
  myfree_movable(rad_vel_clump_z);
  myfree_movable(loc_vel_clump_z);
  myfree_movable(vel_clump_y);
  myfree_movable(rad_vel_clump_y);
  myfree_movable(loc_vel_clump_y);
  myfree_movable(vel_clump_x);
  myfree_movable(rad_vel_clump_x);
  myfree_movable(loc_vel_clump_x);
  myfree_movable(mbe_ratio);
  myfree_movable(tot_mbe);
  myfree_movable(mbe);
  myfree_movable(loc_mbe);
  myfree_movable(enc);
  myfree_movable(mass);
  myfree_movable(loc_mass);
  myfree_movable(rad);
  myfree_movable(r_idx);
  myfree_movable(r);
}
