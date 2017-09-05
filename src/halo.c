
#include "allvars.h"
#include "proto.h"


void halo()
{
  int i, iter;
  double r_test, r2_test, r3_test, r_low, r_high;
  double fac, rho_m, ratio_fac, loc_mvir, mvir, mvir_old, func;
  double loc_vel_cm_x, vel_cm_x, loc_vel_cm_y, vel_cm_y, loc_vel_cm_z, vel_cm_z;
  double loc_angmom_x, angmom_x, loc_angmom_y, angmom_y, loc_angmom_z, angmom_z;
  double angmom, bind, lambda;

  if(!BlockFlag[IO_POS] || !BlockFlag[IO_MASS] || !BlockFlag[IO_ID])
    terminate("Required block not read!");

  box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
    }

  if(All.LengthUnit == 0)
    fac = 1;
  else if(All.LengthUnit == 1)
    fac = 1e-3;
  else if(All.LengthUnit == 2)
    fac = ASTRONOMICAL_UNIT / UNIT_LENGTH;
  else
    terminate("Length unit not implemented!");

  if(!All.ComovingUnits)
    fac *= (1 + All.RedShift);

  fac *= All.HubbleParam;

  for(i = 0; i < 3; i++)
    All.BoxCenter[i] *= fac;

  mpi_printf("X = %g, Y = %g Z = %g\n", All.BoxCenter[0], All.BoxCenter[1], All.BoxCenter[2]);

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

  All.BoxSize *= fac;

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

  rho_m = All.OmegaM * 3 * pow(HUBBLE * All.HubbleParam, 2) / 8 / M_PI / GRAVITY * pow(1 + All.RedShift, 3);

  ratio_fac = 200 * rho_m / SOLAR_MASS * 4 / 3 * M_PI * pow(fac, 3);

  r_low = 0;
  r_high = All.BoxSize;

  r_test = (r_low + r_high) / 2;
  r2_test = r_test * r_test;
  r3_test = r_test * r2_test;

  iter = 0;

  while(r_test != r_low && r_test != r_high)
    {
      mvir_old = mvir;

      for(i = loc_mvir = 0; i < TotNumPart; i++)
	{
	  if(i < NumGas && PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i] < r2_test)
	    loc_mvir += PV.Mass[i];
	}
  
      MPI_Allreduce(&loc_mvir, &mvir, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      func = ratio_fac / mvir * r3_test - 1;

      if(func >= 0)
	r_high = r_test;
      else
	r_low = r_test;
      
      r_test = (r_low + r_high) / 2;
      r2_test = r_test * r_test;
      r3_test = r_test * r2_test;

      if(iter != 0 && dabs(mvir - mvir_old) < 1e-3 * mvir_old)
	break;

      iter++;
    }

  for(i = loc_vel_cm_x = loc_vel_cm_y = loc_vel_cm_z = 0; i < TotNumPart; i++)
    {
      if(i < NumGas && PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      if(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i] < r2_test)
	{
	  loc_vel_cm_x += PV.Mass[i] * PV.VX[i];
	  loc_vel_cm_y += PV.Mass[i] * PV.VY[i];
	  loc_vel_cm_z += PV.Mass[i] * PV.VZ[i];
	}
    }

  MPI_Allreduce(&loc_vel_cm_x, &vel_cm_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_vel_cm_y, &vel_cm_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_vel_cm_z, &vel_cm_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  vel_cm_x /= mvir;
  vel_cm_y /= mvir;
  vel_cm_z /= mvir;

  for(i = loc_angmom_x = loc_angmom_y = loc_angmom_z = 0; i < TotNumPart; i++)
    {
      if(i < NumGas && PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      if(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i] < r2_test)
	{
	  loc_angmom_x += PV.Mass[i] * (PV.Y[i] * (PV.VZ[i] - vel_cm_z) - PV.Z[i] * (PV.VY[i] - vel_cm_y));
	  loc_angmom_y += PV.Mass[i] * (PV.Z[i] * (PV.VX[i] - vel_cm_x) - PV.X[i] * (PV.VZ[i] - vel_cm_z));
	  loc_angmom_z += PV.Mass[i] * (PV.X[i] * (PV.VY[i] - vel_cm_y) - PV.Y[i] * (PV.VX[i] - vel_cm_x));
	}
    }

  MPI_Allreduce(&loc_angmom_x, &angmom_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_angmom_y, &angmom_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_angmom_z, &angmom_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  angmom = sqrt(pow(angmom_x, 2) + pow(angmom_y, 2) + pow(angmom_z, 2));

  angmom *= SOLAR_MASS * UNIT_VELOCITY * fac;

  r_test *= fac;

  mvir *= SOLAR_MASS;

  bind = 3. / 5 * GRAVITY * mvir * mvir / r_test;

  lambda = angmom * sqrt(bind) / GRAVITY / pow(mvir, 5. / 2);

  mpi_printf("Redshift: %g\n", All.RedShift);
  mpi_printf("Virial mass: %g M_sun\n", mvir / SOLAR_MASS);
  mpi_printf("Virial radius: %g pc\n", r_test * 1e3 / UNIT_LENGTH);
  mpi_printf("Spin parameter: %g\n", lambda);
}
