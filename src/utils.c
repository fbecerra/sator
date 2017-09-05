
#include "allvars.h"
#include "proto.h"


void box_center()
{
  int i;
  double fac, sum, glob_sum, *valp, center[3], glob_center[3];

  struct
  {
    double search_val;
    int task;

  } local, global;

  if(!BlockFlag[IO_POS] || !BlockFlag[IO_ID] || !BlockFlag[IO_MASS] || !BlockFlag[IO_RHO])
    terminate_block();

  if(All.FlagCenter == 0)
    {
      if(All.Usage == 0)
	{
	  for(i = center[0] = center[1] = center[2] = sum = 0; i < NumGas; i++)
	    {
	      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
		continue;

	      center[0] += PV.NH[i] * PV.X[i];
	      center[1] += PV.NH[i] * PV.Y[i];
	      center[2] += PV.NH[i] * PV.Z[i];

	      sum += PV.NH[i];
	    }

	  MPI_Allreduce(center, glob_center, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&sum, &glob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  for(i = 0; i < 3; i++)
	    All.BoxCenter[i] = glob_center[i] / glob_sum;
	}
      else
	{
	  local.task = ThisTask;

	  for(i = 0, local.search_val = -DBL_MAX; i < NumGas; i++)
	    {
	      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
		continue;

	      if(PV.NH[i] > local.search_val)
		{
		  local.search_val = PV.NH[i];

		  center[0] = PV.X[i];
		  center[1] = PV.Y[i];
		  center[2] = PV.Z[i];
		}
	    }

	  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	  MPI_Bcast(center, 3, MPI_DOUBLE, global.task, MPI_COMM_WORLD);

	  for(i = 0; i < 3; i++)
	    All.BoxCenter[i] = center[i];
	  /*
	  for(i = 0; i < NumGas; i++)
	    if(PV.X[i] == All.BoxCenter[0] && PV.Y[i] == All.BoxCenter[1] && PV.Z[i] == All.BoxCenter[2])
	      printf("nh = %g\n", PV.NH[i]);
	  */
	}
    }
  else if(All.FlagCenter == 1)
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

      All.BoxCenter[0] = fac * All.CenterX;
      All.BoxCenter[1] = fac * All.CenterY;
      All.BoxCenter[2] = fac * All.CenterZ;

      mpi_printf("fac = %i, x = %f, y = %f, z = %f \n", fac, All.BoxCenter[0], All.BoxCenter[1], All.BoxCenter[2]);
    }
  else if(All.FlagCenter == 2)
    {
      if(!plot_block_present(PLOT_COLLAPSE))
	terminate_block();

      valp = plot_block_val(PLOT_COLLAPSE);

      for(i = center[0] = center[1] = center[2] = sum = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if(PV.NH[i] > 1e12 && PV.NH[i] < 1e13)
	    {
	      center[0] += PV.Mass[i] / PV.NH[i] / (*(valp + i)) * PV.X[i];
	      center[1] += PV.Mass[i] / PV.NH[i] / (*(valp + i)) * PV.Y[i];
	      center[2] += PV.Mass[i] / PV.NH[i] / (*(valp + i)) * PV.Z[i];

	      sum += PV.Mass[i] / PV.NH[i] / (*(valp + i));
	    }
	}

      myfree_movable(PV.Collapse);

      MPI_Allreduce(center, glob_center, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&sum, &glob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < 3; i++)
	All.BoxCenter[i] = glob_center[i] / glob_sum;
    }
  else if(All.FlagCenter == 3)
    {
      if(!plot_block_present(PLOT_COLLAPSE))
	terminate_block();

      valp = plot_block_val(PLOT_COLLAPSE);

      local.task = ThisTask;

      for(i = 0, local.search_val = DBL_MAX; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if((*(valp + i)) < local.search_val)
	    {
	      local.search_val = (*(valp + i));

	      center[0] = PV.X[i];
	      center[1] = PV.Y[i];
	      center[2] = PV.Z[i];
	    }
	}

      //printf("\n", ThisTask, local.search_val, global.search_val, global.task);

      myfree_movable(PV.Collapse);

      MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

      MPI_Bcast(center, 3, MPI_DOUBLE, global.task, MPI_COMM_WORLD);

      for(i = 0; i < 3; i++)
	All.BoxCenter[i] = center[i];
    }
  else if(All.FlagCenter == 4)
    {
      i = TotNumPart - NumSinks;
      mpi_printf("Center ID = %i, x = %f, y = %f, z = %f \n", PV.ID[i], PV.X[i], PV.Y[i], PV.Z[i]);

      All.BoxCenter[0] = PV.X[i];
      All.BoxCenter[1] = PV.Y[i];
      All.BoxCenter[2] = PV.Z[i];
    }
  else
    terminate("Unknown center flag!");
}


void box_rotate()
{
  int i, j, k, *part_flag;
  double r, r_search, fac, search_fac;
  double loc_vel_cm_x, loc_vel_cm_y, loc_vel_cm_z, loc_mass_cm;
  double vel_cm_x, vel_cm_y, vel_cm_z, mass_cm;
  double dvel_x, dvel_y, dvel_z;
  double loc_vec_x, loc_vec_y, loc_vec_z;
  double theta, vec_x, vec_y, vec_z, alpha, beta;
  double rot1[3][3], rot2[3][3], rot12[3][3];
  double x_new, y_new, z_new, vx_new, vy_new, vz_new;

  struct
  {
    double nh_max;
    int task;

  } local, global;

  if(!PV.VX)
    terminate("Velocity not read!");

  if(All.FlagRotate == 1)
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

      if(All.Usage == 0)
	r_search = All.ImgSize / 2;
      else
	r_search = fac * All.BoxSize / 2;

  /* if(0) */
  /*   { */
  /*     search_fac = 100; */
  /*     search_fac2 = pow(search_fac, 2); */

  /*     local.task = ThisTask; */

  /*     for(i = local.nh_max = hsml = 0; i < NumGas; i++) */
  /* 	{ */
  /* 	  if(PV.Mass[i] == 0 && PV.ID[i] == 0) */
  /* 	    continue; */

  /* 	  if(PV.NH[i] > local.nh_max) */
  /* 	    { */
  /* 	      local.nh_max = PV.NH[i]; */
  /* 	      hsml = PV.Hsml[i]; */
  /* 	    } */
  /* 	} */

  /*     MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); */

  /*     MPI_Bcast(&hsml, 1, MPI_DOUBLE, global.task, MPI_COMM_WORLD); */

  /*     hsml2 = pow(hsml, 2); */

  /*     r2_search = search_fac2 * hsml2; */
  /*   } */

      part_flag = mymalloc_movable(&part_flag, "part_flag", NumGas * sizeof(int));

      memset(part_flag, 0, NumGas * sizeof(int));

      for(i = loc_vel_cm_x = loc_vel_cm_y = loc_vel_cm_z = loc_mass_cm = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if(dabs(PV.X[i]) < r_search && dabs(PV.Y[i]) < r_search && dabs(PV.Z[i]) < r_search)
	    {
	      part_flag[i] = 1;

	      loc_vel_cm_x += PV.Mass[i] * PV.VX[i];
	      loc_vel_cm_y += PV.Mass[i] * PV.VY[i];
	      loc_vel_cm_z += PV.Mass[i] * PV.VZ[i];

	      loc_mass_cm += PV.Mass[i];
	    }
	}

      MPI_Allreduce(&loc_vel_cm_x, &vel_cm_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vel_cm_y, &vel_cm_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vel_cm_z, &vel_cm_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(&loc_mass_cm, &mass_cm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      vel_cm_x /= mass_cm;
      vel_cm_y /= mass_cm;
      vel_cm_z /= mass_cm;

      for(i = loc_vec_x = loc_vec_y = loc_vec_z = 0; i < NumGas; i++)
	if(part_flag[i])
	  {
	    dvel_x = PV.VX[i] - vel_cm_x;
	    dvel_y = PV.VY[i] - vel_cm_y;
	    dvel_z = PV.VZ[i] - vel_cm_z;

	    loc_vec_x += PV.Mass[i] * (PV.Y[i] * dvel_z - PV.Z[i] * dvel_y);
	    loc_vec_y += PV.Mass[i] * (PV.Z[i] * dvel_x - PV.X[i] * dvel_z);
	    loc_vec_z += PV.Mass[i] * (PV.X[i] * dvel_y - PV.Y[i] * dvel_x);
	  }

      myfree_movable(part_flag);

      MPI_Allreduce(&loc_vec_x, &vec_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vec_y, &vec_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vec_z, &vec_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      r = sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2));

      vec_x /= r;
      vec_y /= r;
      vec_z /= r;
    }
  else if(All.FlagRotate == 2)
    {
      theta = 2 * M_PI / (All.ImgNumSubSnaps - 1);

      vec_x = 0;
      vec_y = sin(theta);
      vec_z = cos(theta);
    }
  else
    terminate("Unknown FlagRotate!");

  alpha = acos(vec_z / sqrt(pow(vec_y, 2) + pow(vec_z, 2)));

  if(vec_y >= 0)
    alpha *= -1;

  beta = asin(vec_x);

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      rot1[i][j] = rot2[i][j] = rot12[i][j] = 0;

  rot1[0][0] = 1;
  rot1[1][1] = cos(alpha);
  rot1[1][2] = sin(alpha);
  rot1[2][1] = -sin(alpha);
  rot1[2][2] = cos(alpha);

  rot2[0][0] = cos(beta);
  rot2[0][2] = -sin(beta);
  rot2[1][1] = 1;
  rot2[2][0] = sin(beta);
  rot2[2][2] = cos(beta);

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      for(k = 0; k < 3; k++)
	rot12[i][j] += rot2[i][k] * rot1[k][j];

  for(i = 0; i < TotNumPart; i++)
    {
      x_new = rot12[0][0] * PV.X[i] + rot12[0][1] * PV.Y[i] + rot12[0][2] * PV.Z[i];
      y_new = rot12[1][0] * PV.X[i] + rot12[1][1] * PV.Y[i] + rot12[1][2] * PV.Z[i];
      z_new = rot12[2][0] * PV.X[i] + rot12[2][1] * PV.Y[i] + rot12[2][2] * PV.Z[i];

      vx_new = rot12[0][0] * PV.VX[i] + rot12[0][1] * PV.VY[i] + rot12[0][2] * PV.VZ[i];
      vy_new = rot12[1][0] * PV.VX[i] + rot12[1][1] * PV.VY[i] + rot12[1][2] * PV.VZ[i];
      vz_new = rot12[2][0] * PV.VX[i] + rot12[2][1] * PV.VY[i] + rot12[2][2] * PV.VZ[i];

      PV.X[i] = x_new;
      PV.Y[i] = y_new;
      PV.Z[i] = z_new;

      PV.VX[i] = vx_new;
      PV.VY[i] = vy_new;
      PV.VZ[i] = vz_new;
    }
}


double get_H2_cool_rate(int i)
{
  int temp_idx;
  double dtemp, cool_rate;

  temp_idx = imin((int) (log10(PV.Temp[i] / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

  dtemp = PV.Temp[i] - All.TempTable[temp_idx];

  cool_rate = PV.AbH2[i] * PV.NH[i] * (All.CoolTable[temp_idx] + dtemp * All.DCoolTable[temp_idx]);

  return cool_rate;
}


double get_H2_esc_frac(int i)
{
  double esc_frac, gamma, csnd, jeans_length, N_H2_jeans, N_H2_eff, N_H2_LVG;

  esc_frac = 1;

  if(PV.NH[i] > 1e8)
    {
      if(BlockFlag[IO_GAMMA])
	gamma = PV.Gamma[i];
      else
	gamma = GAMMA_ADB;

      csnd = sqrt(gamma * BOLTZMANN * PV.Temp[i] / PROTONMASS);

      jeans_length = csnd * sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

      N_H2_jeans = PV.AbH2[i] * PV.NH[i] / (csnd / jeans_length);

      if(PV.DivVel[i] == 0)
	N_H2_eff = N_H2_jeans;
      else
	{
	  N_H2_LVG = PV.AbH2[i] * PV.NH[i] / dabs(PV.DivVel[i]);

	  N_H2_eff = dmin(N_H2_LVG, N_H2_jeans);
	}

      compute_H2_opacity(PV.Temp[i], N_H2_eff, &esc_frac);
    }

  //esc_frac = dmin(pow(PV.NH[i] / 8e9, -0.45), 1);

  return esc_frac;
}


double get_3b_form_heat(int i)
{
  double abhi, heat_rate;

  abhi = dmax(1 - 2 * PV.AbH2[i] - PV.AbHII[i], 0);

  heat_rate = 4.48 * ELECTRON_VOLT * 5.5e-29 * abhi * abhi * abhi * PV.NH[i] * PV.NH[i] * PV.NH[i] / PV.Temp[i];

  heat_rate += 4.48 * ELECTRON_VOLT * 5.5e-29 / 8 * abhi * abhi * PV.AbH2[i] * PV.NH[i] * PV.NH[i] * PV.NH[i] / PV.Temp[i];

  return heat_rate;
}


void compute_H2_opacity(double temp, double N_H2_eff, double *opac)
{
  double column_min, column_max, logN, logT, diff, dN, dT;
  double opac_tmp[2];
  int i, j;

  column_min = pow(1e1, All.H2OpacColumn[0]);
  column_max = pow(1e1, All.H2OpacColumn[TGCHEM_NUM_H2OP - 1]);

  if(N_H2_eff <= column_min)
    {
      *opac = 1;

      return;
    }
  else if(N_H2_eff >= column_max)
    {
      *opac = 0;

      return;
    }
  else
    {
      logN = log10(N_H2_eff);

      j = (int) (10 * (logN - 17));

      diff = All.H2OpacColumn[j + 1] - All.H2OpacColumn[j];

      dN = (logN - All.H2OpacColumn[j]) / diff;

      logT = log10(temp);

      if(logT <= All.H2OpacTemp[0])
	{
	  i = 1;

	  dT = 0;
	}
      else if(logT >= All.H2OpacTemp[TGCHEM_NUM_H2OP - 1])
	{
	  i = TGCHEM_NUM_H2OP - 1;

	  dT = 0;
	}
      else
	{
	  i = (int) ((logT - 1.5) / 0.03);

	  diff = All.H2OpacTemp[i + 1] - All.H2OpacTemp[i];

	  dT = (logT - All.H2OpacTemp[i]) / diff;
	}

      if(dT > 0)
	{
	  opac_tmp[0] = All.H2Opac[i * TGCHEM_NUM_H2OP + j] + dT * (All.H2Opac[(i + 1) * TGCHEM_NUM_H2OP + j] - All.H2Opac[i * TGCHEM_NUM_H2OP + j]);
	  opac_tmp[1] = All.H2Opac[i * TGCHEM_NUM_H2OP + j + 1] + dT * (All.H2Opac[(i + 1) * TGCHEM_NUM_H2OP + j + 1] - All.H2Opac[i * TGCHEM_NUM_H2OP + j + 1]);
	}
      else
	{
	  opac_tmp[0] = All.H2Opac[i * TGCHEM_NUM_H2OP + j];
	  opac_tmp[1] = All.H2Opac[i * TGCHEM_NUM_H2OP + j + 1];
	}

      *opac = pow(10, opac_tmp[0] + dN * (opac_tmp[1] - opac_tmp[0]));
    }
}


void tgchem_spline_eval(int nval, double *positions, double *values, int nnew, double *new_positions, double *new_values)
{
  int i, j, index;
  double coefficients[4][nval - 1];
  double dpt, pt, a, b, c, d;

  dpt = 0;

  tgchem_spline_coefficients(nval, values, &coefficients[0][0]);

  for(i = 0; i < nnew; i++)
    {
      pt = new_positions[i];

      if(pt < positions[0] || pt > positions[nval - 1])
	{
	  index = -1;
	  new_values[i] = 0;
	}
      else if(pt == positions[nval - 1])
	{
	  index = -1;
	  new_values[i] = values[nval - 1];
	}
      else
	{
	  index = 1;

	  for(j = 0; j < nval - 1; j++)
	    if(pt >= positions[j] && pt < positions[j + 1])
	      {
		index = j;
		dpt = (pt - positions[index]) / (positions[index + 1] - positions[index]);
	      }
	}

      if(index != -1)
	{
	  a = coefficients[0][index];
	  b = coefficients[1][index];
	  c = coefficients[2][index];
	  d = coefficients[3][index];

	  new_values[i] = a + b * dpt + c * dpt * dpt + d * dpt * dpt * dpt;
	}
    }
}


void tgchem_spline_coefficients(int nval, double *values, double *coefficients)
{
  int i;
  double derivatives[nval];

  tgchem_spline_derivatives(nval, values, &derivatives[0]);

  for(i = 0; i < nval - 1; i++)
    {
      *(coefficients + i) = values[i];
      *(coefficients + nval - 1 + i) = derivatives[i];
      *(coefficients + 2 * (nval - 1) + i) = 3 * (values[i + 1] - values[i]) - 2 * derivatives[i] - derivatives[i + 1];
      *(coefficients + 3 * (nval - 1) + i) = 2 * (values[i] - values[i + 1]) + derivatives[i] + derivatives[i + 1];
    }
}


void tgchem_spline_derivatives(int nval, double *values, double *derivatives)
{
  int i, j;
  double Cc[nval][nval + 1];
  double f, sum;

  for(i = 0; i < nval; i++)
    {
      for(j = 0; j < nval + 1; j++)
	Cc[i][j] = 0;

      if(i == 0)
	{
	  Cc[i][i] = 2;
	  Cc[i][i + 1] = 1;
	  Cc[i][nval] = 3 * (values[1] - values[0]);
	}
      else if(i == nval - 1)
	{
	  Cc[i][i - 1] = 1;
	  Cc[i][i] = 2;
	  Cc[i][nval] = 3 * (values[nval - 1] - values[nval - 2]);
	}
      else
	{
	  Cc[i][i - 1] = 1;
	  Cc[i][i] = 4;
	  Cc[i][i + 1] = 1;
	  Cc[i][nval] = 3 * (values[i + 1] - values[i - 1]);
	}
    }

  for(i = 1; i < nval; i++)
    {
      f = Cc[i][i - 1] / Cc[i - 1][i - 1];

      for(j = 0; j < nval + 1; j++)
	Cc[i][j] = Cc[i][j] - f * Cc[i - 1][j];
    }

  derivatives[nval - 1] = Cc[nval - 1][nval] / Cc[nval - 1][nval - 1];

  for(i = nval - 2; i > -1; i--)
    {
      sum = 0;

      for(j = i + 1; j < nval; j++)
	sum += Cc[i][j] * derivatives[j];

      derivatives[i] = (Cc[i][nval] - sum) / Cc[i][i];
    }
}


void endrun()
{
  mpi_printf("Run finished...bye!\n");

  fflush(stdout);

  MPI_Finalize();

  exit(0);
}


void mpi_printf(const char *buf, ...)
{
  if(ThisTask == 0)
    {
      va_list l;

      va_start(l, buf);

      vprintf(buf, l);

      fflush(stdout);

      va_end(l);
    }
}


int myflush(FILE *fstream)
{
  return fflush(fstream);
}


size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
	{
	  if(feof(stream))
	    printf("Read error (fread) on task %d has occured: end of file\n", ThisTask);
	  else
	    printf("Read error (fread) on task %d has occured: %s\n", ThisTask, strerror(errno));

	  fflush(stdout);

	  terminate("Read error!");
	}
    }
  else
    nread = 0;

  return nread;
}


size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("Write error (fwrite) on task %d has occured: %s\n", ThisTask, strerror(errno));

	  fflush(stdout);

	  terminate("write error");
	}
    }
  else
    nwritten = 0;

  return nwritten;
}


double mysort(void *base, size_t nel, size_t width, int (*compar) (const void *, const void *))
{
  double t0, t1;

  t0 = second();

  qsort(base, nel, width, compar);

  t1 = second();

  return timediff(t0, t1);
}


void terminate_args()
{
  terminate("Arguments missing!");
}


void terminate_block()
{
  terminate("Required block not present!");
}


int imin(int a, int b)
{
  if(a < b)
    return a;
  else
    return b;
}


int imax(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}


double dabs(double a)
{
  if(a < 0)
    return -a;
  else
    return a;
}


double dmin(double a, double b)
{
  if(a < b)
    return a;
  else
    return b;
}


double dmax(double a, double b)
{
  if(a > b)
    return a;
  else
    return b;
}


double second()
{
  return MPI_Wtime();
}


double measure_time()
{
  double t, dt;

  t = second();

  dt = t - WallClockTime;

  WallClockTime = t;

  return dt;
}


double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)
    dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;

  return dt;
}
