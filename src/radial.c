
#include "allvars.h"
#include "proto.h"


void radial()
{
  char radial_file[MAX_STRING_LEN], label[MAX_STRING_LEN];
  int i, j, idx, num_blocks, max_string_len, *r_idx, temp_idx;
  int block_geff, block_rates, block_cool, block_collapse, block_sidm_density, block_toomre, block_entro;
  double *r, loc_r_min, r_min, log_r_min, loc_r_max, r_max, log_r_max, *rad, fac, rfac;
  double *loc_vel_cm_x, *vel_cm_x, *loc_vel_cm_y, *vel_cm_y, *loc_vel_cm_z, *vel_cm_z;
  double *loc_vel_rad, *vel_rad, *vel_rot, *loc_vel_turb, *vel_turb, *kep_ratio, *vkep;
  double *loc_angmom_x, *angmom_x, *loc_angmom_y, *angmom_y, *loc_angmom_z, *angmom_z;
  double *frac_rad, *frac_rot, *frac_turb, *loc_csnd, *csnd, *mach;
  double vrad, vradx, vrady, vradz, vrotx, vroty, vrotz, vturbx, vturby, vturbz;
  double *loc_mass, *mass, *loc_logmass, *logmass, *enc, *loc_nh, *nh, *loc_dnh, *dnh, *loc_temp, *temp, *loc_abh2, *abh2;
  double *sigma_mass, *q, *omega_rot;
  double dist, mu, gamma, csound, tcool, theat, tcs, tacc, tff;
  double *loc_mu_rad, *mu_rad, *loc_geff, *geff, *loc_entro, *entro, *loc_escfrac, *escfrac;
  double *loc_pdvrate, *pdvrate, *loc_h2rate, *h2rate;
  double *loc_chemrate, *chemrate, *loc_cool, *cool, *loc_collapse, *collapse;
  double *loc_sidm_mass, *sidm_mass, *loc_sidm_density, *sidm_density;
  double cool_rate, esc_frac, abhi, abe;
  FILE *file;

  max_string_len = MAX_STRING_LEN;

  if(!BlockFlag[IO_POS] || !BlockFlag[IO_MASS] || !BlockFlag[IO_ID]) // || !BlockFlag[IO_VOL])
    terminate("Required block not read!");

  if(BlockFlag[IO_U] && BlockFlag[IO_RHO])
    block_geff = 1;
  else
    block_geff = 0;

  if(BlockFlag[IO_U] && BlockFlag[IO_RHO] && BlockFlag[IO_RATES])
    block_rates = 1;
  else
    block_rates = 0;

//  if(BlockFlag[IO_U] && BlockFlag[IO_RHO] && BlockFlag[IO_CHEM] && BlockFlag[IO_DIVVEL])
    block_cool = 1;
//  else
//    block_cool = 0;

  if(BlockFlag[IO_VEL] && BlockFlag[IO_U] && BlockFlag[IO_RHO])
    block_collapse = 1;
  else
    block_collapse = 0;

  block_toomre = 1;

  if(BlockFlag[IO_SIDM_HSML] && BlockFlag[IO_SIDM_DENSITY])
    block_sidm_density = 1;
  else
    block_sidm_density = 0;

  if(BlockFlag[IO_U] && BlockFlag[IO_RHO])
    block_entro = 1;
  else
    block_entro = 0;

  box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
    }

  if(All.FlagRotate)
    box_rotate();

  num_blocks = 1;

  if(BlockFlag[IO_VEL])
    num_blocks += 6;

  if(BlockFlag[IO_RHO])
    num_blocks += 2;

  if(BlockFlag[IO_U])
    num_blocks++;

  if(BlockFlag[IO_CHEM])
    num_blocks++;

  if(block_geff)
    num_blocks++;

  if(block_rates)
    num_blocks++;

  if(block_cool)
    num_blocks += 2;

  if(block_collapse)
    num_blocks++;

  if(block_toomre)
    num_blocks += 4;

  if(block_sidm_density)
    num_blocks++;

  if(block_entro)
    num_blocks++;

  r = (double *) mymalloc_movable(&r, "r", TotNumPart * sizeof(double));

  for(i = loc_r_max = 0, loc_r_min = DBL_MAX; i < TotNumPart; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      r[i] = sqrt(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i]);

//      if(block_sidm_density && i >= NumGas)
//	loc_r_min = dmin(PV.SIDM_Hsml[i - NumGas], loc_r_min);
//      else
//	loc_r_min = dmin(PV.Hsml[i], loc_r_min);

      loc_r_max = dmax(r[i], loc_r_max);
    }

  MPI_Allreduce(&loc_r_min, &r_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_r_max, &r_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(All.RadMin)
    r_min = All.RadMin;

  log_r_min = log10(r_min);
  log_r_max = log10(r_max);

  rfac = All.RadBins / (log_r_max - log_r_min);

  r_idx = (int *) mymalloc_movable(&r_idx, "r_idx", TotNumPart * sizeof(int));

  for(i = 0; i < TotNumPart; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      if(r[i] == 0)
	r_idx[i] = -1;
      else
	r_idx[i] = rfac * (log10(r[i]) - log_r_min);
    }

  rad = (double *) mymalloc_movable(&rad, "rad", All.RadBins * sizeof(double));

  for(i = 0; i < All.RadBins; i++)
    rad[i] = pow(10., log_r_min + i / rfac);

  loc_mass = (double *) mymalloc_movable(&loc_mass, "loc_mass", All.RadBins * sizeof(double));
  mass = (double *) mymalloc_movable(&mass, "mass", All.RadBins * sizeof(double));
  loc_logmass = (double *) mymalloc_movable(&loc_logmass, "loc_logmass", All.RadBins * sizeof(double));
  sigma_mass = (double *) mymalloc_movable(&sigma_mass, "sigma_mass", All.RadBins * sizeof(double));
  logmass = (double *) mymalloc_movable(&logmass, "logmass", All.RadBins * sizeof(double));
  enc = (double *) mymalloc_movable(&enc, "enc", All.RadBins * sizeof(double));

  memset(loc_mass, 0, All.RadBins * sizeof(double));
  memset(loc_logmass, 0, All.RadBins * sizeof(double));

  if(BlockFlag[IO_VEL])
    {
      loc_vel_cm_x = (double *) mymalloc_movable(&loc_vel_cm_x, "loc_vel_cm_x", All.RadBins * sizeof(double));
      vel_cm_x = (double *) mymalloc_movable(&vel_cm_x, "vel_cm_x", All.RadBins * sizeof(double));
      loc_vel_cm_y = (double *) mymalloc_movable(&loc_vel_cm_y, "loc_vel_cm_y", All.RadBins * sizeof(double));
      vel_cm_y = (double *) mymalloc_movable(&vel_cm_y, "vel_cm_y", All.RadBins * sizeof(double));
      loc_vel_cm_z = (double *) mymalloc_movable(&loc_vel_cm_z, "loc_vel_cm_z", All.RadBins * sizeof(double));
      vel_cm_z = (double *) mymalloc_movable(&vel_cm_z, "vel_cm_z", All.RadBins * sizeof(double));
      loc_vel_rad = (double *) mymalloc_movable(&loc_vel_rad, "loc_vel_rad", All.RadBins * sizeof(double));
      vel_rad = (double *) mymalloc_movable(&vel_rad, "vel_rad", All.RadBins * sizeof(double));
      frac_rad = (double *) mymalloc_movable(&frac_rad, "frac_rad", All.RadBins * sizeof(double));
      loc_angmom_x = (double *) mymalloc_movable(&loc_angmom_x, "loc_angmom_x", All.RadBins * sizeof(double));
      angmom_x = (double *) mymalloc_movable(&angmom_x, "angmom_x", All.RadBins * sizeof(double));
      loc_angmom_y = (double *) mymalloc_movable(&loc_angmom_y, "loc_angmom_y", All.RadBins * sizeof(double));
      angmom_y = (double *) mymalloc_movable(&angmom_y, "angmom_y", All.RadBins * sizeof(double));
      loc_angmom_z = (double *) mymalloc_movable(&loc_angmom_z, "loc_angmom_z", All.RadBins * sizeof(double));
      angmom_z = (double *) mymalloc_movable(&angmom_z, "angmom_z", All.RadBins * sizeof(double));
      vel_rot = (double *) mymalloc_movable(&vel_rot, "vel_rot", All.RadBins * sizeof(double));
      omega_rot = (double *) mymalloc_movable(&omega_rot, "omega_rot", All.RadBins * sizeof(double));
      frac_rot = (double *) mymalloc_movable(&frac_rot, "frac_rot", All.RadBins * sizeof(double));
      kep_ratio = (double *) mymalloc_movable(&kep_ratio, "kep_ratio", All.RadBins * sizeof(double));
      vkep = (double *) mymalloc_movable(&vkep, "vkep", All.RadBins * sizeof(double));
      loc_vel_turb = (double *) mymalloc_movable(&loc_vel_turb, "loc_vel_turb", All.RadBins * sizeof(double));
      vel_turb = (double *) mymalloc_movable(&vel_turb, "vel_turb", All.RadBins * sizeof(double));
      frac_turb = (double *) mymalloc_movable(&frac_turb, "frac_turb", All.RadBins * sizeof(double));
      loc_csnd = (double *) mymalloc_movable(&loc_csnd, "loc_csnd", All.RadBins * sizeof(double));
      csnd = (double *) mymalloc_movable(&csnd, "csnd", All.RadBins * sizeof(double));
      mach = (double *) mymalloc_movable(&mach, "mach", All.RadBins * sizeof(double));

      memset(loc_vel_cm_x, 0, All.RadBins * sizeof(double));
      memset(loc_vel_cm_y, 0, All.RadBins * sizeof(double));
      memset(loc_vel_cm_z, 0, All.RadBins * sizeof(double));
      memset(loc_vel_rad, 0, All.RadBins * sizeof(double));
      memset(loc_angmom_x, 0, All.RadBins * sizeof(double));
      memset(loc_angmom_y, 0, All.RadBins * sizeof(double));
      memset(loc_angmom_z, 0, All.RadBins * sizeof(double));
      memset(loc_vel_turb, 0, All.RadBins * sizeof(double));
      memset(loc_csnd, 0, All.RadBins * sizeof(double));
    }

  if(BlockFlag[IO_RHO])
    {
      loc_nh = (double *) mymalloc_movable(&loc_nh, "loc_nh", All.RadBins * sizeof(double));
      nh = (double *) mymalloc_movable(&nh, "nh", All.RadBins * sizeof(double));
      loc_dnh = (double *) mymalloc_movable(&loc_dnh, "loc_dnh", All.RadBins * sizeof(double));
      dnh = (double *) mymalloc_movable(&dnh, "dnh", All.RadBins * sizeof(double));

      memset(loc_nh, 0, All.RadBins * sizeof(double));
      memset(loc_dnh, 0, All.RadBins * sizeof(double));
    }

  if(BlockFlag[IO_U])
    {
      loc_temp = (double *) mymalloc_movable(&loc_temp, "loc_temp", All.RadBins * sizeof(double));
      temp = (double *) mymalloc_movable(&temp, "temp", All.RadBins * sizeof(double));

      memset(loc_temp, 0, All.RadBins * sizeof(double));
    }

  if(BlockFlag[IO_CHEM])
    {
      loc_abh2 = (double *) mymalloc_movable(&loc_abh2, "loc_abh2", All.RadBins * sizeof(double));
      abh2 = (double *) mymalloc_movable(&abh2, "abh2", All.RadBins * sizeof(double));

      memset(loc_abh2, 0, All.RadBins * sizeof(double));
    }

  if(block_geff)
    {
      loc_mu_rad = (double *) mymalloc_movable(&loc_mu_rad, "loc_mu_rad", All.RadBins * sizeof(double));
      mu_rad = (double *) mymalloc_movable(&mu_rad, "mu_rad", All.RadBins * sizeof(double));
      loc_geff = (double *) mymalloc_movable(&loc_geff, "loc_geff", All.RadBins * sizeof(double));
      geff = (double *) mymalloc_movable(&geff, "geff", All.RadBins * sizeof(double));

      memset(loc_mu_rad, 0, All.RadBins * sizeof(double));
      memset(loc_geff, 0, All.RadBins * sizeof(double));
    }

  if(block_rates)
    {
      loc_pdvrate = (double *) mymalloc_movable(&loc_pdvrate, "loc_pdvrate", All.RadBins * sizeof(double));
      pdvrate = (double *) mymalloc_movable(&pdvrate, "pdvrate", All.RadBins * sizeof(double));

      memset(loc_pdvrate, 0, All.RadBins * sizeof(double));
    }

  if(block_cool)
    {
      loc_escfrac = (double *) mymalloc_movable(&loc_escfrac, "loc_escfrac", All.RadBins * sizeof(double));
      escfrac = (double *) mymalloc_movable(&escfrac, "escfrac", All.RadBins * sizeof(double));
/*      loc_h2rate = (double *) mymalloc_movable(&loc_h2rate, "loc_h2rate", All.RadBins * sizeof(double));
      h2rate = (double *) mymalloc_movable(&h2rate, "h2rate", All.RadBins * sizeof(double));
      loc_chemrate = (double *) mymalloc_movable(&loc_chemrate, "loc_chemrate", All.RadBins * sizeof(double));
      chemrate = (double *) mymalloc_movable(&chemrate, "chemrate", All.RadBins * sizeof(double));
*/
      loc_cool = (double *) mymalloc_movable(&loc_cool, "loc_cool", All.RadBins * sizeof(double));
      cool = (double *) mymalloc_movable(&cool, "cool", All.RadBins * sizeof(double));

      memset(loc_escfrac, 0, All.RadBins * sizeof(double));
/*      memset(loc_h2rate, 0, All.RadBins * sizeof(double));
      memset(loc_chemrate, 0, All.RadBins * sizeof(double));
*/
      memset(loc_cool, 0, All.RadBins * sizeof(double));
    }

  if(block_collapse)
    {
      loc_collapse = (double *) mymalloc_movable(&loc_collapse, "loc_collapse", All.RadBins * sizeof(double));
      collapse = (double *) mymalloc_movable(&collapse, "collapse", All.RadBins * sizeof(double));

      memset(loc_collapse, 0, All.RadBins * sizeof(double));
    }

  if(block_toomre)
      q = (double *) mymalloc_movable(&q, "q", All.RadBins * sizeof(double));

  if(block_sidm_density)
    {
      loc_sidm_mass = (double *) mymalloc_movable(&loc_sidm_mass, "loc_sidm_mass", All.RadBins * sizeof(double));
      sidm_mass = (double *) mymalloc_movable(&sidm_mass, "sidm_mass", All.RadBins * sizeof(double));
      loc_sidm_density = (double *) mymalloc_movable(&loc_sidm_density, "loc_sidm_density", All.RadBins * sizeof(double));
      sidm_density = (double *) mymalloc_movable(&sidm_density, "sidm_density", All.RadBins * sizeof(double));

      memset(loc_sidm_mass, 0, All.RadBins * sizeof(double));
      memset(loc_sidm_density, 0, All.RadBins * sizeof(double));
    }

  if(block_entro)
    {
      loc_entro = (double *) mymalloc_movable(&loc_entro, "loc_entro", All.RadBins * sizeof(double));
      entro = (double *) mymalloc_movable(&entro, "entro", All.RadBins * sizeof(double));

      memset(loc_entro, 0, All.RadBins * sizeof(double));
    }

  for(i = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      j = r_idx[i];

      if(j < 0 || j > All.RadBins - 1)
	continue;

      loc_mass[j] += PV.Mass[i];

      loc_logmass[j] += log10(PV.Mass[i]);

      if(BlockFlag[IO_VEL])
	{
	  loc_vel_cm_x[j] += PV.Mass[i] * PV.VX[i];
	  loc_vel_cm_y[j] += PV.Mass[i] * PV.VY[i];
	  loc_vel_cm_z[j] += PV.Mass[i] * PV.VZ[i];
	}

      if(BlockFlag[IO_RHO])
	loc_nh[j] += PV.Mass[i] * log10(PV.NH[i]);

      if(BlockFlag[IO_U])
	loc_temp[j] += PV.Mass[i] * log10(PV.Temp[i]);

      if(BlockFlag[IO_CHEM])
	loc_abh2[j] += PV.Mass[i] * log10(PV.AbH2[i]);

      if(block_geff)
	loc_mu_rad[j] += PV.Mass[i] * log10(PV.Mu[i]);

      if(block_rates)
	{
	  tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	  theat = PV.U[i] * PV.Rho[i] / dabs(PV.PdVRate[i]);

	  loc_pdvrate[j] += PV.Mass[i] * log10(theat / tff);
	}

      if(block_cool)
	{

          tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

          abe = PV.AbHII[i];

          abhi = dmax(1. - 2. * PV.AbH2[i] - PV.AbHII[i], 0.);

          // HI electronic excitation cooling
          cool_rate = 7.5e-19 * exp(-1.18348e5 / PV.Temp[i]) / (1. + sqrt(PV.Temp[i] / 1e5)) * abhi * abe * PV.NH[i] * PV.NH[i];

          esc_frac = exp(-PV.NH[i] / 1e15);

          loc_escfrac[j] += PV.Mass[i] * log10(esc_frac);

          cool_rate *= esc_frac;

          tcool = PV.U[i] * PV.Rho[i] / cool_rate;

          if(tcool > 0)
            loc_cool[j] += PV.Mass[i] * log10(tcool / tff);
	}

     if(block_entro)
        loc_entro[j] += PV.Mass[i] * log10(PV.Mass[i] * pow(PV.Rho[i], 1 - GAMMA_ADB) * PV.Temp[i]);
         
    }

  if(block_sidm_density)
    for(i = NumGas; i < NumGas + NumPart[1]; i++)
      {
	if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	  continue;

	j = r_idx[i];

	if(j < 0 || j > All.RadBins - 1)
	  continue;

	loc_sidm_mass[j] += PV.Mass[i];
	loc_sidm_density[j] += PV.Mass[i] * log10(PV.SIDM_Density[i - NumGas]);
      }

  MPI_Allreduce(loc_mass, mass, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(loc_logmass, logmass, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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

  for(i = 0; i < All.RadBins; i++)
   {
    if(!mass[i])
      mpi_printf("Warning! Bin %d without mass!\n", i);
   if(i == 0)
     sigma_mass[i] = mass[i] * SOLAR_MASS / M_PI / (fac * rad[i]) / (fac * rad[i]);
   else
     sigma_mass[i] = mass[i] * SOLAR_MASS / M_PI / (fac * rad[i] * fac * rad[i] - fac * rad[i-1] * fac * rad[i-1]);
   }

  if(BlockFlag[IO_VEL])
    {
      MPI_Allreduce(loc_vel_cm_x, vel_cm_x, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_vel_cm_y, vel_cm_y, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_vel_cm_z, vel_cm_z, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

  if(BlockFlag[IO_RHO])
    MPI_Allreduce(loc_nh, nh, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(BlockFlag[IO_U])
    MPI_Allreduce(loc_temp, temp, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(BlockFlag[IO_CHEM])
    MPI_Allreduce(loc_abh2, abh2, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(block_geff)
    MPI_Allreduce(loc_mu_rad, mu_rad, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(block_rates)
    MPI_Allreduce(loc_pdvrate, pdvrate, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(block_cool)
    {
      MPI_Allreduce(loc_escfrac, escfrac, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
/*      MPI_Allreduce(loc_h2rate, h2rate, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_chemrate, chemrate, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
*/
      MPI_Allreduce(loc_cool, cool, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

  if(block_sidm_density)
    {
      MPI_Allreduce(loc_sidm_mass, sidm_mass, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_sidm_density, sidm_density, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
	if(!sidm_mass[i])
	  mpi_printf("Warning! Bin %d without mass!\n", i);
    }

  if(block_entro)
    MPI_Allreduce(loc_entro, entro, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(BlockFlag[IO_RHO])
    {
      for(i = 0; i < All.RadBins; i++)
	nh[i] = pow(10, nh[i] / mass[i]);

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  j = r_idx[i];

	  if(j < 0 || j > All.RadBins - 1)
	    continue;

	  loc_dnh[j] += PV.Mass[i] * pow((PV.NH[i] - nh[j]) / nh[j], 2);
	}

      MPI_Allreduce(loc_dnh, dnh, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
	dnh[i] = log10(sqrt(dnh[i] / mass[i]));
    }

  if(BlockFlag[IO_U])
    for(i = 0; i < All.RadBins; i++)
      temp[i] = pow(10, temp[i] / mass[i]);

  if(block_geff)
    {
      for(i = 0; i < All.RadBins; i++)
	mu_rad[i] = pow(10, mu_rad[i] / mass[i]);

      for(i = 1; i < All.RadBins; i++)
	geff[i] = log10(nh[i] * temp[i] / mu_rad[i] / (nh[i - 1] * temp[i - 1] / mu_rad[i - 1])) / log10(nh[i] / nh[i - 1]);

      geff[0] = geff[1];
    }

  for(i = 0; i < All.RadBins; i++)
    {
      if(i == 0)
	enc[i] = mass[i];
      else
	enc[i] = enc[i - 1] + mass[i];

//      sigma_mass[i] = enc[i] * SOLAR_MASS / M_PI / (fac * rad[i]) / (fac * rad[i]);

      if(BlockFlag[IO_VEL])
	{
	  vel_cm_x[i] /= mass[i];
	  vel_cm_y[i] /= mass[i];
	  vel_cm_z[i] /= mass[i];
	}

      if(BlockFlag[IO_RHO])
	nh[i] = log10(nh[i]);

      if(BlockFlag[IO_U])
	temp[i] = log10(temp[i]);

      if(BlockFlag[IO_CHEM])
	abh2[i] /= mass[i];

      if(block_rates)
	pdvrate[i] /= mass[i];

      if(block_cool)
	{
	  escfrac[i] /= mass[i];
/*	  h2rate[i] /= mass[i];
	  chemrate[i] /= mass[i];
*/
	  cool[i] /= mass[i];
	}

      if(block_sidm_density)
	sidm_density[i] /= sidm_mass[i];

      if(block_entro)
        entro[i] /= mass[i];
    }

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
	}

      MPI_Allreduce(loc_vel_rad, vel_rad, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(loc_angmom_x, angmom_x, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_angmom_y, angmom_y, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_angmom_z, angmom_z, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
	{
	  vel_rad[i] /= mass[i];

//          Accretion rate:  Mdot = -4piR^2 v_rad * rho
//          vel_rad[i] = -4 * M_PI * pow(rad[i] * UNIT_LENGTH / 1e3, 2)  * vel_rad[i] * 1e5 * pow(10, nh[i]) * PROTONMASS / SOLAR_MASS * SEC_PER_YEAR;

	  vel_rot[i] = sqrt(angmom_x[i] * angmom_x[i] + angmom_y[i] * angmom_y[i] + angmom_z[i] * angmom_z[i]) / rad[i] / mass[i];

          omega_rot[i] = vel_rot[i] / rad[i];

	  frac_rot[i] = vel_rot[i] / dabs(vel_rad[i]);

	  vkep[i] = sqrt(GRAVITY * enc[i] * SOLAR_MASS / (fac * rad[i])) / 1e5;

	  kep_ratio[i] = vel_rot[i] / vkep[i];
	}

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  j = r_idx[i];

	  if(j < 0 || j > All.RadBins - 1)
	    continue;

	  vradx = vel_rad[j] * PV.X[i] / rad[j];
	  vrady = vel_rad[j] * PV.Y[i] / rad[j];
	  vradz = vel_rad[j] * PV.Z[i] / rad[j];

	  vrotx = (angmom_y[j] * PV.Z[i] - angmom_z[j] * PV.Y[i]) / rad[j] / rad[j] / mass[j];
	  vroty = (angmom_z[j] * PV.X[i] - angmom_x[j] * PV.Z[i]) / rad[j] / rad[j] / mass[j];
	  vrotz = (angmom_x[j] * PV.Y[i] - angmom_y[j] * PV.X[i]) / rad[j] / rad[j] / mass[j];

	  vturbx = PV.VX[i] - vradx - vrotx;
	  vturby = PV.VY[i] - vrady - vroty;
	  vturbz = PV.VZ[i] - vradz - vrotz;

	  loc_vel_turb[j] += PV.Mass[i] * sqrt(vturbx * vturbx + vturby * vturby + vturbz * vturbz);

	  loc_csnd[j] += PV.Mass[i] * log10(sqrt(PV.Gamma[i] * BOLTZMANN * PV.Temp[i] / PV.Mu[i] / PROTONMASS) / 1e5);
	}

      MPI_Allreduce(loc_vel_turb, vel_turb, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_csnd, csnd, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
	{
	  csnd[i] = pow(10, csnd[i] / mass[i]);

	  frac_rad[i] = dabs(vel_rad[i]) / csnd[i];

	  frac_rot[i] = vel_rot[i] / csnd[i];

	  vel_turb[i] /= mass[i];

	  frac_turb[i] = vel_turb[i] / dabs(vel_rad[i]);

	  mach[i] = vel_turb[i] / csnd[i];

/*          mpi_printf("%f %f %f %f %f %f\n", vel_rad[i], vel_rot[i], vkep[i], vel_turb[i], csnd[i], rad[i]); */
	}
    }

  if(block_toomre)
    for(i = 0; i < All.RadBins; i++)
      q[i] = csnd[i] * 1e5 * omega_rot[i] / SEC_PER_YEAR / M_PI / GRAVITY / sigma_mass[i];

  if(block_collapse)
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

	  dist = sqrt(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i]);

	  tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	  if(BlockFlag[IO_CHEM])
	    mu = PV.Mu[i];
	  else
	    mu = MU_NEUTRAL;

	  if(BlockFlag[IO_GAMMA])
	    gamma = PV.Gamma[i];
	  else
	    gamma = GAMMA_ADB;

	  csound = sqrt(gamma * BOLTZMANN * PV.Temp[i] / mu / PROTONMASS);

	  tcs = fac * dist / csound;

	  vradx = PV.X[i] * (PV.VX[i] - vel_cm_x[j]);
	  vrady = PV.Y[i] * (PV.VY[i] - vel_cm_y[j]);
	  vradz = PV.Z[i] * (PV.VZ[i] - vel_cm_z[j]);

	  vrad = (vradx + vrady + vradz) / dist;

	  tacc = dabs(fac * dist / vrad / 1e5);

	  //loc_collapse[j] += PV.Mass[i] * dmax(tff / tcs, tff / tacc);
	  loc_collapse[j] += PV.Mass[i] * log10(tff / tcs);
	}

      MPI_Allreduce(loc_collapse, collapse, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < All.RadBins; i++)
	collapse[i] /= mass[i];
    }

  for(i = 0; i < All.RadBins; i++)
    {
      rad[i] = log10(rad[i]);
      enc[i] = log10(enc[i]);
      omega_rot[i] = log10(omega_rot[i]);
      sigma_mass[i] = log10(sigma_mass[i]);
    }

  if(ThisTask == 0)
    {
      sprintf(radial_file, "%s/%s/%s_%s.radial", All.Path, All.Base, All.Base, All.SnapNumString);

      if(!(file = fopen(radial_file, "w")))
	terminate("Could not write radial file!");

      fwrite(&num_blocks, sizeof(int), 1, file);
      fwrite(&max_string_len, sizeof(int), 1, file);
      fwrite(&All.RadBins, sizeof(int), 1, file);
      fwrite(&All.LengthUnit, sizeof(int), 1, file);
      fwrite(&All.ComovingUnits, sizeof(int), 1, file);
      fwrite(&log_r_min, sizeof(double), 1, file);
      fwrite(&rad[All.RadBins - 1], sizeof(double), 1, file);
      fwrite(rad, All.RadBins * sizeof(double), 1, file);

      if(BlockFlag[IO_RHO])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "nh");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(nh, All.RadBins * sizeof(double), 1, file);
	}

      memset(label, 0, MAX_STRING_LEN * sizeof(char));
      strcpy(label, "enc");
      fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
      fwrite(enc, All.RadBins * sizeof(double), 1, file);

      if(BlockFlag[IO_U])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "temp");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(temp, All.RadBins * sizeof(double), 1, file);
	}

      if(BlockFlag[IO_CHEM])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "abh2");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(abh2, All.RadBins * sizeof(double), 1, file);
	}

      if(BlockFlag[IO_VEL])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "vrad");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(vel_rad, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "frac_rad");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(frac_rad, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "frac_rot");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(frac_rot, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "kep_ratio");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(kep_ratio, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "frac_turb");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(frac_turb, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "mach");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(mach, All.RadBins * sizeof(double), 1, file);
	}

      if(block_rates)
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "pdvrate");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(pdvrate, All.RadBins * sizeof(double), 1, file);
	}

      if(block_cool)
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "escfrac");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(escfrac, All.RadBins * sizeof(double), 1, file);

/*	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "h2rate");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(h2rate, All.RadBins * sizeof(double), 1, file);

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "chemrate");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(chemrate, All.RadBins * sizeof(double), 1, file);
*/
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "cool");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(cool, All.RadBins * sizeof(double), 1, file);
	}

      if(BlockFlag[IO_CHEM])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "geff");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(geff, All.RadBins * sizeof(double), 1, file);
	}

      if(block_collapse)
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "collapse");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(collapse, All.RadBins * sizeof(double), 1, file);
	}

      if(BlockFlag[IO_RHO])
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "dnh");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(dnh, All.RadBins * sizeof(double), 1, file);
	}

      if(block_sidm_density)
	{
	  memset(label, 0, MAX_STRING_LEN * sizeof(char));
	  strcpy(label, "sidm_density");
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(sidm_density, All.RadBins * sizeof(double), 1, file);
	}

     if(block_toomre)
        {
          memset(label, 0, MAX_STRING_LEN * sizeof(char));
          strcpy(label, "sigma_mass");
          fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
          fwrite(sigma_mass, All.RadBins * sizeof(double), 1, file);

          memset(label, 0, MAX_STRING_LEN * sizeof(char));
          strcpy(label, "csnd");
          fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
          fwrite(csnd, All.RadBins * sizeof(double), 1, file);

          memset(label, 0, MAX_STRING_LEN * sizeof(char));
          strcpy(label, "omega_rot");
          fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
          fwrite(omega_rot, All.RadBins * sizeof(double), 1, file);

          memset(label, 0, MAX_STRING_LEN * sizeof(char));
          strcpy(label, "q");
          fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
          fwrite(q, All.RadBins * sizeof(double), 1, file);
        }

     if(block_entro)
        {
          memset(label, 0, MAX_STRING_LEN * sizeof(char));
          strcpy(label, "entro");
          fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
          fwrite(entro, All.RadBins * sizeof(double), 1, file);
        }

      fclose(file);
    }

  myfree_movable(r);
  myfree_movable(r_idx);
  myfree_movable(rad);
  myfree_movable(loc_mass);
  myfree_movable(mass);
  myfree_movable(sigma_mass);
  myfree_movable(loc_logmass);
  myfree_movable(logmass);
  myfree_movable(enc);

  if(BlockFlag[IO_VEL])
    {
      myfree_movable(loc_vel_cm_x);
      myfree_movable(vel_cm_x);
      myfree_movable(loc_vel_cm_y);
      myfree_movable(vel_cm_y);
      myfree_movable(loc_vel_cm_z);
      myfree_movable(vel_cm_z);
      myfree_movable(loc_vel_rad);
      myfree_movable(vel_rad);
      myfree_movable(frac_rad);
      myfree_movable(loc_angmom_x);
      myfree_movable(angmom_x);
      myfree_movable(loc_angmom_y);
      myfree_movable(angmom_y);
      myfree_movable(loc_angmom_z);
      myfree_movable(angmom_z);
      myfree_movable(vel_rot);
      myfree_movable(omega_rot);
      myfree_movable(frac_rot);
      myfree_movable(kep_ratio);
      myfree_movable(vkep);
      myfree_movable(loc_vel_turb);
      myfree_movable(vel_turb);
      myfree_movable(frac_turb);
      myfree_movable(loc_csnd);
      myfree_movable(csnd);
      myfree_movable(mach);
    }

  if(BlockFlag[IO_RHO])
    {
      myfree_movable(loc_nh);
      myfree_movable(nh);
      myfree_movable(loc_dnh);
      myfree_movable(dnh);
    }

  if(BlockFlag[IO_U])
    {
      myfree_movable(loc_temp);
      myfree_movable(temp);
    }

  if(BlockFlag[IO_CHEM])
    {
      myfree_movable(loc_abh2);
      myfree_movable(abh2);
    }

  if(block_geff)
    {
      myfree_movable(loc_mu_rad);
      myfree_movable(mu_rad);
      myfree_movable(loc_geff);
      myfree_movable(geff);
    }

  if(block_rates)
    {
      myfree_movable(loc_pdvrate);
      myfree_movable(pdvrate);
    }

  if(block_cool)
    {
/*      myfree_movable(loc_h2rate);
      myfree_movable(h2rate);
      myfree_movable(loc_chemrate);
      myfree_movable(chemrate);
*/
      myfree_movable(loc_cool);
      myfree_movable(cool);
      myfree_movable(loc_escfrac);
      myfree_movable(escfrac);
    }

  if(block_collapse)
    {
      myfree_movable(loc_collapse);
      myfree_movable(collapse);
    }

  if(block_toomre)
      myfree_movable(q);

  if(block_sidm_density)
    {
      myfree_movable(loc_sidm_mass);
      myfree_movable(sidm_mass);
      myfree_movable(loc_sidm_density);
      myfree_movable(sidm_density);
    }

  if(block_entro)
    {
      myfree_movable(loc_entro);
      myfree_movable(entro);
    }
}
