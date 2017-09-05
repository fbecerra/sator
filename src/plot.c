
#include "allvars.h"
#include "proto.h"


int plot_block_present(int block)
{
  int val;

  if(block == PLOT_NH && BlockFlag[IO_RHO])
    val = 1;
  else if(block == PLOT_TEMP && BlockFlag[IO_U])
    val = 1;
  else if(block == PLOT_GRAVACC && BlockFlag[IO_GRAVACC])
    val = 1;
  else if(block == PLOT_GRADP && BlockFlag[IO_GRADP])
    val = 1;
#ifdef OLDCHEM
  else if(block == PLOT_H2 && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_HII && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_DII && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_HD && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_HEII && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_HEIII && BlockFlag[IO_CHEM])
    val = 1;
#else
  else if(block == PLOT_HM && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_H2 && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_HII && BlockFlag[IO_CHEM])
    val = 1;
#endif
  else if(block == PLOT_GAMMA && BlockFlag[IO_GAMMA])
    val = 1;
  else if(block == PLOT_ESCFRAC && BlockFlag[IO_RHO])
    val = 1;
  else if(block == PLOT_PDVRATE && BlockFlag[IO_RATES])
    val = 1;
  else if(block == PLOT_H2RATE && BlockFlag[IO_RATES])
    val = 1;
  else if(block == PLOT_CIERATE && BlockFlag[IO_RATES])
    val = 1;
  else if(block == PLOT_CHEMRATE && BlockFlag[IO_RATES])
    val = 1;
  else if(block == PLOT_ALLOWREF && BlockFlag[IO_ALLOWREF])
    val = 1;
  else if(block == PLOT_DIVVEL && BlockFlag[IO_DIVVEL])
    val = 1;
  else if(block == PLOT_COOL && BlockFlag[IO_U] && BlockFlag[IO_RHO] && BlockFlag[IO_CHEM])
    val = 1;
  else if(block == PLOT_COLLAPSE && BlockFlag[IO_POS] && BlockFlag[IO_VEL] && BlockFlag[IO_U] && BlockFlag[IO_RHO])
    val = 1;
  else if(block == PLOT_SIDM_DENSITY && BlockFlag[IO_SIDM_HSML] && BlockFlag[IO_SIDM_DENSITY])
    val = 1;
  else
    val = 0;

  if(val && BlockPlot[block])
    return 1;
  else
    return 0;
}


void plot_block_label(int block, char *label)
{
  if(block == PLOT_NH)
    strcpy(label, "nh");
  else if(block == PLOT_TEMP)
    strcpy(label, "temp");
  else if(block == PLOT_GRAVACC)
    strcpy(label, "gravacc");
  else if(block == PLOT_GRADP)
    strcpy(label, "gradp");
#ifdef OLDCHEM
  else if(block == PLOT_H2)
    strcpy(label, "abh2");
  else if(block == PLOT_HII)
    strcpy(label, "abhii");
  else if(block == PLOT_DII)
    strcpy(label, "abdii");
  else if(block == PLOT_HD)
    strcpy(label, "abhd");
  else if(block == PLOT_HEII)
    strcpy(label, "abheii");
  else if(block == PLOT_HEIII)
    strcpy(label, "abheiii");
#else
  else if(block == PLOT_HM)
    strcpy(label, "abhm");
  else if(block == PLOT_H2)
    strcpy(label, "abh2");
  else if(block == PLOT_HII)
    strcpy(label, "abhii");
#endif
  else if(block == PLOT_GAMMA)
    strcpy(label, "gamma");
  else if(block == PLOT_ESCFRAC)
    strcpy(label, "escfrac");
  else if(block == PLOT_PDVRATE)
    strcpy(label, "pdvrate");
  else if(block == PLOT_H2RATE)
    strcpy(label, "h2rate");
  else if(block == PLOT_CIERATE)
    strcpy(label, "cierate");
  else if(block == PLOT_CHEMRATE)
    strcpy(label, "chemrate");
  else if(block == PLOT_ALLOWREF)
    strcpy(label, "allowref");
  else if(block == PLOT_DIVVEL)
    strcpy(label, "divvel");
  else if(block == PLOT_COOL)
    strcpy(label, "cool");
  else if(block == PLOT_COLLAPSE)
    strcpy(label, "collapse");
  else if(block == PLOT_SIDM_DENSITY)
    strcpy(label, "sidm_density");
  else
    terminate("Unknown block in block_label!");
}


double *plot_block_val(int block)
{
  int i, temp_idx;
  double r, fac, vradx, vrady, vradz, vrad, acc;
  double mu, gamma, csnd, tacc, tcool, tcs, tff, abe, abhi;
  double loc_vel_cm_x, loc_vel_cm_y, loc_vel_cm_z, vel_cm_x, vel_cm_y, vel_cm_z;
  double cool_rate, esc_frac, heat_rate;

  if(block == PLOT_NH)
    return PV.NH;
  else if(block == PLOT_TEMP)
    return PV.Temp;
  else if(block == PLOT_GRAVACC)
    {
      for(i = 0; i < NumGas; i++)
	{
	  acc = 0;

	  acc += PV.GravAccX[i] * PV.GravAccX[i];
	  acc += PV.GravAccY[i] * PV.GravAccY[i];
	  acc += PV.GravAccZ[i] * PV.GravAccZ[i];

	  PV.GravAccX[i] = sqrt(acc);
	}

      return PV.GravAccX;
    }
  else if(block == PLOT_GRADP)
    {
      for(i = 0; i < NumGas; i++)
	{
	  acc = 0;

	  acc += PV.GradPX[i] * PV.GradPX[i];
	  acc += PV.GradPY[i] * PV.GradPY[i];
	  acc += PV.GradPZ[i] * PV.GradPZ[i];

	  PV.GradPX[i] = sqrt(acc);
	}

      return PV.GradPX;
    }
#ifdef OLDCHEM
  else if(block == PLOT_H2)
    return PV.AbH2;
  else if(block == PLOT_HII)
    return PV.AbHII;
  else if(block == PLOT_DII)
    return PV.AbDII;
  else if(block == PLOT_HD)
    return PV.AbHD;
  else if(block == PLOT_HEII)
    return PV.AbHeII;
  else if(block == PLOT_HEIII)
    return PV.AbHeIII;
#else
  else if(block == PLOT_HM)
    return PV.AbHM;
  else if(block == PLOT_H2)
    return PV.AbH2;
  else if(block == PLOT_HII)
    return PV.AbHII;
#endif
  else if(block == PLOT_GAMMA)
    return PV.Gamma;
  else if(block == PLOT_ESCFRAC)
    {
      PV.EscFrac = mymalloc_movable(&PV.EscFrac, "PV.EscFrac", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	PV.EscFrac[i] = exp(-PV.NH[i] / 1e15);

      return PV.EscFrac;
    }
  else if(block == PLOT_PDVRATE)
    {
      PV.PdVCool = mymalloc_movable(&PV.PdVCool, "PV.PdVCool", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.PdVRate[i])
	    {
	      tcool = PV.U[i] * PV.Rho[i] / PV.PdVRate[i];

	      tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	      PV.PdVCool[i] = dabs(tcool / tff);
	    }
	  else
	    PV.PdVCool[i] = 0;
	}

      return PV.PdVCool;
    }
  else if(block == PLOT_H2RATE)
    {
      PV.H2Cool = mymalloc_movable(&PV.H2Cool, "PV.H2Cool", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.H2Rate[i])
	    {
	      tcool = PV.U[i] * PV.Rho[i] / PV.H2Rate[i];

	      tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	      PV.H2Cool[i] = dabs(tcool / tff);
	    }
	  else
	    PV.H2Cool[i] = 0;
	}

      return PV.H2Cool;
    }
  else if(block == PLOT_CIERATE)
    {
      PV.CIECool = mymalloc_movable(&PV.CIECool, "PV.CIECool", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.CIERate[i])
	    {
	      tcool = PV.U[i] * PV.Rho[i] / PV.CIERate[i];

	      tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	      PV.CIECool[i] = dabs(tcool / tff);
	    }
	  else
	    PV.CIECool[i] = 0;
	}

      return PV.CIECool;
    }
  else if(block == PLOT_CHEMRATE)
    {
      PV.ChemCool = mymalloc_movable(&PV.ChemCool, "PV.ChemCool", NumGas * sizeof(double));

      //double abhi;

      for(i = 0; i < NumGas; i++)
	{
	  //abhi = dmax(1 - 2 * PV.AbH2[i] - PV.AbHII[i], 0);

	  //PV.ChemRate[i] = 4.48 * ELECTRON_VOLT * 5.5e-29 / PV.Temp[i] * abhi * abhi * (abhi + PV.H2I[i] / 8) * PV.NH[i] * PV.NH[i] * PV.NH[i];

	  if(PV.ChemRate[i])
	    {
	      tcool = PV.U[i] * PV.Rho[i] / PV.ChemRate[i];

	      tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	      PV.ChemCool[i] = dabs(tcool / tff);
	    }
	  else
	    PV.ChemCool[i] = 0;
	}

      return PV.ChemCool;
    }
  else if(block == PLOT_ALLOWREF)
    {
      PV.AllowRefPlot = mymalloc_movable(&PV.AllowRefPlot, "PV.AllowRefPlot", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	PV.AllowRefPlot[i] = PV.AllowRef[i];

      return PV.AllowRefPlot;
    }
  else if(block == PLOT_DIVVEL)
    {
      for(i = 0; i < NumGas; i++)
	PV.DivVel[i] = dabs(PV.DivVel[i]);

      return PV.DivVel;
    }
  else if(block == PLOT_COOL)
    {
      PV.Cool = mymalloc_movable(&PV.Cool, "PV.Cool", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

          abe = PV.AbHII[i];

          abhi = dmax(1. - 2. * PV.AbH2[i] - PV.AbHII[i], 0.);

          // HI electronic excitation cooling
	  cool_rate = 7.5e-19 * exp(-1.18348e5 / PV.Temp[i]) / (1. + sqrt(PV.Temp[i] / 1e5)) * abhi * abe * PV.NH[i] * PV.NH[i];

	  esc_frac = exp(-PV.NH[i] / 1e15);

	  cool_rate *= esc_frac;

	  tcool = PV.U[i] * PV.Rho[i] / cool_rate;

	  if(tcool > 0 && tcool/tff < 20)
	    PV.Cool[i] = tcool / tff;
	  else
	    PV.Cool[i] = -DBL_MAX;
	}

      return PV.Cool;
    }
  else if(block == PLOT_COLLAPSE)
    {
      PV.Collapse = mymalloc_movable(&PV.Collapse, "PV.Collapse", NumGas * sizeof(double));

      for(i = loc_vel_cm_x = loc_vel_cm_y = loc_vel_cm_z = 0; i < NumGas; i++)
	{
	  loc_vel_cm_x += PV.Mass[i] * PV.VX[i];
	  loc_vel_cm_y += PV.Mass[i] * PV.VY[i];
	  loc_vel_cm_z += PV.Mass[i] * PV.VZ[i];
	}

      MPI_Allreduce(&loc_vel_cm_x, &vel_cm_x, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vel_cm_y, &vel_cm_y, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_vel_cm_z, &vel_cm_z, All.RadBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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
	  r = sqrt(PV.X[i] * PV.X[i] + PV.Y[i] * PV.Y[i] + PV.Z[i] * PV.Z[i]);
	  //r = sqrt((PV.X[i] - 70) * (PV.X[i] - 70) + (PV.Y[i] - 5) * (PV.Y[i] - 5) + PV.Z[i] * PV.Z[i]);

	  if(r == 0)
	    PV.Collapse[i] = -DBL_MAX;
	  else
	    {
	      tff = sqrt(3 * M_PI / 32 / GRAVITY / PV.Rho[i]);

	      if(BlockFlag[IO_CHEM])
		mu = PV.Mu[i];
	      else
		mu = MU_NEUTRAL;

	      if(BlockFlag[IO_GAMMA])
		gamma = PV.Gamma[i];
	      else
		gamma = GAMMA_ADB;

	      csnd = sqrt(gamma * BOLTZMANN * PV.Temp[i] / mu / PROTONMASS);

	      tcs = fac * r / csnd;

	      vradx = PV.X[i] * (PV.VX[i] - vel_cm_x);
	      vrady = PV.Y[i] * (PV.VY[i] - vel_cm_y);
	      vradz = PV.Z[i] * (PV.VZ[i] - vel_cm_z);

	      vrad = (vradx + vrady + vradz) / r;

	      tacc = dabs(fac * r / vrad / 1e5);

	      PV.Collapse[i] = dmax(tff / tcs, tff / tacc);
	      PV.Collapse[i] = tff / tcs;
	    }
	}

      return PV.Collapse;
    }
  else if(block == PLOT_SIDM_DENSITY)
    return PV.SIDM_Density;
  else
    {
      terminate("Unknown block in block_val!");

      return 0;
    }
}
