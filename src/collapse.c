
#include "allvars.h"
#include "proto.h"


void collapse()
{
  char buf[MAX_STRING_LEN];
  int i, temp_idx, num_cases, num_entries, num_iter, case_flag, cool_flag;
  double temp, dtemp, h2e20, h2e31, h2n2, h2n3, f, lambda;
  double int_acc, t_init, nh_init, nh_final, temp_init, eta;
  double fH2_init, n_init, gamma_init, u_init, kdh, nh_thresh;
  double nh_prev, p_prev, escfrac, tau_CIE;
  double t, nh, fH2, fHI, n, u, ff_fac;
  double tcool, tff, ratio;
  double dfH2_dt, dt_fH2;
  double dnh_dt, dt_nh, du_dt, dt_u, dt;
  double p, p_init, gamma, gamma_eff;
  double log_nh, log_temp, log_fH2, log_ratio, log_escfrac;
  double temp_table[TGCHEM_NUM_TEMP], cool_table[TGCHEM_NUM_TEMP], dcool_table[TGCHEM_NUM_TEMP], chem_table[TGCHEM_NUM_TEMP], dchem_table[TGCHEM_NUM_TEMP];
  double temp_H2[TGCHEM_NUM_H2], rate_H2_out[TGCHEM_NUM_TEMP];
  double rate_H2_in[TGCHEM_NUM_H2] = 
    { -25.390, -25.086, -24.791, -24.510, -24.245, -23.999,
      -23.769, -23.552, -23.343, -23.138, -22.934, -22.731,
      -22.527, -22.321, -22.109, -21.884, -21.640, -21.378,
      -21.105, -20.832, -20.567, -20.317, -20.085, -19.870,
      -19.673, -19.493, -19.328, -19.177, -19.040, -18.916,
      -18.805, -18.707, -18.621, -18.547, -18.483, -18.428,
      -18.381, -18.341, -18.307, -18.278, -18.253
    };
  FILE *file;

  num_cases = 5;
  num_entries = 5;
  cool_flag = 0;
  int_acc = 1e-4;

  eta = 1;
  nh_init = 1e7;
  nh_final = 1e15;
  temp_init = 580;
  fH2_init = 1e-3;

  t_init = 0;

  h2e20 = 508.95;
  h2e31 = 852.5;

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    temp_table[i] = pow(10.0, log10(TGCHEM_TEMP_MIN) + i * TGCHEM_LOG_DTEMP);

  for(i = 0; i < TGCHEM_NUM_H2; i++)
    temp_H2[i] = pow(10, 2 + 5e-2 * i);

  tgchem_spline_eval(TGCHEM_NUM_H2, &temp_H2[0], &rate_H2_in[0], TGCHEM_NUM_TEMP, &temp_table[0], &rate_H2_out[0]);

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    {
      temp = temp_table[i];

      if(temp < 5)
	cool_table[i] = 1e-60;
      else if(temp < 1e2)
	{
	  h2n2 = 0.25 * (5 * exp(-h2e20 / temp) / (1 + 5 * exp(-h2e20 / temp)));
	  h2n3 = 0.75 * (7 / 3e0 * exp(-h2e31 / temp) / (1 + 7 / 3e0 * exp(-h2e31 / temp)));
	  f = 2.94e-11 * h2e20 * h2n2 * BOLTZMANN + 4.76e-10 * h2e31 * h2n3 * BOLTZMANN;

	  cool_table[i] = dmax(f, 1e-60);
	}
      else if(temp > 1e4)
	cool_table[i] = pow(10, -18.253);
      else
	cool_table[i] = pow(10, rate_H2_out[i]);

      //mpi_printf("T = %g, Lambda = %g\n", temp_table[i], cool_table[i]);
    }

  for(i = 0; i < TGCHEM_NUM_TEMP - 1; i++)
    dcool_table[i] = (cool_table[i + 1] - cool_table[i]) / (temp_table[i + 1] - temp_table[i]);

  dcool_table[TGCHEM_NUM_TEMP - 1] = dcool_table[TGCHEM_NUM_TEMP - 2];

  n_init = (1 + HE_ABUND - fH2_init) * nh_init;
  //gamma_init = GAMMA_ADB;
  gamma_init = 1 / (fH2_init / (1 - fH2_init + HE_ABUND) + 1 / (GAMMA_ADB - 1)) + 1;
  u_init = 1 / (gamma_init - 1) * n_init * BOLTZMANN * temp_init;
  p_init = n_init * BOLTZMANN * temp_init;

  ff_fac = sqrt(32 / 3 / M_PI * GRAVITY * PROTONMASS / HYDROGEN_MASSFRAC);

  sprintf(buf, "collapse.dat");

  if(!(file = fopen(buf, "w")))
    terminate("Could not open file!\n");

  fwrite(&num_cases, sizeof(int), 1, file);
  fwrite(&num_entries, sizeof(int), 1, file);

  for(case_flag = 0; case_flag < num_cases; case_flag++)
    {
      fseek(file, 4, SEEK_CUR);

      num_iter = 0;

      t = t_init;
      nh = nh_init;
      temp = temp_init;
      fH2 = fH2_init;
      u = u_init;
      p = p_init;

      while(nh < nh_final)
	{
	  nh_prev = nh;
	  p_prev = p;

	  dnh_dt = eta * ff_fac * pow(nh, 3. / 2.);

	  dt_nh = int_acc * nh / dnh_dt;

	  tff = 1 / sqrt(nh) / ff_fac;

	  if(cool_flag == 0)
	    lambda = 1e-33 * fH2 * nh * pow(temp, 4);
	  else
	    {
	      temp_idx = imin((int) (log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

	      dtemp = temp - temp_table[temp_idx];

	      lambda = fH2 * nh * (cool_table[temp_idx] + dtemp * dcool_table[temp_idx]);
	    }

	  kdh = 5.24e-7 * pow(temp, -0.485) * exp(-5.2e4 / temp);

	  //escfrac = dmin(pow(nh / 8e9, -0.45), 1);

	  if(case_flag > 3)
	    {
	      nh_thresh = 4e9;

	      if(nh >= nh_thresh)
		escfrac = 1.45 * (nh / nh_thresh) / (pow(nh / nh_thresh, 1.45) + 0.45);
	      else
		escfrac = 1;
	    }
	  else
	    escfrac = 1;

	  lambda *= escfrac;

	  if(case_flag > 2)
	    lambda -= 4.48 * ELECTRON_VOLT * 5.5e-29 * fHI * fHI * fHI * nh * nh * nh / temp;

	  if(case_flag > 2)
	    lambda -= 4.48 * ELECTRON_VOLT * 5.5e-29 / 8 * fHI * fHI * fH2 * nh * nh * nh / temp;

	  if(case_flag > 6)
	    lambda += 4.48 * ELECTRON_VOLT * kdh * fHI * fH2 * nh * nh;

	  if(case_flag > 4)
	    {
	      tau_CIE = dmax(pow(nh / 1.4e16, 2.8), 1e-5);

	      lambda += 2.289e-49 * pow(temp, 4) * nh * (2 * fH2 * nh) * dmin(1 / ((1 + tau_CIE) * (1 + tau_CIE / 10)), 1);
	    }

	  tcool = u / dabs(lambda);

	  ratio = tcool / tff;

	  if(case_flag == 1)
	    gamma = GAMMA_ADB;
	  else
	    gamma = 1 / (fH2 / (1 - fH2 + HE_ABUND) + 1 / (GAMMA_ADB - 1)) + 1;

	  du_dt = gamma * u * dnh_dt / nh - lambda;

	  dt_u = dabs(int_acc * u / du_dt);

	  dt = dmin(dt_nh, dt_u);

	  if(case_flag > 0)
	    {
	      dfH2_dt = 5.5e-29 * fHI * fHI * fHI * nh * nh / temp;

	      if(case_flag > 5)
		dfH2_dt += 5.5e-29 / 8 * fHI * fHI * fH2 * nh * nh / temp;

	      if(case_flag > 6)
		dfH2_dt -= kdh * fHI * fH2 * nh;

	      dt_fH2 = dabs(int_acc * fH2 / dfH2_dt);

	      dt = dmin(dt_fH2, dt);

	      fH2 += dfH2_dt * dt;

	      if(case_flag == 1)
		fHI = 1;
	      else
		{
		  if(fH2 > 0.5)
		    fH2 = 0.5;

		  fHI = dmax(0, 1 - 2 * fH2);
		}
	    }

	  t += dt;

	  nh = pow(-ff_fac * t / 2 + 1 / sqrt(nh_init) , -2);

	  if(case_flag == 1)
	    n = (1 + HE_ABUND) * nh;
	  else
	    n = (1 + HE_ABUND - fH2) * nh;

	  u += du_dt * dt;

	  temp = (gamma - 1) * u / n / BOLTZMANN;

	  p = n * BOLTZMANN * temp;

	  gamma_eff = log10(p / p_prev) / log10(nh / nh_prev);

	  log_nh = log10(nh);
	  log_temp = log10(temp);
	  log_fH2 = log10(fH2);
	  log_ratio = log10(ratio);
	  log_escfrac = log10(escfrac);

	  fwrite(&log_nh, sizeof(double), 1, file);
	  fwrite(&temp, sizeof(double), 1, file);
	  fwrite(&log_fH2, sizeof(double), 1, file);
	  fwrite(&gamma_eff, sizeof(double), 1, file);
	  fwrite(&log_ratio, sizeof(double), 1, file);

      //printf("%g %g %g %g %g %g\n", p, p_init, nh, nh_init, n, n_init);

      //printf("%g %g %g %g %g %g %g\n", nh, temp, fH2, gamma_eff, dt_fH2, dt_u, dt_nh);

	  //if(case_flag == 6 && num_iter > 150000)
	    //printf("num_iter = %d, dt = %g, nh = %g, temp = %g, gamma = %g, lambda = %g\n", num_iter, dt, nh, temp, GAMMA_ADB * u * dnh_dt / nh, lambda);

	  num_iter++;
	}

      fseek(file, -num_iter * num_entries * 8 - 4, SEEK_CUR);
      fwrite(&num_iter, sizeof(int), 1, file);
      fseek(file, num_iter * num_entries * 8, SEEK_CUR);
    }

  printf("iters = %d, nh = %g, temp = %g, fH2 = %g, tcool / tff = %g\n", num_iter, nh, temp, fH2, ratio);

  fclose(file);
}
