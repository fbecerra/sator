
#include "allvars.h"
#include "proto.h"


void make_image()
{
  char sub_string[MAX_STRING_LEN], outfile[MAX_STRING_LEN], label[MAX_STRING_LEN];
  int pstart, pend;
  int i, j, k, sub, idx, block, dim, num_blocks, max_string_len;
  int min1, max1, min2, max2;
  double x, y, ax1, ax2, ax3, rho2, sum_part, *sum, *glob_sum, img_part, *img, *glob_img, val, *valp, hsml, fac;
  double width, height, size;
  FILE *file;

  if(!BlockFlag[IO_POS] || !BlockFlag[IO_ID] || !BlockFlag[IO_MASS] || !BlockFlag[IO_RHO]) // || !BlockFlag[IO_VOL])
    terminate("Required block not read!");

  if(All.ImgFlagMovie)
    sprintf(All.Ending, "mov");
  else
    sprintf(All.Ending, "img");

  if(All.ImgFlagScreenRatio == 1)
    All.ImgScreenRatio = 4. / 3.;
  else
    All.ImgScreenRatio = 1;

  All.ImgYBins = All.ImgXBins / All.ImgScreenRatio;
  All.ImgHeight = All.ImgWidth / All.ImgScreenRatio;
  All.ImgSize = dmax(All.ImgWidth, All.ImgHeight);

  max_string_len = MAX_STRING_LEN;

  if(All.FlagCenter != 1 || All.Iter == 0)
    box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
 //     if (PV.ID[i] == 1048117401)
 //       idx = i;
    }

  sum = (double *) mymalloc_movable(&sum, "sum", All.ImgXBins * All.ImgYBins * sizeof(double));
  img = (double *) mymalloc_movable(&img, "img", All.ImgXBins * All.ImgYBins * sizeof(double));

  glob_sum = (double *) mymalloc_movable(&glob_sum, "glob_sum", All.ImgXBins * All.ImgYBins * sizeof(double));
  glob_img = (double *) mymalloc_movable(&glob_img, "glob_img", All.ImgXBins * All.ImgYBins * sizeof(double));

  for(block = num_blocks = 0; block < PLOT_NUM_BLOCKS; block++)
    if(plot_block_present(block))
      num_blocks++;

  for(sub = 0; sub < All.ImgNumSubSnaps; sub++)
    {
      if(All.FlagRotate)
	box_rotate();

      width = pow(10, log10(All.ImgWidth) - sub * log10(All.ImgWidth / All.ImgSubWidth) / imax((All.ImgNumSubSnaps - 1), 1));

      height = width / All.ImgScreenRatio;

      size = dmax(width, height);

      if(All.ImgNumSubSnaps > 1)
	sprintf(sub_string, "_%03d", sub);
      else
	sprintf(sub_string, "");

      if(ThisTask == 0)
	{
	  sprintf(outfile, "%s/%s/%s_%s%s.%s", All.Path, All.Base, All.Base, All.SnapNumString, sub_string, All.Ending);
//          sprintf(outfile, "/n/hernquistfs3/fbecerra/nahw1rs2/%s_%s%s.%s", All.Base, All.SnapNumString, sub_string, All.Ending);

	  if(!(file = fopen(outfile, "w")))
	    terminate("Could not write image file!");

	  fwrite(&num_blocks, sizeof(int), 1, file);
	  fwrite(&max_string_len, sizeof(int), 1, file);
	  fwrite(&All.ImgXBins, sizeof(int), 1, file);
	  fwrite(&All.ImgYBins, sizeof(int), 1, file);
	  fwrite(&All.TotNumSinks, sizeof(int), 1, file);
	  fwrite(&width, sizeof(double), 1, file);
	  fwrite(&height, sizeof(double), 1, file);
	  fwrite(&All.Time, sizeof(double), 1, file);
	  fwrite(&All.LengthUnit, sizeof(int), 1, file);
	  fwrite(&All.ComovingUnits, sizeof(int), 1, file);
	}

      if(All.TotNumSinks)
	{
	  Sink = mymalloc_movable(&Sink, "Sink", All.TotNumSinks * sizeof(struct Sink_struct));

	  for(i = 0; i < NumSinks; i++)
	    {
	      idx = TotNumPart - NumSinks + i;

	      Sink[i].Task = 0;

	      Sink[i].X = PV.X[idx];
	      Sink[i].Y = PV.Y[idx];
	      Sink[i].Z = PV.Z[idx];
              mpi_printf("idx=%i, sx=%f, sy=%f, sz=%f, size=%f, ID=%i \n", idx, Sink[i].X, Sink[i].Y, Sink[i].Z, size, PV.ID[idx]);
	    }

	  mpi_distribute_items_to_tasks(Sink, offsetof(struct Sink_struct, Task), &NumSinks, &All.TotNumSinks, sizeof(struct Sink_struct), 0);

	  if(ThisTask == 0)
	    for(i = 0; i < All.TotNumSinks; i++)
	      if(dabs(Sink[i].X) < size / 2 && dabs(Sink[i].Y) < size / 2 && dabs(Sink[i].Z) < size / 2)
		{
		  x = (Sink[i].X + 0.5 * width);
		  y = (Sink[i].Y + 0.5 * height);
                  mpi_printf("sx=%f, sy=%f, sz=%f, x=%f, y=%f \n", Sink[i].X, Sink[i].Y, Sink[i].Z, x, y);

		  fwrite(&x, sizeof(double), 1, file);
		  fwrite(&y, sizeof(double), 1, file);
		}

	  myfree_movable(Sink);
	}

      for(block = 0; block < PLOT_NUM_BLOCKS; block++)
	{
	  if(!plot_block_present(block))
	    continue;

	  memset(label, 0, MAX_STRING_LEN * sizeof(char));

	  plot_block_label(block, label);

	  valp = plot_block_val(block);

	  if(ThisTask == 0)
	    fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);

	  if(block == PLOT_SIDM_DENSITY)
	    {
	      pstart = NumGas;
	      pend = NumGas + NumPart[1];
	    }
	  else
	    {
	      pstart = 0;
	      pend = NumGas;
	    }

	  for(dim = 0; dim < 3; dim++)
	    {
	      if(ThisTask == 0)
		fseek(file, 2 * sizeof(double), SEEK_CUR);

	      memset(sum, 0, All.ImgXBins * All.ImgYBins * sizeof(double));
	      memset(img, 0, All.ImgXBins * All.ImgYBins * sizeof(double));

	      for(i = pstart; i < pend; i++)
		{
		  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
		    continue;

		  if(block == PLOT_SIDM_DENSITY)
		    {
		      val = *(valp + i - NumGas);
		      hsml = PV.SIDM_Hsml[i - NumGas];
		    }
		  else
		    {
		      val = *(valp + i);

                      if(All.LengthUnit == 0)
                        fac = 1;
                      else if(All.LengthUnit == 1)
                        fac = 1e3;
                      else if(All.LengthUnit == 2)
                        fac = UNIT_LENGTH / ASTRONOMICAL_UNIT;
                      else
                        terminate("Length unit not implemented!");

		      hsml = fac * pow(3. / 4 / M_PI, 1. / 3) * pow((PV.Mass[i] * SOLAR_MASS / PV.Rho[i]), 1. / 3) / UNIT_LENGTH; //PV.Hsml[i];
		    }

		  if(val == -DBL_MAX)
		    continue;

		  if(dim == 0)
		    {
		      ax1 = PV.X[i];
		      ax2 = PV.Y[i];
		      ax3 = PV.Z[i];
		    }
		  else if(dim == 1)
		    {
		      ax1 = PV.X[i];
		      ax2 = PV.Z[i];
		      ax3 = PV.Y[i];
		    }
		  else
		    {
		      ax1 = PV.Y[i];
		      ax2 = PV.Z[i];
		      ax3 = PV.X[i];
		    }

		  if(dabs(ax1) < width / 2 && dabs(ax2) < height / 2 && dabs(ax3) < size / 2)
		    {
		      min1 = (ax1 + width / 2 - hsml) / width * All.ImgXBins;
		      min1 = imin(imax(min1, 0), All.ImgXBins - 1);

		      max1 = (ax1 + width / 2 + hsml) / width * All.ImgXBins;
		      max1 = imin(imax(max1, 0), All.ImgXBins - 1);

		      min2 = (ax2 + height / 2 - hsml) / height * All.ImgYBins;
		      min2 = imin(imax(min2, 0), All.ImgYBins - 1);

		      max2 = (ax2 + height / 2 + hsml) / height * All.ImgYBins;
		      max2 = imin(imax(max2, 0), All.ImgYBins - 1);

		      if(block == PLOT_SIDM_DENSITY)
			rho2 = PV.SIDM_Density[i - NumGas] * PV.SIDM_Density[i - NumGas];
		      else
			rho2 = PV.NH[i] * PV.NH[i];

		      sum_part = rho2 * PV.Mass[i];

		      img_part = sum_part * val;

		      for(j = min1; j <= max1; j++)
			for(k = min2; k <= max2; k++)
			  {
			    idx = All.ImgYBins * j + k;

			    sum[idx] += sum_part;
			    img[idx] += img_part;
			  }
		    }
		}

	      MPI_Allreduce(sum, glob_sum, All.ImgXBins * All.ImgYBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	      MPI_Allreduce(img, glob_img, All.ImgXBins * All.ImgYBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	      for(i = 0; i < All.ImgXBins * All.ImgYBins; i++)
		if(glob_img[i] != 0)
		  {
		    glob_img[i] = log10(glob_img[i] / glob_sum[i]);

		    All.ImgMin[block][dim] = dmin(glob_img[i], All.ImgMin[block][dim]);
		    All.ImgMax[block][dim] = dmax(glob_img[i], All.ImgMax[block][dim]);
		  }

	      if(ThisTask == 0)
		fwrite(glob_img, All.ImgXBins * All.ImgYBins * sizeof(double), 1, file);
	    }

	  if(block == PLOT_PDVRATE)
	    myfree_movable(PV.PdVCool);

	  if(block == PLOT_H2RATE)
	    myfree_movable(PV.H2Cool);

	  if(block == PLOT_CIERATE)
	    myfree_movable(PV.CIECool);

	  if(block == PLOT_CHEMRATE)
	    myfree_movable(PV.ChemCool);

	  if(block == PLOT_ESCFRAC)
	    myfree_movable(PV.EscFrac);

	  if(block == PLOT_ALLOWREF)
	    myfree_movable(PV.AllowRefPlot);

	  if(block == PLOT_COOL)
	    myfree_movable(PV.Cool);

	  if(block == PLOT_COLLAPSE)
	    myfree_movable(PV.Collapse);
	}

      if(ThisTask == 0)
	fclose(file);

      if(All.ImgNumSubSnaps > 1)
	mpi_printf("SATOR: Sub %d done!\n", sub);
    }

  myfree_movable(glob_img);
  myfree_movable(glob_sum);
  myfree_movable(img);
  myfree_movable(sum);
}


write_min_max_vals()
{
  char sub_string[MAX_STRING_LEN], outfile[MAX_STRING_LEN];
  int i, sub, block, dim;
  FILE *file;

  if(ThisTask == 0)
    {
      for(i = 0; i < All.NumSnaps; i++)
	{
	  for(sub = 0; sub < All.ImgNumSubSnaps; sub++)
	    {
	      if(All.ImgNumSubSnaps > 1)
		sprintf(sub_string, "_%03d", sub);
	      else
		sprintf(sub_string, "");

	      sprintf(outfile, "%s/%s/%s_%s%s.%s", All.Path, All.Base, All.Base, All.SnapNumString, sub_string, All.Ending);

	      if(!(file = fopen(outfile, "r+")))
		terminate("Could not write image file!");

	      fseek(file, 52 + NumSinksIter[i] * 2 * sizeof(double), SEEK_SET);

	      for(block = 0; block < PLOT_NUM_BLOCKS; block++)
		{
		  if(!plot_block_present(block))
		    continue;

		  fseek(file, MAX_STRING_LEN * sizeof(char), SEEK_CUR);

		  for(dim = 0; dim < 3; dim++)
		    {
		      fwrite(&All.ImgMin[block][dim], sizeof(double), 1, file);
		      fwrite(&All.ImgMax[block][dim], sizeof(double), 1, file);
		      fseek(file, All.ImgXBins * All.ImgYBins * sizeof(double), SEEK_CUR);
		    }
		}

	      fclose(file);
	    }
	}
    }
}
