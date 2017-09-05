
#include "allvars.h"
#include "proto.h"


void pspace()
{
  char pspace_file[MAX_STRING_LEN], label[MAX_STRING_LEN];
  int i, j, block, num_blocks, idx, xidx, yidx, collapse_weight, max_string_len;
  double *valp, *valxp, *valyp, *valx, *weightp;
  double *loc_valy, *valy, *loc_valy_sum, *valy_sum;
  double *loc_val, *val, *loc_val_sum, *val_sum, tot_val;
  double loc_valx_min, valx_min, loc_valx_max, valx_max;
  double loc_valy_min, valy_min, loc_valy_max, valy_max;
  double logvalx, logvaly, xrange, yrange;
  FILE *file;

  collapse_weight = 0;

  max_string_len = MAX_STRING_LEN;

  if(!BlockFlag[IO_ID] || !BlockFlag[IO_MASS] || !BlockFlag[IO_RHO])
    terminate_block();

  box_center();

  for(i = 0; i < TotNumPart; i++)
    {
      PV.X[i] -= All.BoxCenter[0];
      PV.Y[i] -= All.BoxCenter[1];
      PV.Z[i] -= All.BoxCenter[2];
    }

  if(All.FlagRotate)
    box_rotate();

  valxp = PV.NH;

  for(i = 0, loc_valx_min = DBL_MAX, loc_valx_max = 0; i < NumGas; i++)
    {
      if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	continue;

      if(*(valxp + i) <= 0)
	terminate("Something wrong with x value!");

      if(*(valxp + i) < loc_valx_min)
	loc_valx_min = *(valxp + i);

      if(*(valxp + i) > loc_valx_max)
	loc_valx_max = *(valxp + i);
    }

  MPI_Allreduce(&loc_valx_min, &valx_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_valx_max, &valx_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(All.PSpaceXMin)
    valx_min = All.PSpaceXMin;

  if(All.PSpaceXMax)
    valx_max = All.PSpaceXMax;

  valx_min = log10(valx_min);
  valx_max = log10(valx_max);

  xrange = (valx_max - valx_min);

  valx = (double *) mymalloc_movable(&valx, "valx", All.PSpaceBins * sizeof(double));

  loc_valy = (double *) mymalloc_movable(&loc_valy, "loc_valy", All.PSpaceBins * sizeof(double));
  valy = (double *) mymalloc_movable(&valy, "valy", All.PSpaceBins * sizeof(double));

  loc_valy_sum = (double *) mymalloc_movable(&loc_valy_sum, "loc_valy_sum", All.PSpaceBins * sizeof(double));
  valy_sum = (double *) mymalloc_movable(&valy_sum, "valy_sum", All.PSpaceBins * sizeof(double));

  loc_val = (double *) mymalloc_movable(&loc_val, "loc_val", All.PSpaceBins * All.PSpaceBins * sizeof(double));
  val = (double *) mymalloc_movable(&val, "val", All.PSpaceBins * All.PSpaceBins * sizeof(double));

  loc_val_sum = (double *) mymalloc_movable(&loc_val_sum, "loc_val_sum", All.PSpaceBins * All.PSpaceBins * sizeof(double));
  val_sum = (double *) mymalloc_movable(&val_sum, "val_sum", All.PSpaceBins * All.PSpaceBins * sizeof(double));

  for(i = 0; i < All.PSpaceBins; i++)
    valx[i] = valx_min + (i + 0.5) * xrange / All.PSpaceBins;

  for(block = num_blocks = 0; block < PLOT_NUM_BLOCKS; block++)
    if(plot_block_present(block))
      if(block != PLOT_NH)
	num_blocks++;

  if(ThisTask == 0)
    {
      sprintf(pspace_file, "%s/%s/%s_%s.pspace", All.Path, All.Base, All.Base, All.SnapNumString);

      if(!(file = fopen(pspace_file, "w")))
	terminate("Could not write phase space file!");

      memset(label, 0, MAX_STRING_LEN * sizeof(char));
      plot_block_label(PLOT_NH, label);

      fwrite(&num_blocks, sizeof(int), 1, file);
      fwrite(&max_string_len, sizeof(int), 1, file);
      fwrite(&All.PSpaceBins, sizeof(int), 1, file);
      fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
      fwrite(&valx_min, sizeof(double), 1, file);
      fwrite(&valx_max, sizeof(double), 1, file);
      fwrite(valx, All.PSpaceBins * sizeof(double), 1, file);
    }

  for(block = 0; block < PLOT_NUM_BLOCKS; block++)
    {
      if(!plot_block_present(block))
	continue;

      if(block == PLOT_NH)
	continue;

      memset(label, 0, MAX_STRING_LEN * sizeof(char));

      plot_block_label(block, label);

      valyp = plot_block_val(block);

      if(collapse_weight)
	{
	  if(block == PLOT_COLLAPSE)
	    weightp = valyp;
	  else
	    weightp = plot_block_val(PLOT_COLLAPSE);
	}

      for(i = 0, loc_valy_min = DBL_MAX, loc_valy_max = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if(*(valyp + i) < 0 && *(valyp + i) != -DBL_MAX)
	    {
	      mpi_printf("%d %g\n", block, *(valyp + i));
	      terminate("Something wrong with y value!");
	    }

	  if(*(valyp + i) != 0 && *(valyp + i) != -DBL_MAX && *(valyp + i) < loc_valy_min)
	    loc_valy_min = *(valyp + i);

	  if(*(valyp + i) != 0 && *(valyp + i) != -DBL_MAX && *(valyp + i) > loc_valy_max)
	    loc_valy_max = *(valyp + i);
	}

      MPI_Allreduce(&loc_valy_min, &valy_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_valy_max, &valy_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if(All.PSpaceYMin)
	valy_min = All.PSpaceYMin;

      if(All.PSpaceYMax)
	valy_max = All.PSpaceYMax;

      if(block != PLOT_GAMMA)
	{
	  valy_min = log10(valy_min);
	  valy_max = log10(valy_max);
	}

      yrange = (valy_max - valy_min);

      memset(loc_valy, 0, All.PSpaceBins * sizeof(double));
      memset(loc_valy_sum, 0, All.PSpaceBins * sizeof(double));
      memset(loc_val, 0, All.PSpaceBins * All.PSpaceBins * sizeof(double));
      memset(loc_val_sum, 0, All.PSpaceBins * All.PSpaceBins * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  if(PV.Mass[i] == 0 && PV.ID[i] == 0)
	    continue;

	  if(BlockFlag[IO_ALLOWREF])
	    if(!PV.AllowRef[i])
	      continue;

	  if(*(valyp + i) == 0 || *(valyp + i) == -DBL_MAX)
	    continue;

	  logvalx = log10(*(valxp + i));

	  if(block == PLOT_GAMMA)
	    logvaly = *(valyp + i);
	  else
	    logvaly = log10(*(valyp + i));

	  if(logvalx >= valx_min && logvalx <= valx_max)
	    if(logvaly >= valy_min && logvaly <= valy_max)
	      {
		xidx = (logvalx - valx_min) / xrange * All.PSpaceBins;
		xidx = imax(imin(xidx, All.PSpaceBins - 1), 0);

		yidx = (logvaly - valy_min) / yrange * All.PSpaceBins;
		yidx = imax(imin(yidx, All.PSpaceBins - 1), 0);

		loc_valy[xidx] += PV.Mass[i] * logvaly;
		loc_valy_sum[xidx] += PV.Mass[i];

		if(collapse_weight)
		  loc_val[All.PSpaceBins * xidx + yidx] += *(weightp + i);
		else
		  loc_val[All.PSpaceBins * xidx + yidx] += PV.Mass[i];

		loc_val_sum[All.PSpaceBins * xidx + yidx] += 1;
	      }
	}

      MPI_Allreduce(loc_valy, valy, All.PSpaceBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_valy_sum, valy_sum, All.PSpaceBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_val, val, All.PSpaceBins * All.PSpaceBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_val_sum, val_sum, All.PSpaceBins * All.PSpaceBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = tot_val = 0; i < All.PSpaceBins * All.PSpaceBins; i++)
	tot_val += val[i];

      for(i = 0; i < All.PSpaceBins; i++)
	{
	  if(valy_sum[i] == 0)
	    valy[i] = -DBL_MAX;
	  else
	    valy[i] /= valy_sum[i];

	  for(j = 0; j < All.PSpaceBins; j++)
	    {
	      idx = i * All.PSpaceBins + j;

	      if(val_sum[idx] == 0)
		val[idx] = -DBL_MAX;
	      else
		{
		  if(collapse_weight)
		    val[idx] = log10(val[idx] / val_sum[idx]);
		  else
		    val[idx] = log10(val[idx] / tot_val);
		}


	    }
	}

      if(ThisTask == 0)
	{
	  fwrite(label, MAX_STRING_LEN * sizeof(char), 1, file);
	  fwrite(&valy_min, sizeof(double), 1, file);
	  fwrite(&valy_max, sizeof(double), 1, file);
	  fwrite(valy, All.PSpaceBins * sizeof(double), 1, file);
	  fwrite(val, All.PSpaceBins * All.PSpaceBins * sizeof(double), 1, file);
	}

      if(collapse_weight && block != PLOT_COLLAPSE)
	myfree_movable(PV.Collapse);

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

  myfree_movable(val_sum);
  myfree_movable(loc_val_sum);
  myfree_movable(val);
  myfree_movable(loc_val);
  myfree_movable(valy_sum);
  myfree_movable(loc_valy_sum);
  myfree_movable(valy);
  myfree_movable(loc_valy);
  myfree_movable(valx);
}
