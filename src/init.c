
#include "allvars.h"
#include "proto.h"

void finish();
void get_args(int argc, char **argv);
void allocate_mpi_buffers();
void setup_H2();
void init_snap();
void free_H2();
void free_mpi_buffers();


void init(int argc, char ** argv)
{
  WallClockTime = second();

  All.FlagRotate = 0;

  get_args(argc, argv);

  read_params();

  mymalloc_init();

  allocate_mpi_buffers();

/*  setup_H2(); */

  init_snap();
}


void finish()
{
  myfree_movable(NumSinksIter);

/*  free_H2(); */

  free_mpi_buffers();
}


void get_args(int argc, char **argv)
{
  if(argc < 3)
    terminate_args();

  strcpy(All.ParamFile, argv[1]);

  All.Usage = atoi(argv[2]);

  if(All.Usage == 10)
    {
      All.SnapNumStart = 0;
      All.SnapNumEnd = 0;
    }
  else
    {
      if(argc < 5)
	terminate_args();

      strcpy(All.Base, argv[3]);

      All.SnapNumStart = atoi(argv[4]);
    }

  if(All.Usage == 0)
    {
      if(argc < 7)
	terminate_args();

      All.SnapNumEnd = atoi(argv[5]);

      All.ImgWidth = atof(argv[6]);

      if(argc > 7)
	{
	  if(argc < 9)
	    terminate_args();

	  All.ImgSubWidth = atof(argv[7]);
	  All.ImgNumSubSnaps = atoi(argv[8]);
	}
      else
	{
	  All.ImgSubWidth = All.ImgWidth;
	  All.ImgNumSubSnaps = 1;
	}

      if(argc > 9)
	All.FlagRotate = atoi(argv[9]);
    }

  if(All.Usage == 1 || All.Usage == 3)
    {
      if(argc > 5)
	All.SnapNumEnd = atoi(argv[5]);
      else
	All.SnapNumEnd = -1;
    }

  if(All.Usage == 2)
    {
      if(argc < 7)
	terminate_args();

      All.SnapNumEnd = atoi(argv[5]);

      All.RadMin = atof(argv[6]);
    }

  if(All.Usage == 3)
    {
      if(argc < 10)
	terminate_args();

      All.SnapNumEnd = atoi(argv[5]);

      All.MBEMin = atof(argv[6]);

      All.ShiftX = atof(argv[7]);
      All.ShiftY = atof(argv[8]);
      All.ShiftZ = atof(argv[9]);
    }

  if(All.Usage == 8)
    {
      if(argc < 6)
        terminate_args();

      strcpy(All.RemoveBase, argv[5]);
    }

  if(All.Usage == 9)
    {
      if(argc < 7)
	terminate_args();

      strcpy(All.CutBase, argv[5]);

      All.CutSize = atof(argv[6]);
    }

  if(All.SnapNumEnd < All.SnapNumStart)
    All.SnapNumEnd = All.SnapNumStart;

  All.NumSnaps = All.SnapNumEnd - All.SnapNumStart + 1;
}


void allocate_mpi_buffers()
{
  SendCount = mymalloc_movable(&SendCount, "SendCount", NTask * sizeof(int));
  SendOffset = mymalloc_movable(&SendOffset, "SendOffset", NTask * sizeof(int));
  RecvCount = mymalloc_movable(&RecvCount, "RecvCount", NTask * sizeof(int));
  RecvOffset = mymalloc_movable(&RecvOffset, "RecvOffset", NTask * sizeof(int));
}


void setup_H2()
{
  char buf[MAX_STRING_LEN];
  int i, j, count;
  double h2e20, h2e31, h2n2, h2n3, f, temp;
  double temp_H2[TGCHEM_NUM_H2], rate_H2_in[TGCHEM_NUM_H2] = 
    { -25.390, -25.086, -24.791, -24.510, -24.245, -23.999,
      -23.769, -23.552, -23.343, -23.138, -22.934, -22.731,
      -22.527, -22.321, -22.109, -21.884, -21.640, -21.378,
      -21.105, -20.832, -20.567, -20.317, -20.085, -19.870,
      -19.673, -19.493, -19.328, -19.177, -19.040, -18.916,
      -18.805, -18.707, -18.621, -18.547, -18.483, -18.428,
      -18.381, -18.341, -18.307, -18.278, -18.253
    };
  FILE *file;

  All.TempTable = mymalloc_movable(&All.TempTable, "All.TempTable", TGCHEM_NUM_TEMP * sizeof(double));
  All.CoolTable = mymalloc_movable(&All.CoolTable, "All.CoolTable", TGCHEM_NUM_TEMP * sizeof(double));
  All.DCoolTable = mymalloc_movable(&All.DCoolTable, "All.DCoolTable", TGCHEM_NUM_TEMP * sizeof(double));
  All.ChemTable = mymalloc_movable(&All.ChemTable, "All.ChemTable", TGCHEM_NUM_TEMP * sizeof(double));
  All.DChemTable = mymalloc_movable(&All.DChemTable, "All.DChemTable", TGCHEM_NUM_TEMP * sizeof(double));
  All.H2RateOut = mymalloc_movable(&All.H2RateOut, "All.H2RateOut", TGCHEM_NUM_TEMP * sizeof(double));
  All.H2OpacTemp = mymalloc_movable(&All.H2OpacTemp, "All.H2OpacTemp", TGCHEM_NUM_H2OP * sizeof(float));
  All.H2OpacColumn = mymalloc_movable(&All.H2OpacColumn, "All.H2OpacColumn", TGCHEM_NUM_H2OP * sizeof(float));
  All.H2Opac = mymalloc_movable(&All.H2Opac, "All.H2Opac", TGCHEM_NUM_H2OP * TGCHEM_NUM_H2OP * sizeof(float));

  h2e20 = 508.95;
  h2e31 = 852.5;

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    All.TempTable[i] = pow(10.0, log10(TGCHEM_TEMP_MIN) + i * TGCHEM_LOG_DTEMP);

  for(i = 0; i < TGCHEM_NUM_H2; i++)
    temp_H2[i] = pow(10, 2 + 5e-2 * i);

  tgchem_spline_eval(TGCHEM_NUM_H2, temp_H2, rate_H2_in, TGCHEM_NUM_TEMP, All.TempTable, All.H2RateOut);

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    {
      temp = All.TempTable[i];

      if(temp < 5)
	All.CoolTable[i] = 1e-60;
      else if(temp < 1e2)
	{
	  h2n2 = 0.25 * (5 * exp(-h2e20 / temp) / (1 + 5 * exp(-h2e20 / temp)));
	  h2n3 = 0.75 * (7 / 3e0 * exp(-h2e31 / temp) / (1 + 7 / 3e0 * exp(-h2e31 / temp)));
	  f = 2.94e-11 * h2e20 * h2n2 * BOLTZMANN + 4.76e-10 * h2e31 * h2n3 * BOLTZMANN;

	  All.CoolTable[i] = dmax(f, 1e-60);
	}
      else if(temp > 1e4)
	All.CoolTable[i] = pow(10, -18.253);
      else
	All.CoolTable[i] = pow(10, All.H2RateOut[i]);
    }

  for(i = 0; i < TGCHEM_NUM_TEMP - 1; i++)
    All.DCoolTable[i] = (All.CoolTable[i + 1] - All.CoolTable[i]) / (All.TempTable[i + 1] - All.TempTable[i]);

  All.DCoolTable[TGCHEM_NUM_TEMP - 1] = All.DCoolTable[TGCHEM_NUM_TEMP - 2];

  if(ThisTask == 0)
    {
      strcat(buf, "data/tgchem_h2c_table.dat");

      if(!(file = fopen(buf, "r")))
	{
	  printf("Error opening file!\n");
	  abort();
	}

      for(i = 0; i < TGCHEM_NUM_H2OP; i++)
	for(j = 0; j < TGCHEM_NUM_H2OP; j++)
	  fscanf(file, "%g %g %g\n", &All.H2OpacTemp[i], &All.H2OpacColumn[j], &All.H2Opac[i * TGCHEM_NUM_H2OP + j]);

      fclose(file);
    }

  MPI_Bcast(All.H2OpacTemp, TGCHEM_NUM_H2OP, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(All.H2OpacColumn, TGCHEM_NUM_H2OP, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(All.H2Opac, TGCHEM_NUM_H2OP * TGCHEM_NUM_H2OP, MPI_FLOAT, 0, MPI_COMM_WORLD);
}


void init_snap()
{
  int block, dim, read_ok;
  char fname[MAX_STRING_LEN], buf[MAX_STRING_LEN];
  FILE *file;

  sprintf(All.BlockFile, "data/blocks.txt");
  sprintf(All.BlockFileCut, "data/blocks_cut.txt");
#ifdef OLDCHEM
  sprintf(All.BlockFilePlot, "data/blocks_plot_old.txt");
#else
  sprintf(All.BlockFilePlot, "data/blocks_plot.txt");
#endif

  memset(BlockPresent, 0, NUM_BLOCKS * sizeof(int));
  memset(BlockRead, 0, NUM_BLOCKS * sizeof(int));
  memset(BlockFlag, 0, NUM_BLOCKS * sizeof(int));

  memset(BlockPlot, 0, PLOT_NUM_BLOCKS * sizeof(int));

  if(ThisTask == 0)
    {
      if(All.Usage == 8 || All.Usage == 9)
	sprintf(fname, "%s", All.BlockFileCut);
      else
	sprintf(fname, "%s", All.BlockFile);

      if(!(file = fopen(fname, "r")))
	{
	  sprintf(buf, "Can't open blocks file `%s'!", fname);

	  terminate(buf);
	}

      for(block = 0; block < NUM_BLOCKS; block++)
	{
	  read_ok = fscanf(file, "%s %d %d\n", buf, &BlockPresent[block], &BlockRead[block]);

	  if(read_ok < 0)
	    terminate("Blocks missing in blocks file!");
	}

      fclose(file);

      if(!(file = fopen(All.BlockFilePlot, "r")))
	{
	  sprintf(buf, "Can't open plot blocks file `%s'!", All.BlockFilePlot);

	  terminate(buf);
	}

      for(block = 0; block < PLOT_NUM_BLOCKS; block++)
	{
	  read_ok = fscanf(file, "%s %d\n", buf, &BlockPlot[block]);

	  if(read_ok < 0)
	    terminate("Blocks missing in plot blocks file!");
	}

      fclose(file);
    }

  MPI_Bcast(BlockPresent, NUM_BLOCKS * sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(BlockRead, NUM_BLOCKS * sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  for(block = 0; block < NUM_BLOCKS; block++)
    if(BlockPresent[block] && BlockRead[block])
      BlockFlag[block] = 1;

  MPI_Bcast(BlockPlot, PLOT_NUM_BLOCKS * sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  memset(PD, 0, NUM_BLOCKS * sizeof(void));
  memset(&PV, 0, sizeof(struct PartValues));

  NumSinksIter = mymalloc_movable(&NumSinksIter, "NumSinksIter", All.NumSnaps * sizeof(int));

  for(block = 0; block < PLOT_NUM_BLOCKS; block++)
    for(dim = 0; dim < 3; dim++)
      {
	All.ImgMin[block][dim] = DBL_MAX;
	All.ImgMax[block][dim] = -DBL_MAX;
      }
}


void free_H2()
{
  myfree_movable(All.TempTable);
  myfree_movable(All.CoolTable);
  myfree_movable(All.DCoolTable);
  myfree_movable(All.ChemTable);
  myfree_movable(All.DChemTable);
  myfree_movable(All.H2RateOut);
  myfree_movable(All.H2OpacTemp);
  myfree_movable(All.H2OpacColumn);
  myfree_movable(All.H2Opac);
}


void free_mpi_buffers()
{
  myfree_movable(RecvOffset);
  myfree_movable(RecvCount);
  myfree_movable(SendOffset);
  myfree_movable(SendCount);
}
