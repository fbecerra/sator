
#include "allvars.h"
#include "proto.h"


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  init(argc, argv);

  for(All.Iter = 0; All.Iter < All.NumSnaps; All.Iter++)
    {
      if(All.Usage != 10)
	read_snap();

      if(All.Usage == 0)
	{
	  make_image();

	  if(All.Iter == All.NumSnaps - 1)
	    write_min_max_vals();
	}

      if(All.Usage == 1)
	pspace();

      if(All.Usage == 2)
	radial();

      if(All.Usage == 3)
	mbe();

      if(All.Usage == 4)
	halo();

      if(All.Usage == 8)
        remove_part();

      if(All.Usage == 9)
	cut_region();

      if(All.Usage == 10)
	collapse();

      if(All.Usage != 10)
	clear_snap();

      mpi_printf("SATOR: Snap %d done!\n", All.SnapNum);
    }

  finish();

  mpi_printf("SATOR: Took %g sec. Done!\n", measure_time());

  MPI_Finalize();
}
