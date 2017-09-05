
#include "allvars.h"
#include "proto.h"

#define INT 1
#define REAL 2
#define STRING 3


void read_params()
{
  if(ThisTask == 0)
    {
      int i = 0, max_tags, nargs, id[MAX_PARAMS];
      char buf[MAX_STRING_LEN], buf1[MAX_STRING_LEN], buf2[MAX_STRING_LEN];
      char tag[MAX_PARAMS][MAX_STRING_LEN];
      void *addr[MAX_PARAMS];
      FILE *fd;

      strcpy(tag[i], "MaxMemSize");
      addr[i] = &All.MaxMemSize;
      id[i++] = INT;

      strcpy(tag[i], "FileSystem");
      addr[i] = &All.FileSystem;
      id[i++] = INT;

      strcpy(tag[i], "FlagDouble");
      addr[i] = &All.FlagDouble;
      id[i++] = INT;

      strcpy(tag[i], "LengthUnit");
      addr[i] = &All.LengthUnit;
      id[i++] = INT;

      strcpy(tag[i], "ComovingUnits");
      addr[i] = &All.ComovingUnits;
      id[i++] = INT;

      strcpy(tag[i], "FlagCenter");
      addr[i] = &All.FlagCenter;
      id[i++] = INT;

      strcpy(tag[i], "CenterX");
      addr[i] = &All.CenterX;
      id[i++] = REAL;

      strcpy(tag[i], "CenterY");
      addr[i] = &All.CenterY;
      id[i++] = REAL;

      strcpy(tag[i], "CenterZ");
      addr[i] = &All.CenterZ;
      id[i++] = REAL;

      if(All.FlagRotate == 0)
	{
	  strcpy(tag[i], "FlagRotate");
	  addr[i] = &All.FlagRotate;
	  id[i++] = INT;
	}

      if(All.Usage == 0)
	{
	  strcpy(tag[i], "ImgFlagScreenRatio");
	  addr[i] = &All.ImgFlagScreenRatio;
	  id[i++] = INT;

	  strcpy(tag[i], "ImgXBins");
	  addr[i] = &All.ImgXBins;
	  id[i++] = INT;

	  strcpy(tag[i], "ImgFlagMovie");
	  addr[i] = &All.ImgFlagMovie;
	  id[i++] = INT;
	}

      if(All.Usage == 1)
	{
	  strcpy(tag[i], "PSpaceBins");
	  addr[i] = &All.PSpaceBins;
	  id[i++] = INT;

	  strcpy(tag[i], "PSpaceXMin");
	  addr[i] = &All.PSpaceXMin;
	  id[i++] = REAL;

	  strcpy(tag[i], "PSpaceXMax");
	  addr[i] = &All.PSpaceXMax;
	  id[i++] = REAL;

	  strcpy(tag[i], "PSpaceYMin");
	  addr[i] = &All.PSpaceYMin;
	  id[i++] = REAL;

	  strcpy(tag[i], "PSpaceYMax");
	  addr[i] = &All.PSpaceYMax;
	  id[i++] = REAL;
	}

      if(All.Usage == 2)
	{
	  strcpy(tag[i], "RadBins");
	  addr[i] = &All.RadBins;
	  id[i++] = INT;
	}

      if(All.Usage == 3)
	{
	  strcpy(tag[i], "MBEBins");
	  addr[i] = &All.MBEBins;
	  id[i++] = INT;
	}

      if(All.Usage == 9)
	{
	  strcpy(tag[i], "CutCosm");
	  addr[i] = &All.CutCosm;
	  id[i++] = INT;

	  strcpy(tag[i], "CutDM");
	  addr[i] = &All.CutDM;
	  id[i++] = INT;
	}

      max_tags = i;

      if((fd = fopen(All.ParamFile, "r")))
	{
	  while(!feof(fd))
	    {
	      *buf = 0;

	      fgets(buf, MAX_STRING_LEN, fd);

	      nargs = sscanf(buf, "%s%s", buf1, buf2);

	      if(nargs == 1)
		terminate("Missing entry in parameter file!");

	      if(nargs < 1)
		continue;

	      for(i = 0; i < max_tags; i++)
		if(strcmp(buf1, tag[i]) == 0)
		  {
		    switch (id[i])
		      {
		      case INT:
			*((int *) addr[i]) = atoi(buf2);
			break;
		      case REAL:
			*((double *) addr[i]) = atof(buf2);
			break;
		      case STRING:
			strcpy((char *) addr[i], buf2);
			break;
		      }

		    tag[i][0] = 0;
		  }
	    }

	  fclose(fd);
	}
      else
	terminate("Parameter file not found!");

      for(i = 0; i < max_tags; i++)
	if(*tag[i])
	  {
	    printf("Missing entry for '%s' in parameter file\n", tag[i]);
	    terminate("Missing entry in parameter file!");
	  }
    }

  MPI_Bcast(&All, sizeof(struct GlobData), MPI_BYTE, 0, MPI_COMM_WORLD);
}
