
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXBLOCKS 5000
#define MAXCHARS 20
#define MAX_PARAMS 1000
#define MAX_STRING_LEN 100
#define NUM_PART_TYPES 6

#define GAMMA_ADB (5. / 3.)
#define HYDROGEN_MASSFRAC 0.76
#define HE_ABUND ((1. - HYDROGEN_MASSFRAC) / 4. / HYDROGEN_MASSFRAC)
#define DE_ABUND 2.6e-5
#define MU_NEUTRAL (4. / (1. + 3. * HYDROGEN_MASSFRAC))
#define MU_IONIZED (4. / (8. - 5. * (1. - HYDROGEN_MASSFRAC)))
#define PROTONMASS 1.6726e-24
#define BOLTZMANN 1.3806e-16
#define ELECTRON_VOLT 1.60219e-12
#define GRAVITY 6.672e-8
#define SOLAR_MASS 1.989e33
#define HUBBLE 3.2407789e-18
#define SEC_PER_YEAR 3.155e7
#define ASTRONOMICAL_UNIT 1.49598e13
#define SOLAR_RADIUS 6.955e10

#define TGCHEM_NUM_TEMP 1000
#define TGCHEM_NUM_H2 41
#define TGCHEM_NUM_H2OP 101

#define TGCHEM_TEMP_MIN 1e1
#define TGCHEM_TEMP_MAX 1e8
#define TGCHEM_LOG_DTEMP (log10(TGCHEM_TEMP_MAX / TGCHEM_TEMP_MIN) / TGCHEM_NUM_TEMP)

#define TGCHEM_TOL_ABHM 1e-20
#define TGCHEM_TOL_ABH2 1e-20
#define TGCHEM_TOL_ABHII 1e-20

#define UNIT_LENGTH 3.085678e21
#define UNIT_MASS 1.989e43
#define UNIT_VELOCITY 1e5
#define UNIT_TIME (UNIT_LENGTH / UNIT_VELOCITY)
#define UNIT_ENERGY (UNIT_MASS * UNIT_LENGTH * UNIT_LENGTH / UNIT_TIME / UNIT_TIME)

#define myassert(cond) {if(!(cond)) {char termbuf[MAX_STRING_LEN]; sprintf(termbuf, "Assertion failure!\n\task = %d, function %s(), file %s, line %d:\n\t%s\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, #cond); printf("%s", termbuf); myflush(stdout); assert(0);}}

#define mymalloc(x, y) mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define myrealloc(x, y) myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define myrealloc_movable(x, y) myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define myfree(x) myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define myfree_movable(x) myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define terminate(...) {char termbuf1[MAX_STRING_LEN], termbuf2[MAX_STRING_LEN]; sprintf(termbuf1, "Code termination on task %d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s\n%s\n", termbuf1, termbuf2); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1);}

enum io_fields
{
  IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_VOL,
  IO_DELAUNAY,
  IO_GRAVACC,
  IO_GRADP,
  IO_CHEM,
  IO_GAMMA,
  IO_RATES,
  IO_ALLOWREF,
  IO_DIVVEL,
  IO_SIDM_PSUM,
  IO_SIDM_NUMNGB,
  IO_SIDM_NUMTOTALSCATTER,
  IO_SIDM_HSML,
  IO_SIDM_DENSITY,
  IO_SIDM_VELDISP,
  NUM_BLOCKS
};

enum plot_fields
{
  PLOT_NH,
  PLOT_TEMP,
  PLOT_GRAVACC,
  PLOT_GRADP,
#ifdef OLDCHEM
  PLOT_H2,
  PLOT_HII,
  PLOT_DII,
  PLOT_HD,
  PLOT_HEII,
  PLOT_HEIII,
#else
  PLOT_HM,
  PLOT_H2,
  PLOT_HII,
#endif
  PLOT_GAMMA,
  PLOT_ESCFRAC,
  PLOT_PDVRATE,
  PLOT_H2RATE,
  PLOT_CIERATE,
  PLOT_CHEMRATE,
  PLOT_ALLOWREF,
  PLOT_DIVVEL,
  PLOT_COOL,
  PLOT_COLLAPSE,
  PLOT_SIDM_DENSITY,
  PLOT_NUM_BLOCKS
};

extern int ThisTask;
extern int NTask;
extern int PTask;

extern int NumGas;
extern int NumSinks;
extern int NumPart[NUM_PART_TYPES];
extern int TotNumPart;

extern int *NumSinksIter;

extern int SnapListStart[NUM_PART_TYPES];
extern int SnapListEnd[NUM_PART_TYPES];
extern int PartListStart[NUM_PART_TYPES];
extern int PartListEnd[NUM_PART_TYPES];

extern int BlockPresent[NUM_BLOCKS];
extern int BlockRead[NUM_BLOCKS];
extern int BlockFlag[NUM_BLOCKS];

extern int BlockPlot[PLOT_NUM_BLOCKS];

extern int *SendOffset;
extern int *SendCount;
extern int *RecvCount;
extern int *RecvOffset;

extern double WallClockTime;

extern char *PD[NUM_BLOCKS];

extern struct PartList
{
  int type[NUM_PART_TYPES];

} *TaskNumPart, *SnapNumPart;

extern struct GlobData
{
  // Command line arguments

  char ParamFile[MAX_STRING_LEN];
  char Base[MAX_STRING_LEN];
  int Usage;
  int SnapNumStart;
  int SnapNumEnd;

  double ImgWidth;
  double ImgSubWidth;
  int ImgNumSubSnaps;

  char CutBase[MAX_STRING_LEN];
  double CutSize;

  double RadMin;

  double MBEMin;

  // Parameter file

  int MaxMemSize;

  int FileSystem;
  int FlagDouble;
  char BlockFile[MAX_STRING_LEN];
  int LengthUnit;
  int ComovingUnits;

  char BlockFilePlot[MAX_STRING_LEN];
  int FlagCenter;
  double CenterX;
  double CenterY;
  double CenterZ;
  int FlagRotate;

  int ImgFlagScreenRatio;
  int ImgXBins;
  int ImgFlagMovie;

  int PSpaceBins;
  double PSpaceXMin;
  double PSpaceXMax;
  double PSpaceYMin;
  double PSpaceYMax;

  int RadBins;

  int MBEBins;
  double ShiftX;
  double ShiftY;
  double ShiftZ;

  int RemoveCosm;

  char BlockFileCut[MAX_STRING_LEN];
  int CutCosm;
  int CutDM;

  // Other

  int Iter;
  int SnapNum;
  int NumSnaps;
  char Path[MAX_STRING_LEN];
  char SnapNumString[MAX_STRING_LEN];
  char SnapFile[MAX_STRING_LEN];
  char SnapBase[MAX_STRING_LEN];

  double Mass[NUM_PART_TYPES];
  double Time;
  double RedShift;
  int GlobNumPart[NUM_PART_TYPES];
  int NumFiles;
  double BoxSize;
  double OmegaM;
  double OmegaLambda;
  double HubbleParam;

  double BoxCenter[3];

  int TotNumSinks;

  char Ending[MAX_STRING_LEN];
  int ImgYBins;
  double ImgScreenRatio;
  double ImgHeight;
  double ImgSize;
  double ImgMin[PLOT_NUM_BLOCKS][3];
  double ImgMax[PLOT_NUM_BLOCKS][3];

  double CutSizeHalf;

  struct PartList *TaskNumPart;
  struct PartList *SnapNumPart;

  double *TempTable;
  double *CoolTable;
  double *DCoolTable;
  double *ChemTable;
  double *DChemTable;
  double *H2RateOut;
  float *H2OpacTemp;
  float *H2OpacColumn;
  float *H2Opac;

} All;

extern struct PartValues
{
  double *X;
  double *Y;
  double *Z;
  double *VX;
  double *VY;
  double *VZ;
  int *ID;
  double *Mass;
  double *U;
  double *Temp;
  double *Rho;
  double *NH;
  double *Vol;
  double *Hsml;
  double *Delaunay;
  double *GravAccX;
  double *GravAccY;
  double *GravAccZ;
  double *GradPX;
  double *GradPY;
  double *GradPZ;
#ifdef OLDCHEM
  double *AbH2;
  double *AbHII;
  double *AbDII;
  double *AbHD;
  double *AbHeII;
  double *AbHeIII;
#else
  double *AbHM;
  double *AbH2;
  double *AbHII;
#endif
  double *Mu;
  double *Gamma;
  double *EscFrac;
  double *PdVRate;
  double *H2Rate;
  double *CIERate;
  double *ChemRate;
  double *PdVCool;
  double *H2Cool;
  double *CIECool;
  double *ChemCool;
  int *AllowRef;
  double *AllowRefPlot;
  double *DivVel;
  double *Cool;
  double *Collapse;
  double *SIDM_PSum;
  double *SIDM_NumNgb;
  int *SIDM_NumTotalScatter;
  double *SIDM_Hsml;
  double *SIDM_Density;
  double *SIDM_VelDisp;

} PV;


extern struct Sink_struct
{
  int Task;
  double X;
  double Y;
  double Z;

} *Sink;
