
#include "allvars.h"

int ThisTask;
int NTask;
int PTask;

int NumGas;
int NumSinks;
int NumPart[NUM_PART_TYPES];
int TotNumPart;

int *NumSinksIter;

int SnapListStart[NUM_PART_TYPES];
int SnapListEnd[NUM_PART_TYPES];
int PartListStart[NUM_PART_TYPES];
int PartListEnd[NUM_PART_TYPES];

int BlockPresent[NUM_BLOCKS];
int BlockRead[NUM_BLOCKS];
int BlockFlag[NUM_BLOCKS];

int BlockPlot[PLOT_NUM_BLOCKS];

int *SendOffset;
int *SendCount;
int *RecvCount;
int *RecvOffset;

double WallClockTime;

char *PD[NUM_BLOCKS];

struct GlobData All;
struct PartValues PV;
struct Sink_struct *Sink;
