
#include "allvars.h"
#include "proto.h"

static size_t TotBytes;                    /**< The total dimension (in bytes) of dynamic memory available to the current task. */
static size_t AllocatedBytes;
static size_t HighMarkBytes;
static size_t FreeBytes;

static void *Base;                         /**< Base pointer (initial memory address) of the stack. */

static unsigned long Nblocks;              /**< The current number of allocated memory blocks. */

static void **Table;                       /**< Table containing the initial addresses of the allocated memory blocks.*/
static size_t *BlockSize;                  /**< Array containing the size (in bytes) of all the allocated memory blocks. */
static char *MovableFlag;                  /**< Identifies whether a block is movable. */
static void ***BasePointers;               /**< Base pointers containing the initial addresses of movable memory blocks */
static size_t GlobHighMarkBytes = 0;       /**< The maximum number of bytes allocated by all tasks. */
static char *VarName;                      /**< The name of the variable with which the block has been allocated. */
static char *FunctionName;                 /**< The function name that has allocated the memory block. */
static char *FileName;                     /**< The file name where the function that has allocated the block is called. */
static int *LineNumber;                    /**< The line number in FileName where the function that allocated the block has been called. */

void mymalloc_init(void)
{
  size_t n;
  char msg[MAX_STRING_LEN];

  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));
  MovableFlag = (char *) malloc(MAXBLOCKS * sizeof(char));
  BasePointers = (void ***) malloc(MAXBLOCKS * sizeof(void **));
  VarName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber = (int *) malloc(MAXBLOCKS * sizeof(int));

  memset(VarName, 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName, 0, MAXBLOCKS * MAXCHARS);

  n = All.MaxMemSize * ((size_t) 1024 * 1024);

  if(n & 63)
    terminate("Want 64 byte aligned address!");

  if(!(Base = malloc(n)))
    {
      sprintf(msg, "Failed to allocate memory for `Base' (%d Mbytes)\n", All.MaxMemSize);

      terminate(msg);
    }

  TotBytes = FreeBytes = n;

  AllocatedBytes = 0;
  Nblocks = 0;
  HighMarkBytes = 0;
}


void dump_memory_table(void)
{
  int i, cc = 0;
  size_t totBlocksize = 0;
  char *buf;

  buf = (char *) mymalloc("MemoryTable", 200 * (Nblocks + 10));

  cc += sprintf(buf + cc, "-------------------------- Allocated Memory Blocks----------------------------------------\n");
  cc += sprintf(buf + cc, "Task   Nr F             Variable      MBytes   Cumulative         Function/File/Linenumber\n");
  cc += sprintf(buf + cc, "------------------------------------------------------------------------------------------\n");

  for(i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      cc += sprintf(buf + cc, "%4d %5d %d %19s  %10.4f   %10.4f  %s()/%s/%d\n",
	     ThisTask, i, MovableFlag[i], VarName + i * MAXCHARS, BlockSize[i] / (1024.0 * 1024.0),
	     totBlocksize / (1024.0 * 1024.0), FunctionName + i * MAXCHARS, FileName + i * MAXCHARS, LineNumber[i]);
    }

  cc += sprintf(buf + cc, "------------------------------------------------------------------------------------------\n");

  printf("%s", buf);

  fflush(stdout);

  myfree(buf);
}


void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line)
{
  char msg[MAX_STRING_LEN];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      sprintf(msg, "Task %d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask, func, file, line, MAXBLOCKS);

      terminate(msg);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB)\n", ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));

      terminate(msg);
    }

  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 0;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}

void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line)
{
  char msg[MAX_STRING_LEN];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      sprintf(msg, "Task %d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask, func, file, line, MAXBLOCKS);

      terminate(msg);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB)\n", ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));

      terminate(msg);
    }

  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 1;
  BasePointers[Nblocks] = ptr;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}

void myfree_fullinfo(void *p, const char *func, const char *file, int line)
{
  char msg[MAX_STRING_LEN];

  if(Nblocks == 0)
    terminate("No allocated blocks that could be freed!");

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Wrong call of myfree() at %s()/%s/line %d: not the last allocated block!\n", ThisTask, func, file, line);

      terminate(msg);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
  FreeBytes += BlockSize[Nblocks];
}

void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line)
{
  int i, nr;
  char msg[MAX_STRING_LEN];

  if(Nblocks == 0)
    terminate("No allocated blocks that could be freed!");

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Wrong call of myfree_movable() from %s()/%s/line %d - this block has not been allocated!\n", ThisTask, func, file, line);

      terminate(msg);
    }

  if(nr < Nblocks - 1)
    {
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();

	    sprintf(msg, "Task %d: Wrong call of myfree_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n", ThisTask, func, file, line, nr);

	    terminate(msg);
	  }
    }

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  size_t offset = -BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;
      *BasePointers[i] = *BasePointers[i] + offset;
    }

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i - 1] = Table[i];
      BasePointers[i - 1] = BasePointers[i];
      BlockSize[i - 1] = BlockSize[i];
      MovableFlag[i - 1] = MovableFlag[i];

      strncpy(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS - 1);
      LineNumber[i - 1] = LineNumber[i];
    }

  Nblocks -= 1;
}

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  char msg[MAX_STRING_LEN];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    terminate("No allocated blocks that could be reallocated!");

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Wrong call of myrealloc() at %s()/%s/line %d - not the last allocated block!\n", ThisTask, func, file, line);

      terminate(msg);
    }

  AllocatedBytes -= BlockSize[Nblocks - 1];
  FreeBytes += BlockSize[Nblocks - 1];

  if(n > FreeBytes)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Not enough memory in myremalloc(n=%g MB) at %s()/%s/line %d. previous=%g FreeBytes=%g MB\n", ThisTask, n / (1024.0 * 1024.0), func, file, line, BlockSize[Nblocks - 1] / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));

      terminate(msg);
    }

  Table[Nblocks - 1] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  AllocatedBytes += n;
  BlockSize[Nblocks - 1] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}

void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  int i, nr;
  char msg[MAX_STRING_LEN];

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    terminate("No allocated blocks that could be reallocated!");

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: Wrong call of myrealloc_movable() from %s()/%s/line %d - this block has not been allocated!\n", ThisTask, func, file, line);

      terminate(msg);
    }

  if(nr < Nblocks - 1)
    {
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();

	    sprintf(msg, "Task %d: Wrong call of myrealloc_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n", ThisTask, func, file, line, nr);

	    terminate(msg);
	  }
    }

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n > FreeBytes)
    {
      dump_memory_table();

      sprintf(msg, "Task %d: at %s()/%s/line %d: Not enough memory in myremalloc_movable(n=%g MB). previous=%g FreeBytes=%g MB\n", ThisTask, func, file, line, n / (1024.0 * 1024.0), BlockSize[nr] / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));

      terminate(msg);
    }

  size_t offset = n - BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;

      *BasePointers[i] = *BasePointers[i] + offset;
    }

  FreeBytes -= n;
  AllocatedBytes += n;
  BlockSize[nr] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[nr];
}
