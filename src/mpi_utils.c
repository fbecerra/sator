
#include "allvars.h"
#include "proto.h"

static char *SaveData2;

void mpi_exchange_buffers(void *send_buf, int *send_count, int *send_offset, void *recv_buf, int *recv_count, int *recv_offset, int item_size, int commtag, int include_self);
int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
int intpointer_compare(const void *a, const void *b);
void *sort_based_on_field(void *data, int field_offset, int n_items, int item_size);


void mpi_exchange_buffers(void *send_buf, int *send_count, int *send_offset, void *recv_buf, int *recv_count, int *recv_offset, int item_size, int commtag, int include_self)
{
  int ngrp;

  for(ngrp = include_self ? 0 : 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
	    {
	      MPI_Sendrecv((char *) send_buf + send_offset[recvTask] * item_size, send_count[recvTask] * item_size, MPI_BYTE, recvTask, commtag, (char *) recv_buf + recv_offset[recvTask] * item_size, recv_count[recvTask] * item_size, MPI_BYTE, recvTask, commtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }
}


int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical)
{
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  int j, nimport = 0;

  recv_offset[0] = 0;
  send_offset[0] = 0;

  for(j = 0; j < NTask; j++)
    {
      nimport += recv_count[j];

      if(j > 0)
	{
	  send_offset[j] = send_offset[j - 1] + (send_identical ? 0 : send_count[j - 1]);
	  recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
	}
    }
  return nimport;
}


int intpointer_compare(const void *a, const void *b)
{
  if((**(int **) a) < (**(int **) b))
    return -1;

  if((**(int **) a) > (**(int **) b))
    return +1;

  return 0;
}


void *sort_based_on_field(void *data, int field_offset, int n_items, int item_size)
{
  int i;
  char *data2;
  int **perm;

  data2 = mymalloc_movable(&SaveData2, "data2", n_items * item_size);

  SaveData2 = data2;

  perm = mymalloc("perm", n_items * sizeof(*perm));

  for(i = 0; i < n_items; ++i)
    perm[i] = (int *) ((char *) data + i * item_size + field_offset);

  mysort(perm, n_items, sizeof(*perm), intpointer_compare);

  for(i = 0; i < n_items; ++i)
    {
      size_t orig_pos = ((char *) perm[i] - ((char *) data + field_offset)) / item_size;
      myassert(((char *) perm[i] - ((char *) data + field_offset)) % item_size == 0);
      memcpy(data2 + item_size * i, (char *) data + item_size * orig_pos, item_size);
    }

  myfree(perm);

  return (void *) data2;
}


void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size, int commtag)
{
  int i;

  for(i = 0; i < NTask; i++)
    SendCount[i] = 0;

  for(i = 0; i < *n_items; i++)
    {
      int task = *(int *) ((char *) data + i * item_size + task_offset);

      myassert(task >= 0 && task < NTask);

      SendCount[task]++;
    }

  void *data2 = sort_based_on_field(data, task_offset, *n_items, item_size);

  int nimport = mpi_calculate_offsets(SendCount, SendOffset, RecvCount, RecvOffset, 0);

  if(*max_n < nimport)
    {
      data = myrealloc_movable(data, nimport * item_size);

      *max_n = nimport;
    }

  data2 = SaveData2;

  mpi_exchange_buffers(data2, SendCount, SendOffset, data, RecvCount, RecvOffset, item_size, commtag, 1);

  myfree_movable(data2);

  *n_items = nimport;
}
