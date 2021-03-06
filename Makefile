
EXEC = sator
BUILD_DIR = build
SRC_DIR = src

CC = mpicc

OPTIMIZE = -O2 -m64

#OPT += -DOLDCHEM

CFLAGS = $(OPTIMIZE) $(OPT)

SYSTYPE := "$(SYSTYPE)"

ifeq ($(SYSTYPE), "odyssey-opteron")
CC = mpiicc
endif

OBJS = allvars.o collapse.o cut_region.o halo.o init.o image.o main.o mbe.o mpi_utils.o mymalloc.o params.o plot.o pspace.o radial.o remove_part.o snap.o utils.o

INCL = allvars.h proto.h

OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS))
INCL := $(addprefix $(SRC_DIR)/,$(INCL))

LIBS = -lm

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
