#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# don't modify below this line
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
include ../Makefile.common
DISTR_DIR =..
SRC_DIR = src/C
IFLAGS2=$(IFLAGS) -I $(DISTR_DIR)/include
LFLAGS2=$(LFLAGS) -L $(DISTR_DIR)/lib

all: task_control_server

clean:
	cd $(DISTR_DIRC) && rm *.o

task_control_server: 
	$(CC) $(OFLAGS) $(SRC_DIR)/task_control_server.c -o $(OUTDIR)/task_control_server $(IFLAGS2) $(LLFLAGS) $(LFLAGS2) 
