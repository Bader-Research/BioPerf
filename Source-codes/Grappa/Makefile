VERS_MAJ = 1
VERS_MIN = 6

#for MPI, change it to tmcc, and added the line of -l mpi
#then mprun -np grappa -f ...

EXEC  = grappa
SYM   = invdist distmat

OS    = $(shell uname -s | sed -e "s/\-/_/g" -e "s/\.//g")
ZIP   = gzip
MKDIR = mkdir
TAR   = tar
CP    = cp
RM    = rm
LS    = ls

# Select the appropriate choice below for using Concorde
# (Lin-Kernighan TSP solvers)
#USE_CONCORDE = YES
USE_CONCORDE = NO

ifeq ($(USE_CONCORDE),YES)
# Set CONCORDE_ROOT to the location of the compiler concorde libraries
  CONCORDE_ROOT   = $(HOME)/concorde
  CONCORDE_INC    = $(CONCORDE_ROOT)/INCLUDE
  CONCORDE_CFLAGS = -I$(CONCORDE_INC) -DCONCORDE
  CONCORDE_LIB    = $(CONCORDE_ROOT)/concorde.a
else
  CONCORDE_CFLAGS =
  CONCORDE_LIB    =
endif

# Use GNU Multiple Precision Arithmetic Library
#USE_GMP = YES
USE_GMP = NO

ifeq ($(USE_GMP),YES)
# Set GMP_ROOT to the location of GNU Multiple Precision 
  GMP_ROOT   = /usr
  GMP_INC    = $(GMP_ROOT)/include
  GMP_CFLAGS = -I$(GMP_INC) -DGMP
  GMP_LIB    = -L$(GMP_ROOT)/lib -lgmp
else
  GMP_CFLAGS =
  GMP_LIB    =
endif

# GNU C
CC = gcc

# Microsoft C
# CC = cl

# Intel C
# CC = icl

# Portland Group C
# CC = pgcc

# MPI C
#CC = tmcc -DMPBPA

# Default CFLAGS
CFLAGS := -D${OS} $(CONCORDE_CFLAGS) $(GMP_CFLAGS) -DINIT

# TESTING option for CFLAGS
# CFLAGS := -DTESTING $(CFLAGS)

ifeq ($(strip $(CC)), gcc)
CFLAGS := -O3 -Wall -DGCC $(CFLAGS)
endif
ifeq ($(strip $(CC)), pgcc)
CFLAGS := -fast $(CFLAGS)
endif
ifeq ($(strip $(CC)), cl)
CFLAGS := /Ot /Ox $(CFLAGS)
endif
ifeq ($(strip $(CC)), icl)
CFLAGS := -O2 $(CFLAGS)
endif

ifeq ($(strip $(CC)), gcc)
ifeq ($(OS), Linux)
CFLAGS := -mpentiumpro $(CFLAGS)
endif
ifeq ($(OS), SunOS)
#CFLAGS := -mv9 $(CFLAGS)
endif
endif

OBJ_SUFFIX := o
EXE_SUFFIX := 

ifeq ($(OS), CYGWIN_NT_40)
OBJ_SUFFIX := obj
EXE_SUFFIX := exe
endif

SOURCES = convert labeltree main binencode gen_tree inittree circ_order \
	  lk_main bbtsp condense const_tree greedy_median invdist uf \
	  neighborj condense3 getopt3 randomBSD cheaptsp specialtsp \
	  specialinit read_input circular_ordering correction \
	  inversion_median inversion_median_alberto priority_stack hashtable med_util simpleio \
	  vertex_factory growTree bitvector lists \
	  sorting_reversal_median mythread_rwlock all_sorting_reversals
HEADERS = structs.h $(filter-out main.h,$(addsuffix .h,$(SOURCES)))

OBJS = $(addsuffix .${OBJ_SUFFIX},$(SOURCES))

TSTAMP   := $(shell date +"%y%m%d-%H%M")
DISTBASE  = MPBPA
DISTNAME  = $(DISTBASE)-$(TSTAMP)

RELBASE   = GRAPPA
RELNAME   = $(RELBASE)-$(VERS_MAJ).$(VERS_MIN)

LIBS 		:= -lm -lc -lpthread 
# added this part to the end of above "-l mpi"
ifeq ($(OS), SunOS)
LIBS		:= -liberty -lsocket -lposix4 -lnsl ${LIBS}
endif

.PHONY: dist release release_code clean

all: ${EXEC} ${SYM}

${EXEC} :  ${OBJS}
	${CC}  ${CFLAGS}  -o $@  ${OBJS} ${CONCORDE_LIB} ${GMP_LIB} ${LIBS}

%.${OBJ_SUFFIX}:%.c
	$(CC) $(CFLAGS) -c $< -o $@

${SYM}: ${EXEC}
	ln -s $< $@

dist:
	$(RM) -fr $(DISTNAME) $(DISTNAME).tar $(DISTNAME).tar.gz
	$(MKDIR) $(DISTNAME)
	-$(CP) Makefile COPYRIGHT README CHANGES $(addsuffix .c,$(SOURCES)) $(HEADERS) \
	   $(DISTNAME)/.
	$(TAR) cvf $(DISTNAME).tar $(DISTNAME)
	$(ZIP) $(DISTNAME).tar
	$(RM) -fr $(DISTNAME)
	$(LS) -l $(DISTNAME).tar.gz

release:
	$(RM) -fr $(RELNAME) $(RELNAME).tar $(RELNAME).tar.gz
	$(MKDIR) $(RELNAME)
	-$(CP) Makefile COPYRIGHT README CHANGES $(RELNAME)/.
	$(MAKE) $(addprefix $(RELNAME)/,$(addsuffix .c,$(SOURCES)) $(HEADERS))
	$(TAR) cvf $(RELNAME).tar $(RELNAME)
	$(ZIP) $(RELNAME).tar
	$(RM) -fr $(RELNAME)
	$(LS) -l $(RELNAME).tar.gz

$(RELNAME)/%.c : %.c
	-$(CP) $< $(RELNAME)/.

$(RELNAME)/%.h : %.h
	-$(CP) $< $(RELNAME)/.

clean:
	-rm -f *.${OBJ_SUFFIX} *~ ${EXEC} ${EXEC}.${EXE_SUFFIX} ${SYM} \
	       core gmon.out

# DO NOT DELETE

$(OBJS): Makefile
bbtsp.o: structs.h bbtsp.c bbtsp.h
binencode.o: structs.h binencode.c binencode.h
circ_order.o: structs.h invdist.h
circular_ordering.o: structs.h 
condense.o: condense.h condense.c structs.h
const_tree.o: structs.h const_tree.c const_tree.h gen_tree.h
convert.o: structs.h convert.c convert.h binencode.h
gen_tree.o: structs.h gen_tree.c gen_tree.h 
greedy_median.o: structs.h greedy_median.c greedy_median.h
inittree.o: structs.h inittree.c inittree.h convert.h bbtsp.h
inittree.o: binencode.h specialtsp.h cheaptsp.h greedy_median.h
invdist.o: structs.h invdist.h invdist.c uf.h
labeltree.o: structs.h labeltree.c labeltree.h
labeltree.o: bbtsp.h convert.h lk_main.h binencode.h invdist.h
main.o: main.c labeltree.h bbtsp.h inittree.h gen_tree.h binencode.h
main.o: structs.h condense.h const_tree.h uf.h invdist.h neighborj.h
main.o: circ_order.h read_input.h circular_ordering.h
uf.o: uf.c uf.h
neighborj.o: neighborj.h neighborj.c structs.h invdist.h const_tree.h
neighborj.o: labeltree.h inittree.h condense.h
condense3.o: condense3.h condense3.c structs.h
getopt3.o: getopt3.c getopt3.h
randomBSD.o: randomBSD.c randomBSD.h
cheaptsp.o: cheaptsp.h cheaptsp.c structs.h
specialtsp.o: specialtsp.h specialtsp.c structs.h
specialinit.o: specialinit.h specialinit.c structs.h
read_input.o: structs.h read_input.h
correction.o: correction.h

