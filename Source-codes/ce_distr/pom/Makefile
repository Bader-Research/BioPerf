# $Revision: 1.3 $ $Date: 1997/12/06 00:01:56 $
MKSRCS	= ipdb.C \
	  ipdb_df.C \
	  pom.C \
	  miscutil.C \
	  bnastats.C \
	  derive.C \
	  calc.C \
	  jdate.C \
	  pdbutil.C \
	  linkedid.C

HEADER	= pom.h \
	  miscutil.h \
	  config.h \
	  ipdb.h \
	  ipdb_df.h \
	  pdbutil.h \
	  jdate.h \
	  water.h \
	  bnastats.h \
	  derive.h \
	  calc.h \
	  linkedid.h

CC = g++

MKOBJS	= $(MKSRCS:.C=.o)

mkDB: $(HEADER)
	@if [ -f pom.o ]; then rm -f pom.o; fi
	@($(MAKE) $(MKOBJS) "CFLAGS = $(CFLAGS) -DREAD_WRITE -DSGI")
	$(CC) $(DEBUGFLAGS) $(CFLAGS) $(INC) -lm -o $@ $(MKOBJS)
	@rm -f pom.o

clean:
	rm -f *.o *.s *.o\~ *.ii core mkDB

include Deps_ipdb

.C.ln:
	$(LINT) $<

.C.o:
	$(CC) $(DEBUGFLAGS) $(CFLAGS) $(INC) -c $<

.C,v.C:
	$(CO) $<

.h,v.h:
	$(CO) $<

depend: 
	makedepend -f Deps_ipdb  -I  /usr/local/apps/gcc/lib/g++-include $(MKSRCS) > Deps_ipdb
