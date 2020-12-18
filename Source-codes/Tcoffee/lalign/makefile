#
# several parts of the fasta package assume that characters are signed
# if this is not the default for your compiler (e.g. MIPS), be certain to
# use the compiler option that invokes signed char's.
#
#sun (bsd)
CC= cc
# for MIPS: CFLAGS= -O -DUNIX -DBIGMEM -DSFCHAR="':'" -signed 
CFLAGS= -O -DUNIX -DBIGMEM -DSFCHAR="':'"
LFLAGS=  -o
#NRAND= nrand
#IBM RS/6000
NRAND= nrand48
RANFLG= -DRAND32
#HZ=60 for sun, mips, 100 for rs/6000, SGI
HZ=60


PROGS= lalign2list
.c.o:
	$(CC) $(CFLAGS) -c $<
all : $(PROGS)
clean : 
	rm *.o

lalign2list : lalign.o pam.o lsim.o lgetaa.o ndispn.o time.o
	$(CC) lalign.o pam.o lsim.o lgetaa.o ndispn.o time.o $(LFLAGS) lalign2list
	mv lalign2list $(USER_BIN)/lalign2list
lalign.o : lalign.c upam.gbl
pam.o : uascii.gbl upam.gbl
lsim.o : lsim.c zzgmata.gbl upam.gbl
	$(CC)  $(CFLAGS) -c lsim.c
lgetaa.o : nxgetaa.c upam.gbl uascii.gbl  altlib.h
	$(CC) $(CFLAGS) -DNOLIB -c nxgetaa.c
	mv nxgetaa.o lgetaa.o
ndispn.o : upam.gbl
time.o : time.c
	$(CC) $(CFLAGS) -DHZ=$(HZ) -c time.c

