SOURCES =ce.C ../ce/ce_align.C   ../pom/pom.C ../ce/cmp_util.C   ../pom/miscutil.C 
LIBS =  -lm 
#-lstdc++
PROGRAM = ce
MAKEDEP = makedepend -f Deps_ce
OBJECTS= $(SOURCES:.C=.o)
LINTFILES=$(SOURCES:.C=.ln)
INCLUDE_DIR=

CC = g++   ${INCLUDE_DIR} -DREAD_WRITE 
CFLAGS = -DFUNCPROTO  


$(PROGRAM): $(OBJECTS) 
		$(CC)  -o $@ $(OBJECTS) $(LIBS)

$(LINTFILES):	
		$(LINT) $(LINTFILES) 

clean:
	rm -f *.o core ce

clall:
	rm -f *.o core a.out $(PROGRAM)

depend:
	$(MAKEDEP) ${INCLUDE_DIR} $(SOURCES) > Deps_ce

include Deps_ce

ce.o: ce.C 
	${CC} -c ce.C

../ce/ce_draw.o: ../ce/cmp_draw.C
	${CC} -o ../ce/cmp_draw.o  -c  ../ce/cmp_draw.C

../pom/pom.o:  ../pom/pom.C 
	${CC} -o ../pom/pom.o  -c  ../pom/pom.C

../pom/miscutil.o:  ../pom/miscutil.C 
	${CC} -o ../pom/miscutil.o  -c  ../pom/miscutil.C

../ce/cmp_util.o:  ../ce/cmp_util.C 
	${CC} -o ../ce/cmp_util.o  -c  ../ce/cmp_util.C

../ce/ce_align.o:  ../ce/ce_align.C 
	${CC} -o ../ce/ce_align.o  -c  ../ce/ce_align.C

../ce/ce_calc.o:  ../ce/ce_calc.C 
	${CC} -o ../ce/ce_calc.o  -c  ../ce/ce_calc.C

../ce/ce_all.o:  ../ce/ce_all.C 
	${CC} -o ../ce/ce_all.o  -c  ../ce/ce_all.C

../ce/ce_fssp.o:  ../ce/ce_fssp.C 
	${CC} -o ../ce/ce_fssp.o  -c  ../ce/ce_fssp.C

../ce/ce_thr.o:  ../ce/ce_thr.C 
	${CC} -o ../ce/ce_thr.o  -c  ../ce/ce_thr.C

../ce/ce_sp.o:  ../ce/ce_sp.C 
	${CC} -o ../ce/ce_sp.o  -c  ../ce/ce_sp.C

../ce/ce_dom.o:  ../ce/ce_dom.C 
	${CC} -o ../ce/ce_dom.o  -c  ../ce/ce_dom.C

# DO NOT DELETE THIS LINE -- make depend depends on it.
