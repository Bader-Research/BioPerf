FLAGS = -lm -o
CFLAGS = -O2

CC = cc
#CC = gcc #for Suns, OK for other machines too

OBJECT = prd.o prdfali.o prdmul.o prdpred.o prdmem.o prdread.o prdout.o prdmisc.o prdseq.o prdglali.o

predator : $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) $(FLAGS) predator

$(OBJECT) : prd.h prdprot.h

clean:
	rm -f $(OBJECT)

