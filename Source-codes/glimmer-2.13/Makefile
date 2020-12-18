# -*- makefile -*-
# just make headings for each thing and put the headings
# on the line for ALL
PROJECTDIR = .
BINDIR = $(PROJECTDIR)

CC = gcc
CPPC = g++
CFLAGS = -O3 -g -Wall
LDFLAGS = -lm

DEPEND_FILES = *.cc *.c *.h
CLEANABLE_FILES = *.o *~

ALL = adjust anomaly build-icm check codon-usage compare-lists \
       extract generate get-len get-putative glimmer2 long-orfs

.SUFFIXES: .cc .c

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o 

.cc.o:
	$(CPPC) $(CFLAGS) -c $*.cc -o $*.o 

all:    $(ALL)

adjust:	adjust.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ adjust.o delcher.o $(LDFLAGS)

anomaly:	anomaly.o delcher.o gene.o
	$(CPPC) -o $(BINDIR)/$@ anomaly.o delcher.o gene.o $(LDFLAGS)

build-icm:	build-icm.o misc.o
	$(CC) -o $(BINDIR)/$@ build-icm.o misc.o $(LDFLAGS)

check:	check.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ check.o delcher.o $(LDFLAGS)

codon-usage:	codon-usage.o delcher.o gene.o
	$(CPPC) -o $(BINDIR)/$@ codon-usage.o delcher.o gene.o $(LDFLAGS)

compare-lists:	compare-lists.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ compare-lists.o delcher.o $(LDFLAGS)

extract:	extract.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ extract.o delcher.o $(LDFLAGS)

generate:	generate.o delcher.o gene.o
	$(CPPC) -o $(BINDIR)/$@ generate.o delcher.o gene.o $(LDFLAGS)

get-len:	get-len.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ get-len.o delcher.o $(LDFLAGS)

get-putative:	get-putative.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ get-putative.o delcher.o $(LDFLAGS)

glimmer2:	glimmer2.o delcher.o gene.o rnabin.o
	$(CPPC) -o $(BINDIR)/$@ glimmer2.o delcher.o gene.o rnabin.o $(LDFLAGS)

long-orfs:	long-orfs.o delcher.o gene.o
	$(CPPC) -o $(BINDIR)/$@ long-orfs.o delcher.o gene.o $(LDFLAGS)

rnabin:	rnabin.o delcher.o
	$(CPPC) -o $(BINDIR)/$@ rnabin.o delcher.o $(LDFLAGS)

depend:
	makedepend $(DEPEND_FILES)

clean:
	/bin/rm -f $(CLEANABLE_FILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.

