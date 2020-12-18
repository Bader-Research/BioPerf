#
# makefile for fasta3, fasta3_t.  Use makefile.pvm for pvcompxx.
#
# On the version of linux that I use, MININT is defined and causes a warning
# as dropnfa.c is compiled.  Unfortunately, using the Linux definition of
# MININT causes the code to break. Do not change the MININT definition
# in dropnfa.c  (wrp 3/19/1998)
#

CC= gcc -g -O3
#CC=gcc -Wall -pedantic -ansi -g -O
#CC = gcc -g -DDEBUG
#CC= /opt/parasoft/bin.linux2/insure -g -pg -DDEBUG

# EBI uses the following with pgcc, -O3 does not work:
# CC= pgcc -O2 -pipe -mcpu=pentiumpro -march=pentiumpro -fomit-frame-pointer

# this file works for x86 LINUX

# use options below for superfamily validations
#CFLAGS= -DSHOWSIM -DUNIX -DTIMES -DHZ=100 -DSFCHAR="'|'" -c -DMAX_WORKERS=2 -DTHR_EXIT=pthread_exit -DPROGRESS -DSUPERFAMNUM -DUSE_MMAP -D_REENTRANT -DBIG_LIB64 -D_LARGE_FILE_SOURCE -DUSE_FSEEKO -D_FILE_OFFSET_BITS=64 -DHAS_INTTYPES -DSAMP_STATS

# standard options
CFLAGS= -DSHOWSIM -DUNIX -DTIMES -DHZ=100 -DSFCHAR="':'" -c -DMAX_WORKERS=2 -DTHR_EXIT=pthread_exit -DPROGRESS -DFASTA_HOST='"wrpx00.bioch.virginia.edu/fasta/cgi"' -DUSE_MMAP -D_REENTRANT -DHAS_INTTYPES -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_FSEEKO -DSAMP_STATS -DPGM_DOC
# -I/usr/local/include/mysql -DMYSQL_DB 
#
#(for mySQL databases)  (also requires change to Makefile34.common)

LIB_M = -lm
#LIB_M = -L/usr/local/lib/mysql -lmysqlclient -lm
# for mySQL databases

HFLAGS= -o
NFLAGS= -o

# for Linux
THR_SUBS = pthr_subs2
THR_LIBS = -lpthread
THR_CC =

XDIR = /seqprg/bin

DROPNFA_O = drop_nfa.o
DROPGSW_O = dropgsw.o
DROPRSS_O = dropnsw.o
DROPTFA_O = drop_tfa.o

# renamed (fasta33)  programs
include Makefile34m.common
# conventional (fasta3) names
# include Makefile.common

