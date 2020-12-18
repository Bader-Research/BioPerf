//  A. L. Delcher
//
//     File:  ~adelcher/Glimmer/gene.h
//  Version:  1.02  25 Feb 98
//            Remove unused variable  Found
//  Version:  1.03   8 Feb 99
//            Make easier to modify start/stop codons
//
//    Copyright (c) 1997-99 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  Common routines for DNA sequences.



#ifndef  __GENE_H_INCLUDED
#define  __GENE_H_INCLUDED


const unsigned  ATG_MASK = 0x184;
const unsigned  CAA_MASK = 0x211;
const unsigned  CAC_MASK = 0x212;
const unsigned  CAG_MASK = 0x214;
const unsigned  CAT_MASK = 0x218;
const unsigned  CAY_MASK = 0x21a;
const unsigned  CTA_MASK = 0x281;
const unsigned  CTG_MASK = 0x284;
const unsigned  GTG_MASK = 0x484;
const unsigned  RTG_MASK = 0x584;
const unsigned  TAA_MASK = 0x811;
const unsigned  TAG_MASK = 0x814;
const unsigned  TAR_MASK = 0x815;
const unsigned  TCA_MASK = 0x821;
const unsigned  TGA_MASK = 0x841;
const unsigned  TRA_MASK = 0x851;
const unsigned  TTA_MASK = 0x881;
const unsigned  TTG_MASK = 0x884;
const unsigned  TYA_MASK = 0x8a1;
const unsigned  YTA_MASK = 0xa81;
const unsigned  SHIFT_MASK = 0xFF;

const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const int  MAX_LINE = 300;


unsigned  Ch_Mask
    (char);
int  Codon_To_Subscript
    (char * s);
char  Complement
    (char);
char  Filter
    (char);
void  Find_Stop_Codons
    (char [], int, int []);
int  Is_Forward_Start
    (unsigned);
int  Is_Forward_Stop
    (unsigned);
int  Is_Reverse_Start
    (unsigned);
int  Is_Reverse_Stop
    (unsigned);
int  Is_Start
    (char *);
int  Is_Stop
    (char *);
int  Nucleotide_To_Subscript
    (char ch);
int  Read_String
    (FILE *, char * &, long int &, char [], int);
int  Rev_Codon_To_Subscript
    (char * s);
void  Reverse_Complement
    (char [], long int);
int  Rev_Nucleotide_To_Subscript
    (char ch);
char *  Subscript_To_Codon
    (int sub);



#endif
