//
//  Programmer:  A. Delcher
//
//        File:  ~delcher/Glimmer2/get-putative.cc
//
//     Version:  1.00  31 Jul 98
//
//     Version:  1.01   9 Aug 2000
//               Change  gets  to  fgets  to eliminate compiler warning
//
//    Copyright (c) 1997 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program extracts the list of putative genes from Glimmer
//  output.  Input is from  stdin  and output goes to  stdout.
//


#include  "delcher.h"


const long int  MAX_LINE_LEN = 1000;


int  main  ()

  {
   char  S [MAX_LINE_LEN];
   long int  Ct = 0L;

   while  (fgets (S, MAX_LINE_LEN, stdin) != NULL
             && strstr (S, "Putative Genes:") == NULL)
     ;
   while  (fgets (S, MAX_LINE_LEN, stdin) != NULL)
     {
      printf ("%s", S);
      Ct ++;
     }

   fprintf (stderr, "%ld putative genes extracted\n", Ct);

   return  0;
  }