//  Programmer:  A. Delcher
//
//        File:  Glimmer/anomaly.cc
//
//     Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Arthur Delcher, Steven Salzberg,
//    Simon Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program checks for in-frame stop codons
//


#include  "delcher.h"

const long int  MAX_GENE_LEN = 100000;


int  main
    ()

  {
   char  Tag [100], Gene [MAX_GENE_LEN];
   long int  i, Len;

   while  (cin >> Tag >> Gene)
     {
      cout << Tag << ":" << endl;
      Len = strlen (Gene);
      assert (Len < MAX_GENE_LEN);

      for  (i = 0;  i < Len;  i += 1)
        if  (strncmp (Gene + i, "tag", 3) == 0
              || strncmp (Gene + i, "taa", 3) == 0
              || strncmp (Gene + i, "tga", 3) == 0)
            cout << "  stop at offset " << setw (3) << i << endl;
     }

   return  0;
  }
