//
//  Programmer:  A. Delcher
//
//        File:  Glimmer/adjust.cc
//
//     Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Arthur Delcher, Steven Salzberg,
//    Simon Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program removes stop codons from coordinate lists
//


#include  "delcher.h"


int  main  ()

  {
   char  Tag [100];
   long int  Start, Stop;

   while  (cin >> Tag >> Start >> Stop)
     {
      cout << setw (8) << setiosflags (ios :: left) << Tag
           << setw (8) << setiosflags (ios :: right) << Start;
      if  (Start < Stop)
          cout << setw (8) << Stop - 3;
        else
          cout << setw (8) << Stop + 3;
      cout << endl;
     }

   return  0;
  }