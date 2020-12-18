
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa34t25d2 $ - $Id: nrand.c,v 1.1.1.1 1999/10/22 20:56:01 wrp Exp $ */

#include <stdlib.h>
#include <time.h>

irand(n)	/* initialize random number generator */
     int n;
{

  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;

  }
  srand(n);
}

nrand(n)	/* returns a random number between 1 and n where n < 64K) */
     int n;
{
  int rand();
  long rn;

  rn = rand();
#ifdef RAND32
  rn = rn >> 16;
#endif
  rn = rn % n;
  return (int)rn;
}




