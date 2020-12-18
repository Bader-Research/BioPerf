
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa34t25d2 $ - $Id: nrandom.c,v 1.1 2000/03/02 18:04:59 wrp Exp $ */

#include <stdlib.h>
#include <time.h>

void 
irand(n)	/* initialize random number generator */
     int n;
{
  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;
  }
  srandom(n);
}

int
nrand(n)	/* returns a random number between 0 and n-1 where n < 2^24) */
     int n;
{
  int rn;

  rn = random();
  rn = (rn % n);
  return rn;
}

