#include "correction.h"


/* call ede; 
   invdist 	the minimum inversion distance,
   ngene 	the number of genes. 
   The function returns the EDE distance */

int ede(int invdist, int ngene) {
    double ll, tt, kk, pp, dval;
    int newvalue;

    kk=invdist/(ngene+0.0);

    if (kk>=0.999999999999) {	/* the distance correction has singularity at 1 */
	kk=0.999999999999;
    }
    if (kk <= 1-ede_c) return invdist;

    ll=ede_c*kk-ede_a;
    tt=4*ede_a*(1- kk)*kk + ll*ll;
    tt=ll + sqrt(tt);
    pp=tt/(2*(1-kk));
    pp*=ngene;

    dval=pp;
    newvalue = (int)ceil(dval);
    /*if (newvalue-dval > 0) return newvalue-1;*/
    return newvalue;
}


