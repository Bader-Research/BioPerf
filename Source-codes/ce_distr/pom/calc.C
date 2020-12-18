#include <math.h>
#include <iostream.h>
#include "contour.h"
#include "calc.h"

/*****************************************************************************/
/*                                Class Torsion                              */
/*****************************************************************************/

Point Torsion::X3;
Point Torsion::X4;
Point Torsion::X1;
Point Torsion::X2;

/*****
Function:   Torsion::TorsionAngle (was Torsion::TorsionFunction, HW)
Purpose:    This function figures out the torsion angle by doing vector math.
  This function was translated from phipsi.for, a program written in
  FORTRAN IV for the protein data bank.
Parameters: Point &A, Point &B, Point &C, Point &D
Returns:    double
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/

double Torsion::TorsionAngle(IAtom &A, IAtom &B, IAtom &C, IAtom &D) {
	if (A && B && C && D) {
		Point p1((double) A.x(), (double) A.y(), (double) A.z()),
		      p2((double) B.x(), (double) B.y(), (double) B.z()),
		      p3((double) C.x(), (double) C.y(), (double) C.z()),
		      p4((double) D.x(), (double) D.y(), (double) D.z());
		return (TorsionAngle(p1, p2, p3, p4));
	}

	return (999.0);
}
		
double Torsion::TorsionAngle(Point &A, Point &B, Point &C, Point &D) {
  double return_angle;
  int ichk=-1;
  Point zero;

  zero.SetX(0.0);
  zero.SetY(0.0);
  zero.SetZ(0.0);
  Vecdif(A,B,X3);
  DisChk(X3,&ichk);
  if(ichk > 0)
    return(999.0);
  Vecdif(C,B,X4);
  DisChk(X4,&ichk);
  if(ichk > 0)
    return(999.0);
  Cross(X3,X4,X1);
  Vecdif(B,C,X3);
  Vecdif(D,C,X4);
  DisChk(X4,&ichk);
  if(ichk > 0)
    return(999.0);
  Cross(X3,X4,X2);
  return_angle=CommonAngle(X1,zero,X2);
  Cross(X1,X2,X4);
  X1.SetX(Dot(X3,X4));
  if(X1.X() > 0.0)
    return_angle=-return_angle;
  return(return_angle);
}


/*****
Function:   Torsion::Cross
Purpose:    This function takes the cross product of a and b, and puts the
  result in c.  This function was translated from phipsi.for, a program
  written in FORTRAN IV for the protein data bank.  
Parameters: Point &a, Point &b, Point &c
Returns:    (void)
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/
void Torsion::Cross(Point &a, Point &b, Point &c)
{
  c.SetX(a.Y()*b.Z()-b.Y()*a.Z());
  c.SetY(-a.X()*b.Z()+b.X()*a.Z());
  c.SetZ(a.X()*b.Y()-b.X()*a.Y());
}


/*****
Function:   Torsion::Dot
Purpose:    This function returns the dot product of a and b.  This function
  was translated from phipsi.for, a program written in FORTRAN IV for the
  protein data bank.  
Parameters: Point &a, Point &b
Returns:    double
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/
double Torsion::Dot(Point &a, Point &b)
{
  double val=0.0;
  
  val+=a.X()*b.X();
  val+=a.Y()*b.Y();
  val+=a.Z()*b.Z();

  return(val);
}


/*****
Function:   Torsion::Vecdif
Purpose:    This function takes the difference between the vectors a and b,
  and puts the resulting vector in c.  This function was translated from
  phipsi.for, a program written in FORTRAN IV for the protein data bank.  
Parameters: const Point &a, const Point &b, const Point &c
Returns:    double
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/
void Torsion::Vecdif(Point &a, Point &b, Point &c)
{
  c.SetX(a.X()-b.X());
  c.SetY(a.Y()-b.Y());
  c.SetZ(a.Z()-b.Z());
}


/*****
Function:   Torsion::CommonAngle (was Torsion::Angle, HW)
Purpose:    This function will return the angle formed by the three given
  Atoms, centered around the second atom.
Parameters: const Point &a, const Point &b, const Point &c
Returns:    double
Written On: April 28, 1993
Author:     Alex Chin
Added back to PDBLIB V2 by John Biggs  8-9-94
*****/

/*****
Function:   Torsion::Angle
Purpose:    I'm not quite sure what this function does...  This function was
   translated from phipsi.for, a program written in FORTRAN IV for the
   protein data bank.  
Parameters: Point &a, Point &b, Point &c
Returns:    double
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/

double Torsion::CommonAngle(IAtom &A, IAtom &B, IAtom &C) {
	if (A && B && C) {
		Point p1((double) A.x(), (double) A.y(), (double) A.z()),
		      p2((double) B.x(), (double) B.y(), (double) B.z()),
		      p3((double) C.x(), (double) C.y(), (double) C.z());
		return (CommonAngle(p1, p2, p3));
	}

	return (999.0);
}

double Torsion::CommonAngle(Point &a, Point &b, Point &c) {
  double ang;
  double Q;

  Vecdif(a,b,X1);
  Vecdif(c,b,X2);
  
  Q=Amag(X1)*Amag(X2);
  if(Q < .0000001)
    Q=.0000001;
  ang=Dot(X1,X2)/Q;
  if(ang > 1.0)
    ang=1.0;
  if(ang < -1.0)
    ang=-1.0;
  return(57.29577951*acos(ang));
}


/*****
Function:   Torsion::DisChk
Purpose:    This function tells whether the atom should be a dummy atom,
  thereby returning torsion angles of 999 degrees for all calculations using
  this atom.  This function was translated from phipsi.for, a program written
  in FORTRAN IV for the protein data bank.  
Parameters: Point &x, int *ichk
Returns:    (void)
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/
void Torsion::DisChk(Point &x, int *ichk){
  double DSSQ;

  *ichk=0;
  DSSQ=Dot(x,x);
  if((DSSQ < .50) || (DSSQ > 9.00))
    *ichk=1;
  return;
}


/*****
Function:   Torsion::Amag
Purpose:    This function returns the magnitude of the point, that is, the
  square root of the dot product of the point and itself.  This function was
  translated from phipsi.for, a program written in FORTRAN IV for the
  protein data bank.  
Parameters: Point &x
Returns:    double
Written On: May, 1978
Revised On: September, 1979
Author:     L. Andrews, G. Williams, F. Bernstein
Translated: Karthik Sheka
On:         August 9, 1992
*****/
double Torsion::Amag(Point &x)
{
  return(sqrt(Dot(x,x)));
}


/*****************************************************************************/
/*                               Class Calculate                             */
/*****************************************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Function: Calculate::IsPointOnRamachandranContour
   Purpose: checks to see if a point is in contour
Parameters: phi and psi by reference
   Returns: TRUE iff atom in contour
Written On: July 1992
    Author: Paul Liu
----------------------------------------------------------------------------*/
int Calculate::IsPointOnRamachandranContour(const double &phi, const double &psi)
{
  int crossings = 0, i;
  
  for (i = 0; i < contoursize(rcontour1); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour1,
                                contoursize(rcontour1)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour2a); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour2a,
                                contoursize(rcontour2a)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour2b); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour2b,
                                contoursize(rcontour2b)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour2c); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour2c,
                                contoursize(rcontour2c)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour2d); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour2d,
                                contoursize(rcontour2d)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour3); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour3,
                                contoursize(rcontour3)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour4); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour4,
                                contoursize(rcontour4)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour5); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour5,
                                contoursize(rcontour5)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour6); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour6,
                                contoursize(rcontour6)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
  crossings = 0;
 
  for (i = 0; i < contoursize(rcontour7); i++)
    if (IsRayCrossingSegment(i, phi, psi, rcontour7,
                                contoursize(rcontour7)))
      crossings++;
  if (crossings % 2 == 1)
    return TRUE;
 
  return FALSE;
 
} // end Calculate::is_in_contour

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Function: Calculate::is_ray_crossing_segment
   Purpose: checks to see if a ray extended in the negative x direction
            crosses a line segment
Parameters: index to first point in segment
            phi and psi by reference
            extent of polygon
   Returns: TRUE iff the ray crosses the segment
Written On: July 1992
    Author: Paul Liu
----------------------------------------------------------------------------*/
int Calculate::IsRayCrossingSegment(const int &i, const double &phi,
                                          const double &psi,
                                          const torsion rcontour[],
                                          const int extent)
{
double top, t, deltaphi, deltapsi, deltapsinext;
int j, k;              // index to next point in contour and one after that
int isnegative;        // TRUE iff value negative
 
  j = (i + 1) % extent;
 
  deltaphi = rcontour[j].phi - rcontour[i].phi;
  deltapsi = rcontour[j].psi - rcontour[i].psi;
 
  // ignore horizontal lines
  if (deltapsi == 0.0)
    return FALSE;
 
  // don't check anything equal to or below the bottom of the y range
  if (psi <= min(rcontour[i].psi, rcontour[j].psi))
    return FALSE;
 
  top = max(rcontour[i].psi, rcontour[j].psi);
  if (psi > top)
    return FALSE;
 
  // now solve parametric equation P + t(Q - P) for t using y
  t = (psi - rcontour[i].psi) / deltapsi;
 
  // solve for x on segment given t; if this x is to the right of the point,
  // then a ray extended left from the point will not intersect segment
  if ((rcontour[i].phi + t * deltaphi) > phi)
    return FALSE;
 
  if (t == 1.0)
    { // it intersects the end of the line segment; if the next segment
      // forms a /\ or \/ with the current one, ignore
      // therefore, the slope of the next line must be either zero
      // or have the same sign
    k = (j + 1) % extent;
    deltapsinext = rcontour[k].psi - rcontour[j].psi;
    if (deltapsinext == 0.0)
      return TRUE;
    isnegative = (deltapsi < 0.0);
    if (isnegative != (deltapsinext < 0))
      return FALSE;
    }
 
  return TRUE;
 
} // end Calculate::is_ray_crossing_segment


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Function: Calculate::distance
   Purpose: calculates the distance between two atoms
Parameters: two atoms by reference
   Returns: distance
Written On: July 1992
    Author: Paul Liu
----------------------------------------------------------------------------*/
double Calculate::Distance(const Point &a, const Point &b)
{
double dx, dy, dz;
 
  dx = a.X() - b.X();
  dy = a.Y() - b.Y();
  dz = a.Z() - b.Z();
 
  return sqrt(dx*dx + dy*dy + dz*dz);
 
} // end Calculate::distance

 
/*****
Function:   Calculate::Distance
Purpose:    This function returns the distance between two points in
  three-space.
Parameters: double x1, double y1, double z1, double x2, double y2, double z2
Returns:    double
Written On: August 17, 1992
Author:     Karthik Sheka
*****/
double Calculate::Distance(double x1, double y1, double z1,
		    double x2, double y2, double z2)
{
  return(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)));
}


/*****
Function:   Calculate::Distance
Purpose:    This function returns the distance between two points in
  two-space.
Parameters: double x1, double y1, double x2, double y2
Returns:    double
Written On: August 17, 1992
Author:     Karthik Sheka
*****/
double Calculate::Distance(double x1, double y1, double x2, double y2)
{
  return(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

double Calculate::Distance(IAtom &a, IAtom &b) {
	if (a && b) {
		double dx, dy, dz;
		dx = (double) (a.x() - b.x());
		dy = (double) (a.y() - b.y());
		dz = (double) (a.z() - b.z());

		return sqrt(dx*dx + dy*dy + dz*dz);
	}

	return 999.0;
}

/****************************************************************************

    SUBROUTINE RAMREG  -  Determine which part of the Ramachandran plot
                          this residue is in

****************************************************************************/

void Calculate::RamReg(float *phi, float *psi,
		       char **region, int *regno, int *regtyp) {

    static int first = 1;
    static char code[12+1] = "oBbAaLlexyzw";
    static char newcod[24+1] = "XXB b A a L l p ~b~a~l~p";
    static int rtype[12+1] = { 1,4,3,4,3,4,3,3,2,2,2,2 };
    static char mapstr[36*36+1] = "bBBBBBBBBBBbbbxxooooowwwwwoooooxxxxbbBBB"
	    "BBBBBBBBbbbxxooooooooooooooxxxbbbBBBBBBBBBBBbbbxxxxoooooooooo"
	    "ooxxbbbbbBBBBBBBBBBBbbbbxxooooooooooooxxxbbbbBBBBBBBBBBBBb"
	    "bxxxooooooooooooxxxxbbbBBBBBBBBBBBbbbxxxoooooooooooooxxxb" 
	    "bbbBBBBBBBBBbbbbxxxooozzzzzzoooooxxbbbbbBBBBBBBbbbbbbxxzzzzz"
	    "zzzzoooooxxbbbbbbBbbBbbbbbbbbxxzzzzzllzzoooooxxbbbbbbbbbbb"
	    "bbbbxxxxxzzllllzzzoooooxxxbbbbbbbbbbbbxxxxxxxzzllllzzzoooooxxx"
	     "xbbbbbbbbbbbbxxoooozzzlllzzzzooooxxxxxbbbbbbbbbbbxxoooozzll"
	    "llllzzooooxxxyyaaaaaaaaaayyyoooozzllLllzzzooooxxxyaaaaaaaa"
	    "aaayyyoooozzzlLLlzzzooooxxxyaaaaaaAaaaaayyyooozzzlllllzzoooox"
	    "xxyaaaaaAAAAaaayyyyooozzzlllzzzooooxxxyaaaaaAAAAAaaayyyooo"
	    "zzzzlllzzooooxxxyaaaaAAAAAAAaaayyyoozzzlzllzzooooxxxyaaaaa"
	    "AAAAAAAaayyyyozzzzzzzzzooooxxxyyaaaaAAAAAAAAaayyyozzzzzzzzzoo"
	    "ooxxxyyaaaaaAAAAAAAaaayyyoooooooooooooxxxyyyaaaaaAAAAAAAaa"
	    "yyyoooooooooooooxxxyaaaaaaaaAAAAAAaaayyyooooooooooooxxxaay"
	    "aaaaaaaaAAAAaaayyyooooooooooooxxxyyyaaaaaaaaaaaaaaaayyooooooo"
	    "oooooxxxyyyyyaaaaaaaaaaaaayyyooooooooooooxxxoyyyyaaaaaaaaa"
	    "ayyyyyyooooooooooooooooooyyyyyyyyayyyyyyyyoooooooooooooooo" 
	    "oooxxxxbbxxxxxxxxoooooooooooooooooooxxxxxbbxxxxxxxooooooowww"
	    "wwooooooooooxxxxxxbbbbxxxxooooooowwwwwooooooooooxbbbxbbbbb"
	    "xxxxxoooooowwewwwoooooooooxbbbbbbbbbbxxxxoooooowwewwwooooooxxx"
	    "xbbbbbbbbbbbbxxoooooowweewwooooooxxxbbbbbbbbbbbbbxxooooooww"
	    "wwwwooooooxxb";

    /* Local variables */
    static char char__[1];
    static int i__, j, k;
    static float gap;
    static int map[1296+1]	/* was [36][36] */;

/* If this is the first call to this routine, initialise the MAP array */

    if (first) {
	first = 0;
	for (i__ = 1; i__ <= 36; ++i__) {
	    for (j = 1; j <= 36; ++j) {
		*char__ = mapstr[(i__ - 1) * 36 + (j - 1)];
		for (k = 1; k <= 12; ++k) {
		    if (*char__ == code[k - 1]) {
			map[36 - i__ + 1 + j * 36 - 37] = k;
		    }
		}
	    }
	}
	gap = 10.f;
    }

/* Determine which part of the Ramachandran plot the residue is in */

    if (*phi > 180.f || *psi > 180.f) {
	*regno = 1;
	strncpy(*(region), "XX", 2);
	*regtyp = 0;
    } else {
	i__ = (int) ((*psi + 180.f) / gap + 1);
	if (i__ < 1) {
	    i__ = 1;
	}
	if (i__ > 36) {
	    i__ = 36;
	}
	j = (int) ((*phi + 180.f) / gap + 1);
	if (j < 1) {
	    j = 1;
	}
	if (j > 36) {
	    j = 36;
	}
	*regno = map[i__ + j * 36 - 37];
	i__ = (*regno << 1) - 1;
	strncpy(*(region), newcod + (i__ - 1), 2);
	*regtyp = rtype[*regno - 1];
    }
}
