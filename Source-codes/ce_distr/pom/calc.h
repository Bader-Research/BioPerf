/**************************************************************************
Copyright (C) 1992, 1993 by The Trustees of Columbia University
in the City of New York

OO-PDB: Object Oriented Protein Data Bank Project.
Biochemistry and Molecular Biophysics and Computer Science Departments


 Filename: calc.h
 Author:   Karthik Sheka, Weider Chang
 Date:     Oct 18, 1994
 
 modified for use in calc.C by Helge Weissig, Ph.D., June 2, 1997
 
 **************************************************************************/

#include "pom.h"

#ifndef _CALC_H_
#define _CALC_H_

#ifndef CONTOUR_H

typedef struct {
  double phi, psi;
} torsion;
/*Forward declaration for contour.h structure.*/
#endif

#define ABS(X) (X>0)?X:-(X);

class Point {

private:
  double _x, _y, _z;

public:
  Point(double x, double y, double z) { _x = x; _y = y; _z = z; }
  Point() {_x=_y=_z=0.0;}
  double X() const {return _x;}
  double Y() const {return _y;}
  double Z() const {return _z;}
  void SetX(double x){_x=x;}
  void SetY(double y){_y=y;}
  void SetZ(double z){_z=z;}
  void operator +=(const Point& point) { 
    _x+=point.X(); _y+=point.Y(); _z+=point.Z(); 
  }
  void operator -=(const Point& point){ 
    _x-=point.X(); _y-=point.Y(); _z-=point.Z(); 
  }
  void operator *=(const Point& point) {
    _x*=point.X(); _y*=point.Y(); _z*=point.Z();
  } 
  void operator /=(const Point& point) {
    _x/=point.X(); _y/=point.Y(); _z/=point.Z();
  } 
  void operator /=(double d){_x/=d; _y/=d; _z/=d;}
  void operator *=(double d) {_x*=d; _y*=d; _z*=d;}
};

class Torsion{
private:
  static Point X1;
  static Point X2;
  static Point X3;
  static Point X4;
  
  static void Cross(Point &a, Point &b, Point &c);
  static double Dot(Point &a, Point &b);
  static void Vecdif(Point &a, Point &b, Point &c);
  static void DisChk(Point &x, int *ichk);
  static double Amag(Point &x);
  
public:
  static double TorsionAngle(IAtom &a, IAtom &b, IAtom &c, IAtom &d);
  static double TorsionAngle(Point &a, Point &b, Point &c, Point &d);
  static double CommonAngle(IAtom &a, IAtom &b, IAtom &c );
  static double CommonAngle(Point &a, Point &b, Point &c );

};

class Calculate:public Torsion{
private:
  static int IsRayCrossingSegment(const int &i, const double &phi,
				  const double &psi,
				  const torsion rcontour[],
				  const int extent);
public:
  static double Distance(IAtom &a, IAtom &b);
  static double Distance(const Point &a, const Point &b);
  static double Distance(double x1, double y1, double z1,
			 double x2, double y2, double z2);
  static double Distance(double x1, double y1, double x2, double y2);
  static int IsPointOnRamachandranContour(const double &phi, 
					  const double &psi);
  static void RamReg(float *, float *, char **, int *, int*);
};

#endif

/**************************************************************************
    Copyright (C) 1992, 1993 by The Trustees of Columbia University 
    in the City of New York

    OO-PDB: Object Oriented Protein Data Bank Project.
    Biochemistry and Molecular Biophysics and Computer Science Departments


    Filename: misc.h
    Author:   Karthik Sheka
    Date:     February 15, 1993

    partially incorporated into calc.C by Helge Weissig, Ph.D.
    June 2, 1997

**************************************************************************/
 
#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef YES
#define YES TRUE
#endif
#ifndef NO
#define NO FALSE
#endif
#ifndef max
#define max(A, B) ((A) > (B)? (A) : (B))
#endif
#ifndef min
#define min(A, B) ((A) < (B)? (A) : (B))
#endif
