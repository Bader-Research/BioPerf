///////////////////////////////////////////////////////////////////////////
//  Combinatorial Extension Algorithm for Protein Structure Alignment    //
//                                                                       //
//  Authors: I.Shindyalov & P.Bourne                                     //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  1997-1999   The Regents of the University of California
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this CE
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this CE Software into commercial products
or use for commercial purposes should contact the Technology Transfer
Office, University of California, San Diego, 9500 Gilman Drive, La Jolla,
CA 92093-0910, Ph: (619) 534-5815, FAX: (619) 534-7345.
 
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS CE SOFTWARE, EVEN IF THE
UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
 
THE CE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE CE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
 
*/
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <time.h>
#include <dirent.h>
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
 
#include <unistd.h>

//#include "../mmlib/mmlib.h"
#include "../pom/pom.h"
#include "../pom/miscutil.h"

#ifndef H_CMP_UTIL
#define H_CMP_UTIL
#define PI 3.1415926
//typedef unsigned char int1;
typedef char int1;
///////////////////////////////////////////////////////////////////////////////
//                  cmp_util.C                                               //
///////////////////////////////////////////////////////////////////////////////
void readList(char *file, char ***list, int& nList);
void scanList(char *buffer, char ***list, int& nList);
void toUpperCase(char *buffer);
void toLowerCase(char *buffer);
////////////////////////////////////////////////////////////////////
double dpAlign(double **mat, int nSeq1, int nSeq2, int *align_se1, 
	       int *align_se2, int &lAlign, double gapI, double gapE, 
	       int isGlobal1, int isGlobal2);
///////////////////////////////////////////////////////////////////////////////
void alignToSeq(int *align1, int1 *seq1, int nse1, int1 *seq2, int nse2, 
		int1 **seq1a, int1 **seq2a, int& nSeqa, 
		int& seqMatchFrom, int& seqMatchTo, int endsMode=0);
void alignToRef(int *align1, int nse1, int nse2, 
		int **seAlign1, int **seAlign2, int **alignSE1, int **alignSE2,
		int& nSeqa, int endsMode);
void seqCharToInt(char *seq1, char *seq2, int nseq);
void align_se_to_se_align(int* align_se, int lcmp, int nse, int *se_align);
void align_se_to_align_rep(int* align_se, int lcmp, int nse, int *align_rep);
int lali_of_trace(int *trace, int nTrace);
void align_se_to_trace(int* align_se, int lcmp, int4 **trace, int& nTrace);
void align_se_to_trace(int* align_se1, int* align_se2, int lcmp, int4 **trace, 
		       int& nTrace);
void trace_to_align_se(int *trace, int nTrace, int** align_se, int& lcmp);
void sortAlign(int item, int *align_se, int *ent, int lcmp, int ncmp);
///////////////////////////////////////////////////////////////////////////////
int sup_str(XYZ *, XYZ *, int,  double [20]);
int matinv(double [20][20], double [20], int n);
double det3(double r[3][3]);
double calc_rmsd(XYZ *, XYZ *, int,  double [20]);
void rot_mol(XYZ *molA, XYZ *molB, int nAtom, double d_[20]);
XYZ *arrayToXYZ(float *array, int nXYZ);
///////////////////////////////////////////////////////////////////////////////
double dAngle(XYZ &A, XYZ &B, XYZ &C, XYZ &D);
void Cross(XYZ &a, XYZ &b, XYZ &c);
double Dot(XYZ &a, XYZ &b);
void Vecdif(XYZ &a, XYZ &b, XYZ &c);
double Angle(XYZ &a, XYZ &b, XYZ &c);
void DisChk(XYZ &x, int *ichk);
double Amag(XYZ &x);
///////////////////////////////////////////////////////////////////////////////
//                  cmp_tree.C                                               //
///////////////////////////////////////////////////////////////////////////////
void printTree(double **d, char **names, int *top1, int *top2, int n, 
	       char method);
int getTree(double **d, char **tree_pic, int *top0, int *top1, int n_, 
	    char met_set, int *tree_ind);
double nj_f(int i1, int i2, int n, double **d);
///////////////////////////////////////////////////////////////////////////////
//                  cmp_draw.C                                               //
///////////////////////////////////////////////////////////////////////////////
class Draw
{
public:
  static unsigned char sym[256][14];
  char colorTable[256][3];
  char *map;
  int nX, nY;
  int mapSize;
  Draw(int nX_, int nY_);
  ~Draw();
  void outtext(int ix, int iy, unsigned char *buffer, int col);
  void line(int x1, int y1, int x2, int y2, int col);
  void rect(int ix, int iy, int lx, int ly, int col);
  void distr(int ix, int iy, char *title, 
	     int xmin, int xmax, int ymin, int ymax, 
	     int *array, int nx, int ny, int width, int xmarks,
	     int ymarks, int col1, int col2);
  void write(char *file = NULL);
  void setColor(int index, int red, int green, int blue);
};
/////////////////////////////////////////////////////////////////////
//       Kabsch (1976) superposition routines

double kabsch(XYZ *strBuf1, XYZ *strBuf2, int lali_str, double d[20]);
int cubic_roots (double c[4], double s[3]);
void eigen_values (float m[3][3], float values[3], float vectors[3][3]);
void transformpoint (float p_new[3], float m[3][3], float p[3]);
void matrix_multiply (float m[3][3], float a[3][3], float b[3][3]);
void matrix_transpose (float result[3][3], float a[3][3]);
void normalise (float a[3]);

/////////////////////////////////////////////////////////////////////
#endif
