/*
 
Copyright (c)  1995-2000   The Regents of the University of California
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this PDB
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this PDB Software into commercial products
or use for commercial purposes should contact the Technology Transfer
Office, University of California, San Diego, 9500 Gilman Drive, La Jolla,
CA 92093-0910, Ph: (619) 534-5815, FAX: (619) 534-7345.
 
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS PDB SOFTWARE, EVEN IF THE
UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
 
THE PDB SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE PDB SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
 
*/
///////////////////////////////////////////////////////////////////////////
//  POM PDB Parser (C) 1995-2000 I.Shindyalov, H.Weissig, P. Bourne      //
///////////////////////////////////////////////////////////////////////////
#ifndef H_IPDBDF
#define H_IPDBDF

#include "pom.h"
#include "ipdb.h"

#ifndef MINIMAL
	static	const	int	NUM_DF		= 9;
#else
	static	const	int	NUM_DF		= 4;
#endif

class DF {
   public:
	static Tlpdb *tlpdb;

	static int assignDF(char *, int, int, int, int, int);
	static int assignDF(char *, int, int, int);
	static int assignDF(int, int, int);

	static void bfac(int, int, int, int);			// 1
	static void ss_propens(int, int, int, int);		// 2
	static void obsNspr(int, int, int, int);		// 3
	static void collect(int, int, int, int);		// 4
	static void collFDS(int, int, int, int);		// 5
	static void calcChi1(int, int, int, int);		// 6
	static void proteinFeatures(char *, int, int, int, int);// 7
	static void mooseReport(int, int, int, int);		// 8
	static void status_com(int);	        // 9

	static void printFeatures();

   private:
	DF();	// static class
	static char *flt2str(float *);
	static int indexIds(Property *,char *, int4 **);
};

#endif
