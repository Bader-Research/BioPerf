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
#include "pom.h"
#include <stdlib.h>

#ifndef H_BNASTATS
#define H_BNASTATS

#define MAXVALUE 999.0

enum ValueType { DIHEDRAL, COMMON, BONDLENGTH, DOUBLES };

struct ValsNStddev {
        public:
                double value;
                double stddev;
};

class Summary;
class Stats;

class Summary {
	public:
		char *name;		// the name of the value
		double mean, stdDev;	// mean and stddev for small molecule
		double *values;		// array of values summarized
		int num;		// the number of values summarized
		double total;		// the sum of all values
		double min, max;	// overall minimum and maximum
		int minPos, maxPos;	// the position of min and max
		char *minAA, *maxAA;	// the amino acids of min and max
		char minKS, maxKS;	// KS for the above
		int adjust;

		Summary();
		Summary(const char *, int, double, double, ValueType);
		~Summary();
		void reset();
		void add(ISubentity &, double);
		double stddev();
		double getFDS();
		void toString(char *, char **);
	
	private:
		ValueType myType;
};
	
class Stats {
	public:
		static void printDihedrals(IEntity &, char *);
		static void printDihedrals(IEntity &);
		static Summary *sumDihedrals(IEntity &, int &);

		static void printBondAngles(IEntity &, char *);
		static void printBondAngles(IEntity &);
		static Summary *sumBondAngles(IEntity &, int &);

		static void printBondLengths(IEntity &, char *);
		static void printBondLengths(IEntity &);
		static Summary *sumBondLengths(IEntity &, int &);

		static void printResGeo(ISubentity &);
		static double *getResGeo(ISubentity &, int &);

		static ValsNStddev *getSmallResVals(ISubentity &, int &);
		static ValsNStddev *getSmallResVals(char *, char, int &);

		static double getFDS(ISubentity &);
		static double getFDS(IEntity &);
		static Summary getAllFDS(IEntity &);
};

#endif
