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
/* $Revision: 2.1 $ $Date: 1999/08/31 16:21:19 $ */
#include "pom.h"

class LinkedId {
    public:
	int index;
	static DB *db;

	LinkedId(int);
	~LinkedId();
	int getChain(int **);
	static int getChain(int, int **);
	int getObsNSpr(int **);
	static int getObsNSpr(int, int **);
	static void replacers(int **, int);

    private:
	static Property id_com;
	static Property obs_com;
	static Property spr_com;
	static Property relDat;
	static int open;

	static void init();
	static void addChain(int, int **, int *, int **);
	static void addToArray(int **, int, int);
};

class IVersion {
    public:
	static DB *db;

	IVersion();
	IVersion(int);
	~IVersion();
	operator void*() { return((void*) (iId < nIds)); };
	void operator++() { nextSet(); };
	int getChain(int **);

    private:
	int *versions, nVers, iversions, *dbIds, iId, nIds;

	void nextSet();
};
