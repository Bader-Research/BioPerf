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
#include "bnastats.h"
#include "calc.h"

template <class A> A abs(A num) { return (num < 0.0 ? -(num) : num); };

template <class A> void DeleteArray(A ** array, int num) {
	for (int i = 0; i < num; i++)
		if (*(*array+i)) delete [] (*array+i);
	if (num) delete [] (*array);
	(*array) = NULL;
};

Summary::Summary(const char *valueName, int nSE, double meanVal,
			   double stanDev, ValueType type) {
	name = new char[strlen(valueName)+1];
	strcpy(name, valueName);
	mean = meanVal;
	stdDev = stanDev;
	myType = type;

	minAA = new char[4]; maxAA = new char[4];

	values = new double[nSE];
	for (int n = 0; n < nSE; n++) values[n] = 0.0;

	adjust = (type == DIHEDRAL ? 360 : 0);

	reset();
}

Summary::Summary() {
	name = NULL;
	myType = DIHEDRAL;
	values = NULL;
	adjust = 0;
	reset();
}

Summary::~Summary() {
	if (name) delete [] name;
	if (values) delete [] values;
	if (minAA) delete [] minAA;
	if (maxAA) delete [] maxAA;
}

void Summary::reset() {
	num             = 0;
	total           = 0.0;
	min             = -MAXVALUE;
	max             = MAXVALUE;

	if (myType != DIHEDRAL) { max *= -1; min *= -1; }

	minPos = maxPos = -1;
	minKS  = maxKS  = '\0';
}

void Summary::add(ISubentity &iS, double value) {
	if (value < MAXVALUE) {
		values[num++] = value;
		total += value; if (value < 0.0) total += adjust;
		if (adjust) {
			if (value < 0.0 && value > min) {
				min    = value;
				strcpy(minAA, iS.code3());
				minPos = iS() + 1;
				minKS = iS.ks();
			}

			if (value > 0.0 && value < max) {
				max    = value;
				strcpy(maxAA, iS.code3());
				maxPos = iS() + 1;
				maxKS = iS.ks();
			}
		} else {
			if (value < min) {
				min = value;
				strcpy(minAA, iS.code3());
				minPos = iS() + 1;
				minKS = iS.ks();
			}

			if (value > max) {
				max = value;
				strcpy(maxAA, iS.code3());
				maxPos = iS() + 1;
				maxKS = iS.ks();
			}
		}
	}
}

double Summary::stddev() {
	if (num) {
		double ssq = 0.0, diff;
		double ave = total/num;

		for (int i = 0; i < num; i++) {
			diff = values[i] - ave;
			if (values[i] < 0.0) diff += adjust;
			ssq += diff*diff;
		}

		ssq /= num;
		return sqrt(ssq);
	}

	return MAXVALUE;
}

double Summary::getFDS() {
	if (num && stdDev) {
		double ave = total/num;
		ave = (ave < 0.0 && mean > 0.0 ? -1.0*ave : ave);
		double diff = ave - mean;
		diff = (diff > 0.0 ? diff : -1.0*diff);

		return (diff/stdDev);
	}

	return MAXVALUE;
}

void Summary::toString(char *format, char **out) {
	*(*out) = '\0';
	if (num) {
		sprintf((*out), format, name, num, total/num, stddev(), min,
			minAA, minPos, max, maxAA, maxPos);
	}
}

void Stats::printDihedrals(IEntity &iE) {
	char *f = {"%-10s %5.d %7.2f %7.2f  %7.2f %s %3.d  %7.2f %s %3.d\n"};
	printf("name        total   ave   stddev");
	printf("  max (-)   at     min (+)   at\n");

	printDihedrals(iE, f);
}

void Stats::printDihedrals(IEntity &iE, char *format) {
	int count;
	Summary *all = sumDihedrals(iE, count);

	char *out = new char[500];

	for (int n = 0; n < count; n++) {
		all[n].toString(format, &out); printf("%s", out);
	}
	printf("\n");

	delete [] out;
	if (all) delete [] all;
}

Summary *Stats::sumDihedrals(IEntity &iE, int &count) {

	IAtom CA, N, C, CA1, N1, C1, CG, CB;
	double phi = MAXVALUE, psi, omega, phi1, chi1;
	char *amino, ks;
	int i = 1, nSE;

	nSE = iE.NSE();

	count = 11;
	Summary *out = new Summary[count];

	Summary *psiS = new Summary("Psi", nSE, MAXVALUE, 0.0, DIHEDRAL);
	Summary *psiSGLY = new Summary("Psi (GLY)",nSE,MAXVALUE,0.0,DIHEDRAL);
	Summary *psiSHelix = new Summary("Psi helix",nSE,-39.4,11.3,DIHEDRAL);
	psiSHelix->adjust = 0;
	psiSHelix->min = MAXVALUE; psiSHelix->max = -MAXVALUE;

	Summary *phiS = new Summary("Phi", nSE, MAXVALUE, 0.0, DIHEDRAL);
	Summary *phiSGLY = new Summary("Phi (GLY)",nSE,MAXVALUE,0.0,DIHEDRAL);
	Summary *phiSPRO = new Summary("Phi (PRO)", nSE, -65.4, 11.2, DIHEDRAL);
	phiSPRO->adjust = 0;
	phiSPRO->min = MAXVALUE; phiSPRO->max = -MAXVALUE;

	Summary *phiSHelix = new Summary("Phi helix",nSE,-65.3,11.9,DIHEDRAL);
	phiSHelix->adjust = 0;
	phiSHelix->min = MAXVALUE; phiSHelix->max = -MAXVALUE;

	Summary *omS = new Summary("Omega", nSE, 180.0, 5.8, DIHEDRAL);

	Summary *chi1gm = new Summary("Chi1 g(-)",nSE, 64.1, 15.7, DIHEDRAL);
	Summary *chi1trans = new Summary("CHI1 trans",nSE,183.6,16.8,DIHEDRAL);
	Summary *chi1gp = new Summary("Chi1 g(+)", nSE, -66.7, 15.0, DIHEDRAL);
	chi1gp->min = chi1trans->min = chi1gm->min = MAXVALUE;
	chi1gp->max = chi1trans->max = chi1gm->max = -MAXVALUE;
	chi1gp->adjust = chi1trans->adjust = chi1gm->adjust = 0;

	ISubentity iS(iE);
	ISubentity iS0;

	CA = iS.findAtom(" CA ");
	N  = iS.findAtom(" N  ");
	C  = iS.findAtom(" C  ");
	CB = iS.findAtom(" CB ");

	amino = iS.code3();

	if (!strcmp(amino, "VAL") ||
	    !strcmp(amino, "ILE"))
		CG = iS.findAtom(" CG1");

	else if (!strcmp(amino, "THR"))
		CG = iS.findAtom(" OG1");

	else if (!strcmp(amino, "SER"))
		CG = iS.findAtom(" OG ");

	else CG = iS.findAtom(" CG ");

	ks = iS.ks();

	iS0 = iS; ++iS;

	for (; iS; ++iS, ++i) {
		CA1 = iS.findAtom(" CA ");
		N1  = iS.findAtom(" N  ");
		C1  = iS.findAtom(" C  ");

		phi1  = Torsion::TorsionAngle(C, N1, CA1, C1);
		psi   = Torsion::TorsionAngle(N, CA, C, N1);
		omega = Torsion::TorsionAngle(CA, C, N1, CA1);
		chi1  = Torsion::TorsionAngle(CG, CB, CA, N);

		if (!iS.monType() || !iS0.monType()) {
			phi1 = psi = omega = chi1 = MAXVALUE;
		}

		if (ks == 'H') {
			phiSHelix->add(iS0, phi);
			psiSHelix->add(iS0, psi);
		}

		if (strcmp(amino, "GLY")) {
			psiS->add(iS0, psi);
			if (strcmp(amino, "PRO")) phiS->add(iS0, phi);
		}

		if (!strcmp(amino, "PRO")) phiSPRO->add(iS0, phi);

		if (!strcmp(amino, "GLY")) {
			psiSGLY->add(iS0, psi);
			phiSGLY->add(iS0, phi);
		}

		omS->add(iS0, omega);

		if (!strcmp(amino, "GLY") ||
		    !strcmp(amino, "ALA") ||
		    !strcmp(amino, "ASX") ||
		    !strcmp(amino, "PRO") ||
		    !strcmp(amino, "GLX")) chi1 = MAXVALUE;

		if (chi1 < -120.0) chi1 += 360;

		if (chi1 > 120.0) chi1trans->add(iS0, chi1);
		else if (chi1 >= 0.0) chi1gm->add(iS0, chi1);
		else chi1gp->add(iS0, chi1);

		CA = CA1; N = N1; C = C1;
		phi = phi1;
		amino = iS.code3(); ks = iS.ks(); iS0 = iS;

		CB = iS.findAtom(" CB ");
		if (!strcmp(amino, "VAL") ||
		    !strcmp(amino, "ILE"))
			CG = iS.findAtom(" CG1");

		else if (!strcmp(amino, "THR"))
			CG = iS.findAtom(" OG1");

		else if (!strcmp(amino, "SER"))
			CG = iS.findAtom(" OG ");

		else CG = iS.findAtom(" CG ");
	}

	psiS->add(iS, MAXVALUE);
	omS->add(iS, MAXVALUE);

	if (ks == 'H') phiSHelix->add(iS0, phi);

	if (strcmp(amino, "GLY") && strcmp(amino, "PRO")) {
		phiS->add(iS0, phi);
	} else {
		if (!strcmp(amino, "PRO")) {
			phiSPRO->add(iS0, phi);
		} else {
			phiSGLY->add(iS0, phi);
		}
	}

	out[0]  = (*phiS);
	out[1]  = (*phiSPRO);
	out[2]  = (*phiSGLY);
	out[3]  = (*phiSHelix);
	out[4]  = (*psiS);
	out[5]  = (*psiSGLY);
	out[6]  = (*psiSHelix);
	out[7]  = (*omS);
	out[8]  = (*chi1gm);
	out[9]  = (*chi1trans);
	out[10] = (*chi1gp);

	return out;
}

void Stats::printBondAngles(IEntity &iE, char *format) {
	char *out = new char[500];
	int count;

	Summary *all = sumBondAngles(iE, count);
	for (int n = 0; n < count; n++) {
		all[n].toString(format, &out); printf("%s", out);
	}

	printf("\n");
	delete [] out;
	if (all) delete [] all;
}

void Stats::printBondAngles(IEntity &iE) {
	char *f = {"%-15s %5.d %7.2f %7.2f  %7.2f %s %3.d  %7.2f %s %3.d\n"};
	printf("name             total   ave    stddev");
	printf("   min     at       max     at\n");

	printBondAngles(iE, f);
}

Summary *Stats::sumBondAngles(IEntity &iE, int &count) {
	IAtom N, CA, CB, C, O, N1, CA1;
	IAtom *iA0 = new IAtom(); iA0->entities = iE.entities;

	char *amino;
	int i = 1, nSE;

	nSE = iE.NSE();

	count = 20; Summary *out = new Summary[count];

	Summary *N_CA_C = new Summary("N-CA-C", nSE, 111.2, 2.8, COMMON);
	Summary *N_CA_CPro = new Summary("N-CA-C (P)", nSE, 111.8, 2.5, COMMON);
	Summary *N_CA_CGly = new Summary("N-CA-C (G)", nSE, 112.5, 2.9, COMMON);
	Summary *N_CA_CB = new Summary("N-CA-CB", nSE, 110.5, 1.7, COMMON);
	Summary *N_CA_CBAla = new Summary("N-CA-CB (A)",nSE,110.4,1.5,COMMON);
	Summary *N_CA_CBPro = new Summary("N-CA-CB (P)",nSE,103.0,1.1,COMMON);
	Summary *N_CA_CBITV=new Summary("N-CA-CB (I,T,V)",nSE,111.5,1.7,COMMON);
	Summary *CA_C_O = new Summary("CA-C-O", nSE, 120.8, 1.7, COMMON);
	Summary *CA_C_OGly = new Summary("CA-C-O (G)", nSE, 120.8, 2.1, COMMON);
	Summary *CA_C_N1 = new Summary("CA-C-N", nSE, 116.2, 2.0, COMMON);
	Summary *CA_C_N1Pro = new Summary("CA-C-N (P)",nSE,116.9,1.5,COMMON);
	Summary *CA_C_N1Gly = new Summary("CA-C-N (G)",nSE,116.4,2.1,COMMON);
	Summary *CB_CA_C = new Summary("CB-CA-C", nSE, 110.1, 1.9, COMMON);
	Summary *CB_CA_CAla = new Summary("CB-CA-C (A)",nSE,110.5,1.5,COMMON);
	Summary *CB_CA_CITV=new Summary("CB-CA-C (I,T,V)",nSE,109.1,2.2,COMMON);
	Summary *O_C_N1 = new Summary("O-C-N", nSE, 123.0, 1.6, COMMON);
	Summary *O_C_N1Pro = new Summary("O-C-N (P)", nSE, 122.0, 1.4, COMMON);
	Summary *C_N1_CA1 = new Summary("C-N-CA", nSE, 121.7, 1.8, COMMON);
	Summary *C_N1_CA1Pro = new Summary("C-N-CA (P)",nSE,122.6,5.0,COMMON);
	Summary *C_N1_CA1Gly = new Summary("C-N-CA (G)",nSE,120.6,1.7,COMMON);

	ISubentity iS(iE);
	ISubentity iS0;

	CA = iS.findAtom(" CA ");
	CB = iS.findAtom(" CB ");
	O  = iS.findAtom(" O  ");
	N  = iS.findAtom(" N  ");
	C  = iS.findAtom(" C  ");
	amino = iS.code3();

	iS0 = iS; ++iS;

	for (; iS; ++iS, ++i) {
		CA1 = iS.findAtom(" CA ");
		N1  = iS.findAtom(" N  ");

		if (!iS.monType() || !iS0.monType()) {
			CA = CB = O = N = C = CA1 = N1 = (*iA0);
		}

		if (strcmp(amino, "ALA") && strcmp(amino, "ILE") &&
		    strcmp(amino, "THR") && strcmp(amino, "VAL") &&
		    strcmp(amino, "PRO")) {
			O_C_N1->add(iS0, Torsion::CommonAngle(O, C, N1));
			N_CA_CB->add(iS0, Torsion::CommonAngle(N, CA, CB));
			CB_CA_C->add(iS0, Torsion::CommonAngle(CB, CA, C));
		}

		if (!strcmp(amino, "ALA")) {
		   N_CA_CBAla->add(iS0, Torsion::CommonAngle(N, CA, CB));
		   CB_CA_CAla->add(iS0, Torsion::CommonAngle(CB, CA, C));
		}

		if (!strcmp(amino, "ILE") || !strcmp(amino, "THR") ||
		    !strcmp(amino, "VAL")) {
		   N_CA_CBITV->add(iS0, Torsion::CommonAngle(N, CA, CB));
		   CB_CA_CITV->add(iS0, Torsion::CommonAngle(CB, CA, C));
		}

		if (!strcmp(amino, "PRO")) {
		    CA_C_N1Pro->add(iS0, Torsion::CommonAngle(CA, C, N1));
		    N_CA_CPro->add(iS0, Torsion::CommonAngle(N, CA, C));
		    N_CA_CBPro->add(iS0, Torsion::CommonAngle(N, CA, CB));
		    CB_CA_C->add(iS0, Torsion::CommonAngle(CB, CA, C));
		    O_C_N1Pro->add(iS0, Torsion::CommonAngle(O, C, N1));
		}

		if (strcmp(amino, "GLY") && strcmp(amino, "PRO")) {
			CA_C_O->add(iS0, Torsion::CommonAngle(CA, C, O));
			CA_C_N1->add(iS0, Torsion::CommonAngle(CA, C, N1));
			N_CA_C->add(iS0, Torsion::CommonAngle(N, CA, C));
		} else {
		   if (!strcmp(amino, "GLY")) {
		     CA_C_OGly->add(iS0, Torsion::CommonAngle(CA, C, O));
		     CA_C_N1Gly->add(iS0, Torsion::CommonAngle(CA, C, N1));
		     N_CA_CGly->add(iS0, Torsion::CommonAngle(N, CA, C));
		   } else {
		     CA_C_O->add(iS0, Torsion::CommonAngle(CA, C, O));
		   }
		}

		amino = iS.code3();

		if (strcmp(amino, "GLY") && strcmp(amino, "PRO")) {
		   C_N1_CA1->add(iS, Torsion::CommonAngle(C, N1, CA1));
		} else {
		   if (!strcmp(amino, "GLY")) {
		     C_N1_CA1Gly->add(iS,Torsion::CommonAngle(C,N1,CA1));
		   } else {
		     C_N1_CA1Pro->add(iS,Torsion::CommonAngle(C,N1,CA1));
		   }
		}

		CB = iS.findAtom(" CB ");
		O  = iS.findAtom(" O  ");
		C  = iS.findAtom(" C  ");
		N = N1; CA = CA1; iS0 = iS;
	}

	if (strcmp(amino, "ALA") && strcmp(amino, "ILE") &&
	    strcmp(amino, "THR") && strcmp(amino, "VAL") &&
	    strcmp(amino, "PRO")) {
		N_CA_CB->add(iS0, Torsion::CommonAngle(N, CA, CB));
		CB_CA_C->add(iS0, Torsion::CommonAngle(CB, CA, C));
	}

	if (!strcmp(amino, "ALA")) {
	   N_CA_CBAla->add(iS0, Torsion::CommonAngle(N, CA, CB));
	   CB_CA_CAla->add(iS0, Torsion::CommonAngle(CB, CA, C));
	}

	if (!strcmp(amino, "ILE") || !strcmp(amino, "THR") ||
	    !strcmp(amino, "VAL")) {
	   N_CA_CBITV->add(iS0, Torsion::CommonAngle(N, CA, CB));
	   CB_CA_CITV->add(iS0, Torsion::CommonAngle(CB, CA, C));
	}

	if (!strcmp(amino, "PRO")) {
	    N_CA_CPro->add(iS0, Torsion::CommonAngle(N, CA, C));
	    N_CA_CBPro->add(iS0, Torsion::CommonAngle(N, CA, CB));
	}

	if (strcmp(amino, "GLY") && strcmp(amino, "PRO")) {
		CA_C_O->add(iS0, Torsion::CommonAngle(CA, C, O));
		N_CA_C->add(iS0, Torsion::CommonAngle(N, CA, C));
	} else {
	   if (!strcmp(amino, "GLY")) {
	     CA_C_OGly->add(iS0, Torsion::CommonAngle(CA, C, O));
	     N_CA_CGly->add(iS0, Torsion::CommonAngle(N, CA, C));
	   }
	}

	out[0]  = (*N_CA_C);
	out[1]  = (*N_CA_CPro);
	out[2]  = (*N_CA_CGly);
	out[3]  = (*N_CA_CB);
	out[4]  = (*N_CA_CBAla);
	out[5]  = (*N_CA_CBPro);
	out[6]  = (*N_CA_CBITV);
	out[7]  = (*CA_C_O);
	out[8]  = (*CA_C_OGly);
	out[9]  = (*CA_C_N1);
	out[10] = (*CA_C_N1Pro);
	out[11] = (*CA_C_N1Gly);
	out[12] = (*CB_CA_C);
	out[13] = (*CB_CA_CAla);
	out[14] = (*CB_CA_CITV);
	out[15] = (*O_C_N1);
	out[16] = (*O_C_N1Pro);
	out[17] = (*C_N1_CA1);
	out[18] = (*C_N1_CA1Pro);
	out[19] = (*C_N1_CA1Gly);

	delete iA0;

	return out;
}

void Stats::printBondLengths(IEntity &iE) {
	char *f = {"%-13s %5.d %7.2f %7.2f  %7.2f %s %3.d  %7.2f %s %3.d\n"};
	printf("name          total   ave    stddev");
	printf("   min     at       max     at\n");
	printBondLengths(iE, f);
}

void Stats::printBondLengths(IEntity &iE, char *format) {
	int count;
	Summary *all = sumBondLengths(iE, count);

	char *out = new char[500];

	for (int n = 0; n < count; n++) {
		all[n].toString(format, &out); printf("%s", out);
	}
	printf("\n");

	delete [] out;
	if (all) delete [] all;
}

Summary *Stats::sumBondLengths(IEntity &iE, int &count) {
	IAtom C0, N, C, O, CA, CB;

	IAtom *iA0 = new IAtom(); iA0->entities = iE.entities;

	char *amino;
	int i = 1, nSE;

	nSE = iE.NSE();

	count = 11;
	Summary *out = new Summary[count];

	Summary *C0_N = new Summary("C-N", nSE, 1.329, 0.014, BONDLENGTH);
	Summary *C0_NPro = new Summary("C-N (PRO)",nSE,1.341,0.016,BONDLENGTH);
	Summary *C_O = new Summary("C-O", nSE, 1.231, 0.02, BONDLENGTH);
	Summary *CA_C = new Summary("CA-C", nSE, 1.525, 0.021, BONDLENGTH);
	Summary *CA_CGly = new Summary("CA-C (GLY)",nSE,1.516,0.018,BONDLENGTH);
	Summary *CA_CB = new Summary("CA-CB", nSE, 1.53, 0.02, BONDLENGTH);
	Summary *CA_CBAla=new Summary("CA-CB (ALA)",nSE,1.521,0.033,BONDLENGTH);
	Summary *CA_CBITV =
		new Summary("CA-CB (I,T,V)", nSE, 1.54, 0.027, BONDLENGTH);
	Summary *N_CA = new Summary("N-CA", nSE, 1.458, 0.019, BONDLENGTH);
	Summary *N_CAGly = new Summary("N-CA (GLY)",nSE,1.451,0.016,BONDLENGTH);
	Summary *N_CAPro = new Summary("N-CA (PRO)",nSE,1.466,0.015,BONDLENGTH);

	for (ISubentity iS(iE); iS; ++iS, ++i) {
		N  = iS.findAtom (" N  ");
		C  = iS.findAtom (" C  ");
		O  = iS.findAtom (" O  ");
		CA = iS.findAtom (" CA ");
		CB = iS.findAtom (" CB ");

		if (!iS.monType()) { C0 = (*iA0); continue; }

		amino = iS.code3();

		if (i != 1) {
			if (strcmp(amino, "PRO")) {
			    C0_N->add(iS, Calculate::Distance(C0, N));
			} else {
			    C0_NPro->add(iS, Calculate::Distance(C0, N));
			}
		}

		C_O->add(iS, Calculate::Distance(C, O));

		if (strcmp(amino, "GLY")) {
			CA_C->add(iS, Calculate::Distance(CA, C));
			if (strcmp(amino, "PRO")) {
			    N_CA->add(iS, Calculate::Distance(N, CA));
			} else {
			    N_CAPro->add(iS, Calculate::Distance(N, CA));
			}
		} else {
			N_CAGly->add(iS, Calculate::Distance(N, CA));
			CA_CGly->add(iS, Calculate::Distance(CA, C));
		}

		if (strcmp(amino, "ALA") && strcmp(amino, "ILE") &&
		    strcmp(amino, "THR") && strcmp(amino, "VAL")) {
			CA_CB->add(iS, Calculate::Distance(CA, CB));
		}

		if (!strcmp(amino, "ALA")) {
			CA_CBAla->add(iS, Calculate::Distance(CA, CB));
		}

		if (!strcmp(amino, "ILE") || !strcmp(amino, "THR") ||
		    !strcmp(amino, "VAL")) {
			CA_CBITV->add(iS, Calculate::Distance(CA, CB));
		}

		C0 = C;
	}

	out[0]  = (*C0_N);
	out[1]  = (*C0_NPro);
	out[2]  = (*C_O);
	out[3]  = (*CA_C);
	out[4]  = (*CA_CGly);
	out[5]  = (*CA_CB);
	out[6]  = (*CA_CBAla);
	out[7]  = (*CA_CBITV);
	out[8]  = (*N_CA);
	out[9]  = (*N_CAGly);
	out[10] = (*N_CAPro);

	delete iA0;

	return out;
}

void Stats::printResGeo(ISubentity &iS) {

	if (!iS) return;

	int count;
	double *values = getResGeo(iS, count);

	char **names = new char *[count];
	names[0]  = "Phi";
	names[1]  = "Psi";
	names[2]  = "Omega";
	names[3]  = "Chi1 g(-)";
	names[4]  = "Chi1 trans";
	names[5]  = "Chi1 g(+)";
	names[6]  = "N-CA-C";
	names[7]  = "N-CA-CB";
	names[8]  = "CA-C-O";
	names[9]  = "CA-C-N";
	names[10]  = "CB-CA-C";
	names[11]  = "O-C-N";
	names[12]  = "C-N-CA";
	names[13] = "C-N";
	names[14] = "C-O";
	names[15] = "CA-C";
	names[16] = "CA-CB";
	names[17] = "N-CA";

	printf("%s %d\n", iS.code3(), iS());
	for (int n = 0; n < count; n++) {
		if (values[n] == MAXVALUE) continue;
		printf("%-7s:  %7.2f\n", names[n], values[n]);
	}
	printf("\n");

	delete [] values;
}

double *Stats::getResGeo(ISubentity &iS, int &count) {
	count = 18;

	double *values = new double[count];
	char *amino = iS.code3();
	double chi1;

	IAtom N, C, O, CA, CB, CA1, N1, C0, CG;

	N  = iS.findAtom(" N  ");
	C  = iS.findAtom(" C  ");
	O  = iS.findAtom(" O  ");
	CA = iS.findAtom(" CA ");
	CB = iS.findAtom(" CB ");

	if (!strcmp(amino, "VAL") ||
	    !strcmp(amino, "ILE"))
		CG = iS.findAtom(" CG1");
	else if (!strcmp(amino, "THR"))
		CG = iS.findAtom(" OG1");
	else if (!strcmp(amino, "SER"))
		CG = iS.findAtom(" OG ");
	else CG = iS.findAtom(" CG ");

	ISubentity iS0 = iS.neighbor(-1);
	ISubentity iS1 = iS.neighbor(1);

	if (iS0() >= 0) {
		C0 = iS0.findAtom(" C  ");

		values[0] = Torsion::TorsionAngle(C0, N, CA, C);
		values[12] = Torsion::CommonAngle(C0, N, CA);
		values[13] = Calculate::Distance(C0, N);
	} else {
		values[0] = MAXVALUE;
		values[12] = MAXVALUE;
		values[13] = MAXVALUE;
	}

	if (iS1() >= 0) {
		N1  = iS1.findAtom(" N  ");
		CA1 = iS1.findAtom(" CA ");

		values[1] = Torsion::TorsionAngle(N, CA, C, N1);
		values[2] = Torsion::TorsionAngle(CA, C, N1, CA1);
		values[9] = Torsion::CommonAngle(CA, C, N1);
		values[11] = Torsion::CommonAngle(O, C, N1);
	} else {
		values[1] = MAXVALUE;
		values[2] = MAXVALUE;
		values[9] = MAXVALUE;
		values[11] = MAXVALUE;
	}

	values[3] = values[4] = values[5] = MAXVALUE;

	chi1	  = Torsion::TorsionAngle(CG, CB, CA, N);

	if (chi1 < -120.0) chi1 += 360;

	if (chi1 > 120.0) values[4] = chi1;
	else if (chi1 >= 0.0) values[3] = chi1;
	else values[5] = chi1;

	values[6] = Torsion::CommonAngle(N, CA, C);
	values[7] = Torsion::CommonAngle(N, CA, CB);
	values[8] = Torsion::CommonAngle(CA, C, O);
	values[10] = Torsion::CommonAngle(CB, CA, C);

	values[14] = Calculate::Distance(C, O);
	values[15] = Calculate::Distance(CA, C);
	values[16] = Calculate::Distance(CA, CB);
	values[17] = Calculate::Distance(N, CA);

	return values;
}

ValsNStddev *Stats::getSmallResVals(ISubentity &iS, int &count) {
	return getSmallResVals(iS.code3(), iS.ks(), count);
}

ValsNStddev *Stats::getSmallResVals(char *amino, char ks, int &count) {
	// small molecule data taken from procheck (Engh & Huber, 1991)
	count = 18;
	ValsNStddev *vals = new ValsNStddev[count];
	
	for (int n = 0; n < count; n++) {
		vals[n].value = MAXVALUE;
		vals[n].stddev = 0.0;
	}

	if (strcmp(amino, "ALA") && strcmp(amino, "VAL") &&
	    strcmp(amino, "LEU") && strcmp(amino, "ILE") &&
	    strcmp(amino, "CYS") && strcmp(amino, "MET") &&
	    strcmp(amino, "PRO") && strcmp(amino, "PHE") &&
	    strcmp(amino, "TYR") && strcmp(amino, "TRP") &&
	    strcmp(amino, "ASP") && strcmp(amino, "ASN") &&
	    strcmp(amino, "GLU") && strcmp(amino, "GLN") &&
	    strcmp(amino, "HIS") && strcmp(amino, "SER") &&
	    strcmp(amino, "THR") && strcmp(amino, "ARG") &&
	    strcmp(amino, "LYS") && strcmp(amino, "GLY") &&
	    strcmp(amino, "ASX") && strcmp(amino, "GLX")) return vals;

	if (ks == 'H') {
		vals[0].value = -65.3; vals[0].stddev = 11.9;
		vals[1].value = -39.4; vals[1].stddev = 11.3;
	}

	vals[2].value = 180.0; vals[2].stddev = 5.8;
	vals[14].value = 1.231; vals[14].stddev = 0.020;

	if (!strcmp(amino, "PRO")) {
		vals[0].value = -65.4; vals[0].stddev = 11.2;
		vals[13].value = 1.341; vals[13].stddev = 0.016;
		vals[17].value = 1.466; vals[17].stddev = 0.015;
		vals[12].value = 122.60; vals[12].stddev = 5.00;
		vals[9].value = 116.90; vals[9].stddev = 1.50;
		vals[6].value = 111.80; vals[6].stddev = 2.50;
		vals[7].value = 103.00; vals[7].stddev = 1.10;
		vals[11].value = 122.00; vals[11].stddev = 1.40;
	} else {
		vals[13].value = 1.329; vals[13].stddev = 0.014;
		vals[11].value = 123.00; vals[11].stddev = 1.60;
	}

	if (!strcmp(amino, "GLY")) {
		vals[15].value = 1.516; vals[15].stddev = 0.018;
		vals[17].value = 1.451; vals[17].stddev = 0.016;
		vals[12].value = 120.60; vals[12].stddev = 1.70;
		vals[9].value = 116.40; vals[9].stddev = 2.10;
		vals[8].value = 120.80; vals[8].stddev = 2.10;
		vals[6].value = 112.50; vals[6].stddev = 2.90;
	} else {
		vals[15].value = 1.525; vals[15].stddev = 0.021;
		vals[8].value = 120.80; vals[8].stddev = 1.70;
		if (strcmp(amino, "PRO")) {
			vals[17].value = 1.458; vals[17].stddev = 0.019;
			vals[12].value = 121.70; vals[12].stddev = 1.80;
			vals[9].value = 116.20; vals[9].stddev = 2.00;
			vals[6].value = 111.20; vals[6].stddev = 2.80;
		}
	}

	if (!strcmp(amino, "ALA")) {
		vals[16].value = 1.521; vals[16].stddev = 0.033;
		vals[10].value = 110.50; vals[10].stddev = 1.50;
		vals[7].value = 110.40; vals[7].stddev = 1.50;
	} else {
		if (!strcmp(amino, "ILE") || !strcmp(amino, "THR") ||
		    !strcmp(amino, "VAL")) {
			vals[16].value = 1.540; vals[16].stddev = 0.027;
			vals[10].value = 109.10; vals[10].stddev = 2.20;
			vals[7].value = 111.50; vals[7].stddev = 1.70;
		} else {
			if (strcmp(amino, "GLY")) {
			    vals[16].value = 1.530; vals[16].stddev = 0.020;
			    vals[10].value = 110.10; vals[10].stddev = 1.90;
			    vals[7].value = 110.50; vals[7].stddev = 1.70;
			}
		}
	}

        if (strcmp(amino, "GLY") && strcmp(amino, "ALA") &&
            strcmp(amino, "ASX") && strcmp(amino, "PRO") &&
            strcmp(amino, "GLX")) {
                vals[3].value = 64.1; vals[3].stddev = 15.7;
                vals[4].value = 183.6; vals[4].stddev = 16.8;
                vals[5].value = -66.7; vals[5].stddev = 15.0;
        }

	return vals;
}

double Stats::getFDS(ISubentity &iS) {
	int num;
	ValsNStddev *norm = getSmallResVals(iS, num);
	double *vals = getResGeo(iS, num);

	int nMeans = 0;
	double val_, diff, fds = 0.0;

	for (int n = 0; n < num; n++) {
		if (norm[n].stddev && vals[n] != MAXVALUE) {
			val_ = (vals[n] < 0.0 && norm[n].value > 0.0
						? -1.0*vals[n] : vals[n]);
			diff = val_ - norm[n].value;
			diff = (diff > 0.0 ? diff : -1.0*diff);
			fds += diff/norm[n].stddev;
			nMeans++;
		}
	}

	delete [] norm;
	delete [] vals;

	if (nMeans) return fds/nMeans;
	return MAXVALUE;
}

double Stats::getFDS(IEntity &iE) {
	double fds, total = 0.0;
	int num = 0;

	for (ISubentity iS(iE); iS; ++iS) {
		fds = getFDS(iS);
		if (fds == MAXVALUE) continue;
		total += fds;
		num++;
	}

	if (num) return total/num;
	return MAXVALUE;
}

Summary Stats::getAllFDS(IEntity &iE) {

	char *name = new char[100];
	strcpy(name, iE.name());
	strcat(name, " FDS");

	Summary *all = new Summary(name, iE.NSE(), MAXVALUE, 0.0, DOUBLES);
	delete [] name;
	
	for (ISubentity iS(iE); iS; ++iS) all->add(iS, getFDS(iS));

	return (*all);
}
