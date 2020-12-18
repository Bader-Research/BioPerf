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
#include "derive.h"
#include "miscutil.h"

int Derive::numHOH(int index) {

	static Property i_enc("i_enc.com", 1);
	static Property n_enp("n_enp.com", 1);
	static Property n_enc("n_enc.com", 1);
	static Property name_enc("name.enc", 0);

        int num, pos, n, count, i, x, j; char **hets;
        char *name = new char[100]; name[0] = '\0';
 
        pos = *i_enc.item4(index) + *n_enp.item2(index);
        num = *n_enc.item2(index) - *n_enp.item2(index);
 
        hets = new char*[num];
 
        for (n = 0; n < num; n++) {
                hets[n] = name_enc.item1(pos+n, 1);
        }
 
        count = 0;
 
        for (n = 0; n < num; n++) {
                if (!hets[n]) continue;
 
                for (x = 0; x < strlen(hets[n]); x++) {
                        if (*(hets[n]+x) == ':') {
                                strcpy(name, (hets[n]+x+1));
                                break;
                        }
                }
 
                delete [] hets[n]; hets[n] = NULL;
                if (!name || strcmp(name, "HOH")) continue;
 
                count++;
                for (j = n+1; j < num; j++) {
                    if (!hets[j]) continue;
                    if (!strncmp(name, (hets[j]+x+1), strlen(name))) {
                        count++;
                        delete [] hets[j]; hets[j] = NULL;
                    }
                }
 
                name[0] = '\0';
        }
 
        delete [] hets;
        delete [] name;
        return count;
}

int Derive::numAtoms(int index) {

	int atoms = 0, num, first;

	static Property n_enc("n_enc.com", 1);
	static Property i_enc("i_enc.com", 1);
	static Property nxyz("n_xyz.enc", 1);

	num = *n_enc.item2(index);
	first = *i_enc.item4(index);

	for (int i = first; i < first+num; i++) atoms += *nxyz.item4(i);

	return atoms;
}

int Derive::numResidues(int index) {

	int residues = 0, num, first;

	static Property n_enp("n_enp.com", 1);
        static Property i_enp("i_enp.com", 1);
        static Property n_se("n_se.enp", 1);

        num = *n_enp.item2(index);
        first = *i_enp.item4(index);

        for (int i = first; i < first+num; i++)
		residues += (int) *n_se.item2(i);

        return residues;
}

PDBVersionType Derive::getPDBVersionType(int index) {
	static Property obs("obs_dat.com", 1);
	static Property spr("spr_dat.com", 1);

	if (*obs.item4(index) == 0 && *spr.item4(index) == 0) return CURRENT;
	if (*obs.item4(index) != 0) return OBSLTE;
	if (*obs.item4(index) == 0 && *spr.item4(index) != 0) return SPRSDE;
	return UNK;
}

void Derive::matchType(PDBVersionType vType,
		       int *ids,
		       int nIds,
		       int **result,
		       int *nResult) {

	*nResult = 0;
	int *tmp = new int[nIds], nTmp = 0, i;

	for (i = 0; i < nIds; i++) {
		if (Derive::getPDBVersionType(ids[i]) == vType) {
			tmp[nTmp++] = ids[i];
		}
	}

	*result = new int[nTmp];
	for (i = 0; i < nTmp; i ++) *(*result + i) = tmp[i];

	delete [] tmp;
	*nResult = nTmp;
}

char *Derive::getTitle(char *chainId) {
	static Property name_enp("name.enp");
	char chain[7]; chain[6] = '\0';
	int i;

	strncpy(chain, chainId, 6);

	/* upper case only */
	for (i = 0; i < 6; i++)
		if (chain[i] > 96 && chain[i] < 123) chain[i] -= 32;

	i = name_enp.find(chain);

	return getTitle(i);
}

char *Derive::getTitle(int index) {
	static Property name_enp("name.enp");
	static Property i_com_enp("i_com.enp");
	static Property compnd("compnd.com");
	char *com = NULL, *title, *semPtr, *chnPtr;
	char chain[7]; chain[6] = '\0';

	title = NULL;

	if (index >= 0) {

	    int iC = *i_com_enp.item4(index);
	    addText(&com, compnd.item1(iC));
	    strcpy(chain, name_enp.item1(index));

	    /* see if this entry complies with PDB format 2.n and higher */
	    if (!strchr(com, ':') && !strchr(com, ';')) {

		/* it doesn't */
		oldFormat:
		    addText(&com, " - CHAIN ");
		    addText(&com, &chain[5]);
		    title = strdup(com);

	    } else {

		/* it does */
		for (int p = 0; p < strlen(com); p++) {

		    /* find a reference to chains */
		    while (strncmp(&com[p], "CHAIN: ", 7) && com[++p]) {}

		    if (!com[p]) break;

		    semPtr = strchr(&com[p], ';');
		    chnPtr = strchr(&com[p+7], chain[5]);

		    if ((chnPtr && semPtr && chnPtr < semPtr) ||
			!strncmp(&com[p+7], "NULL", 4) ||
			(chnPtr && !semPtr)) {

			/* find a reference to molecule (backwards) */
		        while (strncmp(&com[p], "MOLECULE: ", 10) && p--) {}
			if (strncmp(&com[p], "MOLECULE: ", 10)) goto oldFormat;
			semPtr = strchr(&com[p], ';');

		        if (p) {
			    *semPtr = '\0';
			    title = strdup(&com[p+10]);
			}

			break;
		    }
		}

		if (!title) goto oldFormat;
	     }
	 }

	 delete [] com;
	 return title;
}

/* this subroutine will return the molecule type of a chain as follows:

	P	=> protein
	D	=> DNA
	R	=> RNA
	N	=> nucleic acid (any)
	X	=> protein/nucleic acid complex
	C	=> carbohydrate
	U	=> unknown

*/

char Derive::getType(int index) {
	static Property i_enp("i_enp.com");
	static Property n_enp("n_enp.com");
	static Property type_enp("type.enp");

	int iEnp = *i_enp.item4(index);
	char type = 'U';
	int N, P, D, R, C, U;
	N = P = D = R = C = U = 0;

	for (int n = 0; n < *n_enp.item2(index); n++) {
	    switch(*type_enp.item1(n+iEnp)) {
		case 'N': N = 1; break;
		case 'P': P = 1; break;
		case 'D': D = 1; break;
		case 'R': R = 1; break;
		case 'C': C = 1; break;
		case 'U': U = 1; break;
	    }
	}

	if (N && !(P+C+U))			type = 'N';
	else if (P && !(N+D+R+C+U))		type = 'P';
	else if (D && !(N+P+R+C+U))		type = 'D';
	else if (R && !(N+P+D+C+U))		type = 'R';
	else if (C && !(N+P+D+R+U))		type = 'C';
	else if (U && !(N+P+D+R+C))		type = 'U';
	else if (D*R && !(P+C+U))		type = 'N';
	else if (P)				type = 'X';

	return type;
}
