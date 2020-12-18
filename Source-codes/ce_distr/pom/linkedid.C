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
/* $Revision: 2.2 $ $Date: 2000/03/30 23:17:03 $ */
#include "linkedid.h"
#include "miscutil.h"

DB *LinkedId::db = NULL;

Property LinkedId::id_com;
Property LinkedId::obs_com;
Property LinkedId::spr_com;
Property LinkedId::relDat;

int LinkedId::open = 0;

void LinkedId::init() {
		id_com.open("id.com", 1);
		obs_com.open("obs.com", 1);
		spr_com.open("spr.com", 1);
		relDat.open("date_int.com", 1);

		open = 1;
}

LinkedId::~LinkedId() {
	id_com.close();
	obs_com.close();
	spr_com.close();
	open = 0;
}

LinkedId::LinkedId(int i) {
	index = i;
	if (!open) init();
}

int LinkedId::getChain(int idx, int **chain) {

	if (!open) init();
	
	int *checked = new int[id_com.nItems];
	int *ch = NULL, nChain = 0;
	int sorted = 0, tmp;

	for (int n = 0; n < id_com.nItems; n++) checked[n] = 0;

        addChain(idx, &ch, &nChain, &checked);

	delete [] checked;

	(*chain) = ch;

	while (!sorted) {
	    sorted = 1;
	    for (int n = 0; n < nChain - 1; n++) {
		if (*relDat.item4((*chain)[n]) > *relDat.item4((*chain)[n+1])) {
			tmp = (*chain)[n];
			(*chain)[n] = (*chain)[n+1];
			(*chain)[n+1] = tmp;
			sorted = 0;
		}

		if (*relDat.item4((*chain)[n]) == *relDat.item4((*chain)[n+1])) {
			char *id1 = id_com.item1((*chain)[n], 1);
			char *id2 = id_com.item1((*chain)[n+1], 1);

			if (id1[0] > id2[0]) {
				tmp = (*chain)[n];
				(*chain)[n] = (*chain)[n+1];
				(*chain)[n+1] = tmp;
			}

			delete [] id1;
			delete [] id2;
		}
	    }
	}

	return nChain;
}

int LinkedId::getChain(int **chain) {
	
	int *checked = new int[id_com.nItems];
	for (int n = 0; n < id_com.nItems; n++) checked[n] = 0;

	int *ch = NULL, nChain = 0;
        addChain(index, &ch, &nChain, &checked);

	delete [] checked;

	int sorted = 0, tmp;

	(*chain) = ch;

	while (!sorted) {
	    sorted = 1;
	    for (int n = 0; n < nChain - 1; n++) {
		if (*relDat.item4((*chain)[n]) > *relDat.item4((*chain)[n+1])) {
			tmp = (*chain)[n];
			(*chain)[n] = (*chain)[n+1];
			(*chain)[n+1] = tmp;
			sorted = 0;
		}

		if (*relDat.item4((*chain)[n]) == *relDat.item4((*chain)[n+1])) {
			char *id1 = id_com.item1((*chain)[n], 1);
			char *id2 = id_com.item1((*chain)[n+1], 1);

			if (id1[0] > id2[0]) {
				tmp = (*chain)[n];
				(*chain)[n] = (*chain)[n+1];
				(*chain)[n+1] = tmp;
			}

			delete [] id1;
			delete [] id2;
		}
	    }
	}

	return nChain;
}

void LinkedId::addChain(int i, int **chain, int *nChain, int **check) {
        int n = 0;
        if (!(*chain) || !(*check)[(*chain)[n]]) {
                LinkedId::addToArray(chain, i, (*nChain)); (*nChain)++;
        }
 
        while(n < (*nChain)) {
                if ((*check)[(*chain)[n]]) {
                        n++; continue;
                } else (*check)[(*chain)[n]] = 1;
 
                int *hist = NULL;
                int add = getObsNSpr((*chain)[n], &hist);
                for (int j = 0; j < add; j++) {
                        if ((*check)[hist[j]]) continue;
                        LinkedId::addToArray(chain, hist[j], (*nChain));
                        (*nChain)++;
                }

		if (hist) delete [] hist;
 
                n++; if (n >= (*nChain)) break;
 
                addChain((*chain)[n], chain, nChain, check);
        }
}

int LinkedId::getObsNSpr(int idx, int **history) {
	int nObsIds, *obsIds, nSprIds, *sprIds, *temp, nTemp, n;

        nObsIds = obs_com.getItemSize(idx);
        obsIds = obs_com.item4(idx, 1);
 
        nSprIds = spr_com.getItemSize(idx);
        sprIds = spr_com.item4(idx, 1);
 
        temp = NULL; nTemp = 0; n;
 
        for (n = 0; n < nObsIds; n++) {
                if (obsIds[n] != -1)
			LinkedId::addToArray(&temp, obsIds[n], nTemp++);
        }
 
        for (n = 0; n < nSprIds; n++) {
                if (sprIds[n] != -1)
			LinkedId::addToArray(&temp, sprIds[n], nTemp++);
        }
 
	delete [] sprIds;
	delete [] obsIds;

        (*history) = temp;
        return nTemp;
}

int LinkedId::getObsNSpr(int **history) {
        int nObsIds = obs_com.getItemSize(index);
        int *obsIds = obs_com.item4(index, 1);
 
        int nSprIds = spr_com.getItemSize(index);
        int *sprIds = spr_com.item4(index, 1);
 
        int *temp = NULL, nTemp = 0, n;
 
        for (n = 0; n < nObsIds; n++) {
                if (obsIds[n] != -1) LinkedId::addToArray(&temp, obsIds[n], nTemp++);
        }
 
        for (n = 0; n < nSprIds; n++) {
                if (sprIds[n] != -1) LinkedId::addToArray(&temp, sprIds[n], nTemp++);
        }
 
	if (nSprIds) delete [] sprIds;
	if (nObsIds) delete [] obsIds;

        (*history) = temp;
        return nTemp;
}

void LinkedId::replacers(int **array, int nArray) {
	int *obsIds;

	if (!open) init();

	for (int n = 0; n < nArray; n++) {
		obsIds = obs_com.item4(*(*array + n), 0);
		*(*array + n) = obsIds[0];
	}

	bubbleSort(array, nArray);
}

void LinkedId::addToArray(int ** array, int value, int Narray) {
	int *array_ = new int[Narray + 1];
	for(int i = 0; i < Narray; i++) array_[i] = *(*array + i);
	array_[Narray] = value;
	if (Narray) delete [] (*array);
	(*array) = array_;
}

DB *IVersion::db = NULL;

IVersion::IVersion() {
	Property id_com("id.com", 0);
	nIds = id_com.nItems;
	dbIds = new int[nIds];
	int i;
	for (i = 0; i < nIds; i++) dbIds[i] = 1;

	LinkedId::db = db;

	nVers = LinkedId::getChain(0, &versions);
	for (i = 0; i < nVers; i++) dbIds[versions[i]] = 0;

	iId = 0; iversions = 0;
}

IVersion::~IVersion() {
	delete [] dbIds;
	iId = iversions = -1;

	db = NULL;
}

void IVersion::nextSet() {
	delete [] versions;

	while (iId < nIds && !dbIds[iId]) iId++;

	if (iId >= nIds) return;

	nVers = LinkedId::getChain(iId, &versions);
	for (int i = 0; i < nVers; i++) dbIds[versions[i]] = 0;

	iversions = 0;
}

int IVersion::getChain(int **array) {
	if ((*array)) delete [] (*array);
	(*array) = new int[nVers];
	for (int i = 0; i < nVers; i++) *(*array + i) = versions[i];
	return nVers;
}
