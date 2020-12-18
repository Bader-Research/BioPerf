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
#include "ipdb.h"
#include "ipdb_df.h"
#include "pom.h"
#include "jdate.h"
#include "bnastats.h"
#include "derive.h"
#include "miscutil.h"
#include "linkedid.h"
#include "pdbutil.h"
#include "config.h"

#include <stdio.h>
#include <sys/types.h>
#include "sys/times.h"
#include "sys/param.h"
///////////////////////////////////////////////////////////////////////////
Tlpdb * DF::tlpdb=NULL;
////////////////////////////////////////////////////////////////////
void _addToArray(int4 ** array, int4 value, int Narray);

void _addToArray(int4 ** array, int4 value, int Narray) {
	int4 * array_;
	array_=new int4[Narray+1];
	for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
	array_[Narray]=value;
	if(Narray) delete [] (*array);
	(*array)=array_;
	
};
//////////////////////////////////////////////////////////////////
template <class A> void DeleteArray(A ** array, int Narray) {
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};
//////////////////////////////////////////////////////////////////
extern AminoacidProperty properties[];
///////////////////////////////////////////////////////////////////////////
int DF::assignDF(int nCom1, int nCom2, int log_mode) {
	 assignDF(DBPATH, nCom1, nCom2, 1, NUM_DF, log_mode);
}
///////////////////////////////////////////////////////////////////////////
int DF::assignDF(char *db, int nCom1, int nCom2, int log_mode) {
	assignDF(db, nCom1, nCom2, 1, NUM_DF, log_mode);
}
///////////////////////////////////////////////////////////////////////////
void DF::printFeatures() {

    int num = 1;

    printf("\n");

#ifndef MINIMAL

    printf("%d - generate enp filt bfac and c_a (for QuickPDB)\n", num++);
    printf("%d - create ss propens\n", num++);

#endif

    printf("%d - generate obs.com, spr.com\n", num++);
    printf("%d - calculate FDS collection and property objects\n", num++);
    printf("%d - calculate chi1 property objects\n", num++);

#ifndef MINIMAL

    printf("%d - create df mm, fast prop, pat5x3 for env3 & stat5\n", num++);
    printf("    and propens for env3\n");
    printf("%d - create moose-report properties\n", num++);

#endif
    printf("%d - calculate Collection objects: res.col & misc.col\n", num++);
    printf("%d - generate status.com\n", num++);

    printf("\n");

}
///////////////////////////////////////////////////////////////////////////
int DF::assignDF(char *db_path,
		 int nCom1,
		 int nCom2,
		 int df1,
		 int df2,
		 int log_mode) {

  if (!log_mode) printf("calculating DF %d to %d for compounds %d to %d\n",
                        df1, df2, nCom1, nCom2);

#ifdef MINIMAL

	df1 += 2;
	df2 += 2;

	if (df2 > 7) df2 == 7;

#endif

  DB *db=new DB;
  db->setPath(db_path);
  LinkedId::db = db;
  if(!db->testFile("code3.mon")) {
	printf("No database available in the directory %s\n", db_path);
  }
 
  Monomers monomers;

  Property id("id.com", 1);
  if (nCom2 > id.nItems) nCom2 = id.nItems;
  id.close();

  for(int iDF = df1; iDF <= df2; iDF++) {
	switch (iDF) {
		case 1:
			bfac(nCom1, nCom2, log_mode, df2);
			break;

		case 2:
			ss_propens(nCom1, nCom2, log_mode, df2);
			break;
		
		case 3:
			obsNspr(nCom1, nCom2, log_mode, df2);
			break;
		
		case 4:
			collFDS(nCom1, nCom2, log_mode, df2);
			break;

		case 5:
			calcChi1(nCom1, nCom2, log_mode, df2);
			break;

		case 6:
			proteinFeatures(db_path, nCom1, nCom2, log_mode, df2);
			break;

		case 7:
			mooseReport(nCom1, nCom2, log_mode, df2);
			break;

		case 8:
			collect(nCom1, nCom2, log_mode, df2);
			break;

		case 9:
			status_com(log_mode);
			break;

		default:
			break;
	}
  }

}
///////////////////////////////////////////////////////////////////////////
int DF::indexIds(Property *id_com, char *ids, int4 **index) {
	int nIds = strlen(ids)/5;

	if (!nIds) return nIds;

	char **idArray = NULL, *iId;
	(*index) = new int4[nIds];

	char *ids_ = new char[nIds*5+1];
	strcpy(ids_, ids);

	int n;

	for (n = 0; n < nIds; n++) (*index)[n] = -1;

	BufferToStringArray(nIds, &ids_, &idArray);

	for (n = 0; n < nIds; n++) {
		for (int ii = 0; ii < (*id_com).nItems; ii++) {
			iId = (*id_com).item1(ii, 0);
			if (!strcmp(iId, idArray[n])) {
				(*index)[n] = ii;
				break;
			}
		}
	}

	DeleteArray(&idArray, nIds);

	return nIds;
}
				

void DF::bfac(int nCom1, int nCom2, int log_mode, int df2) {

  if(!P_BFAC_FLT_ENP && !P_C_A_ENP && !P_N_ENP_COM) return;

  Property bfac_flt_enp;
  Property c_a_enp;

  Property n_enp_com;

  if(P_BFAC_FLT_ENP) bfac_flt_enp.open("bfac_flt.enp", 0);
  if(P_C_A_ENP) c_a_enp.open("c_a.enp", 0);
  if(P_N_ENP_COM) n_enp_com.open("n_enp.com", 0);


  int iProp, i;
  if(!nCom1)
    {
      if(P_BFAC_FLT_ENP) bfac_flt_enp.clear();
      if(P_C_A_ENP) c_a_enp.clear();
    }

  Entities entities1; 

  int nEnt, iEnt; int2 *nProp; flt4 **prop;
  int nAv; double pAv; flt4 **ca;
  IAtom ca_atom;

  for(int iC=nCom1; iC<nCom2; iC++)
    {
      if(!log_mode)
	printf("DF 1/%d, COM %d (%d/%d)\n", df2, iC,
					    iC-nCom1+1,
					    (nCom2-nCom1));
      entities1.addCom(iC, 0);
      if(P_N_ENP_COM) nEnt=*n_enp_com.item2(iC);

      if(nEnt > 0) {

	nProp=new int2 [nEnt];
	prop=new flt4* [nEnt];
	ca=new flt4* [nEnt];
	iEnt=0;
	for(IEntity iE1(&entities1); iEnt<nEnt; ++iE1, iEnt++)
	  {
	    nProp[iEnt]=iE1.NSE();
	    prop[iEnt]=new flt4 [iE1.NSE()];
	    ca[iEnt]=new flt4 [iE1.NSE()*3];
	    for(i=0; i<iE1.NSE(); i++) 
	      {
		*(prop[iEnt]+i)=-1.0;
		*(ca[iEnt]+i*3)=2e10;
		*(ca[iEnt]+i*3+1)=2e10;
		*(ca[iEnt]+i*3+2)=2e10;
	      }
	    
	    i=0;
	    for(ISubentity iS1(iE1); iS1; ++iS1, i++)
	      {
		ca_atom=iS1.findAtom(" CA ");
		if(ca_atom)
		  {
		    *(ca[iEnt]+i*3)=ca_atom.x();
		    *(ca[iEnt]+i*3+1)=ca_atom.y();
		    *(ca[iEnt]+i*3+2)=ca_atom.z();
		  }
	      }
	    
	    i=0;
	    for(ISubentity iS11(iE1); iS11; ++iS11)
	      {
		if(!iS11) continue;
		
		for(IAtom iA(iS11); iA; ++iA, i++)
		  {
		    if(iA.bfac()<=0.0) goto tmf_no;
		  }
	      }
	    
	    i=0;
	    for(ISubentity iS12(iE1); iS12; ++iS12, i++)
	      {
		if(!iS12) continue;
		if(!iS12.monType()) continue;
		
		pAv=0.0; nAv=0;
		for(IAtom iA(iS12); iA; ++iA)
		  {
		    pAv+=iA.bfac(); nAv++;
		  }
		
		if(!nAv) continue;
		*(prop[iEnt]+i)=pAv/nAv;
	      }
	  tmf_no: 
	    if(P_BFAC_FLT_ENP) bfac_flt_enp.addItem(prop[iEnt], iE1.NSE());
	    if(P_C_A_ENP) c_a_enp.addItem(ca[iEnt], iE1.NSE()*3);
	  }
	
	delete [] nProp;
	DeleteArray(&prop, nEnt);
	DeleteArray(&ca, nEnt);
      }
    }
}

void DF::ss_propens(int nCom1, int nCom2, int log_mode, int df2) {
  /*
  int4 tEnt1; int2 tEnt2;
  int1 *seq, seq_; char *ks;

  if(!MOOSE_PATTERNS_GROUP) return;

  int i, j, iCom, iEnt, nEnt, aa_code_; int4 nSeq; 
  
  int4 *aa_freq, nAA_freq;

  Property i_enp_com("i_enp.com");
  Property ks_pp3("ks.pp3");
  Property iseq_enp("iseq.enp", 1);
  Property k_s_enp("k_s.enp", 1);
  Property nse_enp("n_se.enp");
  
  tEnt1=*i_enp_com.item4(nCom1);
  tEnt2=*i_enp_com.item4(nCom2);

  if(!nCom1)
    {
      aa_freq=new int4 [60];
      for(i=0; i<20; i++) for(j=0; j<3; j++) aa_freq[i*3+j]=0;
    }
  else
    {
      aa_freq=ks_pp3.item4(0, 1);
    }

  for(i=tEnt1; i<tEnt2; i++)
    {
      if(!log_mode)
	printf("DF 2/%d, COM %d (%d/%d)\n", df2, i,
					    i-nCom1+1,
					    (nCom2-nCom1));

      nSeq = *nse_enp.item2(i);
      seq=iseq_enp.item1(i);
      ks=k_s_enp.item1(i);
      
      for(j=0; j<nSeq; j++)
	{
	  aa_code_=seq[j];
	  if (aa_code_ < 0 || aa_code_ > 19) continue;
	  
	  if(ks[j] == 'H' || ks[j] == 'G' || ks[j] == 'I') 
	    {
	      aa_freq[aa_code_*3]++; goto ks_found;
	    }
	  if(ks[j]=='E') 
	    {
	      aa_freq[aa_code_*3+1]++; goto ks_found;
	    }
	  aa_freq[aa_code_*3+2]++;
	  
	ks_found: ;
	}
    }

  ks_pp3.clear();
  ks_pp3.addItem(aa_freq, 60);

  delete [] aa_freq;
  */
}


void DF::obsNspr(int nCom1, int nCom2, int log_mode, int df2) {
    int4 *index, nIds = 0;
    char *ids;

    if(!P_OBS_COM && !P_SPR_COM && !P_OBS_IDS_COM && !P_SPR_IDS_COM &&
       !P_OBS_DAT_COM && !P_SPR_DAT_COM) return;

    Property id_com("id.com", 1);

    Property obs_com;
    Property spr_com;
    Property obs_ids_com;
    Property spr_ids_com;

    if(P_OBS_COM) obs_com.open("obs.com", 1);
    if(P_SPR_COM) spr_com.open("spr.com", 1);
    if(P_OBS_IDS_COM) obs_ids_com.open("obs_ids.com", 1);
    if(P_SPR_IDS_COM) spr_ids_com.open("spr_ids.com", 1);

    int4 *justOne = new int4[1]; justOne[0] = -1;

    int4 iC;

    for(iC = nCom1; iC < nCom2; iC++) {

	if (!log_mode)
		printf("DF 3/%d, COM %d (%d/%d)\n", df2, iC,
						    iC-nCom1+1,
						    (nCom2-nCom1));

	/* why reset a possibly fixed db????
	if (!nCom1) {
		if(P_OBS_COM) obs_com.setItem(iC, justOne, 1);
		if(P_SPR_COM) spr_com.setItem(iC, justOne, 1);
	}
	*/

	ids = obs_ids_com.item1(iC, 0);
	nIds = indexIds(&id_com, ids, &index);

	if (nIds) {
		if(P_OBS_COM) obs_com.setItem(iC, index, nIds);
		delete [] index;
	}

	ids = spr_ids_com.item1(iC, 0);
	nIds = indexIds(&id_com, ids, &index);

	if (nIds) {
		if(P_SPR_COM) spr_com.setItem(iC, index, nIds);
		delete [] index;
	}
    }

    delete [] justOne;

    Property obsdat;
    Property sprdat;
    Property current_com;

    if(P_OBS_DAT_COM) obsdat.open("obs_dat.com", 1);
    if(P_SPR_DAT_COM) sprdat.open("spr_dat.com", 1);

    int4 *supers, nSupers, flag, i, j;

    int **obsIdIdx = new int *[spr_com.nItems];
    int *nObsIdIdx = new int[spr_com.nItems];
    for (iC = 0; iC < spr_com.nItems; iC++) {
	nObsIdIdx[iC] = obs_com.getItemSize(iC);
	if (nObsIdIdx[iC]) obsIdIdx[iC] = obs_com.item4(iC, 1);
    }

    /* go through all SPRSDE records to match OBSLTE records */
    for (iC = 0; iC < spr_com.nItems; iC++) {

	/* get my obsoleted entries */
	index = spr_com.item4(iC, 0);
	nIds = spr_com.getItemSize(iC);

	/* check to see that the obsoleted entries refer back to me */
	for (i = 0; i < nIds; i++) {
		if (index[i] == -1) continue;

		/* get the superseeding entries */
		supers = obs_com.item4(index[i], 1);
		nSupers = obs_com.getItemSize(index[i]);

		/* fix a missing OBSLTE record for a single entry */
		if (nSupers == 1 &&
		    supers[0] == -1 &&
		    !*obsdat.item4(index[i])) {

			obsIdIdx[index[i]][0] = (int4) iC;
			obsdat.setItem(index[i], (int4)*sprdat.item4(iC));

			delete [] supers;
			continue;

		}

		/* now we have some and need to check if I'm in */
		flag = 1;
		for (j = 0; j < nSupers; j++) {
			/* there I am! Hold your fire! */
			if (supers[j] == iC) { flag = 0; break; }
		}

		if (flag)  {
		    /* add me to the superseeding entries */
		    _addToArray(&obsIdIdx[index[i]], iC, nObsIdIdx[index[i]]);
		    nObsIdIdx[index[i]]++;
		}

		delete [] supers;
	}
    }

    if(P_OBS_COM) {
	obs_com.clear();
	for (iC = 0; iC < spr_com.nItems; iC++) {
		obs_com.addItem(obsIdIdx[iC], nObsIdIdx[iC]);
		if (nObsIdIdx[iC]) delete [] obsIdIdx[iC];
	}
	obs_com.close();
    }

    delete [] obsIdIdx;
    delete [] nObsIdIdx;
			
    if(P_SPR_DAT_COM) sprdat.close();
    if(P_OBS_DAT_COM) obsdat.close();
    if(P_SPR_COM) spr_com.close();
    if(P_OBS_IDS_COM) obs_ids_com.close();
    if(P_SPR_IDS_COM) spr_ids_com.close();

    if(P_SPR_COM) spr_com.open("spr.com", 1);
    if(P_OBS_DAT_COM) obsdat.open("obs_dat.com", 1);

    int numIds = 0, *versions, alive, n, *currents = new int[spr_com.nItems];

    /* set the current entry to either myself, my last in line or -1 for
       a deleted entry */

    for (iC = 0; iC < spr_com.nItems; iC++) {

	if (!*obsdat.item4(iC)) currents[iC] = iC;
	else {
		alive = 0; numIds = 0;
		numIds = LinkedId::getChain(iC, &versions);

		for (n = 0; n < numIds; n++) 
			/* use the first entry to keep my line alive */
			if (!*obsdat.item4(versions[n])) {
				currents[iC] = versions[n];
				alive = 1;
				break;
			}

		if (!alive) currents[iC] = -1;
		if (numIds) { delete [] versions; numIds = 0; }
	}
    }

    if(P_CURRENT_COM) current_com.open("current.com", 1);
    current_com.clear();
    for (iC = 0; iC < spr_com.nItems; iC++)
	current_com.addItem((int4) currents[iC]);
    
    delete [] currents;

    if(P_OBS_DAT_COM) obsdat.close();
    if(P_SPR_COM) spr_com.close();
    if(P_CURRENT_COM) current_com.close();

}

/*
     the method below populates misc.col, a "collective" Collection object
     the indices of misc.col refer to the following information:

     index	what		quality		incl.	obs only
						NMR

	 0	unit cell	unchanged	yes	yes
	 1			changed		
	 2			n/d (NMR)
	 3	# of atoms	unchanged	yes	yes
	 4			increased
	 5			decreased
	 6	Z-value		unchanged	no	yes
	 7			changed
	 8	Spacegroup	unchanged	no	yes
	 9			changed
	10	non X-ray structures		yes	no
	11	X-ray structures		no	no
	12	deleted structure		yes	yes
	13-19	# of replacing structures	yes	yes
	20	replaced by >7 structures	yes	yes
	21-27	# of replaced structures	yes	no
	28	replaced >7 structures		yes	no
	29	replaced			yes	yes
	30	current				yes	no!
	31	# of HOH	unchanged	yes	yes
	32			increased
	33			decreased
	34	has only OBSLTE record		yes	I guess :)
	35	has SPRSDE and OBSLTE record	yes	yes
	36	has only SPRSDE record		yes	no
	37	has neither			yes	no!
	38      resolution      unchanged       no      yes
	39                      increased
        40                      decreased

	below are currents only!

	41	x-ray structures		no	no!
	42	nmr structures			yes	no!
	43	theoretic. mdl.			no	no!

	molecule type

	44	protein structures		yes	no!
	45	protein/dna complexes		yes	no!
	46	protein/rna complexes		yes	no!
	47	protein/na complexes		yes	no!
	48	dna structures			yes	no!
	49	rna structures			yes	no!
	50	dna/rna structures		yes	no!
	51	carbohydrates			yes	no!
	52	other				yes	no!


     notes: 1) a qualitative entry is the index of the replaced structure.
	       (i.e. if entry #0 is replaced by #1 which has more atoms,
	       #0 will be stored in misc.col #4)
	    2) a ! in the 'obs only' field means that obs are excluded
	    3) if an entry is replaced by multiple entries, only the first
	       of the replacing entries is considered.
	    4) changes in #atoms for NMR only reflect 1st model
*/
	
void DF::collect(int nCom1, int nCom2, int log_mode, int df2) {
	int4 iC, *obsIds, nObsIds, ii, index, nChains, nonXrayFlag = 0, atoms = 0,
	    tmp = 0, zval, nSprIds, *sprIds, year, nWater = 0;

	if(!OBS_COL_GROUP) return;
	
	int *dates;

	double res;
	char *string, *space;

	double lRes = 0.5;
	double hRes = 6.0;
	int nCRes = (int) ((hRes - lRes)*20);

	int lCha = 0;
	int hCha = 25;
	int nCChains = (hCha - lCha)*2;

	Property id("id.com", 1);
	Property expdta("expdta.com", 1);
	Property res_com("res.com", 1);
	Property obs_com("obs.com", 1);
	Property spr_com("spr.com", 1);
	Property obs_dat("obs_dat.com", 1);
	Property current_com("current.com", 1);
	Property spr_dat("spr_dat.com", 1);
	Property dep_dat("date_int.com", 1);
	Property rel_dat("reldat.com", 1);
	Property n_enp_com("n_enp.com", 1);
	Property unit("unitcell.com", 1);
	Property zval_com("zval.com", 1);
	Property spcgrp("spcgrp.com", 1);

	Collection res_col;
	Collection chains_col;
	Collection misc_col;
	Collection depyear;
	Collection relyear;
	Collection obsyear;
	Collection spryear;

	#ifdef ALPHA
		int now = time(0);
	#else
		long int now = time(0);
	#endif

	tm *tmstruc = localtime(&now);
	year = (*tmstruc).tm_year - 72;

	// we always regenrate the full objects!
	res_col.create("res.col", D_INT4, nCRes);
	chains_col.create("chains.col", D_INT4, nCChains);
	misc_col.create("misc.col", D_INT4, 53);

	depyear.create("depyear.col", D_INT4, year+1);
	relyear.create("relyear.col", D_INT4, year+1);
	obsyear.create("obsyear.col", D_INT4, year+1);
	spryear.create("spryear.col", D_INT4, year+1);

	Property i_enp("i_enp.com");
	Property type_enp("type.enp");

	for (iC = 0; iC < nCom2; iC++) {

		if (!log_mode)
			printf("DF 8/%d, COM %d (%d/%d)\n", df2, iC,
							    iC+1,
							    nCom2);


		dates = JDate::jdate(*spr_dat.item4(iC));
		year = dates[2] - 1972;
		if (year >= 0 && year < spryear.getCollectionSize())
						spryear.add(year, iC);

		delete [] dates;

		dates = JDate::jdate(*rel_dat.item4(iC));
		year = dates[2] - 1972;
		if (year >= 0 && year < relyear.getCollectionSize())
						relyear.add(year, iC);

		delete [] dates;

		dates = JDate::jdate(*dep_dat.item4(iC));
		year = dates[2] - 1972;
		if (year >= 0 && year < depyear.getCollectionSize())
						depyear.add(year, iC);

		delete [] dates;

		dates = JDate::jdate(*obs_dat.item4(iC));
		year = dates[2] - 1972;
		if (year >= 0 && year < obsyear.getCollectionSize())
						obsyear.add(year, iC);
		
		delete [] dates;

		nObsIds = obs_com.getItemSize(iC);
		obsIds = obs_com.item4((int4)iC, 0);

		nSprIds = spr_com.getItemSize(iC);
		sprIds = spr_com.item4((int4)iC, 0);

		// # of superseded entries

		if (sprIds[0] != -1) {
			index = nSprIds + 20;
			if (index < 28) misc_col.add(index, iC);
			else misc_col.add(28, iC);
		}

		string = flt2str(unit.itemf(iC,0));

		// replaced/deleted/current?

		if (nObsIds == 1 && obsIds[0] == -1) {
			if (*obs_dat.item4(iC)) {	// deleted
				misc_col.add(12, iC);
				misc_col.add(29, iC);
			} else {
				if (*current_com.item4(iC) != -1)
					misc_col.add(30, iC);	// current
			}
		} else misc_col.add(29, iC);		// replaced

		/* has only OBSLTE record */
		if (*obs_dat.item4(iC) && !*spr_dat.item4(iC))
					misc_col.add(34, iC);

		/* has only SPRSDE record */
		if (*spr_dat.item4(iC) && !*obs_dat.item4(iC))
					misc_col.add(36, iC);

		/* has both OBSLTE and SPRSDE records */
		if (*spr_dat.item4(iC) && *obs_dat.item4(iC))
					misc_col.add(35, iC);

		/* has neither */
		if (!*spr_dat.item4(iC) && !*obs_dat.item4(iC))
					misc_col.add(37, iC);

		// # of replacing entries

		index = 12 + nObsIds;
		if (index < 20) misc_col.add(index, iC);
		else misc_col.add(20, iC);

		nonXrayFlag = 0;
		if (*expdta.item1(iC) != 9) nonXrayFlag = 1;

		if (nonXrayFlag) {
			misc_col.add(10, iC);	// model and nmr or other

			// model
			if (*expdta.item1(iC) == 8 && obsIds[0] < 0)
							misc_col.add(43, iC);

			// nmr
			if (*expdta.item1(iC) == 4 && obsIds[0] < 0)
							misc_col.add(42, iC);
		} else {
			misc_col.add(11, iC);	// x-ray
			if (obsIds[0] < 0) misc_col.add(41, iC);
		}

		// # structure type

	        char type = Derive::getType(iC);

		// mixed structures need more work
		if (type == 'X' || type == 'N') {

		    int pflag, dflag, rflag;
		    pflag = dflag = rflag = 0;

		    int iEnp = *i_enp.item4(iC);

		    for (int n = 0; n < *n_enp_com.item2(iC); n++) {
			switch(*type_enp.item1(n+iEnp)) {
			    case 'P': pflag = 1; break;
			    case 'D': dflag = 1; break;
			    case 'R': rflag = 1; break;
			    case ' ': type = 'O'; break;
			}
		    }

		    if (type != 'O') {
			    if (pflag*dflag*rflag) type = 'A';	// protein/na
			    else if (pflag*dflag)  type = 'B';  // protein/dna
			    else if (pflag*rflag)  type = 'F';  // protein/rna
			    else if (rflag*dflag)  type = 'N';  // dna/rna
			    else if (rflag)	   type = 'R';  // rna
			    else if (dflag)	   type = 'D';  // dna
		    }

		}

		switch(type) {
		    case 'P': misc_col.add(44, iC); break;
		    case 'B': misc_col.add(45, iC); break;
		    case 'F': misc_col.add(46, iC); break;
		    case 'A': misc_col.add(47, iC); break;
		    case 'D': misc_col.add(48, iC); break;
		    case 'R': misc_col.add(49, iC); break;
		    case 'N': misc_col.add(50, iC); break;
		    case 'C': misc_col.add(51, iC); break;
		    case 'O': misc_col.add(52, iC); break;
		}

		// continue only for replaced structures

		if (obsIds[0] < 0) {
			if (string) delete [] string;
			continue;
		}

		if (nonXrayFlag) misc_col.add(2, iC);

		// resolution change

		res = (double) *res_com.itemf(iC);
		if (res >= lRes && res <= hRes &&
		    *res_com.itemf(obsIds[0]) >= lRes &&
		    *res_com.itemf(obsIds[0]) <= hRes) {

			index = (int) (*res_com.itemf(obsIds[0])*1000
				       - res*1000)/100 + nCRes/2;

			res_col.add(index, iC);
		}

		// number of chains

		if (res > *res_com.itemf(obsIds[0])) misc_col.add(39, iC);
		else if (res < *res_com.itemf(obsIds[0])) misc_col.add(40, iC);
		     else misc_col.add(38, iC);

		nChains = *n_enp_com.item2(iC);

		if (nChains >= lCha && nChains <= hCha &&
		    *n_enp_com.item2(obsIds[0]) >= lCha &&
		    *n_enp_com.item2(obsIds[0]) <= hCha) {
			
			index = *n_enp_com.item2(obsIds[0]) - nChains + 
				hCha - lCha;
			
			chains_col.add(index, iC);
		}

		// number of atoms

		atoms = Derive::numAtoms(iC);
		tmp = Derive::numAtoms(obsIds[0]);
	
		if (tmp > atoms) misc_col.add(4, iC);
		else if (tmp < atoms) misc_col.add(5, iC);
		     else misc_col.add(3, iC);

		// change in # of HOH

		// fprintf(stderr, "calculating HOH for %d\n", iC);
		nWater = Derive::numHOH(iC);
		tmp = Derive::numHOH(obsIds[0]);

		if (tmp > nWater) misc_col.add(32, iC);
		else if (tmp < nWater) misc_col.add(33, iC);
		     else misc_col.add(31, iC);
	    
		if (!nonXrayFlag) {
		    // unit cell changes

		    char *temp = flt2str(unit.itemf(obsIds[0],0));

		    if (!strcmp(string, temp)) {
			misc_col.add(0, iC);
		    } else {
			misc_col.add(1, iC);
		    }

		    delete [] temp;

		    // Z-value changes

		    zval = (int) *zval_com.item2(iC);
		    tmp = (int) *zval_com.item2(obsIds[0]);

		    if (tmp == zval) misc_col.add(6, iC);
		    else misc_col.add(7, iC);

		    // space group changes

		    space = spcgrp.item1(iC, 0);

		    if (!strcmp(space, spcgrp.item1(obsIds[0], 0)))
			misc_col.add(8, iC);
		    else misc_col.add(9, iC);
		}

		if (string) delete [] string;
	}

	res_col.save();
	spryear.save();
	obsyear.save();
	depyear.save();
	relyear.save();
}


char *DF::flt2str(float *f) {
        char *out = new char[43];
        sprintf(out, "%.2f %.2f %.2f %.2f %.2f %.2f", f[0],
                                                      f[1],
                                                      f[2],
                                                      f[3],
                                                      f[4],
                                                      f[5]);
	// printf("%s\n", out);
        return out;
}

void DF::collFDS(int nCom1, int nCom2, int log_mode, int df2) {

  if(!FDS_GROUP) return;

	Monomers mons;
	Entities ents;

	Property n_enp("n_enp.com", 1);
	Property i_enp("i_enp.com", 1);
	Property res("res.com", 1);

	Collection fds_col;
	Property fds_enp;
	Property fds_com;

	int iEnt, nEnt, entIn, iCol, max = 30, iCom;
	double fds;

	if(!nCom1) {
		fds_col.create("fds.col", D_INT4, max+1);
		fds_enp.create("fds.enp", D_FLT4, 1, 0);
		fds_com.create("fds.com", D_FLT4, 1, 0);
	} else {
		fds_col.open("fds.col");
		fds_enp.open("fds.enp", 0);
		fds_com.open("fds.com", 0);
	}

	int numFDS, numEnts = 0;
	double sumFDS;

	for(iCom = nCom1; iCom < nCom2; iCom++) {
		if (!log_mode)
			printf("DF 4/%d, COM %d (%d/%d)\n", df2, iCom,
							    iCom-nCom1+1,
							    nCom2-nCom1);
		iEnt = 0;
		ents.addCom(iCom, 0);
		nEnt = *n_enp.item2(iCom);
		entIn = *i_enp.item4(iCom);
		sumFDS = 0; numFDS = 0;

		for (IEntity iE(&ents); iEnt < nEnt; ++iE, ++iEnt, entIn++) {
			fds = Stats::getFDS(iE);

			fds_enp.addItem((float) fds);

			if (fds == MAXVALUE) continue;
			numFDS++; sumFDS += fds;

		}

		if (numFDS && res.itemf(iCom)) {

			iCol = (int) (sumFDS/(0.5*numFDS));
			if (iCol > max) iCol = max;

			fds_col.add(iCol, iCom);
			fds_com.addItem((float) sumFDS/numFDS);

		} else {

			fds_com.addItem((float) 0.);
		}
	}

	fds_col.save();
}

void DF::calcChi1(int nCom1, int nCom2, int log_mode, int df2) {

  if(!P_CHI1_ENP && !P_CHI1_COM) return;

        Monomers mons;
        Entities ents;

        Property id("id.com", 1);
        Property n_enp("n_enp.com", 1);
        Property i_enp("i_enp.com", 1);

        Property chi1_enp;
        Property chi1_com;

        int iEnt, nEnt, entIn;

	if(!nCom1) {
		if(P_CHI1_ENP) chi1_enp.create("chi1.enp", D_FLT4, 6, 0);
		if(P_CHI1_COM) chi1_com.create("chi1.com", D_FLT4, 1, 0);
	} else {
		if(P_CHI1_ENP) chi1_enp.open("chi1.enp", 0);
		if(P_CHI1_COM) chi1_com.open("chi1.com", 0);
	}

        Summary *statsSum;
        int count, total;
        float *array = new float[8];

	for(int iCom = nCom1; iCom < nCom2; iCom++) {
		if (!log_mode)
			printf("DF 5/%d, COM %d (%d/%d)\n", df2, iCom,
							    iCom-nCom1+1,
							    nCom2-nCom1);
		iEnt = 0;
		ents.addCom(iCom, 0);
		nEnt = *n_enp.item2(iCom);
		entIn = *i_enp.item4(iCom);

                total = 0;
                for (int i = 0; i < 7; i++) array[i] = 0.;

		for (IEntity iE(&ents); iEnt < nEnt; ++iE, ++iEnt, entIn++) {
                    statsSum = Stats::sumDihedrals(iE, count);

                    if (statsSum[8].num) {
                        array[0] = (float) statsSum[8].total/statsSum[8].num;
                        array[1] = (float) statsSum[8].stddev();
                    }

                    if (statsSum[9].num) {
                        array[2] = (float) statsSum[9].total/statsSum[9].num;
                        array[3] = (float) statsSum[9].stddev();
                    }

                    if (statsSum[10].num) {
                        array[4] = (float) statsSum[10].total/statsSum[10].num;
                        array[5] = (float) statsSum[10].stddev();
                    }

		    array[6] += (float) (array[1]*statsSum[8].num +
				        array[3]*statsSum[9].num +
				        array[5]*statsSum[10].num);

		    total += statsSum[8].num + statsSum[9].num +
			     statsSum[10].num;

		    if(P_CHI1_ENP) chi1_enp.addItem(array);
		    delete [] statsSum;
                }

                if (total) {
                        array[6] /= total;
                }

                if(P_CHI1_COM) chi1_com.addItem(array[6]);
	}

	delete [] array;
}

void DF::proteinFeatures(char *db,
			 int nCom1,
			 int nCom2,
			 int log_mode,
			 int df2) {

  if(!MOOSE_PATTERNS_GROUP) return;


	  char *prefTab[] = {"exp", "pol", "bfac", "spol", "shyd", "svol", 
			     "sexp", "siso"};

	  char *nameTab[] = {"Exposure", 
			     "Polarity",
			     "Temperature factor",
			     "Static Polarity", 
			     "Static Hydrophobicity",
			     "Static Volume", 
			     "Static Exposure",
			     "Static Isoelectric Point"};
	  
	  int iProp, iPT0 = 0, i, j, k, l, seq_, nStates = 6,
	      iMin, iMax, fSize, ind, iS1, iS2, aa_code_, nSum,
	      nnHits1 = 240,  nnPatterns2 = (int4) pow(nStates, 5), 
	      nnHits2 = 6*nnPatterns2;

	  int1 *seq;

	  int2 *eHits1 = new int2[nnHits1], *eHits2, iEnt;

	  int4 tEnt1, tEnt2, *aa_freq, *aa_freq8, *pDistr0_8, nDistr0_8,
	       sePropAvI[160], nsePropAv, nStates8; int4 nProp;

	  int1 *sePropAvC;
	  
	  char buffer[100], *comEnt_, **comEnt, *propPac, file_buf[20];

	  flt4 *states8, *prop;

	  double *states = new double[nStates],
		 *sePropAvD = new double[160], 
		 clSize, pMin, pMax, pSpan, dProp, sum;

	  Property i_enp_com("i_enp.com");
	  Property iseq_flt_enp("iseq_flt.enp", 0);
	  Property prop8_ave("prop8.ave", 0);
	  Property prop8_sts("prop8.sts", 0);
	  Property prop_tmp;
	  Property glob_1b;

	  Collection imm;
	  Collection ip6;
	  Property pp6;
	  Property pp8;
	  
	  tEnt1 = *i_enp_com.item4(nCom1);
	  tEnt2 = *i_enp_com.item4(nCom2);
	  
//	  printf("%d entities found\n", tEnt);

	  for (i = 0; i<160; i++) {
	      sePropAvI[i] = 0;
	      sePropAvD[i] = 0.0;
	  }

	  if (!nCom1) {
	      pDistr0_8 = new int4 [8000];
	      states8 = new flt4 [8*(nStates+1)];
	      for (i = 0; i < 8000; i++) pDistr0_8[i] = 0;
	      sePropAvC = new int1 [160];
	      aa_freq = new int4 [20*nStates];
	      aa_freq8 = new int4 [160];
	  } else {
	      states8 = prop8_sts.itemf(0);
	      sePropAvC = prop8_ave.item1(0);
	  }

	  for (int iPT = iPT0; iPT < 8; iPT++) {
	      eHits2 = new int2 [nnHits2];

	      clSize = iPT > 1  ?  5.0  :  2.0;

	      strcpy(file_buf, prefTab[iPT]);
	      strcat(file_buf, ".imm");

	      if (nCom1) imm.open(file_buf);
	      else imm.create(file_buf, D_INT2, nnHits1);

	      strcpy(file_buf, prefTab[iPT]);
	      strcat(file_buf, "_1b.enp");
	      glob_1b.open(file_buf);

	      if (!nCom1) glob_1b.clear();

	      strcpy(file_buf, prefTab[iPT]);
	      strcat(file_buf, ".ip6");

	      if (nCom1) ip6.open(file_buf);
	      else ip6.create(file_buf, D_INT2, nnHits2);

	      if (iPT < 3) {
		  strcpy(file_buf, prefTab[iPT]);
		  strcat(file_buf, ".pp6");
		  pp6.open(file_buf);

		  strcpy(file_buf, prefTab[iPT]);
		  strcat(file_buf, ".pp8");
		  pp8.open(file_buf);
	      }

	      if (!nCom1) {
		  int4 *pDistr0 = pDistr0_8+1000*iPT;
		  pMin = 100000.0; pMax = 0.0; 

		  if (iPT < 3) {
		      for (i = 0; i < 20; i++) for (j = 0; j < nStates; j++) 
			aa_freq[i*nStates+j] = 0;
		      for (i = 0; i < 20; i++) for (j = 0; j < 8; j++) 
			aa_freq8[i*8+j] = 0;
		  }

		  for (iEnt = tEnt1; iEnt < tEnt2; iEnt++) {
		      if (!log_mode)
			printf("DF 6/%d, COM %d (%d/%d)\n", df2, iEnt,
							    iEnt-tEnt1+1,
							    tEnt2-tEnt1);
		      
		      if (iPT == 0) {
			  prop_tmp.open("exp_flt.enp");
			  prop = prop_tmp.itemfn(iEnt, nProp);
		      }

		      if (iPT == 1) {
			  prop_tmp.open("pol_flt.enp");
			  prop = prop_tmp.itemfn(iEnt, nProp);
		      }

		      if (iPT == 2) {
			  prop_tmp.open("bfac_flt.enp");
			  prop = prop_tmp.itemfn(iEnt, nProp);
		      }

		      if (iPT > 2) {

			  seq = iseq_flt_enp.item1n(iEnt, nProp);
			  prop = new flt4 [nProp];

			  for (i = 0; i < nProp; i++) {
			      prop[i] = -1.0;
			      if (seq[i] > 19) continue;
			      aa_code_ = seq[i];

			      switch(iPT) {
				case 3:
				 dProp = properties[aa_code_].polarity;
				 break;

				case 4:
			         dProp = properties[aa_code_].hydrophobicity;
				 break;

				case 5:
			         dProp = properties[aa_code_].volume;
				 break;

				case 6:
			         dProp = properties[aa_code_].meanExposure;
				 break;

				case 7:
			         dProp = properties[aa_code_].isoelectricPoint;
				 break;
			      }

			      prop[i] = dProp;

			  }
		      }
		      
		      for (iS1 = 0; iS1 < nProp; iS1++) {
			  dProp = prop[iS1];
			  if (dProp < 0.0) continue;
			  if (pMin > dProp) pMin = dProp;
			  if (pMax < dProp) pMax = dProp;
			  iProp = (int) (dProp*10.0);
			  if (iProp > 999) iProp = 999;
			  pDistr0[iProp]++;
		      }

		      if (iPT > 2) delete [] prop;
		  }
		  
		  nSum = 0;
		  for (i = 0; i < 1000; i++) nSum += pDistr0[i];
		  
		  int pSum = 0, iState = 0;
		  for (i = 0; i < 1000; i++) {
		      pSum += pDistr0[i];
		      if (((double)pSum)/nSum >= (iState+1.0)/nStates) {
			  states[iState] = (i+1)*0.1;
			  iState++;
			  if (iState == nStates-1) break;
		      }
		  }
		  
		  for (i = 0; i < nStates-1; i++)
			states8[iPT*(nStates+1)+i+1] = states[i];

		  states8[iPT*(nStates+1)] = pMin;
		  states8[iPT*(nStates+1) + nStates] = pMax; 
		  pSpan = pMax-pMin;     
		  

//	          printf("States :  ");
//	          for (i = 0; i < nStates+1; i++) 
//		      printf("%.2f ", states8[iPT*(nStates+1)+i]*0.01);
//	          printf("\n");

	      } else {		// if (nCom1)

		  if (iPT < 3) {
		      aa_freq = pp6.item4(0, 1);
		      aa_freq8 = pp8.item4(0, 1);
		  }
		  
		  pMin = states8[iPT*(nStates+1)]; 
		  pMax = states8[iPT*(nStates+1)+nStates]; 
		  pSpan = pMax-pMin;

		  for (i = 0; i < nStates-1; i++) 
		    states[i] = states8[iPT*(nStates+1)+i+1];
	      }

	      states[nStates-1] = 100.0;
	      
	      for (iEnt = tEnt1; iEnt < tEnt2; iEnt++) {
		  if (!log_mode)
			printf("DF 7/%d, COM %d (%d/%d)\n", df2, iEnt,
							    iEnt-tEnt1+1,
							    tEnt2-tEnt1);
		  
		  if (iPT == 0) {
		      prop_tmp.open("exp_flt.enp");
		      prop = prop_tmp.itemf(iEnt);
		  }

		  if (iPT == 1) {
		      prop_tmp.open("pol_flt.enp");
		      prop = prop_tmp.itemf(iEnt);
		  }

		  if (iPT == 2) {
		      prop_tmp.open("bfac_flt.enp");
		      prop = prop_tmp.itemf(iEnt);
		  }
		  
		  seq = iseq_flt_enp.item1n(iEnt, nProp);

		  if (iPT > 2) {
		      prop = new flt4 [nProp];
		      for (i = 0; i < nProp; i++) {
			  prop[i] = -1.0;
			  if (seq[i] > 19) continue;
			  aa_code_ = seq[i];
			  
		          switch(iPT) {
			     case 3:
				 dProp = properties[aa_code_].polarity;
				 break;

			     case 4:
			         dProp = properties[aa_code_].hydrophobicity;
				 break;

			     case 5:
			         dProp = properties[aa_code_].volume;
				 break;

			     case 6:
			         dProp = properties[aa_code_].meanExposure;
				 break;

			     case 7:
			         dProp = properties[aa_code_].isoelectricPoint;
				 break;
			  }

			  prop[i] = dProp;
		      }
		  }
	  
		  
		  for (i = 0; i < nnHits1; i++) eHits1[i] = 0;
		  for (i = 0; i < nnHits2; i++) eHits2[i] = 0;

		  for (iS1 = 0; iS1 < nProp; iS1++) {
		      if (prop[iS1] < 0.0) continue;
		      dProp = prop[iS1];
		      i = seq[iS1];

		      if (i < 20) {
			  if (iPT < 3)
			    for (k = 0; k < nStates; k++)
			      if (dProp <= states[k]) {
				  aa_freq[i*nStates+k]++;
				  break;
			      }

			  sePropAvD[iPT*20+i] += prop[iS1];
			  sePropAvI[iPT*20+i]++;
		      }
		      
		      for (j = 0; j < 6; j++) {
			  fSize = 5*(j+1);
			  if (iS1+fSize >= nProp) break;
			  
			  iMin = 0; iMax = 19;
			  for (k = 0; k < fSize; k++) {
			      dProp = prop[iS1+k];
			      if (dProp < 0.0) goto exit301;
			      iProp =  (int) (dProp/clSize);
			      if (iProp > iMin) iMin = iProp;
			      if (iProp < iMax) iMax = iProp;
			      if (iMax > 19) iMax = 19;
			      if (iMin > 19) iMin = 19;
			      if (iMax < 0 || iMin < 0) {
				  printf("Property value error at %d %d\n",
					 iS2, iProp);
				  exit(0);
			      }
			  }

			  for (k = iMin; k < 20; k++) {
			      eHits1[k*12+j] = 1;
			  }

			  for (k = 0; k <= iMax; k++) {
			      eHits1[k*12+j+6] = 1;
			  }

		        exit301 :  
			  
			  ind = 0;
			  for (i = 0; i < 5; i++) {
			      sum = 0.0;
			      for (k = 0; k < j+1; k++) {
				  dProp = prop[iS1+i*(j+1)+k];
				  if (dProp < 0.0) goto exit302;
				  sum += dProp;
			      }

			      sum /= (j+1);
			      for (k = 0; k < nStates; k++)
				if (sum < states[k]) {
				    ind = ind*nStates+k;
				    break;
				}
			  } 

			  if (ind < 0 || ind >= nnPatterns2) {
			      printf("Index %d\n", ind); exit(0);
			  }

			  eHits2[nnPatterns2*j+ind] = 1;
			  
		        exit302 :  ;
		      }
		  }
		  
		  for (i = 0; i < nnHits1; i++) {
		      if (eHits1[i]) imm.add(i, iEnt);
		  }
		  
		  for (i = 0; i < nnHits2; i++) {
		      if (eHits2[i]) ip6.add(i, iEnt);
		  }
		  
		  int1 *propPac = new int1 [nProp];

		  for (iS1 = 0; iS1 < nProp; iS1++) {
		      propPac[iS1] = 0;
		      iProp =  (int) prop[iS1];
		      if (iProp < 0) continue;
		      iProp =  (int) ((iProp-pMin)/pSpan*255+1);
		      if (iProp >= 256) iProp = 255;
		      if (iProp < 1) iProp = 1;

		      propPac[iS1] = iProp;

		      if (iPT < 3) {
			  j = iProp/32;
			  i = seq[iS1];
			  if (i >= 0 && i < 20 && j >= 0 && j < 8)
				aa_freq8[i*8+j]++;
		      }
		  }
		  
		  glob_1b.addItem(propPac, nProp);

		  delete [] propPac;

		  if (iPT > 2) delete [] prop;
	      }
	      
	      if (!log_mode)
		printf("Write global property %d :  length = %d\n", 
		       iPT, glob_1b.getPropertySize());
	      
	      sprintf(buffer, "%s/table_%s_mm.html", db, prefTab[iPT]);

	      FILE *tmp = fopen(buffer, "w");
//	      fprintf(tmp, "<TITLE>%s sites in proteins</TITLE>\n", 
//		      nameTab[iPT]);
//	      fprintf(tmp, "<H2>%s sites in proteins</H2>\n", nameTab[iPT]);
	      for (k = 0; k < 2; k++) {

		  fprintf(tmp, 
			  "<HR> <P><H3>Property values %s threshold</H3>\n", 
			  k == 0 ? "below" : "above");

		  for (int p = 0; p < 2; p++) {
		      fprintf(tmp, 
	      "<TABLE BORDER CELLPADDING=5>\n<TR ALIGN=CENTER><TH></TH> \n");

		      for (i = p*10; i < p*10+10; i++) {
			  fprintf(tmp, "<TH> %d </TH>", 
				  (int)(i*clSize+clSize));
		      }

		      fprintf(tmp, "</TR>\n");
		      for (j = 0; j < 6; j++) {
			  fprintf(tmp, "<TR ALIGN=CENTER>\n<TH> %d </TH>\n",
				  5*(j+1));
			  for (i = p*10; i < p*10+10; i++)
			    fprintf(tmp, "<TD>%d</TD>", 
				    ip6.getClassSize(i*12+j+k*6));
			  fprintf(tmp, "</TR>\n");
		      }

		      fprintf(tmp, "</TABLE><P>\n");
		  }
	      }

	      fclose(tmp);
	      
	      imm.save();
	      ip6.save();
	      
	      if (iPT < 3) {
//		  pp6.clear();
		  strcpy(file_buf, prefTab[iPT]);
		  strcat(file_buf, ".pp6");
		  pp6.create(file_buf, D_INT4, 120, 1);
		  pp6.addItem(aa_freq, nStates*20);

//		  pp8.clear();
		  strcpy(file_buf, prefTab[iPT]);
		  strcat(file_buf, ".pp8");
		  pp8.create(file_buf, D_INT4, 160, 1);
		  pp8.addItem(aa_freq8, 160);

		  if (nCom1) {
		      delete [] aa_freq;
		      delete [] aa_freq8;
		  }

	      }

	      delete [] eHits2;

	      for (i = 0; i < 20; i++) 
		sePropAvC[iPT*20+i] =
		   (int1) (sePropAvI[iPT*20+i] > 0  ? 
		           ((sePropAvD[iPT*20+i]/sePropAvI[iPT*20+i]-pMin)
			    /pSpan*254+1)
		           :  0);

	  }

	  delete [] eHits1;
	  delete [] states;
	  delete [] sePropAvD;

	  if (!nCom1) {
//	      prop8_ave.clear();
	      prop8_ave.create("prop8.ave", D_INT1, 160, 1);
	      prop8_ave.addItem(sePropAvC);

//	      prop8_sts.clear();
	      prop8_sts.create("prop8.sts", D_FLT4, 56, 1);
	      prop8_sts.addItem(states8);
	  
	      delete [] pDistr0_8;
	  }
}

void DF::mooseReport(int nCom1, int nCom2, int log_mode, int df2) {


  if(!MOOSE_REPORT_GROUP) return;

	int i, j, k;
	double x, entMolWeight;
	int ii, ssAlpha, ssBeta, nAlpha, nBeta;
	char nEntType[3], entType;
	char *header, *author, *compnd, ec[40], *buf_;
	char *ks;

	int iCom, iEnp, iEnc, iSE;

	int4 nSeq; int2 *seq;

	char *seName;

	Property i_enp_com("i_enp.com");
	Property i_enc_enp("i_enc.enp");
	Property se_enc("se.enc");
	//Property k_s_enp("k_s.enp");
	Property code3_mon("code3.mon");

	Property ec_com("ec.com");
	Property compnd_com("compnd.com");

	Property type_enp("type.enp", 0);
	Property mw_enp("mw.enp", 0);
	Property alpha_c_enp("alpha_c.enp", 0);
	Property beta_c_enp("beta_c.enp", 0);
	Property alpha_n_enp("alpha_n.enp", 0);
	Property beta_n_enp("beta_n.enp", 0);
	Property ss_seg_enp("ss_seg.enp", 0);

	if (!nCom1) {
		ec_com.clear();
		type_enp.clear();
		mw_enp.clear();
		alpha_c_enp.clear();
		beta_c_enp.clear();
		alpha_n_enp.clear();
		beta_n_enp.clear();
		ss_seg_enp.clear();
	}

	for (iCom = nCom1; iCom < nCom2; iCom++) {
		if (!log_mode)
		printf("DF 7/%d, COM %d (%d/%d)\n", df2, iCom,
						    iCom-nCom1+1,
						    nCom2-nCom1);

		compnd = compnd_com.item1(iCom);

		ec[0] = '\0';
		if (strlen(compnd) > 10)
		    for (j = 0; j < strlen(compnd)-5; j++)
			if (!strncmp(&compnd[j], "E.C.", 4) ||
			    !strncmp(&compnd[j], "EC:", 3)) {
			    buf_ = &compnd[j+4];
			    for (i = 0; i < strlen(buf_); i++)
			      if (buf_[i] == ' ' || buf_[i] == ')' ||
				  buf_[i] == ';') break;
			    strncpy(ec, buf_, i); ec[i] = '\0';
			}

		ec_com.addItem(ec);

		nEntType[0] = ' '; nEntType[1] = ' '; nEntType[2] = ' '; 

		for (iEnp = *i_enp_com.item4(iCom);
		    iEnp < *i_enp_com.item4(iCom+1); 
		    iEnp++) {
		  entMolWeight = 0;
		  entType = ' '; ssAlpha = 0; ssBeta = 0;
		  
		  iEnc = *i_enc_enp.item4(iEnp);
		  
		  seq = se_enc.item2n(iEnc, nSeq);
		  //ks = k_s_enp.item1(iEnp);

		  for (iSE = 0; iSE < nSeq; iSE++) {
		      seName = code3_mon.item1(seq[iSE]);

		      for (i = 0; i < 20; i++)
			if (!strcmp(properties[i].code3, seName)) {
			    entMolWeight += properties[i].molWeight;
			    if (entType == ' ') entType = 'P';
			    break;
			}

		      if (entType == 'P') continue;

		      if (!strcmp("  A", seName)) {
			  entMolWeight += 347.22;
		          if (entType == ' ') entType = 'D';
		      } else if (!strcmp("  T", seName)) {
			  entMolWeight += 322.21;
		          if (entType == ' ') entType = 'D';
		      } else if (!strcmp("  C", seName)) {
			  entMolWeight += 323.20;
		          if (entType == ' ') entType = 'D';
		      } else if (!strcmp("  G", seName)) {
			  entMolWeight += 363.22;
		          if (entType == ' ') entType = 'D';
		      } else if (!strcmp("  U", seName)) {
			  entMolWeight += 324.18;
			  entType = 'R';
		      } else if (!strncmp(" +", seName, 2)) {
			  entMolWeight = 0;
			  entType = 'N';
			  break;
		      } else if (!strcmp("UNK", seName)) {
			  entMolWeight = 0;
			  entType = 'U';
			  break;
		      }
		  }

		  if (entType == ' ') entType = 'C';
		  /*
		  nAlpha = 0; nBeta = 0;
		  
		  int iKS = 1, k_, *dkss; char *dks;
		  for (k = 0; k < nSeq; k++)
			if (ks[k] != 'H' && ks[k] != 'E') ks[k] = 'L';
		  for (k = 0; k < nSeq-1; k++) 
		    if (ks[k] != ks[k+1]) iKS++;
		  dks = new char [iKS+1];
		  dkss = new int [iKS]; 
		  iKS = 0; k_ = 0;
		  for (k = 0; k < nSeq; k++)
		    {
		      if (ks[k] == 'H') ssAlpha++;
		      if (ks[k] == 'E') ssBeta++;
		    }
		  for (k = 0; k < nSeq-1; k++) 
		    {
		      if (ks[k] != ks[k+1]) 
			{
			  if (ks[k] == 'H') nAlpha++;
			  if (ks[k] == 'E') nBeta++;
			  dks[iKS] = ks[k];
			  dkss[iKS] = k-k_;
			  k_ = k;
			  iKS++;
			}
		    }
		  
		  if (ks[nSeq-1] == 'H') nAlpha++;
		  if (ks[nSeq-1] == 'E') nBeta++;
		  dks[iKS] = ks[nSeq-1]; dks[iKS+1] = '\0';
		  dkss[iKS] = nSeq-k_; iKS++;
		  */
		  mw_enp.addItem((int4)entMolWeight);
		  type_enp.addItem((int1)entType);

		  //alpha_c_enp.addItem((flt4)(ssAlpha*100.0/nSeq));
		  //beta_c_enp.addItem((flt4)(ssBeta*100.0/nSeq));
		  //alpha_n_enp.addItem((int2)nAlpha);
		  //beta_n_enp.addItem((int2)nBeta);

		  //ss_seg_enp.addItem(dks);

		  //delete [] dks;
		  //delete [] dkss;
		}

	}
}
///////////////////////////////////////////////////////////////////////////
void DF::status_com(int log_mode) {
  
  if(!P_STATUS_COM) return;

  // masks: 0x0 - ok; 0x1 - no distr file; 0x2 - obsolete; 0x4 - updated;

  Property id_com("id.com"), file_com("file.com"), 
    status_com("status.com"), obs_dat_com("obs_dat.com");

  status_com.clear();

  int nCom=id_com.getObjectSize();
  
  for(int ic=0; ic<nCom; ic++) {
    int1 stat=0x0;
    
    char *id_com_ic=id_com.item1(ic);
    for(int jc=ic+1; jc<nCom; jc++) 
      if(!strcmp(id_com_ic, id_com.item1(jc))) {
	stat|=0x4;
	break;
      }
    
    char *cpath_file=file_com.item1(ic);
    char *cpath=cpath_file;
    char *cfile=strrchr(cpath_file, (int)(*PATH_SEPARATOR));
    if(cfile==NULL) {
      cpath="."; cfile=cpath_file;
    }
    else {
      *cfile='\0'; cfile++;
    }
    
    if(!DB::testPathFile(cpath, cfile)) stat|=0x1;
    
    if(*obs_dat_com.item4(ic)!=0) stat|=0x2;
    
    status_com.addItem(stat);
  }
  
}
//////////////////////////////////////////////////////////////////
