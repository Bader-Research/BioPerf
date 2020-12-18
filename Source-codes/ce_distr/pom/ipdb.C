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
// 10/02/98 I.Shindyalov The following modifications to improve 
// chain parsing have been made or copied from older version:
//    1. Treating N- C- terminal caps was introduced 
//       (in Molecule::matchAtomSeqres()).
//    2. Single residue mismatch after first 10 residues was allowed
//       (in Molecule::atchAtomSeqres()).
//    3. Parsing residues only based on resSeq but not resName to fix
//       HET groups insertion has been introduced (in Molecule::getMonomer()).
//    4. Nucl. acids 1-letter code assignment and monomer type 2 assignment
//       have been implemented (in Connectivity::search()).
//    5. Only nucl. acids, proteins with at least one canonical monomer and
//       number of monomers more than 1 are assigned to be polymer chains
//       (in Molecule::sortChains()).
//    6. Assignment of canonical monomer is changed to name-based and by 
//       presence of " CA " atom (in Connectivity::search()).
#include "ipdb.h"
#include "ipdb_df.h"
#include "pdbutil.h"
#include "config.h"
#include "miscutil.h"
#include <stdio.h>
#include <errno.h>
int errno;
#ifndef T3E
#include <sgtty.h>
#endif

//#include "perl_embedded.h"

//////////////////////////////////////////////////////////////////
//
// mkDB                                  - dialog mode
//
// mkDB update main_db distr tmp_db      - update db 
//
// mkDB scratch /scratch/s1/moose_output  
//      moose_pdb.pdb [1USR][is_add]     - create/add scratch db at  
//                                         /scratch/s1/moose_output from 
//                                         pdb-file moose_pdb.pdb with 
//                                         pdb-code 1USR
// 
//////////////////////////////////////////////////////////////////
template <class A> void AddToArray2(A ** array,
				    A value1,
				    A value2,
				    int Narray) {
  A * array_;
  array_=new A[Narray+2];
  for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
  array_[Narray]=value1;
  array_[Narray+1]=value2;
  if(Narray) delete [] (*array);
  (*array)=array_;
};
//////////////////////////////////////////////////////////////////
template <class A> void AddToArrayMULT(A ** array,
				       double * values,
				       int Narray,
				       int Nvalues,
				       int mult) {
  A * array_; int i;
  array_=new A[Narray+Nvalues];
  for(i=0; i<Narray; i++) array_[i]=*(*array+i);
  for(i=0; i<Nvalues; i++) array_[Narray+i]=(A)(values[i]*mult);
  if(Narray) delete [] (*array);
  (*array)=array_;
};
///////////////////////////////////////////////////////////////////////////
template <class A> void InsertItems(A **array,
				    int n,
				    A *vals,
				    int nVals,
				    int *nObj,
				    int diff) {

   (*nObj) += diff; A *newArray = new A[(*nObj)];
   int i, old, nw = 0;

   for (old = 0; old < n; old++, nw++) newArray[nw] = *(*array + old);

   for (i = 0; i < nVals; i++, nw++) newArray[nw] = vals[i];

   old = n + nVals - diff;
   for (; nw < (*nObj); nw++, old++) newArray[nw] = *(*array + old);

   if (nVals) delete [] (*array);
   (*array) = newArray;
};
////////////////////////////////////////////////////////////////////
template <class A> void AddToArray(A ** array,
				   A value,
				   int Narray,
				   int cacheSize = DEFAULT_CACHE) {

	if (leftOver(Narray, cacheSize)) {
		*(*array + Narray) = value;
	} else {
		A * array_;
		array_=new A[Narray+cacheSize+1];
		    
		for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
		array_[Narray]=value;
		if(Narray) delete [] (*array);
		(*array)=array_;
	}
};
///////////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
template <class A> void AddToArray(A ** array,
				   A value,
				   int2 Narray,
				   int cacheSize = DEFAULT_CACHE) {

	if (leftOver(Narray, cacheSize)) {
		*(*array + Narray) = value;
	} else {
		A * array_;
		array_=new A[Narray+cacheSize+1];
		    
		for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
		array_[Narray]=value;
		if(Narray) delete [] (*array);
		(*array)=array_;
	}
};
#endif
///////////////////////////////////////////////////////////////////////////
template <class A> void AddToArrayN(A ** array,
				    A * values,
				    int Narray,
				    int Nvalues,
				    int cacheSize = DEFAULT_CACHE) {
	if (Nvalues <= 0) return;

	if (leftOver(Narray, cacheSize) &&
	    Nvalues < (cacheSize - leftOver(Narray, cacheSize))) {
	    for (int i = 0; i < Nvalues; i++) *(*array+Narray+i) = values[i];
	} else {
	    A * array_; int i;
	    array_=new A[Narray+Nvalues+cacheSize];
	    for(i=0; i<Narray; i++) array_[i]=*(*array+i);
	    for(i=0; i<Nvalues; i++) array_[Narray+i]=values[i];
	    if(Narray) delete [] (*array);
	    (*array)=array_;
	}
};
//////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
template <class A> void AddToArrayN(A ** array,
				    A * values,
				    int2 Narray,
				    int2 Nvalues,
				    int cacheSize = DEFAULT_CACHE) {
	if (Nvalues <= 0) return;

	if (leftOver(Narray, cacheSize) &&
	    Nvalues < (cacheSize - leftOver(Narray, cacheSize))) {
	    for (int i = 0; i < Nvalues; i++) *(*array+Narray+i) = values[i];
	} else {
	    A * array_; int i;
	    array_=new A[Narray+Nvalues+cacheSize];
	    for(i=0; i<Narray; i++) array_[i]=*(*array+i);
	    for(i=0; i<Nvalues; i++) array_[Narray+i]=values[i];
	    if(Narray) delete [] (*array);
	    (*array)=array_;
	}
};
#endif
//////////////////////////////////////////////////////////////////
template <class A> void ExtendArray(A ** array,
				    int NarrayOld,
				    int NarrayNew) {
	A * array_;
	array_=new A[NarrayNew];
	for(int i=0; i<NarrayOld; i++) array_[i]=*(*array+i);
	if(NarrayOld) delete [] (*array);
	(*array)=array_;
};
//////////////////////////////////////////////////////////////////
template <class A> void ExtendArrayValue(A ** array,
					 int NarrayOld,
					 int NarrayNew,
					 A value) {
	int i;
	A * array_;
	array_=new A[NarrayNew];
	for(i=0; i<NarrayOld; i++) array_[i]=*(*array+i);
	for(i=NarrayOld; i<NarrayNew; i++) array_[i]=value;
	if(NarrayOld) delete [] (*array);
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
#ifndef INT2_OFF
template <class A> void DeleteArray(A ** array, int2 Narray) {
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};
#endif
//////////////////////////////////////////////////////////////////
Monomer::Monomer()
{
  nAtom=0;
  atomNames=NULL;
  pdbNum[0]='\0';
  nBond=0;
  bond=NULL;
  prevA=0;
  nextA=0;
  bondDist=2.1;
  pdbBreak=0;
  bondDist2=bondDist*bondDist;
}
//////////////////////////////////////////////////////////////////
void Monomer::makeBonds(Monomer * prevM, Monomer * nextM)
{
  int2 i, j; flt4 * xyzI, * xyzJ;
  for(i=0, xyzI=xyz; i<nAtom-1; i++, xyzI+=3)
    for(j=i+1, xyzJ=xyzI+3; j<nAtom; j++, xyzJ+=3)
      if(checkBond(xyzI, xyzJ))
        {
	  AddToArray2(&bond, (int2)(i+1), (int2)(j+1), nBond*2);
	  nBond++;
	}
  checkTriangles();

  if(!strcmp(name, "HOH")) return;

  if(prevM)
    {
      if(!strcmp(prevM->name, "HOH")) goto break_1;
      if((i=prevM->findAtom(" C  "))>=0 && (j=findAtom(" N  "))>=0)
	if(checkBond(prevM->xyz+i*3, xyz+j*3))
	  {
	    prevA=j+1;
	    goto break_1;
	  }

      for(i=0, xyzI=xyz; i<nAtom; i++, xyzI+=3)
	for(j=0, xyzJ=prevM->xyz; j<prevM->nAtom; j++, xyzJ+=3)
	  if(checkBond(xyzI, xyzJ))
	    {
	      prevA=i+1;
	      goto break_1;
	    }
    }

 break_1:
  if(nextM)
    {
      if(!strcmp(nextM->name, "HOH")) return;
      if((i=findAtom(" C  "))>=0 && (j=nextM->findAtom(" N  "))>=0)
	if(checkBond(xyz+i*3, nextM->xyz+j*3))
	  {
	    nextA=i+1;
	    return;
	  }

      for(i=0, xyzI=xyz; i<nAtom; i++, xyzI+=3)
	for(j=0, xyzJ=nextM->xyz; j<nextM->nAtom; j++, xyzJ+=3)
	  if(checkBond(xyzI, xyzJ))
	    {
	      nextA=i+1;
	      return;
	    }
    }
}
//////////////////////////////////////////////////////////////////
int Monomer::checkBond(flt4 * xyzI, flt4 * xyzJ)
{
  double dx, dy, dz;
  if((dx=fabs((double)*xyzI-(double)*xyzJ))>bondDist) return(0);
  if((dy=fabs((double)*(xyzI+1)-(double)*(xyzJ+1)))>bondDist) return(0);
  if((dz=fabs((double)*(xyzI+2)-(double)*(xyzJ+2)))>bondDist) return(0);
  if(dx*dx+dy*dy+dz*dz>bondDist2) return(0);
  return(1);
}
//////////////////////////////////////////////////////////////////
double Monomer::dist2(flt4 * xyzI, flt4 * xyzJ)
{
  double dx, dy, dz;
  dx=fabs((double)*xyzI-(double)*xyzJ);
  dy=fabs((double)*(xyzI+1)-(double)*(xyzJ+1));
  dz=fabs((double)*(xyzI+2)-(double)*(xyzJ+2));
  return(dx*dx+dy*dy+dz*dz);
}
//////////////////////////////////////////////////////////////////
void Monomer::removeBond(int iBond)
{
  int2 *bond_, i, j;
  bond_=new int2 [nBond*2-2];
  i=0;
  for(j=0; j<nBond; j++)
    {
      if(j==iBond) continue;
      bond_[i*2]=bond[j*2];
      bond_[i*2+1]=bond[j*2+1];
      i++;
    }
  delete [] bond;
  bond=bond_; nBond--;
}
//////////////////////////////////////////////////////////////////
void Monomer::checkTriangles()
{
// Checking for possible triangles
  int i, j, k, iAtom, jAtom, iAtom_, jAtom_;
  double Dk, Di, Dj;

  for(k=0; k<nBond*2-4; k+=2)
    {
      iAtom=bond[k];
      jAtom=bond[k+1];
      for(i=k+2; i<nBond*2; i++)
        {
	  if(bond[i]==iAtom)
	    {
	      iAtom_=bond[i%2?i-1:i+1];
	      for(j=k+2; j<nBond*2; j++)
	        {
	          if(bond[j]==jAtom)
	            {
	              jAtom_=bond[j%2?j-1:j+1];
	              if(iAtom_==jAtom_)
		        {
			  Dk=dist2(xyz+iAtom*3-3, xyz+jAtom*3-3);
			  Di=dist2(xyz+iAtom*3-3, xyz+iAtom_*3-3);
			  Dj=dist2(xyz+jAtom*3-3, xyz+jAtom_*3-3);
		          if(Dk>Di && Dk>Dj)
		            {
			      removeBond(k/2);
			    }
			  else
			    {
			      if(Di>Dj)
			        {
			          removeBond(i/2);
			        }
			      else
			        {
			          removeBond(j/2);
                                }
			    }
                        }
		    }
		}
            }
	}
    }
}
//////////////////////////////////////////////////////////////////
int Monomer::findAtom(char * atom)
{
  for(int i=0; i<nAtom; i++)
    if(!strncmp(atomNames+i*4, atom, 4)) return(i);
  return(-1);
}
//////////////////////////////////////////////////////////////////
Monomer::~Monomer()
{
  if(nAtom)
    {
      delete [] atomNames;
      delete [] xyz;
      delete [] tmf;
    }
  if(nBond) delete [] bond;
}
//////////////////////////////////////////////////////////////////
Molecule::Molecule(Connectivity * con_, Tlpdb * tlpdb_)
{
  con=con_;
  tlpdb=tlpdb_;
  db=tlpdb->db;
  con->db=db;
  nChain=0;
  nPolyChain=0;

#ifndef MINIMAL
  header=NULL;
  compnd=NULL;
  title=NULL;
  source=NULL;
  author=NULL;
  jrnl=NULL;
#endif

  seqNchain=0;
  
  strcpy(id, "XXXX");
  resolution=-1;

#ifndef MINIMAL
  expdta=9;
  expdta_txt = NULL;
  ssBond = NULL;
  site = NULL;
  ndbMap = NULL;
#endif

  relDate = 0;
  dateInt = 0;

  obsDate = 0;
  obsCount = 0;
  sprDate = 0;
  sprCount = 0;

  spcgrp = NULL;

  cell = new float[6];
  Zval = -1;
  rvalue = 0.;

  energyCutoff=-0.5;
  distanceCutoff=5.2;

  
  if(P_N_XYZ_ENC) n_xyz_enc.open("n_xyz.enc", 0);
  if(P_XYZ_ENC) xyz_enc.open("xyz.enc", 0);
  if(P_SE_XYZ_ENC) se_xyz_enc.open("se_xyz.enc", 0);
  if(P_XYZ_SE_ENC) xyz_se_enc.open("xyz_se.enc", 0);

  if(P_NAME_ENC) name_enc.open("name.enc", 0);
  if(P_NAME_ENP) name_enp.open("name.enp", 0);
  if(P_N_SE_ENC) n_se_enc.open("n_se.enc", 0);
  if(P_N_SE_ENP) n_se_enp.open("n_se.enp", 0);
  if(P_SE_ENC) se_enc.open("se.enc", 0);
  if(P_N_ENP_COM) n_enp_com.open("n_enp.com", 0);
  if(P_N_ENC_COM) n_enc_com.open("n_enc.com", 0);
  if(P_I_ENP_COM) i_enp_com.open("i_enp.com", 0);
  if(P_I_ENC_COM) i_enc_com.open("i_enc.com", 0);

  if(P_I_COM_ENC) i_com_enc.open("i_com.enc", 0);
  if(P_I_COM_ENP) i_com_enp.open("i_com.enp", 0);

  if(P_I_ENC_ENP) i_enc_enp.open("i_enc.enp", 0);
  if(P_I_ENP_ENC) i_enp_enc.open("i_enp.enc", 0);

#ifndef MINIMAL
  if(P_SEN_PDB_ENC) sen_pdb_enc.open("sen_pdb.enc", 0);
#endif

  if(P_BFAC_ENC) bfac_enc.open("bfac.enc", 0);

#ifndef MINIMAL
  if(P_EXP_ENP) exp_enp.open("exp.enp", 0);
  if(P_POL_ENP) pol_enp.open("pol.enp", 0);
#endif

  if(P_K_S_ENP) k_s_enp.open("k_s.enp", 0);
  if(P_SEQ_ENP) seq_enp.open("seq.enp", 0);

  if(P_RES_COM) res_com.open("res.com", 0);

#ifndef MINIMAL
  if(P_EXPDTA_COM) expdta_com.open("expdta.com", 0);
  if(P_EXPDTA_TXT_COM) expdta_txt_com.open("expdta_txt.com", 0);
  if(P_EXP_FLT_ENP) exp_flt_enp.open("exp_flt.enp", 0);
  if(P_POL_FLT_ENP) pol_flt_enp.open("pol_flt.enp", 0);
  if(P_SEQ_FLT_ENP) seq_flt_enp.open("seq_flt.enp", 0);

  if(P_ISEQ_FLT_ENP) iseq_flt_enp.open("iseq_flt.enp", 0);
  if(P_ISEQ_ENP) iseq_enp.open("iseq.enp", 0);
#endif

  if(P_SE_TYPE_ENP) se_type_enp.open("se_type.enp", 0);

  if(P_RELDAT_COM) rel_dat_com.open("reldat.com", 0);
  if(P_CURRENT_COM) current_com.open("current.com", 0);

  if(P_OBS_DAT_COM) obs_date_com.open("obs_dat.com", 0);
  if(P_OBS_IDS_COM) obs_ids_com.open("obs_ids.com", 0);
  if(P_OBS_COM) obs_com.open("obs.com", 0);
  if(P_SPR_COM) spr_date_com.open("spr_dat.com", 0);
  if(P_SPR_IDS_COM) spr_ids_com.open("spr_ids.com", 0);
  if(P_SPR_COM) spr_com.open("spr.com", 0);

  if(P_UNITCELL_COM) unit_cell_com.open("unitcell.com", 0);
  if(P_SPCGRP_COM) space_grp_com.open("spcgrp.com", 0);
  if(P_ZVAL_COM) z_val_com.open("zval.com", 0);
  if(P_RVAL_COM) rval_com.open("rval.com", 0);

  if(P_ID_COM) id_com.open("id.com", 0);
  if(P_FILE_COM) file_com.open("file.com", 0);
  if(P_DATE_TEX_COM) date_tex_com.open("date_tex.com", 0);
  if(P_DATE_INT_COM) date_int_com.open("date_int.com", 0);

#ifndef MINIMAL
  if(P_TITLE_COM) title_com.open("title.com", 0);
  if(P_COMPND_COM) compnd_com.open("compnd.com", 0);
  if(P_SOURCE_COM) source_com.open("source.com", 0);
  if(P_DATE_TEX_COM) date_tex_com.open("date_tex.com", 0);
  if(P_DATE_INT_COM) date_int_com.open("date_int.com", 0);
  if(P_HEADER_COM) header_com.open("header.com", 0);
  if(P_AUTH_COM) auth_com.open("auth.com", 0);
  if(P_JRNL_COM) jrnl_com.open("jrnl.com", 0);
  if(P_SITE_COM) site_com.open("site.com", 0);
  if(P_SSBOND_COM) ssbond_com.open("ssbond.com", 0);
  if(P_NDBMAP_COM) ndbmap_com.open("ndbmap.com", 0);
#endif
}
//////////////////////////////////////////////////////////////////
int Molecule::makeMolecule(char * fileName)
{
  file=fileName;

  int iXYZ, i, j;
  bufferStatus=0;
  modelStatus=0;
  chainStatus=0;

  if (obsCount) DeleteArray(&obsIds, obsCount); obsCount = 0;
  if (sprCount) DeleteArray(&sprIds, sprCount); sprCount = 0;

  obsDate = 0;
  sprDate = 0;
  relDate = 0;
  dateInt = 0;

  if(nChain)
    {
      for(i=0; i<nChain; i++)
	DeleteArray(&seNum[i], nMon[i]);
      delete [] seNum;
      delete [] nMon;
      DeleteArray(&mon, nChain);
      delete [] comNXYZ;
      DeleteArray(&comXYZ, nChain);
      DeleteArray(&comTmf, nChain);
      delete [] chainCode;

#ifndef MINIMAL
      if(header) delete [] header; header=NULL;
      if(compnd) delete [] compnd; compnd=NULL;
      if(title) delete [] title; title=NULL;
      if(source) delete [] source; source=NULL;
      if(author) delete [] author; author=NULL;
      if(ssBond) delete [] ssBond; ssBond = NULL;
      if(site) delete [] site; site = NULL;
      if(ndbMap) delete [] ndbMap; ndbMap = NULL;
#endif

      if(spcgrp) delete [] spcgrp; spcgrp = NULL;

      DeleteArray(&ks, nChain);

#ifndef MINIMAL
      if(jrnl) delete [] jrnl; jrnl=NULL;
      DeleteArray(&exp, nChain);
      DeleteArray(&pol, nChain);
#endif

      DeleteArray(&monXYZ, nChain);

    }
  nChain=0;
  nPolyChain=0;
  resolution=-1;

#ifndef MINIMAL
  expdta=9;
  if (expdta_txt) delete [] expdta_txt;
  expdta_txt=NULL;
#endif

  if(seqNchain)
    {
      for(i=0; i<seqNchain; i++) DeleteArray(&seq[i], seqL[i]);
      delete [] seq;
      delete [] seqL;
      delete [] seqChainCode;
    }
  seqNchain=0;
  if(!(pdbFile=fopen(fileName,"r")))
    {
      fprintf(stderr, "****** Can't open %s: ", fileName);
      perror("");
      return(errno);
    }

  Monomer * M=NULL, * prevM=NULL, * nextM=NULL;
  int addChainStatus;

  M=getMonomer();
  while(M)
    {
      nextM=getMonomer();
      M->makeBonds(prevM, nextM);

      addChainStatus=0;
      if(!nChain) addChainStatus=1;
      if(prevM)
	{
	  if(prevM->pdbChain!=M->pdbChain) addChainStatus=1;
	}
      if(!M->prevA) addChainStatus=1;
      if(M->pdbBreak) addChainStatus=1;
      if(addChainStatus) addChain(M);

      M->index=con->search(M);
      addMonomer(M);

      if(prevM) delete prevM;
      prevM=M;
      M=nextM;
    }

  fclose(pdbFile);

  if(prevM) delete prevM;
  if(nextM) delete nextM;
  if(M) delete M;

  matchAtomSeqres();
  sortChains();

  monXYZ=new int2* [nChain];
  for(i=0; i<nChain; i++)
    {
      monXYZ[i]=new int2 [nMon[i]+1];
      iXYZ=0;
      for(j=0; j<nMon[i]; j++)
	{
	 *(monXYZ[i]+j)=iXYZ;
	 iXYZ+=*con->monNatom.item2(*(mon[i]+j));
	}
      *(monXYZ[i]+nMon[i])=iXYZ;
    }

  ks=new char * [nChain];

#ifndef MINIMAL
  exp=new flt4 * [nChain];
  pol=new flt4 * [nChain];
#endif

  int iE, iSE, NSE;
  for(iE=0; iE<nChain; iE++)
    {
      NSE=nMon[iE];
      ks[iE]=new char [NSE+1];
      *(ks[iE]+NSE)='\0';

#ifndef MINIMAL
      exp[iE]=new flt4 [NSE];
      pol[iE]=new flt4 [NSE];
#endif

      for(iSE=0; iSE<NSE; iSE++)
	{

#ifndef MINIMAL
	  *(*(exp+iE)+iSE)=0.0;
	  *(*(pol+iE)+iSE)=0.0;
#endif

	  *(*(ks+iE)+iSE)=' ';
	}
    }
  if(tlpdb->cboxes[5] && P_K_S_ENP) assignKS();

#ifndef MINIMAL
  if(tlpdb->cboxes[6] && (P_EXP_ENP || P_EXP_FLT_ENP || P_POL_ENP || 
  P_POL_FLT_ENP || P_EXP_1B_ENP || P_POL_1B_ENP ||
  P_EXP_IP6 || P_POL_IP6)) assignEnv();
#endif

  if(P_SEQ_FLT_ENP || P_ISEQ_ENP || P_ISEQ_FLT_ENP || P_EXP_FLT_ENP || 
     P_POL_FLT_ENP || P_SE_TYPE_ENP) assignFiltProp();
      
#ifdef PURIFY
	if (purify_clear_leaks()) {
		purify_new_inuse();
		fprintf(stderr, "purify found leaks in makeMolecule!\n");
	}
#endif

  return(0);
}
//////////////////////////////////////////////////////////////////
void Molecule::writeDB()
{
  int err, i, j, k;
  
  if(!db->testFile("id.com"))
    {
      printf("Molecule::writeDB - db->testFile error\n"); exit(0);
    }


  if(P_UNITCELL_COM) unit_cell_com.addItem(cell);
  for (i = 0; i < 6; i++) cell[i] = 0;

  if (spcgrp) compressText(&spcgrp);

  if(P_SPCGRP_COM) space_grp_com.addItem(spcgrp);

  if(P_RELDAT_COM) rel_dat_com.addItem(relDate);

  if(P_ZVAL_COM) z_val_com.addItem(Zval); Zval = -1;
  if(P_RVAL_COM) rval_com.addItem(rvalue); rvalue = 0.;

  if (P_CURRENT_COM) current_com.addItem((int4) -1);

  if(P_OBS_COM) obs_com.addItem((int4) -1);
  if(P_OBS_DAT_COM) obs_date_com.addItem(obsDate);
  if (obsCount) {
	char *obsIdsTmp;
	StringArrayToBuffer(obsCount, &obsIdsTmp, &obsIds);
	if(P_OBS_IDS_COM) obs_ids_com.addItem(obsIdsTmp);
	delete [] obsIdsTmp;
  } else {
	if(P_OBS_IDS_COM) obs_ids_com.addItem("\0");
  }

  if(P_SPR_COM) spr_com.addItem((int4) -1);
  if(P_SPR_DAT_COM) spr_date_com.addItem(sprDate);
  if (sprCount) {
	char *sprIdsTmp;
	StringArrayToBuffer(sprCount, &sprIdsTmp, &sprIds);
	if(P_SPR_IDS_COM) spr_ids_com.addItem(sprIdsTmp);
	delete [] sprIdsTmp;
  } else {
	if(P_SPR_IDS_COM) spr_ids_com.addItem("\0");
  }

  for(i=0; i<nChain; i++) 
    if(P_N_XYZ_ENC) n_xyz_enc.addItem(comNXYZ[i]);
  for(i=0; i<nChain; i++) 
    if(P_XYZ_ENC) xyz_enc.addItem(comXYZ[i], comNXYZ[i]*3);
      
  if(P_ID_COM) id_com.addItem(id);
  if(P_FILE_COM) file_com.addItem(file);

  if(P_DATE_TEX_COM) date_tex_com.addItem(date);
  if(P_DATE_INT_COM) date_int_com.addItem(dateInt);

#ifndef MINIMAL
  compressText(&title);
  if(P_TITLE_COM) title_com.addItem(title);

  compressText(&compnd);
  if(P_COMPND_COM) compnd_com.addItem(compnd);

  compressText(&source);
  if(P_SOURCE_COM) source_com.addItem(source);

  if(P_HEADER_COM) header_com.addItem(header);

  compressText(&author);
  if(P_AUTH_COM) auth_com.addItem(author);

  compressText(&jrnl);
  if(P_JRNL_COM) jrnl_com.addItem(jrnl);

  if(P_SITE_COM) site_com.addItem(site);
  if(P_SSBOND_COM) ssbond_com.addItem(ssBond);

  compressText(&ndbMap);
  if(P_NDBMAP_COM) ndbmap_com.addItem(ndbMap);
#endif

  int4 nCom=n_enp_com.getObjectSize();
  int4 nEnp=n_se_enp.getObjectSize();
  int4 nEnc=n_se_enc.getObjectSize();

  //int2 nPolyChain=0;

  for(i=0; i<nChain; i++) 
    if(P_N_SE_ENC) n_se_enc.addItem(nMon[i]);
  for(i=0; i<nChain; i++) 
    if(P_SE_ENC) se_enc.addItem(mon[i], nMon[i]);
  for(i=0; i<nChain; i++) 
    if(P_I_COM_ENC) i_com_enc.addItem(nCom);

  for(i=0; i<nChain; i++) 
    if(P_I_ENP_ENC) i_enp_enc.addItem((int4)(i<nPolyChain?nEnp+i:-1));
  /*
  for(i=0; i<nChain; i++) {
    if(nMon[i]>1) {
      nPolyChain=i+1;
      //if(P_N_SE_ENP) n_se_enp.addItem(nMon[i]);
    }
  }
  */
  for(i=0; i<nPolyChain; i++) 
    if(P_N_SE_ENP) n_se_enp.addItem(nMon[i]);

  for(i=0; i<nPolyChain; i++) 
    if(P_I_COM_ENP) i_com_enp.addItem(nCom);

  for(i=0; i<nPolyChain; i++) 
    if(P_I_ENC_ENP) i_enc_enp.addItem((int4)(nEnc+i));
     
  if(P_I_ENP_COM) i_enp_com.addItem(n_se_enp.getObjectSize());
  if(P_I_ENC_COM) i_enc_com.addItem(n_se_enc.getObjectSize());

  if(P_N_ENC_COM) n_enc_com.addItem(nChain);
  if(P_N_ENP_COM) n_enp_com.addItem(nPolyChain);

  for(i=0; i<nChain; i++)
    {
      if(P_BFAC_ENC) bfac_enc.addItem(comTmf[i], comNXYZ[i]);
    }

#ifndef MINIMAL
  for(i=0; i<nPolyChain; i++) 
    if(P_EXP_ENP) exp_enp.addItem(exp[i], nMon[i]);
  for(i=0; i<nPolyChain; i++) 
    if(P_POL_ENP) pol_enp.addItem(pol[i], nMon[i]);
#endif

  for(i=0; i<nPolyChain; i++) 
    if(P_K_S_ENP) k_s_enp.addItem(ks[i]);

  char *comEnt_, **comEnt=new char* [nChain], *comSeq;
  int4 *comSExyz, *comXYZse;

  for(i=0; i<nChain; i++)
    {
      comEnt[i]=new char [14];

      if(i < nPolyChain)
	{
	  if(chainCode[i]==' ')
	    {
	      sprintf(comEnt[i],"%s:_", id);
	    }
	  else
	    {
	      sprintf(comEnt[i],"%s:%c", id, chainCode[i]);
	    }
	}
      else
	{
	  sprintf(comEnt[i],"%s:%s", id, con->monCode3.item1(*(mon[i])));
	}

      int nSameEnt=0, entNamLen=strlen(comEnt[i]);
      for(j=0; j<i; j++)
	if(!strncmp(comEnt[j], comEnt[i], entNamLen)) nSameEnt++;
      
      if(nSameEnt) sprintf(comEnt[i],"%s:%d", comEnt[i], nSameEnt);
      
      if(P_NAME_ENC) name_enc.addItem(comEnt[i]);
      
      if(i < nPolyChain)
	{
	  if(P_NAME_ENP) name_enp.addItem(comEnt[i]);

	  comSeq=new char[nMon[i]+1]; 
	  *(comSeq+nMon[i])='\0';
	  for(j=0; j<nMon[i]; j++)
	    {
	      *(comSeq+j)=*con->monCode1.item1(*(mon[i]+j));
	    }
	  
	  if(P_SEQ_ENP) seq_enp.addItem(comSeq);
	  delete [] comSeq;
	}
      
      comSExyz=new int4 [nMon[i]+1];
      comXYZse=new int4 [comNXYZ[i]];
      int iXYZ=0, nAtom;
      for(j=0; j<nMon[i]; j++)
	{
	  comSExyz[j]=iXYZ;
	  nAtom=*con->monNatom.item2(*(mon[i]+j));
	  for(k=iXYZ; k<iXYZ+nAtom; k++) comXYZse[k]=j;
	  iXYZ+=nAtom;
	}
      comSExyz[nMon[i]]=iXYZ;

      if(P_SE_XYZ_ENC) se_xyz_enc.addItem(comSExyz, nMon[i]+1);
      if(P_XYZ_SE_ENC) xyz_se_enc.addItem(comXYZse, comNXYZ[i]);

      delete [] comSExyz;
      delete [] comXYZse;

#ifndef MINIMAL
      char *sen_tmp; 
      StringArrayToBuffer(nMon[i], &sen_tmp, &seNum[i]);
      if(P_SEN_PDB_ENC) sen_pdb_enc.addItem(sen_tmp);
      delete [] sen_tmp;
#endif

    }

  DeleteArray(&comEnt, nChain);

  if(P_RES_COM) res_com.addItem(resolution);

#ifndef MINIMAL
  if(P_EXPDTA_COM) expdta_com.addItem(expdta);
  if(P_EXPDTA_TXT_COM) {
	compressText(&expdta_txt);
	if (!expdta_txt) expdta_txt_com.addItem("X-RAY DIFFRACTION");
	else expdta_txt_com.addItem(expdta_txt);
  }
#endif
      
#ifdef PURIFY
	if (purify_clear_leaks()) {
		purify_new_inuse();
		fprintf(stderr, "purify found leaks in writeDB!\n");
	}
#endif

/*
  FILE * CONN;

  CONN=fopen("db_mol","a");

  int iEnt, i;

  fprintf(CONN, "\n\n---------- %s ----------\n", id);
  fprintf(CONN, "\n %s \n", compnd);
  fprintf(CONN, "\n %s \n", source);

  for(iEnt=0; iEnt<nChain; iEnt++)
    {
      fprintf(CONN, "\n---------- Entity #%d ----------\n",iEnt);
      for(i=0; i<nMon[iEnt]; i++)
        {
	  if(i%20==0)fprintf(CONN, "\n");
	  fprintf(CONN, "%4d", *(mon[iEnt]+i));
	}

      fprintf(CONN, "\n---------- property 1 ----------\n");
      for(i=0; i<seqsL; i++)
        {
	  if(i%20==0)fprintf(CONN, "\n");
	  fprintf(CONN, "%4d",prop1[i]/10);
	}

      fprintf(CONN, "\n---------- coords %d----------------------\n",
	comNXYZ[iEnt]);

      for(i=0; i<comNXYZ[iEnt]*3; i++)
        {
	  if(i%9==0)fprintf(CONN, "\n");
	  if(i%3==0)fprintf(CONN, "   ");
	  fprintf(CONN, "%.3f",*(comXYZ[iEnt]+i));
	}
    }

  fclose(CONN);
*/

 // test for extra poly-chains and missing pdb chains
	/*
   for(int ic1 = 0; ic1 < nPolyChain - 1; ic1++)
     for(int ic2 = ic1 + 1; ic2 < nPolyChain; ic2++) {
       if(chainCode[ic1] == chainCode[ic2]) 
	 printf("chain test: in %s '%c' used again\n", id, chainCode[ic1]);
     }

   for(int ic1 = 0; ic1 < nPolyChain - 1; ic1++) {
     int isInSeqres = 0;
     for(int ic2 = 0; ic2 < seqNchain; ic2++)
       if(chainCode[ic1] == seqChainCode[ic2]) {
	 isInSeqres = 1; break;
       }
     if(!isInSeqres) 
       printf("chain test: in %s '%c' not in SEQRES\n", id, chainCode[ic1]);
   }

   for(int ic2 = 0; ic2 < seqNchain; ic2++) {
     int isInAtom = 0;
     for(int ic1 = 0; ic1 < nPolyChain; ic1++) 
       if(chainCode[ic1] == seqChainCode[ic2]) {
	 isInAtom = 1; break;
       }
     if(!isInAtom) 
       printf("chain test: in %s '%c' is not a polymer\n", id, seqChainCode[ic2]);
   }
	*/
}
//////////////////////////////////////////////////////////////////
void Molecule::addChain(Monomer * M)
{
  AddToArray(&nMon, (int2)0, nChain);
  AddToArray(&mon, (int2 *)NULL, nChain);
  AddToArray(&seNum, (char **)NULL, nChain);
  AddToArray(&comNXYZ, (int4)0, nChain);
  AddToArray(&comXYZ, (flt4 *)NULL, nChain);
  AddToArray(&comTmf, (flt4 *)NULL, nChain);
  AddToArray(&chainCode, M->pdbChain, nChain);
  nChain++;
}
//////////////////////////////////////////////////////////////////
void Molecule::addMonomer(Monomer * M)
{
  char buf[6], *seNum_; int status;
  status=sscanf(M->pdbNum, "%s", buf);
  if(status!=0 && status!=EOF)
    {
      seNum_=new char [strlen(buf)+1];
      strcpy(seNum_, buf);
    }
  else
    {
      seNum_=new char [1]; seNum_[0]='\0';
    }
  AddToArray(&seNum[nChain-1], seNum_, nMon[nChain-1]);
  AddToArray(&mon[nChain-1], (int2)M->index, nMon[nChain-1]);
  AddToArrayN(&comXYZ[nChain-1], M->xyz, comNXYZ[nChain-1]*3, M->nAtom*3);
  AddToArrayN(&comTmf[nChain-1], M->tmf, comNXYZ[nChain-1], (int) M->nAtom);
  nMon[nChain-1]++;
  comNXYZ[nChain-1]+=M->nAtom;
}
//////////////////////////////////////////////////////////////////
char *expdtaValues[]={"ELECTRON", 
		      "FIBER DI", 
		      "FLUORESC",
		      "NEUTRON",
		      "NMR",
		      "NMR SOLI",
		      "NMR SOLU",
		      "SYNCHROT",
		      "THEORETI",
		      "X-RAY DI",
		      NULL};
//////////////////////////////////////////////////////////////////
Monomer * Molecule::getMonomer()
{
  Monomer * M=NULL;
  flt4 xyz_[3], tmf_;
  char atomNames_[5];
  char tmp[10];
  int missingSome = 0, l;
  int monomerStatus=0;
  char *rem3 = NULL;
 read_1:
  if(modelStatus==2) return(NULL);
  if(!bufferStatus)
    if(!fgets(buffer, 100, pdbFile))
      {
	return(M);
      }
    else
      {
	if(buffer[strlen(buffer)-1]<32) buffer[strlen(buffer)-1]='\0';
      }
  bufferStatus=0;
  if(!strncmp("ENDMDL", buffer, 6))
    {
      if(!tlpdb->cboxes[3]) modelStatus=2;
      return(M);
    }
  if(!strncmp("MODEL ", buffer, 6) && !tlpdb->cboxes[3])
    {
      // only increment modelStatus if we had a non-empty model
      if (nChain) modelStatus++;
      if(modelStatus==2) return(M);
    }
  if(!strncmp("HEADER", buffer, 6)) addHeader(buffer);

#ifndef MINIMAL
  if(!strncmp("COMPND", buffer, 6)) addCompnd(buffer);
  if(!strncmp("SOURCE", buffer, 6)) addSource(buffer);
  if(!strncmp("AUTHOR", buffer, 6)) addAuthor(buffer);
  if(!strncmp("JRNL  ", buffer, 6)) addJrnl(buffer);
#endif

  if(!strncmp("SEQRES", buffer, 6)) addSeqres(buffer);
  if(!strncmp("CRYST1", buffer, 6)) addCryst(buffer);
  if(!strncmp("MODEL ", buffer, 6)) chainStatus=1;

  if (!strncmp("OBSLTE", buffer, 6)) addObslte(buffer);
  if (!strncmp("SPRSDE", buffer, 6)) addSprsde(buffer);
  if (!strncmp("REVDAT   1", buffer, 10)) addRelDate(buffer);

  if (!strncmp("SSBOND", buffer, 6)) {
	buffer[16] =
	buffer[21] =
	buffer[30] =
	buffer[35] = '\0';
	if (buffer[15] != ' ') addText(&ssBond, &buffer[15]);
	sprintf(tmp,"%d", atoi(&buffer[17]));
	addText(&ssBond, tmp);
	addText(&ssBond, ":");
	if (buffer[29] != ' ') addText(&ssBond, &buffer[29]);
	sprintf(tmp,"%d", atoi(&buffer[31]));
	addText(&ssBond, tmp);
	addText(&ssBond, ";");
  }

  if(!strncmp("DBREF", buffer, 5) && !ndbMap) {
	if (!strncmp("NDB", &buffer[26], 3)) {
		buffer[42] = '\0';
		addText(&ndbMap, &buffer[29]);
	}
  }
	
  if(!strncmp("SITE", buffer, 4)) {
	if (strlen(buffer) < 80)
	    for (int j = strlen(buffer); j < 79; j++)
		buffer[j] = ' ';

	buffer[14] =
	buffer[17] =
	buffer[21] =
	buffer[32] =
	buffer[43] =
	buffer[54] =
	buffer[61] = '\0';

	if (!missingSome) {
	    addText(&site, &buffer[11]);               // name
	    addText(&site, "=");
	    missingSome = atoi(&buffer[15]);            // numRes
	}

	for (int j = 0; j < 4; j++, missingSome--) {
	    if (!missingSome) break;

	    addText(&site, &buffer[18+j*11]);          // resName
	    if (buffer[22+j*11] != ' ') {
		l = strlen(site);
		addText(&site, " ");
		site[l] = buffer[22+j*11];           // chainId
		site[l+1] = '\0';
	    }

	    sprintf(tmp, "%d", atoi(&buffer[23+j*11]));
	    addText(&site, tmp);                       // seqNum

	    if (buffer[27+j*11] != ' ') {
		l = strlen(site);
		addText(&site, " ");
		site[l] = buffer[27+j*11];           // insCode ???
		site[l+1] = '\0';
	    }

	    if (missingSome != 1) addText(&site, ":");
	}

	if (!missingSome) addText(&site, ";");
  }

#ifndef MINIMAL
  if (!strncmp("TITLE ", buffer, 6)) {
                    buffer[71]='\0';
                    addText(&title, &buffer[10]);
  }

  if(!strncmp("EXPDTA", buffer, 6)) {
      buffer[71]='\0';
      addText(&expdta_txt, &buffer[10]);
      for(int i=0; expdtaValues[i]; i++)
	  if(!strncmp(expdtaValues[i], &buffer[10], strlen(expdtaValues[i]))) {
	      expdta=i; break;
	  }
  }
#endif

  if(!strncmp("REMARK   2", buffer, 10) && resolution==-1)
    {
      for(int i=11; i<60; i++)
	{
	  if(!strncmp("RESOLUTION", &buffer[i], 10))
	    {
	      resolution=atof(&buffer[i+11]);
	      break;
	    }
        }
    }
  /*  
  if (!strncmp("REMARK   3", buffer, 10)) {

      buffer[70] = '\0';
      addText(&rem3, &buffer[10]);

  } else if (rem3) {

	SV *tmp_text = newSV(0);
	AV *rval;
	Perl::init();
	sv_setpv(tmp_text, rem3);
	Perl::substitute(&tmp_text, "s/[ ]+/ /mg");
	Perl::substitute(&tmp_text, "s/^ //mg");

	if (Perl::matches(tmp_text, "/R VALUE.*?:\\s*([0-9\\.]+)/", &rval)) {
		rvalue = (flt4) SvNV(*av_fetch(rval, 0, FALSE));
	}

	if (rvalue <= 0.) {
	    Perl::substitute(&tmp_text, "s/\\n/ /mg");
	    if (Perl::matches(tmp_text,
			      "/R[- ]?VALUE I?S?\\s+([0-9\\.]+)/",
			      &rval))
		rvalue = (flt4) SvNV(*av_fetch(rval, 0, FALSE));

	    else if (Perl::matches(tmp_text,
			      "/R VALUE\\s+([0-9\\.]+)/",
			      &rval))
		rvalue = (flt4) SvNV(*av_fetch(rval, 0, FALSE));

	    else if (Perl::matches(tmp_text,
			      "/R VALUE\\s+\\d\\s+([0-9\\.]+)/",
			      &rval))
		rvalue = (flt4) SvNV(*av_fetch(rval, 0, FALSE));
	}

	if (rvalue < 0. || rvalue >= 1.) rvalue = 0.;
	delete [] rem3; rem3 = NULL;
	SvREFCNT_dec(tmp_text);

  }
  */
  if(!strncmp("ATOM  ", buffer, 6) ||
     !strncmp("HETATM", buffer, 6))
    {
      if(strlen(buffer)<55) goto read_1;
      if(buffer[16]!=' ' && buffer[16]!='A' && buffer[16]!='1') goto read_1;
      if(buffer[13]=='H') goto read_1;
      if(buffer[13]=='Q') goto read_1;
      if(!strncmp("EXC", &buffer[17], 3)) goto read_1;
      if(!monomerStatus)
	{
	  M=new Monomer;
	  strncpy(M->name, &buffer[17], 3);
	  M->name[3]='\0';
	  M->pdbChain=buffer[21];
	  M->pdbNum[5]='\0';
	  strncpy(M->pdbNum, &buffer[22], 5);
	  monomerStatus=1;
	}
      else
	{
	  if(//strncmp(M->name, &buffer[17], 3) ||
	    M->pdbChain!=buffer[21] ||
	    strncmp(M->pdbNum, &buffer[22], 5)) goto exit_1;
	}
      buffer[66]='\0'; tmf_=atof(&buffer[60]);
      buffer[54]='\0'; xyz_[2]=atof(&buffer[46]);
      buffer[46]='\0'; xyz_[1]=atof(&buffer[38]);
      buffer[38]='\0'; xyz_[0]=atof(&buffer[30]);
      atomNames_[4]='\0';
      if(buffer[12]=='\0') buffer[12]=' ';
      if(buffer[13]=='\0') buffer[13]=' ';
      if(buffer[14]=='\0') buffer[14]=' ';
      if(buffer[15]=='\0') buffer[15]=' ';
      strncpy(atomNames_, &buffer[12], 4);
      addText(&M->atomNames, atomNames_);
      AddToArrayN(&M->xyz, xyz_, M->nAtom*3, 3);
      AddToArray(&M->tmf, tmf_, M->nAtom);
      M->nAtom++;
    }
  goto read_1;
 exit_1:
  bufferStatus=1;
      
#ifdef PURIFY
	if (purify_clear_leaks()) {
		purify_new_inuse();
		fprintf(stderr, "purify found leaks in getMonomer!\n");
	}
#endif

  return(M);
}
//////////////////////////////////////////////////////////////////
void Molecule::addHeader(char * buffer)
{
  if(tlpdb->scratch_mode==0) {
    buffer[66]='\0';
    strcpy(id, &buffer[62]);
  }
  buffer[59]='\0';
  strcpy(date, &buffer[50]);

  dateInt = toJday(date);
  // fprintf(stderr, "%s -> %d\n", date, dateInt);

#ifndef MINIMAL
  buffer[50]='\0';
  addText(&header, &buffer[10]);
#endif

}
//////////////////////////////////////////////////////////////////
void Molecule::addCryst(char *buffer) {

	// grab the unit cell params
	char *tmp = new char[10];
	tmp[9] = '\0'; int count = 0;

	int n;

	// first the dims a, b, c
	for (n = 0; n < 3; n++) {
		strncpy(tmp, &buffer[6+n*9], 9);
		cell[count++] = atof(tmp);
	}

	tmp[7] = '\0';

	// now the angles alpha, beta, gamma
	for (n = 0; n < 3; n++) {
		strncpy(tmp, &buffer[33+n*7], 7);
		cell[count++] = atof(tmp);
	}

	// now we get the space group
	char tmp2 = buffer[66];
	buffer[66] = '\0';
	addText(&spcgrp, &buffer[55]);
	buffer[66] = tmp2;

	// and last not least the Z value
	buffer[70] = '\0';
	Zval = (int2) atoi(&buffer[66]);
	delete [] tmp;
}
//////////////////////////////////////////////////////////////////
void Molecule::addObslte(char *buffer) {

	if (!obsDate) {
		buffer[20] = '\0';
		obsDate = toJday(&buffer[11]);
	}

	char *tmpId = new char[5];
	for (int n = 0; n < 8; n++) {
		if (buffer[31+n*5] > 48) {
			buffer[35+n*5] = '\0';
			strcpy(tmpId, &buffer[31+n*5]);
			AddBuffer(&obsIds, tmpId, obsCount);
			obsCount++;
		} else {
			break;
		}
	}
	delete [] tmpId;
}
//////////////////////////////////////////////////////////////////
void Molecule::addSprsde(char *buffer) {

	if (!sprDate) {
		buffer[20] = '\0';
		sprDate = toJday(&buffer[11]);
	}

	char *tmpId = new char[5];
	for (int n = 0; n < 8; n++) {
		if (buffer[31+n*5] > 48) {
			buffer[35+n*5] = '\0';
			strcpy(tmpId, &buffer[31+n*5]);
			AddBuffer(&sprIds, tmpId, sprCount);
			sprCount++;
		} else {
			break;
		}
	}
	delete [] tmpId;
}
//////////////////////////////////////////////////////////////////
void Molecule::addRelDate(char *buffer) {
	buffer[22] = '\0';
	relDate = toJday(&buffer[13]);
}
//////////////////////////////////////////////////////////////////
void Molecule::addCompnd(char * buffer)
{
  buffer[71]='\0';
  addText(&compnd, &buffer[10]);
}
//////////////////////////////////////////////////////////////////
void Molecule::addSource(char * buffer)
{
  buffer[71]='\0';
  addText(&source, &buffer[10]);
}
//////////////////////////////////////////////////////////////////
void Molecule::addAuthor(char * buffer)
{
  buffer[71]='\0';
  addText(&author, &buffer[10]);
}
//////////////////////////////////////////////////////////////////
void Molecule::addJrnl(char * buffer)
{
  buffer[71]='\0';
  addText(&jrnl, &buffer[19]);
}
//////////////////////////////////////////////////////////////////
void Molecule::addSeqres(char * buffer)
{
  int i, addChainStatus=0;
  char * seq_;
  if(!seqNchain) addChainStatus=1;
  if(seqNchain)
    if(seqChainCode[seqNchain-1]!=buffer[11]) addChainStatus=1;
  if(addChainStatus)
    {
      AddToArray(&seq, (char **)NULL, seqNchain);
      AddToArray(&seqL, 0, seqNchain);
      AddToArray(&seqChainCode, buffer[11], seqNchain);
      seqNchain++;
    }

  for(i=0; i<13; i++)
    {
      buffer[10] = '\0';
      buffer[17] = '\0';
      if(atoi(&buffer[8]) == 0) {
	for(int is = 0; is < atoi(&buffer[13]); is++) {
	  seq_=new char [4];
	  strcpy(seq_, "UNK");
	  AddToArray(&seq[seqNchain-1], seq_, seqL[seqNchain-1]);
	  seqL[seqNchain-1]++;
	}
	break;
      }

      if (strlen(&buffer[19+4*i]) < 3) break;

      if(strncmp("   ", &buffer[19+4*i], 3) &&
	 strncmp("EXC", &buffer[19+4*i], 3))
	{
	  seq_=new char [4];
	  seq_[3]='\0';
	  strncpy(seq_, &buffer[19+4*i], 3);
	  AddToArray(&seq[seqNchain-1], seq_, seqL[seqNchain-1]);
	  seqL[seqNchain-1]++;
      
#ifdef PURIFY
	if (purify_clear_leaks()) {
		purify_new_inuse();
		fprintf(stderr, "purify found leaks in addSeqres!\n");
		fprintf(stderr, "addChainStatus: %d\n", addChainStatus);
		fprintf(stderr, "seqNchain: %d\n", seqNchain);
		fprintf(stderr, "buffer:\n%s\n", buffer);
		fprintf(stderr, "seqChainCode[%d]: \"%c\"\n",
						seqNchain-1,
						seqChainCode[seqNchain-1]);
		fprintf(stderr, "seq position: %d\n", i);
		fprintf(stderr, "seq length: %d\n", seqL[seqNchain-1]);
	}
#endif

        }
      
    }

}
//////////////////////////////////////////////////////////////////
/*
void Molecule::matchAtomSeqres()
{
  int i, j, k1, k2, *seqMatch, matchIndex, matchStatus, matchChain,
    *chainStatus, nChain_, firstChain; int2 index; 

  int2 *mon_, nMon_, NXYZ_; flt4 *xyz_, *tmf_; char **seNum_;
  int2 **_mon_, *_nMon_, *_NXYZ_; flt4 **_tmf_, **_xyz_;
  char *_chainCode_, ***_seNum_;  char *_seNum;
  Monomer *M=new Monomer;

  if(nChain) chainStatus=new int [nChain];
  for(j=0; j<nChain; j++) chainStatus[j]=0;

  for(i=0; i<seqNchain; i++)
    {
      seqMatch=new int [seqL[i]];
      matchIndex=0;
      matchStatus=0;

      nMon_=0; NXYZ_=0;
      for(j=0; j<seqL[i]; j++) seqMatch[j]=-1;
      for(j=0; j<nChain; j++)
	{
	  if(chainStatus[j]) continue;
	  if(seqChainCode[i]!=chainCode[j]) continue;
	  for(k1=matchIndex; k1<seqL[i] && k1+nMon[j]<=seqL[i]; k1++)
            {
	      for(k2=0; k2<nMon[j]; k2++)
		if(strcmp(*(seq[i]+k1+k2), con->monCode3.item1(*(mon[j]+k2))))
		  goto break_1;
	      matchStatus++;
	      for(k2=0; k2<nMon[j]; k2++) seqMatch[k1+k2]=j;
	      matchIndex=k1+nMon[j];
	      chainStatus[j]=1;
	      break;
	    break_1: ;
            }
        }

      if(matchStatus)
	{
	  firstChain=-1;
	  for(j=0; j<seqL[i]; )
	    if(seqMatch[j]==-1)
	      {
		strcpy(M->name, *(seq[i]+j));
		index=con->search(M);
		AddToArray(&mon_, index, nMon_);
		_seNum=new char [1]; _seNum[0]='\0';
		AddToArray(&seNum_, _seNum, nMon_);
		nMon_++; j++;
	      }
	    else
	      {
		matchChain=seqMatch[j];
		if(firstChain==-1) firstChain=matchChain;
		AddToArrayN(&mon_, mon[matchChain], nMon_,
		  nMon[matchChain]);
		AddToArrayN(&seNum_, seNum[matchChain], nMon_,
		  nMon[matchChain]);
		AddToArrayN(&xyz_, comXYZ[matchChain], NXYZ_*3,
		  comNXYZ[matchChain]*3);
		AddToArrayN(&tmf_, comTmf[matchChain], NXYZ_,
		  comNXYZ[matchChain]);
		NXYZ_+=comNXYZ[matchChain];
		nMon_+=nMon[matchChain]; j+=nMon[matchChain];
	      }
	  for(j=0; j<nChain; j++)
	    if(chainStatus[j]==1)
	      {
		delete [] mon[j];
		delete [] seNum[j];
		delete [] comXYZ[j];
		delete [] comTmf[j];
		chainStatus[j]=2;
	      }
	  chainStatus[firstChain]=3;
	  mon[firstChain]=mon_;
	  seNum[firstChain]=seNum_;
	  comXYZ[firstChain]=xyz_;
	  comTmf[firstChain]=tmf_;
	  nMon[firstChain]=nMon_;
	  comNXYZ[firstChain]=NXYZ_;
	}
      delete [] seqMatch;
    }

  nChain_=0;
  for(j=0; j<nChain; j++)
    if(chainStatus[j]==0 || chainStatus[j]==3) nChain_++;
  _mon_=new int2 * [nChain_];
  _seNum_=new char ** [nChain_];
  _xyz_=new flt4 * [nChain_];
  _tmf_=new flt4 * [nChain_];
  _nMon_=new int2 [nChain_];
  _NXYZ_=new int2 [nChain_];
  _chainCode_=new char [nChain_];

  nChain_=0;
  for(j=0; j<nChain; j++)
    {
      if(chainStatus[j]==0 || chainStatus[j]==3)
	{
	  _mon_[nChain_]=mon[j];
	  _seNum_[nChain_]=seNum[j];
	  _xyz_[nChain_]=comXYZ[j];
	  _tmf_[nChain_]=comTmf[j];
	  _nMon_[nChain_]=nMon[j];
	  _NXYZ_[nChain_]=comNXYZ[j];
	  _chainCode_[nChain_]=chainCode[j];
	  nChain_++;
	}
    }
  if(nChain)
    {
      delete [] mon;
      delete [] seNum;
      delete [] comXYZ;
      delete [] comTmf;
      delete [] comNXYZ;
      delete [] nMon;
      delete [] chainCode;
      mon=_mon_;
      seNum=_seNum_;
      comXYZ=_xyz_;
      comTmf=_tmf_;
      nMon=_nMon_;
      comNXYZ=_NXYZ_;
      chainCode=_chainCode_;
      nChain=nChain_;
    }

  if(nChain) delete [] chainStatus;
  delete M;
}
*/
//////////////////////////////////////////////////////////////////
void Molecule::matchAtomSeqres()
{
  int i, j, k1, k2, *seqMatch, matchIndex, matchStatus, matchChain,
    *chainStatus, nChain_, firstChain; int2 index; 
    int nTermHet=0, cTermHet=0, nMismatch;

  int2 *mon_, nMon_; int NXYZ_; flt4 *xyz_, *tmf_; char **seNum_;
  int2 **_mon_, *_nMon_; flt4 **_tmf_, **_xyz_;
  int4 *_NXYZ_;
  char *_chainCode_, ***_seNum_;  char *_seNum;
  Monomer *M=new Monomer;

  if(nChain) chainStatus=new int [nChain];
  for(j=0; j<nChain; j++) chainStatus[j]=0;

  for(i=0; i<seqNchain; i++)
    {
      seqMatch=new int [seqL[i]+2]; // allow single terminal het-groups
      matchIndex=0;
      matchStatus=0;
      nTermHet=0;
      cTermHet=0;

      nMon_=0; NXYZ_=0;
      for(j=0; j<seqL[i]+2; j++) seqMatch[j]=-1;
      for(j=0; j<nChain; j++)
	{
	  if(chainStatus[j]) continue;
	  if(seqChainCode[i]!=chainCode[j]) continue;

	  // special case: handling single n-terminal het group which 
	  // doesn't exist in SEQRES
	  if(matchStatus==0 && nMon[j]>10
	     && strcmp(*(seq[i]), con->monCode3.item1(*(mon[j])))) 
	     {
	       for(k2=0; k2<nMon[j] && k2<10; k2++)
		 if(strcmp(*(seq[i]+k2), con->monCode3.item1(*(mon[j]+1+k2)))) 
		   goto break_n_term;
	       seqMatch[matchStatus]=j; matchStatus++; nTermHet=1;
	     }
	break_n_term:
	  //

	  for(k1=matchIndex; k1<seqL[i] && k1+nMon[j]<=seqL[i]+1; k1++)
            {
	      // special case: handling single c-terminal het group which 
	      // doesn't exist in SEQRES
	      cTermHet=0;
	      if(k1+nMon[j]-nTermHet==seqL[i]+1) 
		{
		  cTermHet=1;
		}
	      //
	      nMismatch=0;
	      for(k2=0; k2<nMon[j]-nTermHet-cTermHet; k2++)
		if(strcmp(*(seq[i]+k1+k2), 
			  con->monCode3.item1(*(mon[j]+k2+nTermHet))))
		  {
		    if(k2<10 || nMismatch>1) goto break_1;
		    nMismatch++;
		  }
	      matchStatus++;
	      for(k2=0; k2<nMon[j]-nTermHet-cTermHet; k2++) 
		seqMatch[k1+k2+nTermHet]=j;
	      matchIndex=k1+nMon[j]+nTermHet-cTermHet;
	      chainStatus[j]=1;
	      // special case (cont): handling single c-terminal het group
	      // which doesn't exist in SEQRES
	      if(cTermHet==1)
		{
		  seqMatch[matchIndex]=j; matchIndex++;
		}
	      //
	      break;
	    break_1: ;
	      cTermHet=0;
            }
        }

      if(matchStatus)
	{
	  firstChain=-1;
	  for(j=0; j<seqL[i]+nTermHet+cTermHet; )
	    if(seqMatch[j]==-1)
	      {
		strcpy(M->name, *(seq[i]+j-nTermHet));
		index=con->search(M);
		AddToArray(&mon_, index, nMon_);
		_seNum=new char [1]; _seNum[0]='\0';
		AddToArray(&seNum_, _seNum, nMon_);
		nMon_++; j++;
	      }
	    else
	      {
		matchChain=seqMatch[j];
		if(firstChain==-1) firstChain=matchChain;
		AddToArrayN(&mon_, mon[matchChain], nMon_,
		  nMon[matchChain]);
		AddToArrayN(&seNum_, seNum[matchChain], nMon_,
		  nMon[matchChain]);
		AddToArrayN(&xyz_, comXYZ[matchChain], NXYZ_*3,
		  comNXYZ[matchChain]*3);
		AddToArrayN(&tmf_, comTmf[matchChain], NXYZ_,
		  comNXYZ[matchChain]);
		NXYZ_+=comNXYZ[matchChain];
		nMon_+=nMon[matchChain]; j+=nMon[matchChain];
	      }
	  for(j=0; j<nChain; j++)
	    if(chainStatus[j]==1)
	      {
		delete [] mon[j];
		delete [] seNum[j];
		delete [] comXYZ[j];
		delete [] comTmf[j];
		chainStatus[j]=2;
	      }
	  chainStatus[firstChain]=3;
	  mon[firstChain]=mon_;
	  seNum[firstChain]=seNum_;
	  comXYZ[firstChain]=xyz_;
	  comTmf[firstChain]=tmf_;
	  nMon[firstChain]=nMon_;
	  comNXYZ[firstChain]=NXYZ_;
	}
      delete [] seqMatch;
    }

  nChain_=0;
  for(j=0; j<nChain; j++)
    if(chainStatus[j]==0 || chainStatus[j]==3) nChain_++;
  _mon_=new int2 * [nChain_];
  _seNum_=new char ** [nChain_];
  _xyz_=new flt4 * [nChain_];
  _tmf_=new flt4 * [nChain_];
  _nMon_=new int2 [nChain_];
  _NXYZ_=new int4 [nChain_];
  _chainCode_=new char [nChain_];

  nChain_=0;
  for(j=0; j<nChain; j++)
    {
      if(chainStatus[j]==0 || chainStatus[j]==3)
	{
	  _mon_[nChain_]=mon[j];
	  _seNum_[nChain_]=seNum[j];
	  _xyz_[nChain_]=comXYZ[j];
	  _tmf_[nChain_]=comTmf[j];
	  _nMon_[nChain_]=nMon[j];
	  _NXYZ_[nChain_]=comNXYZ[j];
	  _chainCode_[nChain_]=chainCode[j];
	  nChain_++;
	}
    }
  if(nChain)
    {
      delete [] mon;
      delete [] seNum;
      delete [] comXYZ;
      delete [] comTmf;
      delete [] comNXYZ;
      delete [] nMon;
      delete [] chainCode;
      mon=_mon_;
      seNum=_seNum_;
      comXYZ=_xyz_;
      comTmf=_tmf_;
      nMon=_nMon_;
      comNXYZ=_NXYZ_;
      chainCode=_chainCode_;
      nChain=nChain_;
    }

  if(nChain) delete [] chainStatus;
  delete M;
}
//////////////////////////////////////////////////////////////////
void Molecule::sortChains()
{
  if(!nChain) return;

  int j, nChain_;
  int2 **_mon_, *_nMon_;
  int4 *_NXYZ_;
  char *_chainCode_, ***_seNum_;  
  flt4 **_tmf_, **_xyz_;

  _mon_=new int2 * [nChain];
  _seNum_=new char ** [nChain];
  _xyz_=new flt4 * [nChain];
  _tmf_=new flt4 * [nChain];
  _nMon_=new int2 [nChain];
  _NXYZ_=new int4 [nChain];
  _chainCode_=new char [nChain];

  nChain_=0;

  int *useChain = new int [nChain];
  for(j=0; j<nChain; j++) useChain[j] = 0;
  
  for(j=0; j<nChain; j++)
    {
      int nCanon = 0;
      if(nMon[j] > 1) 
	for(int is = 0; is < nMon[j]; is++) 
	  if(*con->monType.item1(*(mon[j]+is)) != 0) {
	    nCanon++;
	    if(nCanon > 1) break; 
	  }
      
      if(nCanon > 1)
	{
	  _mon_[nChain_]=mon[j];
	  _seNum_[nChain_]=seNum[j];
	  _xyz_[nChain_]=comXYZ[j];
	  _tmf_[nChain_]=comTmf[j];
	  _nMon_[nChain_]=nMon[j];
	  _NXYZ_[nChain_]=comNXYZ[j];
	  _chainCode_[nChain_]=chainCode[j];
	  nChain_++;
	  useChain[j] = 1;
        }
    }

  nPolyChain = nChain_;

  for(j=0; j<nChain; j++)
    {
      if(!useChain[j])
	{
	  _mon_[nChain_]=mon[j];
	  _seNum_[nChain_]=seNum[j];
	  _xyz_[nChain_]=comXYZ[j];
	  _tmf_[nChain_]=comTmf[j];
	  _nMon_[nChain_]=nMon[j];
	  _NXYZ_[nChain_]=comNXYZ[j];
	  _chainCode_[nChain_]=chainCode[j];
	  nChain_++;
        }
     }
  delete [] useChain;

  if(nChain)
    {
      delete [] mon;
      delete [] seNum;
      delete [] comXYZ;
      delete [] comTmf;
      delete [] comNXYZ;
      delete [] nMon;
      delete [] chainCode;
      mon=_mon_;
      seNum=_seNum_;
      comXYZ=_xyz_;
//      comOcc=_occ_;
      comTmf=_tmf_;
      nMon=_nMon_;
      comNXYZ=_NXYZ_;
      chainCode=_chainCode_;
    }

// Small piece for assigning numbers if there is no pdb number
 char buf[10], *_seNum;
 int ins, i;
 for(j=0; j<nChain; j++)
    {
      ins=1;
      for(i=0; i<nMon[j]; i++)
	{
	  if(strlen(*(seNum[j]+i))==0)
	    {
	      sprintf(buf, "#%d", ins); ins++;
	      _seNum=new char [strlen(buf)+1];
	      strcpy(_seNum, buf);
	      delete [] *(seNum[j]+i);
	      *(seNum[j]+i)=_seNum;
	    }
	      
	}
    }
}
//////////////////////////////////////////////////////////////////
/*
int Molecule::assignKS()
{
  int i, j, iC;
  int NSE, k, strandStatus, breakStatus;
  char *seq, *turn3, *turn4, *turn5, *pbridge, *abridge, *bend,
  *helix3, *helix4, *helix5, *summary;
  int monCA1, monCA2, monCA3, CA1, CA2, CA3, aCA1, aCA2, aCA3;
  XYZ ca1, ca2, ca3;
//  FILE * fks=fopen("ks","wa");

  for(iC=0; iC<nChain; iC++)
    {
      NSE=nMon[iC];
      if(NSE<2) continue;
      seq=new char [NSE];
      turn3=new char [NSE];
      turn4=new char [NSE];
      turn5=new char [NSE];
      pbridge=new char [NSE];
      abridge=new char [NSE];
      bend=new char [NSE];
      helix3=new char [NSE];
      helix4=new char [NSE];
      helix5=new char [NSE];
      summary=new char [NSE];
      if(!seq || !turn3 || !turn4 || !turn5 || !pbridge || !abridge ||
	!bend || !helix3 || !helix4 ||
	!helix5 || !summary) return(1);
      for(NSE=0; NSE<nMon[iC]; NSE++)
	{
	  seq[NSE]=*con->monCode1.item1(*(mon[iC]+NSE));
	  turn3[NSE]=' ';
	  turn4[NSE]=' ';
	  turn5[NSE]=' ';
	  pbridge[NSE]=' ';
	  abridge[NSE]=' ';
	  bend[NSE]=' ';
	  helix3[NSE]=' ';
	  helix4[NSE]=' ';
	  helix5[NSE]=' ';
	  summary[NSE]=' ';
	}

      NSE=nMon[iC];
      for(i=0; i<NSE; i++)
	{
	  if(i<NSE-3)
	    {
	      if(HBondingEnergy(iC, i, "O", i+3, "N"))
		for(k=0; k<4; k++) turn3[i+k]='3';
	    }
	  if(i<NSE-4)
	    {
	      if(HBondingEnergy(iC, i, "O", i+4, "N"))
		for(k=0; k<5; k++) turn4[i+k]='4';
	    }
	  if(i<NSE-5)
	    {
	      if(HBondingEnergy(iC, i, "O", i+5, "N"))
		for(k=0; k<6; k++) turn5[i+k]='5';
	    }

	  if(i<NSE-2)
	    {
	      for(j=0; j<i-2; j++)
		{
		  if((HBondingEnergy(iC, i, "O", j+1, "N") &&
		      HBondingEnergy(iC, i+2, "N", j+1, "O"))
		     ||
		     (HBondingEnergy(iC, i+1, "N", j, "O") &&
		      HBondingEnergy(iC, i+1, "O", j+2, "N")))
		    {
		      pbridge[i+1]='x';
		      pbridge[j+1]='x';
		    }
		  
		  if((HBondingEnergy(iC, i, "O", j+2, "N") &&
		      HBondingEnergy(iC, i+2, "N", j, "O"))
		     ||
		     (HBondingEnergy(iC, i+1, "N", j+1, "O") &&
		      HBondingEnergy(iC, i+1, "O", j+1, "N")))
		    {
		      abridge[i+1]='X';
		      abridge[j+1]='X';
		    }
		  
		}
	    }

	  if(i<NSE-4)
	    {
	      monCA1=*(mon[iC]+i);
	      monCA2=*(mon[iC]+i+2);
	      monCA3=*(mon[iC]+i+4);
	      CA1=findAtom(monCA1, " CA ");
	      CA2=findAtom(monCA2, " CA ");
	      CA3=findAtom(monCA3, " CA ");
	      if(CA1>=0 && CA2>=0 && CA3>=0)
		{
	          aCA1=*(monXYZ[iC]+i);
	          aCA2=*(monXYZ[iC]+i+2);
		  aCA3=*(monXYZ[iC]+i+4);
		  CA1+=aCA1; CA2+=aCA2; CA3+=aCA3;
		  ca1.X=*(comXYZ[iC]+CA1*3);
		  ca1.Y=*(comXYZ[iC]+CA1*3+1);
		  ca1.Z=*(comXYZ[iC]+CA1*3+2);
		  ca2.X=*(comXYZ[iC]+CA2*3);
		  ca2.Y=*(comXYZ[iC]+CA2*3+1);
		  ca2.Z=*(comXYZ[iC]+CA2*3+2);
		  ca3.X=*(comXYZ[iC]+CA3*3);
		  ca3.Y=*(comXYZ[iC]+CA3*3+1);
		  ca3.Z=*(comXYZ[iC]+CA3*3+2);
		  if(ca1.dist(ca2)+ca2.dist(ca3)!=0.0)
		    if(ca1.dist(ca3)/(ca1.dist(ca2)+ca2.dist(ca3))<0.82)
		      bend[i+2]='S';
		}
	    }
	}

      for(i=0; i<NSE; i++)
	{
	  if(i<NSE-4)
	    {
	      breakStatus=0;
	      for(j=0; j<5; j++)
		if(turn3[i+j]==' ') breakStatus=1;
	      if(!breakStatus)
		for(j=0; j<3; j++) helix3[i+1+j]='G';
	    }

	  if(i<NSE-5)
	    {
	      breakStatus=0;
	      for(j=0; j<6; j++)
		if(turn4[i+j]==' ') breakStatus=1;
	      if(!breakStatus)
	        for(j=0; j<4; j++) helix4[i+1+j]='H';
	    }

	  if(i<NSE-6)
	    {
	      breakStatus=0;
	      for(j=0; j<7; j++)
		if(turn5[i+j]==' ') breakStatus=1;
	      if(!breakStatus)
		for(j=0; j<5; j++) helix5[i+1+j]='I';
	    }
	}

      for(i=0; i<NSE; i++)
	{
	  summary[i]=bend[i];
	  if(turn3[i]!=' ' || turn4[i]!=' ' || turn5[i]!=' ')
	    summary[i]='T';
	  if(helix5[i]!=' ') summary[i]=helix5[i];
	  if(helix3[i]!=' ') summary[i]=helix3[i];
	  if(pbridge[i]!=' ') summary[i]=pbridge[i];
	  if(abridge[i]!=' ') summary[i]=abridge[i];
	  if(helix4[i]!=' ') summary[i]=helix4[i];
	}

      for(i=0; i<NSE; i++)
	{
	  strandStatus=0;
	  if(summary[i]=='x' || summary[i]=='X')
	    {
	      if(i>0)
		if(summary[i-1]=='E') strandStatus=1;
	      if(i<NSE-1)
		if(summary[i+1]=='x' || summary[i+1]=='X') strandStatus=1;
	      summary[i]=strandStatus?'E':'B';
	    }
	  *(*(ks+iC)+i)=summary[i];
	}
*/
/*
      fprintf(fks, "\n");
      for(i=0; i<NSE; i+=50)
	{
	  fprintf(fks, "\n");

	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", bend[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", pbridge[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", abridge[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn5[j]);
	    }

	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn4[j]);
	    }

	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn3[j]);
	    }

	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", summary[j]);
	    }

	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", seq[j]);
	    }

	}
*/
      /*
      delete [] seq;
      delete [] turn3;
      delete [] turn4;
      delete [] turn5;
      delete [] pbridge;
      delete [] abridge;
      delete [] bend;
      delete [] helix3;
      delete [] helix4;
      delete [] helix5;
      delete [] summary;
    }
//   fprintf(fks, "\n");
//   fclose(fks);
  return(0);
}
*/
///////////////////////////////////////////////////////////////////////////
int Molecule::assignKS()
{
  int iC, i, j;
  int NSE, k, strandStatus, breakStatus;
  char *seq, *turn3, *turn4, *turn5, *pbridge, *abridge, *bend,
  *helix3, *helix4, *helix5, *summary;
  int monCA1, monCA2, monCA3, CA1, CA2, CA3, aCA1, aCA2, aCA3;
  XYZ ca1, ca2, ca3;
  //FILE * fks=fopen("ks","a");
  
  for(iC=0; iC<nChain; iC++)
    {
      NSE=nMon[iC];
      if(NSE<2) continue;
      seq=new char [NSE];
      turn3=new char [NSE];
      turn4=new char [NSE];
      turn5=new char [NSE];
      pbridge=new char [NSE];
      abridge=new char [NSE];
      bend=new char [NSE];
      helix3=new char [NSE];
      helix4=new char [NSE];
      helix5=new char [NSE];
      summary=new char [NSE];
      if(!seq || !turn3 || !turn4 || !turn5 || !pbridge || !abridge ||
	 !bend || !helix3 || !helix4 ||
	 !helix5 || !summary) return(1);
      for(NSE=0; NSE<nMon[iC]; NSE++)
	{
	  seq[NSE]=*con->monCode1.item1(*(mon[iC]+NSE));
	  turn3[NSE]=' ';
	  turn4[NSE]=' ';
	  turn5[NSE]=' ';
	  pbridge[NSE]=' ';
	  abridge[NSE]=' ';
	  bend[NSE]=' ';
	  helix3[NSE]=' ';
	  helix4[NSE]=' ';
	  helix5[NSE]=' ';
	  summary[NSE]=' ';
	}
      
      NSE=nMon[iC];
      for(i=1; i<NSE; i++)
	{
	  
	  if(i<NSE-3)
	    {
	      if(HBondingEnergy(iC, i, "O", i+3, "N")) {
		if(turn3[i] == ' ' || turn3[i] == '3') turn3[i] =  '>';
		else if(turn3[i] == '<') turn3[i] = 'X';
		if(turn3[i+1] == ' ') turn3[i+1] = '3';
		if(turn3[i+2] == ' ') turn3[i+2] = '3';
		turn3[i+3] = '<';
	      }
	    }
	  if(i<NSE-4)
	    {
	      if(HBondingEnergy(iC, i, "O", i+4, "N")) {
		if(turn4[i] == ' ' || turn4[i] == '4') turn4[i]  ='>';
		else if(turn4[i] == '<') turn4[i] = 'X';
		if(turn4[i+1] == ' ') turn4[i+1]='4';
		if(turn4[i+2] == ' ') turn4[i+2]='4';
		if(turn4[i+3] == ' ') turn4[i+3]='4';
		turn4[i+4]='<';
	      }
	    }
	  if(i<NSE-5)
	    {
	      if(HBondingEnergy(iC, i, "O", i+5, "N")) {
		if(turn5[i] == ' ' || turn5[i] == '5') turn5[i]  ='>';
		else if(turn5[i] == '<') turn5[i] = 'X';
		if(turn5[i+1] == ' ') turn5[i+1]='5';
		if(turn5[i+2] == ' ') turn5[i+2]='5';
		if(turn5[i+3] == ' ') turn5[i+3]='5';
		if(turn5[i+4] == ' ') turn5[i+4]='5';
		turn5[i+5]='<';
	      }
	    }
	  
	  if(i<NSE-2)
	    {
	      for(j=0; j<i-2; j++)
		{
		  if((HBondingEnergy(iC, i, "O", j+1, "N") &&
		      HBondingEnergy(iC, i+2, "N", j+1, "O"))
		     ||
		     (HBondingEnergy(iC, i+1, "N", j, "O") &&
		      HBondingEnergy(iC, i+1, "O", j+2, "N")))
		    {
		      pbridge[i+1]='x';
		      pbridge[j+1]='x';
		    }
		  
		  if((HBondingEnergy(iC, i, "O", j+2, "N") &&
		      HBondingEnergy(iC, i+2, "N", j, "O"))
		     ||
		     (HBondingEnergy(iC, i+1, "N", j+1, "O") &&
		      HBondingEnergy(iC, i+1, "O", j+1, "N")))
		    {
		      abridge[i+1]='X';
		      abridge[j+1]='X';
		    }
		  
		}
	    }
	  
	  if(i<NSE-4)
	    {
	      monCA1=*(mon[iC]+i);
	      monCA2=*(mon[iC]+i+2);
	      monCA3=*(mon[iC]+i+4);
	      CA1=findAtom(monCA1, " CA ");
	      CA2=findAtom(monCA2, " CA ");
	      CA3=findAtom(monCA3, " CA ");
	      if(CA1>=0 && CA2>=0 && CA3>=0)
		{
	          aCA1=*(monXYZ[iC]+i);
	          aCA2=*(monXYZ[iC]+i+2);
		  aCA3=*(monXYZ[iC]+i+4);
		  CA1+=aCA1; CA2+=aCA2; CA3+=aCA3;
		  ca1.X=*(comXYZ[iC]+CA1*3);
		  ca1.Y=*(comXYZ[iC]+CA1*3+1);
		  ca1.Z=*(comXYZ[iC]+CA1*3+2);
		  ca2.X=*(comXYZ[iC]+CA2*3);
		  ca2.Y=*(comXYZ[iC]+CA2*3+1);
		  ca2.Z=*(comXYZ[iC]+CA2*3+2);
		  ca3.X=*(comXYZ[iC]+CA3*3);
		  ca3.Y=*(comXYZ[iC]+CA3*3+1);
		  ca3.Z=*(comXYZ[iC]+CA3*3+2);
		  if(ca1.dist(ca2)+ca2.dist(ca3)!=0.0)
		    if(ca1.dist(ca3)/(ca1.dist(ca2)+ca2.dist(ca3))<0.82)
		      bend[i+2]='S';
		}
	    }
	}
      
      int hBegin, hEnd = -1;
      for(i=0; i<NSE; i++)
	if(i<NSE-4 && i > hEnd)
	  {
	    hBegin = -1; hEnd = -1;
	    if((turn3[i] == '>' || turn3[i] == 'X') && 
	       (turn3[i+1] == '>' || turn3[i+1] == 'X')) hBegin = i+1;
	    if(hBegin != -1)
	      for(j = i+3; j < NSE-2; j++) {
		if(turn3[j+1] == ' ') break;
		if(hEnd != -1 && turn3[j+1] == '>') break;
		if((turn3[j] == '<' || turn3[j] == 'X') 
		   && (turn3[j+1] == '<' || turn3[j+1] == 'X')) hEnd = j;
	      }
	    if(hBegin != -1 && hEnd != -1) {
	      for(j = hBegin; j <= hEnd; j++) helix3[j]='G';
	    }
	  }

      hEnd = -1;
      for(i=0; i<NSE; i++)
	if(i<NSE-5 && i > hEnd)
	  {
	    hBegin = -1; hEnd = -1;
	    if((turn4[i] == '>' || turn4[i] == 'X') && 
	       (turn4[i+1] == '>' || turn4[i+1] == 'X')) hBegin = i+1;
	    if(hBegin != -1)
	      for(j = i+4; j < NSE-2; j++) {
		if(turn4[j+1] == ' ') break;
		if(hEnd != -1 && turn4[j+1] == '>') break;
		if((turn4[j] == '<' || turn4[j] == 'X') 
		   && (turn4[j+1] == '<' || turn4[j+1] == 'X')) hEnd = j;
	      }
	    if(hBegin != -1 && hEnd != -1) {
	      for(j = hBegin; j <= hEnd; j++) helix4[j]='H';
	    }
	  }
      
      hEnd = -1;
      for(i=0; i<NSE; i++)
	if(i<NSE-6 && i > hEnd)
	  {
	    hBegin = -1; hEnd = -1;
	    if((turn5[i] == '>' || turn5[i] == 'X') && 
	       (turn5[i+1] == '>' || turn5[i+1] == 'X')) hBegin = i+1;
	    if(hBegin != -1)
	      for(j = i+5; j < NSE-2; j++) {
		if(turn5[j+1] == ' ') break;
		if(hEnd != -1 && turn5[j+1] == '>') break;
		if((turn5[j] == '<' || turn5[j] == 'X') 
		   && (turn5[j+1] == '<' || turn5[j+1] == 'X')) hEnd = j;
	      }
	    if(hBegin != -1 && hEnd != -1) {
	      for(j = hBegin; j <= hEnd; j++) helix5[j]='I';
	    }
	  }
      
      
      for(i=0; i<NSE; i++)
	{
	  summary[i]=bend[i];

	 
	  if(turn3[i]!=' ') {
	    int nb_flag = 1;
	    if(i != 0) 
	      if(turn3[i-1] == ' ') nb_flag = 0;
	    if(i != NSE-1)
	      if(turn3[i+1] == ' ') nb_flag = 0;
	    if(nb_flag) summary[i]='T';
	  }
	  if(turn4[i]!=' ') {
	    int nb_flag = 1;
	    if(i != 0) 
	      if(turn4[i-1] == ' ') nb_flag = 0;
	    if(i != NSE-1)
	      if(turn4[i+1] == ' ') nb_flag = 0;
	    if(nb_flag) summary[i]='T';
	  }
	  if(turn5[i]!=' ') {
	    int nb_flag = 1;
	    if(i != 0) 
	      if(turn5[i-1] == ' ') nb_flag = 0;
	    if(i != NSE-1)
	      if(turn5[i+1] == ' ') nb_flag = 0;
	    if(nb_flag) summary[i]='T';
	  }
	  if(helix5[i]!=' ') summary[i]=helix5[i];
	  if(helix3[i]!=' ') summary[i]=helix3[i];
	  if(pbridge[i]!=' ') summary[i]=pbridge[i];
	  if(abridge[i]!=' ') summary[i]=abridge[i];
	  if(helix4[i]!=' ') summary[i]=helix4[i];
	}
      
      for(i=0; i<NSE; i++)
	{
	  strandStatus=0;
	  if(summary[i]=='x' || summary[i]=='X')
	    {
	      if(i>0)
		if(summary[i-1]=='E') strandStatus=1;
	      if(i<NSE-1)
		if(summary[i+1]=='x' || summary[i+1]=='X') strandStatus=1;
	      summary[i]=strandStatus?'E':'B';
	    }
	}

      hBegin = -1; hEnd = -1;
      for(i=0; i<NSE; i++) {
	if(summary[i] == 'H') {
	  if(hBegin == -1) hBegin = i;
	  hEnd = i;
	}
	else {
	  if(hBegin != -1 && hEnd != -1) {
	    if(hEnd - hBegin < 3)
	      for(j = hBegin; j <= hEnd; j++)
		summary[j] = 'T';
	    hBegin = -1; hEnd = -1;
	  }
	}
      }

      hBegin = -1; hEnd = -1;
      for(i=0; i<NSE; i++) {
	if(summary[i] == 'G') {
	  if(hBegin == -1) hBegin = i;
	  hEnd = i;
	}
	else {
	  if(hBegin != -1 && hEnd != -1) {
	    if(hEnd - hBegin < 2)
	      for(j = hBegin; j <= hEnd; j++)
		summary[j] = 'T';
	    hBegin = -1; hEnd = -1;
	  }
	}
      }

      hBegin = -1; hEnd = -1;
      for(i=0; i<NSE; i++) {
	if(summary[i] == 'I') {
	  if(hBegin == -1) hBegin = i;
	  hEnd = i;
	}
	else {
	  if(hBegin != -1 && hEnd != -1) {
	    if(hEnd - hBegin < 4)
	      for(j = hBegin; j <= hEnd; j++)
		summary[j] = 'T';
	    hBegin = -1; hEnd = -1;
	  }
	}
      }

      for(i=0; i<NSE; i++)
	*(*(ks+iC)+i)=summary[i];


	/*
      fprintf(fks, "\n\n%4.4s:%c\n", id, chainCode[iC]);
      fprintf(fks, "\n");
      for(i=0; i<NSE; i+=50)
	{
	  fprintf(fks, "\n");
	  
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", bend[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", pbridge[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", abridge[j]);
	    }
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn5[j]);
	    }
	  
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn4[j]);
	    }
	  
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", turn3[j]);
	    }
	  
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", summary[j]);
	    }
	  
	  fprintf(fks, "\n");
	  for(j=i; j<((NSE<i+50)?NSE:i+50); j++)
	    {
	      if(j%10==0) fprintf(fks, " ");
	      fprintf(fks, "%c", seq[j]);
	    }
	  
	}
	*/
      delete [] seq;
      delete [] turn3;
      delete [] turn4;
      delete [] turn5;
      delete [] pbridge;
      delete [] abridge;
      delete [] bend;
      delete [] helix3;
      delete [] helix4;
      delete [] helix5;
      delete [] summary;
    }
  //fprintf(fks, "\n");
  //fclose(fks);
  return(0);
}
///////////////////////////////////////////////////////////////////////////
#ifndef MINIMAL
int Molecule::assignEnv()
{
#define NN_ARRAY_SIZE 1000
  int i, i1, i2, j, j1, j2, j_, iE, iSE, jSE, jAtom, nAtom;
  int NSE, k, nn, iW, nnArray[NN_ARRAY_SIZE];
  char nnTypes[NN_ARRAY_SIZE], nnRange[NN_ARRAY_SIZE];
  char atomType;
  XYZ a1, a2;
  flt4 aX, aY, aZ, bX, bY, bZ, wX, wY, wZ,
    *waterX=new flt4[nWaterData],
    *waterY=new flt4[nWaterData],
    *waterZ=new flt4[nWaterData];
  double dX, dY, dZ;
  flt4 rAtom=1.8, rWater=1.4;
  flt4 rCont=(rAtom+rWater)*2.0, rWaterCont=rAtom+rWater;
  double dCont2=((double)rCont)*((double)rCont),
    dWaterCont2=((double)rWaterCont)*((double)rWaterCont);
  long nW, nP;

  for(iW=0; iW<nWaterData; iW++)
    {
      waterX[iW]=(waterData[iW*3]*rWaterCont);
      waterY[iW]=(waterData[iW*3+1]*rWaterCont);
      waterZ[iW]=(waterData[iW*3+2]*rWaterCont);
    }

  for(iE=0; iE<nChain; iE++)
    {
      NSE=nMon[iE];

      for(iSE=0; iSE<NSE; iSE++)
	{
	  *(*(exp+iE)+iSE)=-1.0;
	  *(*(pol+iE)+iSE)=-1.0;
	  nW=0;
	  nP=0;
	  i1=*(monXYZ[iE]+iSE);
	  i2=*(monXYZ[iE]+iSE+1);
	  for(i=i1; i<i2; i++)
	    {
	      aX=*(comXYZ[iE]+i*3);
	      aY=*(comXYZ[iE]+i*3+1);
	      aZ=*(comXYZ[iE]+i*3+2);

              nn=0;

	      for(jSE=0; jSE<NSE; jSE++)
		{
		  j1=*(monXYZ[iE]+jSE);
		  j2=*(monXYZ[iE]+jSE+1);
		  for(j=j1; j<j2; j++)
		    {
		      if(i==j) continue;
		      bX=*(comXYZ[iE]+j*3);
		      if(fabs(aX-bX)>rCont) continue;
		      bY=*(comXYZ[iE]+j*3+1);
		      if(fabs(aY-bY)>rCont) continue;
		      bZ=*(comXYZ[iE]+j*3+2);
		      if(fabs(aZ-bZ)>rCont) continue;
		      dX=aX-bX; dY=aY-bY; dZ=aZ-bZ;
		      if(dX*dX+dY*dY+dZ*dZ>dCont2) continue;
		      if(nn>=NN_ARRAY_SIZE) continue;
		      atomType=*(con->monAtom.item1(*(mon[iE]+jSE))+
				 (j-j1)*4+1);
		      nnArray[nn]=j; 
		      nnTypes[nn]=(atomType=='N' || atomType=='O')?1:0; 
		      nnRange[nn]= (jSE==iSE || (abs(jSE-iSE)==1 && j-j1<4))?
			0:1;
		      nn++;
		    }
		}

	      for(iW=0; iW<nWaterData; iW++)
		{
		  wX=aX+waterX[iW];
		  wY=aY+waterY[iW];
		  wZ=aZ+waterZ[iW];
		  for(j_=0; j_<nn; j_++)
		    {
		      j=nnArray[j_];
		      bX=*(comXYZ[iE]+j*3);
		      if(fabs(wX-bX)>rWaterCont) continue;
		      bY=*(comXYZ[iE]+j*3+1);
		      if(fabs(wY-bY)>rWaterCont) continue;
		      bZ=*(comXYZ[iE]+j*3+2);
		      if(fabs(wZ-bZ)>rWaterCont) continue;
		      dX=wX-bX; dY=wY-bY; dZ=wZ-bZ;
		      if(dX*dX+dY*dY+dZ*dZ>dWaterCont2) continue;
                      goto cont1;
		    }
		  nW++; continue;
                 cont1: if(nnTypes[j_] && nnRange[j_]) nP++;
		}
	    }
	  *(*(exp+iE)+iSE)=(i2>i1?((double)nW)/
			    ((double)((i2-i1)*nWaterData))*100.0:-1);
	  *(*(pol+iE)+iSE)=(i2>i1?((double)nP)/
			    ((double)((i2-i1)*nWaterData))*100.0:-1);
	}
    }
  delete [] waterX;
  delete [] waterY;
  delete [] waterZ;
  return(0);
}
#endif
///////////////////////////////////////////////////////////////////////////
extern int aa_decode[25];
extern char aa_code[20];
///////////////////////////////////////////////////////////////////////////
int Molecule::assignFiltProp()
{
  int nEnt = nPolyChain, iE1, iS1, monomer;

//  Property exp_flt_enp("exp_flt.enp", 0, db);
//  Property pol_flt_enp("pol_flt.enp", 0, db);
//  Property seq_flt_enp("seq_flt.enp", 0, db);
//  Property iseq_flt_enp("iseq_flt.enp", 0, db);
//  Property iseq_enp("iseq.enp", 0, db);
//  Property se_type_enp("se_type.enp", 0, db);


#ifndef MINIMAL

  char **prop0=new char* [nEnt];
  char **prop3=new char* [nEnt];
  flt4 **prop1=new flt4* [nEnt];
  flt4 **prop2=new flt4* [nEnt];
  int1 **prop_is=new int1* [nEnt];
  int1 **prop_isf=new int1* [nEnt];

#endif

  int1 **se_type=new int1* [nEnt];
  for(iE1=0; iE1<nEnt; ++iE1)
    {

#ifndef MINIMAL

      prop0[iE1]=new char [nMon[iE1]+1]; *(prop0[iE1]+nMon[iE1])='\0';
      prop3[iE1]=new char [nMon[iE1]+1]; *(prop3[iE1]+nMon[iE1])='\0';
      prop_is[iE1]=new int1 [nMon[iE1]];
      prop_isf[iE1]=new int1 [nMon[iE1]]; 
      prop1[iE1]=new flt4 [nMon[iE1]];
      prop2[iE1]=new flt4 [nMon[iE1]];

#endif
      se_type[iE1]=new int1 [nMon[iE1]];
      char seq_, aa_code_;
      for(int i=0; i<nMon[iE1]; i++) 
	{
	  monomer=*(mon[iE1]+i);
	  *(se_type[iE1]+i)=*con->monType.item1(monomer);

#ifndef MINIMAL

	  seq_=*con->monCode1.item1(monomer);
	  *(prop0[iE1]+i)=*con->monNatom.item2(monomer)<4?'X':seq_;
	  *(prop3[iE1]+i)=seq_;

	  int monType_ = *con->monType.item1(monomer);

	  aa_code_=20;
	  if(seq_ != 'X' && monType_ == 1)
	    {
	      seq_-=65;
	      if(seq_<0 || seq_>24)
		{
		  printf("\nMolecule::assignFiltProp - seq error %d\n", seq_); 
		  exit(0);
		}
	      aa_code_=aa_decode[seq_]-1;
	      if(aa_code_<0 || aa_code_>19)
		{
		  printf("\nMolecule::assignFiltProp - seq error* %d (%c)\n", 
			 aa_code_, seq_); 
		  exit(0);
		}
	    }
	  *(prop_is[iE1]+i)=aa_code_;
	  *(prop_isf[iE1]+i)=*con->monNatom.item2(monomer)<4?20:aa_code_;
	  *(prop1[iE1]+i)=(*(prop0[iE1]+i)=='X')?-1.0:*(exp[iE1]+i);
	  *(prop2[iE1]+i)=(*(prop0[iE1]+i)=='X')?-1.0:*(pol[iE1]+i);

#endif

	}
	  
#ifndef MINIMAL

      if(P_SEQ_FLT_ENP) seq_flt_enp.addItem(prop0[iE1]);
      if(P_ISEQ_ENP) iseq_enp.addItem(prop_is[iE1], nMon[iE1]);
      if(P_ISEQ_FLT_ENP) iseq_flt_enp.addItem(prop_isf[iE1], nMon[iE1]);
      if(P_EXP_FLT_ENP) exp_flt_enp.addItem(prop1[iE1], nMon[iE1]);
      if(P_POL_FLT_ENP) pol_flt_enp.addItem(prop2[iE1], nMon[iE1]);

#endif

      if(P_SE_TYPE_ENP) se_type_enp.addItem(se_type[iE1], nMon[iE1]);
    }

  DeleteArray(&se_type, nEnt);

#ifndef MINIMAL

  DeleteArray(&prop0, nEnt);
  DeleteArray(&prop1, nEnt);
  DeleteArray(&prop2, nEnt);
  DeleteArray(&prop3, nEnt);
  DeleteArray(&prop_is, nEnt);
  DeleteArray(&prop_isf, nEnt);

#endif


  return(0);
}
///////////////////////////////////////////////////////////////////////////
int Molecule::HBondingEnergy(int iC, int iS1, char * atom1,
  int iS2, char * atom2)
{
  int iMon1=*(mon[iC]+iS1), iMon2=*(mon[iC]+iS2), iMonP,
    C, O, N, CA, CN, iA1=*(monXYZ[iC]+iS1), iA2=*(monXYZ[iC]+iS2), iAP;
  if(!strncasecmp(atom1,"O",1))
    {
      if((O=findAtom(iMon1, " O  "))<0) return(0);
      if((N=findAtom(iMon2, " N  "))<0) return(0);
      if(fabs(*(comXYZ[iC]+O*3)-*(comXYZ[iC]+N*3))>distanceCutoff) return(0);
      if(fabs(*(comXYZ[iC]+O*3+1)-*(comXYZ[iC]+N*3+1))>distanceCutoff) return(0);
      if(fabs(*(comXYZ[iC]+O*3+2)-*(comXYZ[iC]+N*3+2))>distanceCutoff) return(0);
      if((C=findAtom(iMon1, " C  "))<0) return(0);
      if((CA=findAtom(iMon2, " CA "))<0) return(0);
      if(iS2<1) return(0);
      iMonP=*(mon[iC]+iS2-1);
      iAP=*(monXYZ[iC]+iS2-1);
      if((CN=findAtom(iMonP, " C  "))<0) return(0);
      C+=iA1; O+=iA1; N+=iA2; CA+=iA2; CN+=iAP;
    }
  else
    {
      if((O=findAtom(iMon2, " O  "))<0) return(0);
      if((N=findAtom(iMon1, " N  "))<0) return(0);
      if(fabs(*(comXYZ[iC]+O*3)-*(comXYZ[iC]+N*3))>distanceCutoff) return(0);
      if(fabs(*(comXYZ[iC]+O*3+1)-*(comXYZ[iC]+N*3+1))>distanceCutoff) return(0);
      if(fabs(*(comXYZ[iC]+O*3+2)-*(comXYZ[iC]+N*3+2))>distanceCutoff) return(0);
      if((C=findAtom(iMon2, " C  "))<0) return(0);
      if((CA=findAtom(iMon1, " CA "))<0) return(0);
      if(iS1<1) return(0);
      iMonP=*(mon[iC]+iS1-1);
      iAP=*(monXYZ[iC]+iS1-1);
      if((CN=findAtom(iMonP, " C  "))<0) return(0);
      C+=iA2; O+=iA2; N+=iA1; CA+=iA1; CN+=iAP;
   }
  XYZ c, o, n, ca, h, cca, cn;
  c.X=*(comXYZ[iC]+C*3); c.Y=*(comXYZ[iC]+C*3+1);
  c.Z=*(comXYZ[iC]+C*3+2);
  o.X=*(comXYZ[iC]+O*3); o.Y=*(comXYZ[iC]+O*3+1);
  o.Z=*(comXYZ[iC]+O*3+2);
  n.X=*(comXYZ[iC]+N*3); n.Y=*(comXYZ[iC]+N*3+1);
  n.Z=*(comXYZ[iC]+N*3+2);
  ca.X=*(comXYZ[iC]+CA*3); ca.Y=*(comXYZ[iC]+CA*3+1);
  ca.Z=*(comXYZ[iC]+CA*3+2);
  cn.X=*(comXYZ[iC]+CN*3); cn.Y=*(comXYZ[iC]+CN*3+1);
  cn.Z=*(comXYZ[iC]+CN*3+2);

  cca=ca; cca+=cn; cca/=2.0; h=n; h-=cca; h/=0.69; h+=n;

  double e=0.42*0.20*332.0*(1.0/o.dist(n)+1.0/c.dist(h)-1.0/o.dist(h)
			   -1.0/c.dist(n));
//  printf("\n");
//  printf("C: %.2f, %.2f, %.2f\n", C->d_x(), C->d_y(), C->d_z());
//  printf("O: %.2f, %.2f, %.2f\n", O->d_x(), O->d_y(), O->d_z());
//  printf("N: %.2f, %.2f, %.2f\n", N->d_x(), N->d_y(), N->d_z());
//  printf("CA: %.2f, %.2f, %.2f\n", CA->d_x(), CA->d_y(), CA->d_z());
//  printf("CN: %.2f, %.2f, %.2f\n", CN->d_x(), CN->d_y(), CN->d_z());
//  printf("H: %.2f, %.2f, %.2f\n", h.X, h.Y, h.Z);

//  printf("HB: %d%s-%s %d%s-%s %.1f, %.1f, %.1f, %.1f %.3f\n",
//	SE1->_seqNum, SE1->_name, atom1, 
//	SE2->_seqNum, SE2->_name, atom2, 
//	o.dist(n),  c.dist(h), o.dist(h), c.dist(n), e);
  return(e<energyCutoff);
}
//////////////////////////////////////////////////////////////////
int Molecule::findAtom(int iMon, char * atom)
{
  char *atom_tmp=con->monAtom.item1(iMon);
  int nAtom_tmp=*con->monNatom.item2(iMon);
  for(int i=0; i<nAtom_tmp; i++)
    if(!strncmp(atom, atom_tmp+i*4, 4)) return(i);
  return(-1);
}
//////////////////////////////////////////////////////////////////
Connectivity::Connectivity()
{
  db=NULL;
}
//////////////////////////////////////////////////////////////////
int Connectivity::read()
{
  int i, err; int2 dummyI, monPrevNext[2];
  char * dummyC;

  if(P_CODE3_MON) monCode3.open("code3.mon");
  if(P_CODE1_MON) monCode1.open("code1.mon");
  if(P_N_BOND_MON) monNbond.open("n_bond.mon");
  if(P_BOND_MON) monBond.open("bond.mon");
  if(P_N_ATOM_MON) monNatom.open("n_atom.mon");
  if(P_ATOM_MON) monAtom.open("atom.mon");
  if(P_PREV_MON) monPrev.open("prev.mon");
  if(P_NEXT_MON) monNext.open("next.mon");
  if(P_TYPE_MON) monType.open("type.mon");
  
  return(0);
}
//////////////////////////////////////////////////////////////////
int Connectivity::search(Monomer *M)
{
  int i, j; int1 monType_; int4 pos;
  char *monCode3_, monCode1_, *monAtom_, *atomNames_; 
  int2 *monBond_;
  int nMon=monCode3.getObjectSize();
  int next_search=0;

  while((i=monCode3.find(M->name, 0, pos, next_search))!=-1)
    {
      next_search=1;
      if(M->nAtom!=*monNatom.item2(i)) continue;
      if(M->nBond!=*monNbond.item2(i)) continue;
      if(M->prevA-1!=*monPrev.item2(i) || M->nextA-1!=*monNext.item2(i)) 
	continue;
      atomNames_=monAtom.item1(i);
      if(M->nAtom) if(strcmp(M->atomNames, atomNames_)) goto break_1;
      monBond_=monBond.item2(i);
      for(j=0; j<*monNbond.item2(i)*2; j++)
	if(M->bond[j]!=monBond_[j]) goto break_1;
      return(i);
     break_1: ;
    }

  if(P_CODE3_MON) monCode3.addItem(M->name);
  if(P_N_ATOM_MON) monNatom.addItem(M->nAtom);
  if(P_N_BOND_MON) monNbond.addItem(M->nBond);
  if(P_ATOM_MON) monAtom.addItem(M->atomNames);
  if(P_BOND_MON) monBond.addItem(M->bond, (M->nBond*2));

  monCode3_=M->name;

  monCode1_ = 'X';
  monType_ = 1;
  if(!strcmp(monCode3_,"ALA")) monCode1_='A';
  if(!strcmp(monCode3_,"VAL")) monCode1_='V';
  if(!strcmp(monCode3_,"LEU")) monCode1_='L';
  if(!strcmp(monCode3_,"ILE")) monCode1_='I';
  if(!strcmp(monCode3_,"CYS")) monCode1_='C';
  if(!strcmp(monCode3_,"MET")) monCode1_='M';
  if(!strcmp(monCode3_,"PRO")) monCode1_='P';
  if(!strcmp(monCode3_,"PHE")) monCode1_='F';
  if(!strcmp(monCode3_,"TYR")) monCode1_='Y';
  if(!strcmp(monCode3_,"TRP")) monCode1_='W';
  if(!strcmp(monCode3_,"ASP")) monCode1_='D';
  if(!strcmp(monCode3_,"ASN")) monCode1_='N';
  if(!strcmp(monCode3_,"GLU")) monCode1_='E';
  if(!strcmp(monCode3_,"GLN")) monCode1_='Q';
  if(!strcmp(monCode3_,"HIS")) monCode1_='H';
  if(!strcmp(monCode3_,"SER")) monCode1_='S';
  if(!strcmp(monCode3_,"THR")) monCode1_='T';
  if(!strcmp(monCode3_,"ARG")) monCode1_='R';
  if(!strcmp(monCode3_,"LYS")) monCode1_='K';
  if(!strcmp(monCode3_,"GLY")) monCode1_='G';

  if(!strcmp(monCode3_,"  A") || !strcmp(monCode3_," +A")) {
    monCode1_='A'; monType_ = 2;
  }
  
  if(!strcmp(monCode3_,"  T") || !strcmp(monCode3_," +T")) {
    monCode1_='T'; monType_ = 2;
  }
  
  if(!strcmp(monCode3_,"  C") || !strcmp(monCode3_," +C")) {
    monCode1_='C'; monType_ = 2;
  }
  
  if(!strcmp(monCode3_,"  G") || !strcmp(monCode3_," +G")) {
    monCode1_='G'; monType_ = 2;
  }
  
  if(!strcmp(monCode3_,"  U") || !strcmp(monCode3_," +U")) {
    monCode1_='U'; monType_ = 2;
  }

  if(monCode1_ == 'X' && M->findAtom(" CA ") == -1) monType_ = 0;

  if(P_CODE1_MON) monCode1.addItem((int1)monCode1_);
  if(P_TYPE_MON) monType.addItem((int1)monType_);

  if(P_PREV_MON) monPrev.addItem((int2)(M->prevA-1));
  if(P_NEXT_MON) monNext.addItem((int2)(M->nextA-1));

  return(monCode3.getObjectSize()-1);
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
int main(int, char**, char **env=NULL);
main(int argc, char **argv, char **env)
{

  Tlpdb tlpdb;

  if(argc<2)
    {
      tlpdb.dialog();
    }
  else
    {
      if(!strcmp(argv[1], "update")) tlpdb.update(argc, argv);
      if(!strcmp(argv[1], "scratch")) tlpdb.scratch(argc, argv);
    }
}
//////////////////////////////////////////////////////////////////
/*
//            old update using log-file
//
void Tlpdb::update(int argc, char **argv) {
  if(argc<5) exit(0);
  strcpy(path, argv[4]);
  log_mode=1;

  char buffer[100], pdb_file[100];
  cboxes[1]=1;
  db->setPath(argv[2]);
  strcpy(db_path, db->path[0]);
  sprintf(log_file, "%sipdb_log", db_path);
  
  // printf("Log: %s\n", log_file);
  
  FILE *flog=fopen(log_file, "a");
  if(!flog) exit(0);
  fclose(flog);

  flog=fopen(argv[3], "r");
  if(!flog) exit(0);
  
  // printf("Opened %s\n", MIRRORLOG);
  
  char *time_buf, *time_buf_, *file_buf, *file_ext;
  if(argc<6)
    {
      time_buf=new char[9];
      time_t clock=time(NULL);
      tm *tm_=localtime(&clock);
      strftime(time_buf, 9, "%D", tm_);
    }
  else
    {
      time_buf=argv[5];
    }
  
  // printf("Doing update for %s\n", time_buf);
  // the following was modified to scan a special log file generated
  // by the mirror script
  
  while(fgets(buffer, 100, flog))
    {
      // if(strncmp("Added", buffer, 5)) continue;
      time_buf_=strrchr(buffer, ' ')+1;
      if(!time_buf_) continue;
      if(!strncmp(time_buf, time_buf_, 8))
	{
	  file_ext=time_buf_-4;
	  
	  if(!strncmp("ent", file_ext, 3)) {
	    sprintf(pdb_file, "%s/%0.14s", argv[4], buffer+4); 
	    AddBuffer(&loadList, pdb_file, nLoadList); nLoadList++;
	    continue;
	  }
	  
	  if(!strncmp("pdb", file_ext, 3)) {
	    sprintf(pdb_file, "%s/%0.8s", argv[4], buffer+4);
	    AddBuffer(&loadList, pdb_file, nLoadList); nLoadList++;
	  }
	  
	  // printf("%s\n", pdb_file);
	}
    }
  
  if(nLoadList) 
    {
      if(!con) con=new Connectivity;
      if(!mol) mol=new Molecule(con, this);
      db=mol->db;
      
      
      Property id_com("id.com");
      int nCom1, nCom2;
      nCom1=id_com.getObjectSize();
      id_com.close();
      load();
      id_com.open("id.com");
      nCom2=id_com.getObjectSize();
      DF::assignDF(db_path, nCom1, nCom2, log_mode);
      
      if(log_mode) {
	char buf[100];
	FILE *LOG;
	if((LOG=fopen(log_file,"a")) == NULL) {
	  printf("Error opening logfile!\n");
	  exit(0);
	}
	time_t clock=time(NULL);
	tm *tm_=localtime(&clock);
	strftime(buf, 100, "%D %T", tm_);
	fprintf(LOG, "%s finished loading and df for %d-%d\n",
		buf, nCom1+1, nCom2);
	fclose(LOG);
      }
      
    }
}
*/
//////////////////////////////////////////////////////////////////
void Tlpdb::update(int argc, char **argv) {
  if(argc<5) exit(0);

  db->setPath(argv[2]);
  Property lock_pid_misc("lock_pid.misc");
  int pid = *lock_pid_misc.item4(0);

  if(pid != -1) {
    //if(getpgid(pid) == -1) pid = -1;
  }
  if(pid == -1) {
    pid = getpid();
    lock_pid_misc.setItem(0, pid);
    lock_pid_misc.saveHeader();
  }
  else {
    printf("DB locked by %d\n", pid);
    exit(0);
  }


  strcpy(path, argv[3]);
  log_mode=1;

  char buffer[100], pdb_file[100];
  cboxes[1]=1;
  db->setPath(argv[2]);
  strcpy(db_path, db->path[0]);
  
  // check with db
  
  int searchResult=db->testFile("id.com");
  
  if(searchResult)
    {
      
      Property id_com("id.com");
      nDBList=id_com.getObjectSize();
      dbList=new char* [nDBList];
      for(int i=0; i<nDBList; i++) dbList[i]=id_com.item1(i, 1);
    }
  else
    {
	fprintf(stderr, "Could not check with db at %s - exiting!\n",
		db_path);
	exit(0);
    }
  

  cboxes[4] = 1; // check with db!
  int total = findFiles(path, "pdb*.ent", 0);

  loadList = new char* [total]; 
  nLoadList=0; total = 0;

  findFiles(path, files, 1);

  fprintf(stderr, "%d entries selected\n", nLoadList);

  
  if(nLoadList) 
    {
      
      char command[100];
      sprintf(command, "rm -rf %s", argv[4]); system(command);
      sprintf(command, "cp -pfr %s %s", argv[2], argv[4]); system(command);

      db->setPath(argv[4]);
      strcpy(db_path, db->path[0]);
      sprintf(log_file, "%sipdb_log", db_path);
  
      if(!con) con=new Connectivity;
      if(!mol) mol=new Molecule(con, this);
  
      
      Property id_com("id.com");
      int nCom1, nCom2;
      nCom1=id_com.getObjectSize();
      id_com.close();
      load();
      id_com.open("id.com");
      nCom2=id_com.getObjectSize();
      DF::assignDF(db_path, nCom1, nCom2, log_mode);
      
      if(log_mode) {
	char buf[100];
	FILE *LOG;
	if((LOG=fopen(log_file,"a")) == NULL) {
	  printf("Error opening logfile!\n");
	  exit(0);
	}
	time_t clock=time(NULL);
	tm *tm_=localtime(&clock);
	strftime(buf, 100, "%D %T", tm_);
	fprintf(LOG, "%s finished loading and df for %d-%d\n",
		buf, nCom1+1, nCom2);
	fclose(LOG);
      }

      // delete mol;
      // delete con;

      id_com.close();
      
      sprintf(command, "mv %s %s.bak", argv[2], argv[2]); system(command);
      sprintf(command, "ln -s %s %s", argv[4], argv[2]); system(command);
      sprintf(command, "cp -pf %s/* %s.bak", argv[4], argv[2]); system(command);
      sprintf(command, "rm -f %s", argv[2]); system(command);
      sprintf(command, "mv %s.bak %s", argv[2], argv[2]); system(command);

    }



}
//////////////////////////////////////////////////////////////////
void Tlpdb::scratch(int argc, char **argv)
{
  if(argc<4) exit(0);
  log_mode=1;
  scratch_mode=1;

  char buffer[100], pdb_file[100];
  cboxes[1]=0;
  if(argc>=6) cboxes[1]=1;
  
  db->setPath(argv[2]);
  strcpy(db_path, db->path[0]);
  sprintf(log_file, "%s/ipdb_log", db->path[0]);
  
  AddBuffer(&loadList, argv[3], nLoadList); nLoadList++;

  strcpy(scratch_id, "USR1");
  if(argc>=5) {
    strncpy(scratch_id, argv[4], 4);
    scratch_id[4]='\0';
  }

  if(nLoadList) 
    {
      load();
      Property id_com("id.com");
      int nCom2=id_com.getObjectSize();
      DF::assignDF(db_path, 0, nCom2, 1);
    }
}
//////////////////////////////////////////////////////////////////
void Tlpdb::dialog() {
  char mode, *buffer = new char[100], line[100];

  log_mode = 0;

  printf("Command (scan = s, load = l, calc df = d, quit = q)? [q] "); 
  scanf("%s", buffer); mode=buffer[0];

  if(mode=='s') {
      printf("Generate list into which file? [%s]\t", list1); 

      getInput(&list1);

      printf("Check with obsolete release (y/n)? [y]\t");
      isObs=(getYesNo('n') ? 0 : 1);

      if(isObs) {
	strcpy(path, OBSDIR);
	printf("obsolete directory? [%s]\t", path); getInput(&path);
	
	strcpy(files, "pdb*.ent");
	printf("obsolete frame? [%s]\t", files); getInput(&files);
      }

      printf("Check with current release (y/n)? [y]\t");

      isCur = (getYesNo('y') ? 1 : 0);

      if (isCur) {
	  strcpy(currentPath, CURRDIR);
	  printf("Path for current release? [%s]\t", currentPath);
	  getInput(&currentPath);
	  
	  strcpy(currFiles, "pdb*.ent");
	  printf("pdb file frame? [%s]\t\t", currFiles);
	  getInput(&currFiles);
      }
      
      printf("Include subdirectories (y/n)? [n]\t"); 

      cboxes[0] = (getYesNo('n') ? 0 : 1);

      printf("Check with second list (y/n)? [n]\t"); 

      cboxes[2] = (getYesNo('n') ? 0 : 1);

      if(cboxes[2]) {
	  printf("File with second list? [%s]\t", list2); 
	  getInput(&list2);
      }

      printf("Check with db (y/n)? [n]\t\t"); 

      cboxes[4] = (getYesNo('n') ? 0 : 1);

      if(cboxes[4]) {
	  strcpy(db_path, DBPATH);
	  printf("Path for db? [%s]\t", db_path); 
	  getInput(&db_path);
      }

      // printf("got %s", currentPath); exit(0);

      printf("Sort files (y/n)? [y]\t\t\t");

      sort_flag = (getYesNo('y') ? 1 : 0);

      scan();
  }

  if(mode=='l') {
      printf("Add to db (y/n)? [n]\t"); 

      cboxes[1] = (getYesNo('n') ? 0 : 1);

      printf("All nmr models (y/n)? [n]\t"); 

      cboxes[3] = (getYesNo('n') ? 0 : 1);

      strcpy(db_path, DBPATH);
      printf("Path for db? [%s]\t", db_path);
      getInput(&db_path);

      getlist:
      printf("File with list? [%s]\t", list1);
      getInput(&list1);

      if (makeLoadListF()) load();
      else goto getlist;
  }

  if(mode=='d') {
      printf("Path for db? [%s] ", db_path);
      getInput(&db_path);

      db->setPath(db_path);

      int nCom1, nCom2, ndf1, ndf2;
      Property id_com("id.com",0);
      nCom2 = id_com.getObjectSize();
      nCom1 = 0;

      printf("%s contains %d compounds\n\n", db_path, nCom2);
      
      printf("Compounds to process? [0,%d]\t", nCom2); 
      getInput(&nCom1, &nCom2);

      ndf1 = 1; ndf2 = NUM_DF;
      DF::printFeatures();
      printf("Features to calculate? [1,%d]\t", NUM_DF); 
      getInput(&ndf1, &ndf2);

      printf("Processing\tCompounds %d to %d\n\t\tFeatures %d to %d\n",
							      nCom1, nCom2,
							      ndf1, ndf2);

      DF::assignDF(db_path, nCom1, nCom2, ndf1, ndf2, log_mode);

      if(log_mode) {
	char buf[100];
	FILE *LOG;
	if((LOG=fopen(log_file,"a")) == NULL) {
	      printf("Error opening logfile!\n");
	      exit(0);
	}
	time_t clock=time(NULL);
	tm *tm_=localtime(&clock);
	strftime(buf, 100, "%D %T", tm_);
	fprintf(LOG, "%s finished loading and df for %d-%d\n",
			 buf, nCom1+1, nCom2);
	fclose(LOG);
      }

  }

  delete [] buffer;
  return;
}
//////////////////////////////////////////////////////////////////
int Tlpdb::getYesNo(char defC) {
	char result = '\0';
	fflush(stdin);
#ifdef T3E 
#define NO_IOCTL 1
#endif

#ifdef SGI
#define NO_IOCTL 1
#endif

#ifndef NO_IOCTL

	struct sgttyb t;

	ioctl( 0, TIOCGETP, &t );

	t.sg_flags &= ~ECHO; ioctl( 0, TIOCSETP, &t );
#else
	result=getchar();
#endif

	while (result != 'y' && result != 'n' && result != '\n')
		result = (char) getchar();

#ifndef NO_IOCTL
	t.sg_flags |= ECHO; ioctl( 0, TIOCSETP, &t );
#endif

	if (result == '\n' || result == defC) {
		printf("%c\n", defC);
		return 1;
	}

	printf("%c\n", result);
	return 0;
}
//////////////////////////////////////////////////////////////////
void Tlpdb::getInput(int *arg1, int *arg2) {
	fflush(stdin);
	char *buffer = new char[255];
	fgets(buffer, 255, stdin);
	if (buffer[0] != '\n') {
		sscanf(buffer, "%d,%d", arg1, arg2);
	}
	delete [] buffer;
}
//////////////////////////////////////////////////////////////////
void Tlpdb::getInput(char **arg) {
	fflush(stdin);
	char *buffer = new char[255];
	fgets(buffer, 255, stdin);
	if (buffer[0] != '\n') {
		sscanf(buffer, "%s\n", (*arg));
	}
	delete [] buffer;
}
//////////////////////////////////////////////////////////////////
Tlpdb::~Tlpdb() {
	delete [] currentPath;
	delete [] files;
	delete [] list2;
	delete [] db_path;
	delete [] log_file;
	delete [] currFiles;
	delete [] path;
}
//////////////////////////////////////////////////////////////////
Tlpdb::Tlpdb()
{
  DF::tlpdb=this;
  con=NULL;
  mol=NULL;
  db=new DB;
  log_mode=0;
  scratch_mode=0;

  currentPath	= new char[100];
  files		= new char[100];
  list1		= new char[100];
  list2		= new char[100];
  db_path	= new char[100];
  log_file	= new char[100];
  currFiles	= new char[100];
  path		= new char[100];

  strcpy(path, OBSDIR);
  strcpy(files, "pdb*.ent");
  strcpy(list1, "list1");
  strcpy(list2, "list2");
  strcpy(db_path, DBPATH);

  cboxes[0]=1; cboxes[1]=0; cboxes[2]=0; cboxes[3]=0; cboxes[4]=0;
  cboxes[5]=1;cboxes[6]=1;

  loadList=NULL; nLoadList=0;
  checkList=NULL; nCheckList=0;
  dbList=NULL; nDBList=0;

//  interruptStatus=0;
}
//////////////////////////////////////////////////////////////////
void Tlpdb::scan()
{
  int i;
  FILE *listF;
  char buffer[100], *buffer_;      
  
  if(nLoadList)
    {
      DeleteArray(&loadList, nLoadList);
      nLoadList=0;
    }
  
  if(nCheckList)
    {
      DeleteArray(&checkList, nCheckList);
      nCheckList=0;
    }
  
  if(nDBList)
    {
      DeleteArray(&dbList, nDBList);
      nDBList=0;
    }

  if(cboxes[2])
    {
      if((listF=fopen(list2, "r"))!=NULL)
	{
	  while(fgets(buffer, 100, listF))
	    {
	      if(strlen(buffer)>0)
		{
		  buffer[strlen(buffer)-1]='\0';
		  AddBuffer(&checkList, buffer, nCheckList);
		  nCheckList++;
		}
	    }
	  fclose(listF);
	}
    }
  
  if(cboxes[4]) {		// check with db
      db = new DB;
      db->setPath(db_path);

      if(!con) con=new Connectivity;
      if(!mol) mol=new Molecule(con, this);
      
      int searchResult=db->testFile("id.com");
      
      if(searchResult)
        {
	  
	  Property id_com("id.com");
	  nDBList=id_com.getObjectSize();
          dbList=new char* [nDBList];
          for(i=0; i<nDBList; i++)
            {
	      dbList[i]=id_com.item1(i, 1);
	    }
        }
      else fprintf(stderr,
		   "Could not find id.com in %s - no db checking!\n",
		   db_path);
  }

  int total=0;

  if(isObs) total += findFiles(path, files, 0);
  if(isCur) total += findFiles(currentPath, currFiles, 0);

  loadList = new char* [total]; 
  nLoadList=0; total = 0;

  if(isObs) findFiles(path, files, 1);
  if(isCur) findFiles(currentPath, currFiles, 1);

  printf("%d entries selected\n", nLoadList);

  int j, jMin, cmpCode, pos = strlen(loadList[0])-7; char *tmp;

  if (sort_flag) {
      printf("Sorting ...\n");

      for (i = 0; i < nLoadList-1; i++) {
	  jMin = i;
	  for(j = i; j < nLoadList; j++) {
	      cmpCode = strcmp(loadList[jMin]+pos, loadList[j]+pos);
	      if (cmpCode == 0) 
		cmpCode = strcmp(loadList[jMin]+pos-1, loadList[j]+pos-1);
		if (cmpCode > 0) jMin=j;
	  }

	  if (i != jMin) {
	      tmp = loadList[jMin];
	      loadList[jMin] = loadList[i];
	      loadList[i] = tmp;
	  }
      }
  }

  if((listF=fopen(list1, "w"))!=NULL) {
      for(i=0; i<nLoadList; i++)
	{
	  fprintf(listF, "%s\n", loadList[i]);
	}
      fclose(listF);
  }
}
//////////////////////////////////////////////////////////////////
int Tlpdb::findFiles(char *path_, char *searchString, int mode) {
  int i, j;
  int searchStatus;
  DIR *handleDir; dirent *d; 
  char *searchDir; struct stat stat_buffer;
  int nFiles = 0;
  
  if(cboxes[0]) {		// include subdirs
    handleDir=opendir(path_);
    if(!handleDir) return nFiles; 

    //if (mode) printf("Searching for subdirs in %s\n", path_);
    
    while((d=readdir(handleDir))) {

      if(!strcmp(d->d_name, ".") || !strcmp(d->d_name, ".."))
	continue;
      
      searchDir=NULL;
      addText(&searchDir, path_);
      addText(&searchDir, PATH_SEPARATOR);
      addText(&searchDir, d->d_name);
      
      if(stat(searchDir, &stat_buffer)==-1) {
	printf("stat failed on path '%s'\n", searchDir);
	exit(0);
      }
      
      if(stat_buffer.st_mode & S_IFDIR)
	nFiles += findFiles(searchDir, searchString, mode);
      delete [] searchDir;
    }
    
    closedir(handleDir);
  }
  
  handleDir=opendir(path_);
  if(!handleDir) return nFiles; 

  //if (mode) printf("Searching entries in %s...", path_);
  
  while((d=readdir(handleDir))) {
    if(file_cmp(d->d_name, searchString)) continue;
    
    if(cboxes[2]) {		// check with second list
      for(i=0; i<nCheckList; i++)
	for(j=0; j<strlen(d->d_name)-strlen(checkList[i])+1; j++)
	  if (!strncasecmp(checkList[i],
			   d->d_name+j,
			   strlen(checkList[i]))) goto checkOK;
      continue;
    }
    
  checkOK:
    if(cboxes[4]) {		// check with db
      
      for(i=0; i<nDBList; i++)
	for(j=0; j<strlen(d->d_name)-3; j++)
	  if(!strncasecmp(dbList[i], d->d_name+j, 4))
	    goto checkCONT;
      
      goto checkOKOK;
      
    checkCONT: continue;
      
    }
    
  checkOKOK:
    
    if(mode) {
      loadList[nLoadList]=NULL;
      addText(&loadList[nLoadList], path_);
      if(path_[strlen(path_)-1]!=*PATH_SEPARATOR)
	addText(&loadList[nLoadList], PATH_SEPARATOR);
      addText(&loadList[nLoadList], d->d_name);
      nLoadList++;
    }
    
    nFiles++;
  }
  
  closedir(handleDir); 
  //if (mode) printf(" found %d entries\n", nFiles);
  return nFiles;
}
//////////////////////////////////////////////////////////////////
int Tlpdb::makeLoadListF() {
  FILE *loadListF=fopen(list1,"r");
  char buffer[100];

  if(!loadListF)
    {
      printf("can't open %s\n", list1); return 0;
    }
  
  if(nLoadList) DeleteArray(&loadList, nLoadList);
  nLoadList=0;

  while(fscanf(loadListF, "%s", buffer)!=EOF)
    {
      AddBuffer(&loadList, buffer, nLoadList); nLoadList++;
    }

  if(!nLoadList) 
    {
      printf("Nothing to load\n");
      return 0;
    }

  return 1;
}
//////////////////////////////////////////////////////////////////
void Tlpdb::load()
{
  int i, handle; 

#ifdef ALPHA
  int startTime = 0, totalTime = 0, thisTime = 0;
#else
  long int startTime = 0, totalTime = 0, thisTime = 0;
#endif
  
  DIR *handleDir; dirent *d; 
//  struct date date_; struct time time_;
  FILE *fhandle;
  
  db->setPath(db_path);

  if(!cboxes[1])
    {
      if(db->testFile("code3.mon") && !scratch_mode) {
	  printf
	    ("Database already exist at this location. Overwrite (y/n)? [n] "); 
	  if (getYesNo('n')) return;
      }

      createProperties();

      if(P_LOCK_PID_MISC) {
	Property lock_pid_misc("lock_pid.misc");
	lock_pid_misc.addItem(-1);
      }

    }

  Property lock_pid_misc("lock_pid.misc");
  int pid = *lock_pid_misc.item4(0);

  if(pid != -1) {
    //if(getpgid(pid) == -1 || pid == getpid()) pid = -1;
  }
  if(pid == -1) {
    pid = getpid();
    lock_pid_misc.setItem(0, pid);
    lock_pid_misc.saveHeader();
  }
  else {
    printf("DB locked by %d\n", pid);
    exit(0);
  }


  if(!con) con=new Connectivity;
  if(!mol) mol=new Molecule(con, this);
  db=mol->db;

  if(scratch_mode) {
    strcpy(mol->id, scratch_id);
  }

  char buf[100];
  char *fileList_;
  int iFile=1, err; long space;

  con->read();

  Property id("id.com", 0);

  for(i=0; i<nLoadList; i++)
    {
      fileList_=loadList[i];
      
      FILE *LOG;

      if(log_mode)
	{
	  LOG=fopen(log_file,"a");
	  time_t clock=time(NULL);
	  tm *tm_=localtime(&clock);
	  strftime(buf, 100, "%D %T", tm_);
	  fprintf(LOG, "%s Loading %s", buf, fileList_);
	}
      else
	{
	  startTime = time(0);
	  printf("Processing [%d] %s", iFile, fileList_);
	  fflush(stdout);
	}

      if((err=mol->makeMolecule(fileList_))) {
	  fprintf(stderr, "Error %d in making [%d]: %s\n",
						 err, iFile, fileList_);
	}

	if (log_mode) {
		if (err) fprintf(LOG, " - failed!");
		fprintf(LOG, "\n");
		fclose(LOG);
	}

	if (err) continue;

      
#ifdef PURIFY
	if (purify_clear_leaks()) {
		purify_new_inuse();
		fprintf(stderr, "purify found leaks!\n");
	}
#endif

      if(!log_mode) {
	  printf(" * WRITING DB *");
	  fflush(stdout);
      }

      mol->writeDB();

      if(!log_mode) {
	thisTime = time(0) - startTime;
	totalTime += thisTime;
	printf(" %2ds (%ds)\n", thisTime, totalTime);
	thisTime = 0;
      }

      iFile++;
/*
   if(interruptStatus)
   {
   interruptStatus=0;
   break;
   }
*/
    }

  delete mol;
  delete con;
/*
  printf("Finish?");
  char c[10];
  scanf("%s", c);
*/
  lock_pid_misc.setItem(0, -1);
  lock_pid_misc.saveHeader();

}
//////////////////////////////////////////////////////////////////
CreatePropertiesList createPropertiesList[]={
//           load|size|type|prop/coll|on/off
"code3.mon",   1, 4, D_INT1, 0, P_CODE3_MON, 
"code1.mon",   1, 1, D_INT1, 0, P_CODE1_MON,
"n_bond.mon",  1, 1, D_INT2, 0, P_N_BOND_MON,
"bond.mon",    1, 0, D_INT2, 0, P_BOND_MON,
"n_atom.mon",  1, 1, D_INT2, 0, P_N_ATOM_MON,
"atom.mon",    1, 0, D_INT1, 0, P_ATOM_MON,
"prev.mon",    1, 1, D_INT2, 0, P_PREV_MON,
"next.mon",    1, 1, D_INT2, 0, P_NEXT_MON,
"type.mon",    1, 1, D_INT1, 0, P_TYPE_MON,

"id.com",      1, 5, D_INT1, 0, P_ID_COM,
"file.com",    1, 0, D_INT1, 0, P_FILE_COM,
"status.com",  1, 1, D_INT1, 0, P_STATUS_COM,

#ifndef MINIMAL
"title.com",   0, 0, D_INT1, 0, P_TITLE_COM,
"compnd.com",  0, 0, D_INT1, 0, P_COMPND_COM,
"source.com",  0, 0, D_INT1, 0, P_SOURCE_COM,
#endif

"date_tex.com",0, 0, D_INT1, 0, P_DATE_TEX_COM,
"date_int.com",1, 1, D_INT4, 0, P_DATE_INT_COM,

#ifndef MINIMAL
"header.com",  0, 0, D_INT1, 0, P_HEADER_COM,
"auth.com",    0, 0, D_INT1, 0, P_AUTH_COM,
"jrnl.com",    0, 0, D_INT1, 0, P_JRNL_COM,
"expdta.com",  1, 1, D_INT1, 0, P_EXPDTA_COM,
"expdta_txt.com",  1, 0, D_INT1, 0, P_EXPDTA_TXT_COM,
"ec.com",      0, 0, D_INT1, 0, P_EC_COM,
"ssbond.com",  1, 0, D_INT1, 0, P_SSBOND_COM,
"site.com",    0, 0, D_INT1, 0, P_SITE_COM,
"ndbmap.com",  1, 0, D_INT1, 0, P_NDBMAP_COM,
#endif

"res.com",     1, 1, D_FLT4, 0, P_RES_COM,
"n_enp.com",   1, 1, D_INT2, 0, P_N_ENP_COM,
"i_enp.com",   1, 1, D_INT4, 0, P_I_ENP_COM,
"n_enc.com",   1, 1, D_INT2, 0, P_N_ENC_COM,
"i_enc.com",   1, 1, D_INT4, 0, P_I_ENC_COM,

"reldat.com",  1, 1, D_INT4, 0, P_RELDAT_COM,
"obs_ids.com", 1, 0, D_INT1, 0, P_OBS_IDS_COM,
"obs.com",     1, 0, D_INT4, 0, P_OBS_COM,
"current.com", 1, 0, D_INT4, 0, P_CURRENT_COM,
"obs_dat.com", 1, 1, D_INT4, 0, P_OBS_DAT_COM,
"spr_ids.com", 1, 0, D_INT1, 0, P_SPR_IDS_COM,
"spr.com",     1, 0, D_INT4, 0, P_SPR_COM,
"spr_dat.com", 1, 1, D_INT4, 0, P_SPR_DAT_COM,

"unitcell.com",0, 6, D_FLT4, 0, P_UNITCELL_COM,
"spcgrp.com",  0, 0, D_INT1, 0, P_SPCGRP_COM,
"zval.com",    1, 1, D_INT2, 0, P_ZVAL_COM,
"rval.com",    1, 1, D_FLT4, 0, P_RVAL_COM,

"name.enc",    1, 0, D_INT1, 0, P_NAME_ENC,
"i_com.enc",   1, 1, D_INT4, 0, P_I_COM_ENC,
"i_enp.enc",   1, 1, D_INT4, 0, P_I_ENP_ENC,
"n_se.enc",    0, 1, D_INT2, 0, P_N_SE_ENC,
"se.enc",      0, 0, D_INT2, 0, P_SE_ENC,

#ifndef MINIMAL
"sen_pdb.enc", 0, 0, D_INT1, 0, P_SEN_PDB_ENC,
#endif

"n_xyz.enc",   0, 1, D_INT4, 0, P_N_XYZ_ENC,
"xyz.enc",     0, 0, D_FLT4, 0, P_XYZ_ENC,
"bfac.enc",    0, 0, D_FLT4, 0, P_BFAC_ENC,
"se_xyz.enc",  0, 0, D_INT4, 0, P_SE_XYZ_ENC,
"xyz_se.enc",  0, 0, D_INT4, 0, P_XYZ_SE_ENC,

"name.enp",    1, 0, D_INT1, 0, P_NAME_ENP,
"i_com.enp",   1, 1, D_INT4, 0, P_I_COM_ENP,
"i_enc.enp",   1, 1, D_INT4, 0, P_I_ENC_ENP,
"type.enp",    1, 1, D_INT1, 0, P_TYPE_ENP,

#ifndef MINIMAL
"mw.enp",      1, 1, D_INT4, 0, P_MW_ENP,
#endif

"n_se.enp",    1, 1, D_INT2, 0, P_N_SE_ENP,

#ifndef MINIMAL
"alpha_c.enp", 1, 1, D_FLT4, 0, P_ALPHA_C_ENP,
"beta_c.enp",  1, 1, D_FLT4, 0, P_BETA_C_ENP,
"alpha_n.enp", 1, 1, D_INT2, 0, P_ALPHA_N_ENP,
"beta_n.enp",  1, 1, D_INT2, 0, P_BETA_N_ENP,
"ss_seg.enp",  0, 0, D_INT1, 0, P_SS_SEG_ENP,
#endif

"seq.enp",     1, 0, D_INT1, 0, P_SEQ_ENP,

#ifndef MINIMAL
"seq_flt.enp", 1, 0, D_INT1, 0, P_SEQ_FLT_ENP,
"iseq.enp",    1, 0, D_INT1, 0, P_ISEQ_ENP,
"iseq_flt.enp",1, 0, D_INT1, 0, P_ISEQ_FLT_ENP,
"c_a.enp",     0, 0, D_FLT4, 0, P_C_A_ENP,
#endif

"k_s.enp",     1, 0, D_INT1, 0, P_K_S_ENP,
"fds.enp",     0, 1, D_FLT4, 0, P_FDS_ENP,

"chi1.enp",    0, 6, D_FLT4, 0, P_CHI1_ENP,
"chi1.com",    1, 1, D_FLT4, 0, P_CHI1_COM,
"fds.com",     1, 1, D_FLT4, 0, P_FDS_COM,

"se_type.enp", 0, 0, D_INT1, 0, P_SE_TYPE_ENP,

#ifndef MINIMAL
"exp.enp",     1, 0, D_FLT4, 0, P_EXP_ENP,
"exp_flt.enp", 1, 0, D_FLT4, 0, P_EXP_FLT_ENP,
"pol.enp",     1, 0, D_FLT4, 0, P_POL_ENP,
"pol_flt.enp", 1, 0, D_FLT4, 0, P_POL_FLT_ENP,
"bfac_flt.enp",1, 0, D_FLT4, 0, P_BFAC_FLT_ENP,

"exp_1b.enp",  1, 0, D_INT1, 0, P_EXP_1B_ENP,
"pol_1b.enp",  1, 0, D_INT1, 0, P_POL_1B_ENP,
"bfac_1b.enp", 1, 0, D_INT1, 0, P_BFAC_1B_ENP,
"sexp_1b.enp", 1, 0, D_INT1, 0, P_SEXP_1B_ENP,
"spol_1b.enp", 1, 0, D_INT1, 0, P_SPOL_1B_ENP,
"shyd_1b.enp", 1, 0, D_INT1, 0, P_SHYD_1B_ENP,
"svol_1b.enp", 1, 0, D_INT1, 0, P_SVOL_1B_ENP,
"siso_1b.enp", 1, 0, D_INT1, 0, P_SISO_1B_ENP,

"exp.pp6",     1, 120, D_INT4, 0, P_EXP_PP6,
"pol.pp6",     1, 120, D_INT4, 0, P_POL_PP6,
"bfac.pp6",    1, 120, D_INT4, 0, P_BFAC_PP6,
"exp.pp8",     1, 160, D_INT4, 0, P_EXP_PP8,
"pol.pp8",     1, 160, D_INT4, 0, P_POL_PP8,
"bfac.pp8",    1, 160, D_INT4, 0, P_BFAC_PP8,

"ks.pp3",      1, 60, D_INT4, 0, P_KS_PP3,
"prop8.ave",   1, 160, D_INT1, 0, P_PROP8_AVE,
"prop8.sts",   1, 56, D_FLT4, 0, P_PROP8_STS,
#endif

"lock_pid.misc", 1, 1, D_INT4, 0, P_LOCK_PID_MISC,

NULL, 0, 0, D_FLT4, 0, 0};
//////////////////////////////////////////////////////////////////
void Tlpdb::createProperties()
{
  Property tmp_p;
  Collection tmp_c;

  for(int i=0; createPropertiesList[i].property; i++)
    if(createPropertiesList[i].isOn) {
      if(createPropertiesList[i].dataClass==0)
	tmp_p.create(createPropertiesList[i].property, 
		   createPropertiesList[i].type,
		   createPropertiesList[i].size,
		   createPropertiesList[i].load_mode);
      if(createPropertiesList[i].dataClass==1)
	tmp_c.create(createPropertiesList[i].property, 
		   createPropertiesList[i].type,
		   createPropertiesList[i].size);
    }
  
  Property i_enp_com("i_enp.com"); i_enp_com.addItem((int4)0);
  Property i_enc_com("i_enc.com"); i_enc_com.addItem((int4)0);
}
//////////////////////////////////////////////////////////////////
