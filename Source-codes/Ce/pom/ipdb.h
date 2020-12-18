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
// typedef int clockid_t;
#ifndef H_IPDB
#define H_IPDB

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/param.h>
#include <time.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pom.h"
#include "water.h"
//////////////////////////////////////////////////////////////////
int file_cmp(char *, char *);
void compressText(char **);
double ran_();
//////////////////////////////////////////////////////////////////
class Connectivity;
class Molecule;
class TlpdbMessage;
//////////////////////////////////////////////////////////////////
class Tlpdb 
{
public:
  Connectivity * con;
  Molecule * mol;
  DB *db;

  char *path, *currentPath, *files, *list1, *list2, *db_path,
       *log_file, *currFiles;

  int nLoadList, nCheckList, nDBList;
  char **loadList, **checkList, **dbList;
  int cboxes[7], isObs, isCur;
  int sort_flag;
  int log_mode, scratch_mode;
  char scratch_id[5];

  void update(int argc, char **argv);
  void scratch(int argc, char **argv);
  void dialog();

  Tlpdb();
  ~Tlpdb();
  void scan();
  int makeLoadListF();
  void load();
  int findFiles(char *, char *, int);
  void HashTable();
  void createProperties();

  int getYesNo(char);
  void getInput(int *, int *);
  void getInput(char **);

};
//////////////////////////////////////////////////////////////////
struct CreatePropertiesList
{
  char *property;
  int load_mode;
  int4 size;
  DataType type;
  int dataClass;
  int isOn;
};
//////////////////////////////////////////////////////////////////
class Monomer
{
 public:
  double bondDist, bondDist2;
  char name[4];
  int2 nAtom;
  flt4 *xyz, *tmf;
//  flt4 *occ;
  char *atomNames;
  int index;
  char pdbChain;
  char pdbNum[6];
  int pdbBreak;
  int2 nBond;
  int2 *bond;
  int prevA, nextA;

  Monomer();
  void makeBonds(Monomer *, Monomer *);
  int checkBond(flt4 *, flt4 *);
  void removeBond(int);
  double dist2(flt4 *, flt4 *);
  void checkTriangles();
  int findAtom(char *);
  ~Monomer();
};
//////////////////////////////////////////////////////////////////
class Molecule
{
 public:
  Tlpdb * tlpdb;
  Connectivity * con;
  FILE * pdbFile;
  DB *db;
  int2 nChain;
  int2 nPolyChain;
  int2 * nMon;
  int2 ** mon;
  char **ks;
  flt4 **exp;
  flt4 **pol;
  int2 **monXYZ;
  int4 * comNXYZ;
  flt4 ** comXYZ;
//  flt4 ** comOcc;
  flt4 ** comTmf;
  char id[5], date[10], * header, * compnd, * title, * source, * author,
       * jrnl, *file, *site, *ssBond, *ndbMap;

  char buffer[100];
  int bufferStatus;
  int modelStatus;
  flt4 resolution;
  int1 expdta;
  char *expdta_txt;

  char **obsIds, **sprIds;
  int4 dateInt, relDate, obsDate, sprDate;
  int obsCount, sprCount;

  flt4 *cell;
  char *spcgrp;
  int2 Zval;
  flt4 rvalue;

  char *chainCode;
  char ***seNum;
  int chainStatus;
  char *** seq;
  int seqNchain;
  int * seqL;
  char * seqChainCode;

  double energyCutoff;
  flt4 distanceCutoff;

  Property n_xyz_enc;
  Property xyz_enc;
  Property se_xyz_enc;
  Property xyz_se_enc;

  Property id_com;
  Property file_com;
  Property compnd_com;
  Property title_com;
  Property source_com;
  Property date_tex_com;
  Property date_int_com;
  Property header_com;
  Property auth_com;
  Property jrnl_com;
  Property name_enc;
  Property name_enp;
  Property n_se_enc;
  Property n_se_enp;
  Property se_enc;
  Property n_enp_com;
  Property n_enc_com;
  Property i_enp_com;
  Property i_enc_com;

  Property site_com, ssbond_com, ndbmap_com;

  Property rel_dat_com;

  Property obs_ids_com;
  Property obs_com;
  Property obs_date_com;
  Property spr_ids_com;
  Property spr_com;
  Property spr_date_com;
  Property current_com;

  Property unit_cell_com;
  Property space_grp_com;
  Property z_val_com;

  Property i_com_enc;
  Property i_com_enp;

  Property i_enc_enp;
  Property i_enp_enc;

  Property sen_pdb_enc;

  Property bfac_enc;

  Property exp_enp;
  Property pol_enp;
  Property k_s_enp;
  Property seq_enp;

  Property res_com;
  Property rval_com;
  Property expdta_com;
  Property expdta_txt_com;

  Property exp_flt_enp;
  Property pol_flt_enp;
  Property seq_flt_enp;
  Property iseq_flt_enp;
  Property iseq_enp;
  Property se_type_enp;

  Molecule(Connectivity *, Tlpdb *);
  int makeMolecule(char *);
  void writeDB();
  void addChain(Monomer *);
  void addMonomer(Monomer *);
  void addHeader(char *);
  void addObslte(char *);
  void addRelDate(char *);
  void addSprsde(char *);
  void addCompnd(char *);
  void addSource(char *);
  void addAuthor(char *);
  void addJrnl(char *);
  void addSeqres(char *);
  void addCryst(char *);
  void matchAtomSeqres();
  void sortChains();
  int assignKS();
  int assignEnv();
  int assignFiltProp();
  int HBondingEnergy(int, int, char *, int, char *);
  int findAtom(int, char *);
  Monomer * getMonomer();
};
//////////////////////////////////////////////////////////////////
class Connectivity
{
 public:
  DB *db;
  Property monCode3;
  Property monCode1;
  Property monNbond;
  Property monBond;
  Property monNatom;
  Property monAtom;
  Property monPrev;
  Property monNext;
  Property monType;
  
  Connectivity();
  int read();
  int search(Monomer *);
};
//////////////////////////////////////////////////////////////////
#endif
