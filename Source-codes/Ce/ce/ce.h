///////////////////////////////////////////////////////////////////////////
//  Combinatorial Extension Algorithm for Protein Structure Alignment    //
//                                                                       //
//  Authors: I.Shindyalov & P.Bourne                                     //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  1997-2000   The Regents of the University of California
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
#include "../ce/cmp_util.h"

#ifndef H_CE
#define H_CE
////////////////////////////////////////////////////////////////////
//                        ce.C                                    //
////////////////////////////////////////////////////////////////////
class DomRep;

class CE {
 public:
  // data
  int winSize;
  double expThr;

  int isDmapOpen;
  Property dmap5, se_dmap5, n_dmap5, exp_flt;
  int4 *dmap, *se_dmap, *n_dmap;
  flt4 *env_exp;


  CE() {
    isDmapOpen=0;
  }
  //
  //  Commands: 
  //
  //  ce sets (ce_calc.C):
  //
  //  purify_entity_set 
  //  check_list moose_db_path list_file [n_print_length]
  //  diff_list moose_db_path [y_print_list] [y_print_rep]
  //  ran_list e_list_file n_ran - random list for ce run 
  //
  //  ce statistics (ce_calc.C):
  //
  //  calc_dist_rmsd moose_db_path hs_list_file
  //  analyze_dist_rmsd moose_db_path hs_list_file frag_file
  //  calc_rmsd_stat moose_db_path hs_list_file
  //
  //  ce queries (ce_align.C):
  //
  //  align_ent moose_db_path ent_name1 ent_name2
  //  align_ent_cgi moose_db_path scratch_path ent_name1 ent_name2 sim_level
  //  align_sets moose_db_path set1_name set2_name
  //  align_sets_1_1 moose_db_path set1_name set2_name
  //  align_set moose_db_path set_name
  //  a_to_s_cgi moose_db_path scratch_path ent_name1 ent_name2 < alignment
  //
  //  ce all-to-all (ce_all.C):
  //
  //  all_df moose_db_path all_to_all_db_path
  //  all moose_db_path all_to_all_db_path
  //  all_print moose_db_path all_to_all_db_path
  //  all_print_aln moose_db_path all_to_all_db_path ent_name1 ent_name2 
  //
 public:
  
  
  // queries (ce.C)

  void align_ent(int argc, char **argv);
  void scratch_align_ent(int argc, char **argv);
  void align_ent_cgi(int argc, char **argv);
  void a_to_s_cgi(int argc, char **argv);
  void align_sets(int argc, char **argv);
  void align_sets_1_1(int argc, char **argv);
  void align_set(int argc, char **argv);
  
  // calculations (ce_calc.C)
  
  void purify_entity_set(int argc, char **argv);
  void check_list(int argc, char **argv);
  void ran_list(int argc, char **argv);
  void diff_list(int argc, char **argv);
  void calc_dist_rmsd(int argc, char **argv);
  void analyze_dist_rmsd(int argc, char **argv);
  void calc_rmsd_stat(int argc, char **argv);
  void calc_dist_rmsd_old(int argc, char **argv);
  
  // all-to-all (ce_all.C)
  void all_df(int argc, char **argv);
  void all_print(int argc, char **argv);
  void all_print_aln(int argc, char **argv);
  void all_update(int argc, char **argv);
  void all_fix_rep(int argc, char **argv);
  void ce_db_diff(int argc, char **argv);
  void all_fix_diff(int argc, char **argv);

  // ce fssp (ce_fssp.C)
  void ce_fssp_make_db(int argc, char **argv);
  void ce_fssp_aln(int argc, char **argv);
  void ce_fssp_cmp_db(int argc, char **argv);
  void ce_fssp_print(int argc, char **argv);


  // swiss-prot (ce_sp.C)
  void make_db_sp(int argc, char **argv);
  void make_sp_nb(int argc, char **argv);
  void sp_cont_stat(int argc, char **argv);

  // domains (ce_dom.C)
  void ce_dom_stat(int argc, char **argv);
  void ce_dom_assign(int argc, char **argv);
  void get_nb(int irep, int icl, 
	      DomRep *dr, int *rep_tab, int i_rec, int n_rec, 
	      int irep2_, int *list_);
  
};
double vecDist(XYZ *o1, XYZ *e1, XYZ *o2, XYZ *e2);
////////////////////////////////////////////////////////////////////
//                  ce_align.C                                    //
////////////////////////////////////////////////////////////////////
double ce_1(char *name1, char *name2, XYZ *ca1, XYZ *ca2, int nse1, 
	    int nse2, int *align_se1, 
	    int *align_se2, int &lcmp, int winSize,
	    double rmsdThr, double rmsdThrJoin_, double d_[20], int isPrint=1);
double ce_a1(char *name1, char *name2, XYZ *ca1, XYZ *ca2, int nse1, 
	     int nse2, int *align_se1, 
	     int *align_se2, int &lcmp, double rmsdThrJoin_, int isPrint=1);
void origin3(XYZ *p, double *d);
void rot(double m1[3][3], double m2[3][3]);
XYZ rotXYZ(XYZ p, double r[3][3]);
XYZ transXYZ(XYZ p, double d[20]);
/////////////////////////////////////////////////////////////////////
int ce_all_ext(char *ata_path);
void ce_all_mv(char *path, char *file);
void all(char *report_file);
void all_fix_report(int argc, char **argv);
void ce_all_verify(int argc, char **argv);
/////////////////////////////////////////////////////////////////////
double zScore(int winSize, int nTrace, double score);
double zGaps(int winSize, int nTrace, int nGaps);
double zStrAlign(int winSize, int nTrace, double score, int nGap);
double zToP(double z);
double pToZ(double p);
double zByZ(double z1, double z2);
double trace_rmsd(int *trace1, int *trace2, int nTrace, int winSize,
		  XYZ *ca1, XYZ *ca2);
double zStrAlign(int4 *align_se, int lcmp, XYZ *ca1, XYZ *ca2);
/////////////////////////////////////////////////////////////////////
void parseFsspId(char *id);
void cmpDbLoad(DB *db1, char *dbCmp, char mask, double zThr, int *isRep1,
	       int **enpToRep1_);
void cmpDbPrintDiff(DB *db1, char *dbCmp, char mask, double zThr, int *isRep2);
double dpAlignFrag(double **mat, int nSeq1, int nSeq2, int *align_se1, 
		   int *align_se2, int &lAlign, double gapI, 
		   double gap(int i), 
		   int isGlobal1, int isGlobal2, int fragSize);
double gap(int i);
/////////////////////////////////////////////////////////////////////
double aac_dist(int *, int *);
double dipep_dist(int *, int *);
/////////////////////////////////////////////////////////////////////
class ScoreTable {
 public:
  double pScore;
  double pRmsd;
  double pSim;
  int    pLali;
  int    pLgap;

  double cScore;
  double cRmsd;
  int ialn;
  int *align_se_nb;
  int lcmp_nb;
  int isExcluded;

  int nali0;
  int nali3;
  int lali_12;
  int lali_nb;
  int nPosDist;

  int irank;

  int lcmp;
  int *align_se1, *align_se2;

  ScoreTable() {
    isExcluded = 0;
    ialn = -1;
  }
};
/////////////////////////////////////////////////////////////////////
class EntityList {
 public:
  int iEnp;
  int iEnp_n;

  int nse;
  char *seq;
  double maxScore;
  double pscoreAv;

  XYZ *ca;
  char *ks;
  int *ksi;
  double ksc[3];
  float *exp;

  EntityList() {
    iEnp = -1; 
    iEnp_n = -1;
    pscoreAv = 1.0;
  }
};
/////////////////////////////////////////////////////////////////////
class DomCluster;
/////////////////////////////////////////////////////////////////////
class DomNeighbor {
 public:
  int iRep1, iRep2;
  int iEnp1, iEnp2;
  int iEnp1_m, iEnp2_m;
  int ialn;
  int *align_rep, nse_rep;
  int isSel;
  int p11, p12, p21, p22;
  int lali;
  double z;

  DomNeighbor() {
     iRep1 = -1; iRep2 = -1; 
     iEnp1 = -1; iEnp2 = -1; isSel = 0; ialn = -1;
     p11 = -1; p12 = -1; p21 = -1; p22 = -1;
     lali = 0; z = 0;
     align_rep = NULL; nse_rep = 0;
  };
  ~DomNeighbor() {
    if(align_rep) delete [] align_rep;
  };
  DomNeighbor(const DomNeighbor& dn){deepCopy(dn);};
  DomNeighbor& operator=(const DomNeighbor& dn) {deepCopy(dn);};

  void deepCopy(const DomNeighbor& dn) {
    iRep1 = dn.iRep1; iRep2 = dn.iRep2; iEnp1 = dn.iEnp1; iEnp2 = dn.iEnp2; 
    iEnp1_m = dn.iEnp1_m; iEnp2_m = dn.iEnp2_m; ialn = dn.ialn;
    nse_rep = dn.nse_rep;
    if(dn.align_rep) {
      align_rep = new int [nse_rep];
      for(int is = 0; is < nse_rep; is++) align_rep[is] = dn.align_rep[is];
    }
    isSel = dn.isSel; 
    p11 = dn.p11; p12 = dn.p12; p21 = dn.p21; p22 = dn.p22;
    lali = dn.lali; z = dn.z;
  };
};
/////////////////////////////////////////////////////////////////////
class DomCluster {
 public:
  int nDn;
  DomNeighbor **dn;
  int isSel;
  
  DomCluster() {
    nDn = 0; isSel = 0; 
  }
  
  DomCluster(const DomCluster& dc){deepCopy(dc);};
  DomCluster& operator=(const DomCluster& dc) {deepCopy(dc);};

  void deepCopy(const DomCluster& dc) {
    nDn = dc.nDn; 
    if(nDn > 0) {
      dn = new DomNeighbor* [nDn];
      for(int idn = 0; idn < nDn; idn++) dn[idn] = dc.dn[idn];
    }
    isSel = 0;
  };

  void add(DomNeighbor *dn_);
};
/////////////////////////////////////////////////////////////////////
class DomRep {
 public:
  int nDn;
  int nDc;
  DomNeighbor *dn;
  DomCluster *dc;
  int *seq_use, nse;

  DomRep() {
    nDn = 0; nDc = 0;
  };

  void add(DomNeighbor& dn_);
  void add(DomCluster& dc_);
};
/////////////////////////////////////////////////////////////////////
#endif
