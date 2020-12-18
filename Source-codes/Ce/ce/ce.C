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
#include "../ce/ce.h"
////////////////////////////////////////////////////////////////////
extern char aa[22];
double ran_();
////////////////////////////////////////////////////////////////////
template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
  *array=new A* [Narray1];
  for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};
double vecDist(XYZ *o1, XYZ *e1, XYZ *o2, XYZ *e2);
////////////////////////////////////////////////////////////////////
main(int argc, char **argv) {
  if(argc<2) exit(0);

  CE ce;
  if(!strcmp(argv[1], "-")) {
    ce.scratch_align_ent(argc, argv);
  }
}
////////////////////////////////////////////////////////////////////
void CE::scratch_align_ent(int argc, char **argv) {
  // ce scratch_align_ent pdb_file1 chain_id1 pdb_file2 chain_id2 scratch
  //           1             2           3         4         5       6
 
  printf("\nStructure Alignment Calculator, version 1.02, last modified: Jun 15, 2001.\n\n");
  
  
  if(argc < 7) return;

  char *db_tmp_path, *mkdb_command, *ent1, *ent2;

  setText(&db_tmp_path, argv[6]);
  
  char *pdb_file1 = argv[2];
  char *pdb_file2 = argv[4];

  char *chain_id1 = argv[3];
  char *chain_id2 = argv[5];

      
  char *user_file1, *entity_id1;
  setText(&user_file1, pdb_file1);

  setText(&mkdb_command, "pom/mkDB");
  addText(&mkdb_command, " scratch ");
  addText(&mkdb_command, db_tmp_path);
  addText(&mkdb_command, " ");
  addText(&mkdb_command, user_file1);

  //printf("%s<BR>\n", mkdb_command);
  system(mkdb_command);

  setText(&ent1, "USR1:");
  addText(&ent1, chain_id1);

      
  char *user_file2, *entity_id2;
  setText(&user_file2, pdb_file2);


  setText(&mkdb_command, "pom/mkDB");
  addText(&mkdb_command, " scratch ");
  addText(&mkdb_command, db_tmp_path);
  addText(&mkdb_command, " ");
  addText(&mkdb_command, user_file2);
  addText(&mkdb_command, " USR2 add");

  system(mkdb_command);

  setText(&ent2, "USR2:");
  addText(&ent2, chain_id2);

  toUpperCase(ent1);
  toUpperCase(ent2);

  char *com1, *com2;

  setText(&com1, ent1);
  com1[4] = '\0';

  setText(&com2, ent2);
  com2[4] = '\0';

  DB db;

  db.setPath(db_tmp_path);
  
  Property name_enp("name.enp"), nse_enp("n_se.enp"), ca_enp("c_a.enp"),
    seq_enp("seq.enp"), code3_mon("code3.mon"), se_enc("se.enc"), 
    i_enc_enp("i_enc.enp"), id_com("id.com"),
    comp_enp("compnd.com"), i_com_enp("i_com.enp"), i_enp_com("i_enp.com");

  double d_[20];

  int pos;

  
  int iEnp11, iEnp12, iEnp21, iEnp22;

  int iCom1 = id_com.find(com1);

  if(*chain_id1 != '-') {
    iEnp11 = name_enp.find(ent1); iEnp12 = iEnp11 + 1;
  }
  else {
    iEnp11 = *i_enp_com.item4(iCom1); iEnp12 = *i_enp_com.item4(iCom1+1);
  }
  

  int iCom2 = id_com.find(com2);

  if(*chain_id2 != '-') {
    iEnp21 = name_enp.find(ent2); iEnp22 = iEnp21 + 1;
  }
  else {
    iEnp21 = *i_enp_com.item4(iCom2); iEnp22 = *i_enp_com.item4(iCom2+1);
  }
  
  char *name1, *name2;

  setText(&name1, pdb_file1);
  addText(&name1, ":");
  addText(&name1, chain_id1);
  
  setText(&name2, pdb_file2);
  addText(&name2, ":");
  addText(&name2, chain_id2);


  if(iEnp11 == -1) printf("Chain %s not found\n", name1);
  if(iEnp12 == -1) printf("Chain %s not found\n", name2);

  if(iEnp11==-1 || iEnp12==-1) exit(0);

  int isPrint = 2;
  
  for(int ienp1 = iEnp11; ienp1 < iEnp12; ienp1++) 
    for(int ienp2 = iEnp21; ienp2 < iEnp22; ienp2++) {
      

      setText(&name1, pdb_file1);
      addText(&name1, ":");
      addText(&name1, name_enp.item1(ienp1) + 5);
  
      setText(&name2, pdb_file2);
      addText(&name2, ":");
      addText(&name2, name_enp.item1(ienp2) + 5);

      int nse1=*nse_enp.item2(ienp1);
      int nse2=*nse_enp.item2(ienp2);


      char *seq1=seq_enp.item1(ienp1, 1);
      char *seq2=seq_enp.item1(ienp2, 1);
      int2 *se1=se_enc.item2(*i_enc_enp.item4(ienp1), 1);
      int2 *se2=se_enc.item2(*i_enc_enp.item4(ienp2), 1);


      int *align_se1=new int [nse1+nse2];
      int *align_se2=new int [nse1+nse2];
      int lcmp;

      XYZ *ca1 = arrayToXYZ(ca_enp.itemf(ienp1), nse1);
      XYZ *ca2 = arrayToXYZ(ca_enp.itemf(ienp2), nse2);
  

      double z;
  

      z=ce_1(name1, name2, ca1, ca2, nse1, nse2, align_se1, align_se2, lcmp, 
	     8, 3.0, 4.0, d_, isPrint);

      isPrint = 1;

      if(lcmp>0) {
	int lsim = 0, lali = 0;
	for(int l=0; l<lcmp; l++) 
	  if(align_se1[l] != -1 && align_se2[l] != -1) {
	    if(seq1[align_se1[l]] == seq2[align_se2[l]]) lsim++;
	    lali++;
	  }
	printf("Sequence identities = %.1f%%", lsim*100.0/lali);


	int lstep = 70;
    
	for(int l=0; l<lcmp; l+=lstep) {
	  printf("\n");
	  for(int ie=0; ie<2; ie++) {
	    int ienp = (ie == 0 ? ienp1 : ienp2);

	    char *seq = (ie == 0 ? seq1 : seq2);

	    int *align_se = (ie == 0 ? align_se1 : align_se2);

	    printf("\n%8.8s ", (ie == 0 ? "Chain 1:" : "Chain 2:"));
	    int ip=-1;
	    for(int l_=l; l_<l+lstep && l_<lcmp; l_++) 
	      if(align_se[l_]!=-1) {
		ip=align_se[l_]+1; break;
	      }
	    if(ip!=-1) printf("%4d ", ip);
	    else printf("     ");


	    for(int l_=l; l_<l+lstep && l_<lcmp; l_++) 
	      printf("%c", align_se[l_]==-1?'-':(seq[align_se[l_]]));
	    
	  }
	}
	printf("\n");
      


	printf("\n     X2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
	       d_[0], d_[1], d_[2], d_[9]);
	printf("     Y2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
	       d_[3], d_[4], d_[5], d_[10]);
	printf("     Z2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
	       d_[6], d_[7], d_[8], d_[11]);

      }

      if(nse1 > 0) {
	delete [] ca1;
      }
      
      if(nse2 > 0) {
	delete [] ca2;
      }
      
      if(nse1+nse2 > 0) {
	delete [] align_se1;
	delete [] align_se2;
      }

    }
}
///////////////////////////////////////////////////////////////////////////
