///////////////////////////////////////////////////////////////////////////
//  Combinatorial Extension Algorithm for Protein Structure Alignment    //
//                                                                       //
//  Authors: I.Shindyalov & P.Bourne                                     //
///////////////////////////////////////////////////////////////////////////
/*
 
Copyright (c)  1997-1999   The Regents of the University of California
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
#include "../ce/cmp_util.h"
//#include <sys/ddi.h>
///////////////////////////////////////////////////////////////////////////
char aa[23]="AVLICMPFYWDNEQHSTRKG*-";
///////////////////////////////////////////////////////////////////////////
template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
  *array=new A* [Narray1];
  for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};
///////////////////////////////////////////////////////////////////////////
template <class A> void DeleteArray(A *** array, int Narray)
{
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};
////////////////////////////////////////////////////////////////////
void readList(char *file, char ***list, int& nList)
{
  char buffer[100];
  nList=0;
  FILE *flist=fopen(file, "r");
  if(flist==NULL)
    {
      printf("%s opening problem\n", file); exit(0);
    }

  while(fscanf(flist, "%s", buffer)!=EOF) 
    {
      AddBuffer(list, buffer, nList); nList++;
    }
  fclose(flist);
  //printf("%d entities read from %s\n", nList, file);
}
////////////////////////////////////////////////////////////////////
void scanList(char *buffer, char ***list, int& nList)
{
  nList=0;
  for(char *tok=strtok(buffer, ", "); tok; tok=strtok(NULL, ", "))
    {
      AddBuffer(list, tok, nList); nList++;
    }
}
///////////////////////////////////////////////////////////////////////////
void toUpperCase(char *buffer) {
  for(int i=0; i<strlen(buffer); i++) 
    if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
///////////////////////////////////////////////////////////////////////////
void toLowerCase(char *buffer) {
  for(int i=0; i<strlen(buffer); i++) 
    if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
///////////////////////////////////////////////////////////////////////////
/*
int gonnet_aa[]={4, 16, 15, 14, 0, 13, 3, 17, 18, 19, 7, 6, 8, 9, 10, 1, 2, 11,
		12, 5};
double gonnet[]={16.7, 5.3, 4.7, 2.1, 5.7, 3.2, 3.4, 2.0, 2.2, 2.8, 3.9, 3.0, 
		2.4, 4.3, 4.1, 3.7, 5.2, 4.4, 4.7, 4.2,
		7.4, 6.7, 5.6, 6.3, 5.6, 6.1, 5.7, 5.4, 5.4, 5.0, 5.0, 5.3,
		3.8, 3.4, 3.1, 4.2, 2.4, 3.3, 1.9,
		7.7, 5.3, 5.8, 4.1, 5.7, 5.2, 5.1, 5.2, 4.9, 5.0, 5.3, 4.6,
		4.6, 3.9, 5.2, 3.0, 3.3, 1.7,
		12.8, 5.5, 3.6, 4.3, 4.5, 4.7, 5.0, 4.1, 4.3, 4.6, 2.8, 2.6,
		2.9, 3.4, 1.4, 2.1, 0.2, 
		7.6, 5.7, 4.9, 4.9, 5.2, 5.0, 4.4, 4.6, 4.8, 4.5, 4.4, 4.0,
		5.3, 2.9, 3.0, 1.6,
		11.8, 5.6, 5.3, 4.4, 4.2, 3.8, 4.2, 4.1, 1.7, 0.7, 0.8, 1.9,
		0.0, 1.2, 1.2,
		9.0, 7.4, 6.1, 5.9, 6.4, 5.5, 6.0, 3.0, 2.4, 2.2, 3.0, 2.1,
		3.8, 1.6,
		9.9, 7.9, 6.1, 5.6, 4.9, 5.7, 2.2, 1.4, 1.2, 2.3, 0.7, 2.4,
		0.0,
		8.8, 6.9, 5.6, 5.6, 6.4, 3.2, 2.5, 2.4, 3.3, 1.3, 2.5, 0.9,
		7.9, 6.4, 6.7, 6.7, 4.2, 3.3, 3.6, 3.7, 2.6, 3.5, 2.5,
		11.2, 5.8, 5.8, 3.9, 3.0, 3.3, 3.2, 5.1, 7.4, 4.4,
		9.9, 7.9, 3.5, 2.8, 3.0, 3.2, 2.0, 3.4, 3.6,
		8.4, 3.8, 3.1, 3.1, 3.5, 1.9, 3.1, 1.7,
		9.5, 7.7, 8.0, 6.8, 6.8, 5.0, 4.2,
		9.2, 8.0, 8.3, 6.2, 4.5, 3.4,
		9.2, 7.0, 7.2, 5.2, 4.5,
		8.6, 5.3, 4.1, 2.6,
		12.2, 10.3, 8.8,
		13.0, 9.3,
		19.4};
main()
{
  int k1, k2;
  for(int i=0; i<20; i++)
    {
      printf("\n");
      for(int j=0; j<20; j++)
	{
	  k1=19-(gonnet_aa[i]<gonnet_aa[j]?gonnet_aa[i]:gonnet_aa[j]);
	  k2=19-(gonnet_aa[i]>gonnet_aa[j]?gonnet_aa[i]:gonnet_aa[j]);
	  printf("%4.1f, ", gonnet[210-k1*(k1+1)/2-k2-1]);
	}
    }
  printf("\n");
}
*/
///////////////////////////////////////////////////////////////////////////
// 4hhb:a and 4hhb:b
/*
char *seq1="VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"; 
char *seq2="VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"; 

char *seq1="GGGAAAAAAAVVVLLLLLLLIIIIIIIGGG";
char *seq2="CCCAAAAAAALLLLLLLWWWIIIIIIICCC"; 

*/

///////////////////////////////////////////////////////////////////////////
/*
main() {
  int ise1, ise2, nse1=strlen(seq1), nse2=strlen(seq2), i;
  double **mat=new double* [nse1], score;
  int *align_se1=new int [nse1+nse2], *align_se2=new int [nse1+nse2];
  for(ise1=0; ise1<nse1; ise1++) mat[ise1]=new double [nse2];
  for(ise1=0; ise1<nse1; ise1++)
    for(ise2=0; ise2<nse2; ise2++)
      *(mat[ise1]+ise2)=(seq1[ise1]==seq2[ise2]?1.0:0.0);

  int lcmp;
  score=dpAlign(mat, nse1, nse2, align_se1, align_se2, lcmp,
		1.0, 0.1, 0, 0);

  
  printf("sum=%.1f nse1=%d nse2=%d lcmp=%d\n", score, nse1, nse2, lcmp);
 
  for(i=0; i<lcmp; i++) printf("%3d", align_se1[i]);
  printf("\n");
  for(i=0; i<lcmp; i++) printf("%3d", align_se2[i]);
  printf("\n");
  for(i=0; i<lcmp; i++) printf("%c", align_se1[i]!=-1?seq1[align_se1[i]]:'-');
  printf("\n");
  for(i=0; i<lcmp; i++) printf("%c", align_se2[i]!=-1?seq2[align_se2[i]]:'-');
  printf("\n");
}
*/
///////////////////////////////////////////////////////////////////////////
double dpAlign(double **mat, int nSeq1, int nSeq2, int *align_se1, 
	       int *align_se2, int &lAlign, double gapI, double gapE, 
	       int isGlobal1, int isGlobal2) {
  int i, j, is, js, iMax, jMax, k, ge=(gapE!=0.0?1:0);
  double sum, sum_ret, sum_brk;
  char **brk_flg=new char* [nSeq1];
  for(i=0; i<nSeq1; i++) brk_flg[i]=new char [nSeq2];
  /*  
  for(i=0; i<nSeq1; i++)
   {
     printf("\n");
     for(j=0; j<nSeq2; j++)
       {
	 printf("%4d", (int)(*(mat[i]+j)*10));
       }
   }
 printf("\n\n\n");
 */
  if(!ge)
    {
      for(i=nSeq1-1; i>=0; i--)
        for(j=nSeq2-1; j>=0; j--)
          {
	    *(brk_flg[i]+j)=0;
	    if(j<nSeq2-1 && i<nSeq1-1) 
	      {
		sum=*(mat[i+1]+j+1);
	      }
	    else
	      {
		sum=0.0;
		if((isGlobal1 && i!=nSeq1-1) || (isGlobal2 && j!=nSeq2-1)) 
		  sum=-gapI;
	      }
	    if(j+1<nSeq2)
	      for(k=i+2; k<nSeq1; k++)
		{
		  if(*(mat[k]+j+1)-gapI>sum)
		    sum=*(mat[k]+j+1)-gapI;
		}
	    if(i+1<nSeq1)
	      for(k=j+2; k<nSeq2; k++)
		{
		  if(*(mat[i+1]+k)-gapI>sum)
		    sum=*(mat[i+1]+k)-gapI;
		}
	    sum+=*(mat[i]+j);
	    sum_brk=(isGlobal1?-gapI:0.0)+(isGlobal2?-gapI:0.0);
	    if(sum<sum_brk) 
	      {
		sum=sum_brk;
		*(brk_flg[i]+j)=1;
	      }
	    *(mat[i]+j)=sum;
	 }
    }
  else
    {
      for(i=nSeq1-1; i>=0; i--)
        for(j=nSeq2-1; j>=0; j--)
          {
	    *(brk_flg[i]+j)=0;
	    if(j<nSeq2-1 && i<nSeq1-1) 
	      {
		sum=*(mat[i+1]+j+1);
	      }
	    else 
	      {
		sum=0.0;
		if(isGlobal1 && i!=nSeq1-1) sum=-gapI-gapE*(nSeq1-i-1);
		if(isGlobal2 && j!=nSeq2-1) sum=-gapI-gapE*(nSeq2-j-1);
	      }
	    if(j+1<nSeq2)
	      for(k=i+2; k<nSeq1; k++)
		if(*(mat[k]+j+1)-gapI-gapE*(k-i-1)>sum)
		  sum=*(mat[k]+j+1)-gapI-gapE*(k-i-1);
	    if(i+1<nSeq1)
	      for(k=j+2; k<nSeq2; k++)
		if(*(mat[i+1]+k)-gapI-gapE*(k-j-1)>sum)
		  sum=*(mat[i+1]+k)-gapI-gapE*(k-j-1);
	    sum+=*(mat[i]+j);
	    sum_brk=(isGlobal1?-gapI-gapE*(nSeq1-1-i):0.0)+
	      (isGlobal2?-gapI-gapE*(nSeq2-1-j):0.0);
	    if(sum<sum_brk) 
	      {
		sum=sum_brk;
		*(brk_flg[i]+j)=1;
	      }
	    *(mat[i]+j)=sum;
	 }
    }
  /*  
 for(i=0; i<nSeq1; i++)
   {
     printf("\n");
     for(j=0; j<nSeq2; j++)
       {
	 printf("%4d", (int)(*(mat[i]+j)*10));
       }
   }
 printf("\n\n\n");
 for(i=0; i<nSeq1; i++)
   {
     printf("\n");
     for(j=0; j<nSeq2; j++)
       {
	 printf("%4d", (int)(*(brk_flg[i]+j)));
       }
   }
 // exit(0);
 */

  is=0; js=0; lAlign=0;
// no nc-end penalty - begin 
  sum_ret=*(mat[0]);
  for(i=0; i<nSeq1; i++)
    for(j=0; j<nSeq2; j++)
      {
	if(i==0 && j==0) continue;
	sum=*(mat[i]+j);
	if(isGlobal1) sum+=-gapI-gapE*i;
	if(isGlobal2) sum+=-gapI-gapE*j;
	if(sum>sum_ret) 
	  {
	    sum_ret=sum;
	    is=i; js=j;
	  }
      }
  //for(k=0; k<is; k++) align1[k]=-1;
  //for(k=0; k<js; k++) align2[k]=-1;
// no nc-end penalty - end 

  for(i=is, j=js; i<nSeq1 && j<nSeq2; i++, j++)
    {
      iMax=i; jMax=j;
      sum=*(mat[i]+j);
      if(!ge)
        {
          for(k=i+1; k<nSeq1; k++)
	    if(*(mat[k]+j)-gapI>sum)
	      {
	        iMax=k; jMax=j;
		sum=*(mat[k]+j)-gapI;
	      }
      
          for(k=j+1; k<nSeq2; k++)
	    if(*(mat[i]+k)-gapI>sum)
	      {
	        iMax=i; jMax=k;
		sum=*(mat[i]+k)-gapI;
	      }
        }
      else
        {
	  for(k=i+1; k<nSeq1; k++)
	    if(*(mat[k]+j)-gapI-gapE*(k-i)>sum)
	      {
	        iMax=k; jMax=j;
		sum=*(mat[k]+j)-gapI-gapE*(k-i);
	      }
      
          for(k=j+1; k<nSeq2; k++)
	    if(*(mat[i]+k)-gapI-gapE*(k-j)>sum)
	      {
	        iMax=i; jMax=k;
		sum=*(mat[i]+k)-gapI-gapE*(k-j);
	      }
	}

      //printf("%d %d\n", iMax, jMax);
      for(k=i; k<iMax; k++, i++) {
	align_se1[lAlign]=k;
	align_se2[lAlign]=-1; lAlign++;
      }

      for(k=j; k<jMax; k++, j++) {
	align_se1[lAlign]=-1;
	align_se2[lAlign]=k; lAlign++;
      }

      align_se1[lAlign]=iMax;
      align_se2[lAlign]=jMax; lAlign++;
      if(*(brk_flg[i]+j)==1) break;
    }

  for(i=0; i<nSeq1; i++) delete [] brk_flg[i];
  delete [] brk_flg;

  return(sum_ret);
}
///////////////////////////////////////////////////////////////////////////
double aaem[20][20]={ 
 7.6,  5.3,  4.0,  4.4,  5.7,  4.5,  5.5,  2.9,  3.0,  1.6,  4.9,  4.9,  5.2,  5.0,  4.4,  6.3,  5.8,  4.6,  4.8,  5.7, 
 5.3,  8.6,  7.0,  8.3,  5.2,  6.8,  3.4,  5.3,  4.1,  2.6,  2.3,  3.0,  3.3,  3.7,  3.2,  4.2,  5.2,  3.2,  3.5,  1.9, 
 4.0,  7.0,  9.2,  8.0,  3.7,  8.0,  2.9,  7.2,  5.2,  4.5,  1.2,  2.2,  2.4,  3.6,  3.3,  3.1,  3.9,  3.0,  3.1,  0.8, 
 4.4,  8.3,  8.0,  9.2,  4.1,  7.7,  2.6,  6.2,  4.5,  3.4,  1.4,  2.4,  2.5,  3.3,  3.0,  3.4,  4.6,  2.8,  3.1,  0.7, 
 5.7,  5.2,  3.7,  4.1, 16.7,  4.3,  2.1,  4.4,  4.7,  4.2,  2.0,  3.4,  2.2,  2.8,  3.9,  5.3,  4.7,  3.0,  2.4,  3.2, 
 4.5,  6.8,  8.0,  7.7,  4.3,  9.5,  2.8,  6.8,  5.0,  4.2,  2.2,  3.0,  3.2,  4.2,  3.9,  3.8,  4.6,  3.5,  3.8,  1.7, 
 5.5,  3.4,  2.9,  2.6,  2.1,  2.8, 12.8,  1.4,  2.1,  0.2,  4.5,  4.3,  4.7,  5.0,  4.1,  5.6,  5.3,  4.3,  4.6,  3.6, 
 2.9,  5.3,  7.2,  6.2,  4.4,  6.8,  1.4, 12.2, 10.3,  8.8,  0.7,  2.1,  1.3,  2.6,  5.1,  2.4,  3.0,  2.0,  1.9,  0.0, 
 3.0,  4.1,  5.2,  4.5,  4.7,  5.0,  2.1, 10.3, 13.0,  9.3,  2.4,  3.8,  2.5,  3.5,  7.4,  3.3,  3.3,  3.4,  3.1,  1.2, 
 1.6,  2.6,  4.5,  3.4,  4.2,  4.2,  0.2,  8.8,  9.3, 19.4,  0.0,  1.6,  0.9,  2.5,  4.4,  1.9,  1.7,  3.6,  1.7,  1.2, 
 4.9,  2.3,  1.2,  1.4,  2.0,  2.2,  4.5,  0.7,  2.4,  0.0,  9.9,  7.4,  7.9,  6.1,  5.6,  5.7,  5.2,  4.9,  5.7,  5.3, 
 4.9,  3.0,  2.2,  2.4,  3.4,  3.0,  4.3,  2.1,  3.8,  1.6,  7.4,  9.0,  6.1,  5.9,  6.4,  6.1,  5.7,  5.5,  6.0,  5.6, 
 5.2,  3.3,  2.4,  2.5,  2.2,  3.2,  4.7,  1.3,  2.5,  0.9,  7.9,  6.1,  8.8,  6.9,  5.6,  5.4,  5.1,  5.6,  6.4,  4.4, 
 5.0,  3.7,  3.6,  3.3,  2.8,  4.2,  5.0,  2.6,  3.5,  2.5,  6.1,  5.9,  6.9,  7.9,  6.4,  5.4,  5.2,  6.7,  6.7,  4.2, 
 4.4,  3.2,  3.3,  3.0,  3.9,  3.9,  4.1,  5.1,  7.4,  4.4,  5.6,  6.4,  5.6,  6.4, 11.2,  5.0,  4.9,  5.8,  5.8,  3.8, 
 6.3,  4.2,  3.1,  3.4,  5.3,  3.8,  5.6,  2.4,  3.3,  1.9,  5.7,  6.1,  5.4,  5.4,  5.0,  7.4,  6.7,  5.0,  5.3,  5.6, 
 5.8,  5.2,  3.9,  4.6,  4.7,  4.6,  5.3,  3.0,  3.3,  1.7,  5.2,  5.7,  5.1,  5.2,  4.9,  6.7,  7.7,  5.0,  5.3,  4.1, 
 4.6,  3.2,  3.0,  2.8,  3.0,  3.5,  4.3,  2.0,  3.4,  3.6,  4.9,  5.5,  5.6,  6.7,  5.8,  5.0,  5.0,  9.9,  7.9,  4.2, 
 4.8,  3.5,  3.1,  3.1,  2.4,  3.8,  4.6,  1.9,  3.1,  1.7,  5.7,  6.0,  6.4,  6.7,  5.8,  5.3,  5.3,  7.9,  8.4,  4.1, 
 5.7,  1.9,  0.8,  0.7,  3.2,  1.7,  3.6,  0.0,  1.2,  1.2,  5.3,  5.6,  4.4,  4.2,  3.8,  5.6,  4.1,  4.2,  4.1, 11.8};
///////////////////////////////////////////////////////////////////////////
void alignToSeq(int *align1, int1 *seq1, int nse1, int1 *seq2, int nse2, 
		int1 **seq1a, int1 **seq2a, int& nSeqa, 
		int& seq1MatchFrom, int& seq1MatchTo, int endsMode)
{
  int i;
  *seq1a=new int1 [nse1+nse2];
  *seq2a=new int1 [nse1+nse2];
  seq1MatchFrom=-1; seq1MatchTo=-1;
  int seq2MatchFrom, seq2MatchTo;
  int ise1, nSeqaPre=0;
  
  for(ise1=0; ise1<nse1; ise1++)
    {
      if(align1[ise1]>=0 &&  seq1MatchFrom==-1) 
	seq1MatchFrom=ise1;
      if(align1[ise1]>=0) seq1MatchTo=ise1;
    }
  seq2MatchFrom=align1[seq1MatchFrom]; seq2MatchTo=align1[seq1MatchTo];

  if(endsMode==1)
    {
      for(ise1=0; ise1<seq1MatchFrom; ise1++)
	{
	  *(*seq1a+nSeqaPre)=seq1[ise1]; 
	  *(*seq2a+nSeqaPre)=21; 
	  nSeqaPre++;
	}
      for(ise1=0; ise1<seq2MatchFrom; ise1++)
	{
	  *(*seq1a+nSeqaPre)=21; 
	  *(*seq2a+nSeqaPre)=seq2[ise1]; 
	  nSeqaPre++;
	}
    }

  nSeqa=nSeqaPre;

  for(ise1=seq1MatchFrom; ise1<=seq1MatchTo; ise1++)
    {
      if(ise1>0)
	if(align1[ise1-1]!=-1 && align1[ise1]-align1[ise1-1]>1)
	  {
	    for(i=0; i<align1[ise1]-align1[ise1-1]-1; i++) 
	      {
		*(*seq1a+nSeqa)=21; nSeqa++;
	      }
	  }
      *(*seq1a+nSeqa)=seq1[ise1]; nSeqa++;
    }
  nSeqa=nSeqaPre;
  for(ise1=seq1MatchFrom; ise1<=seq1MatchTo; ise1++)
    {
      if(ise1>0)
	if(align1[ise1-1]!=-1 &&align1[ise1]-align1[ise1-1]>1)
	  {
	    for(i=0; i<align1[ise1]-align1[ise1-1]-1; i++)
	      {
		*(*seq2a+nSeqa)=seq2[align1[ise1-1]+i+1]; nSeqa++;
	      }
	  }
      if(align1[ise1]>=0)
	{
	  *(*seq2a+nSeqa)=seq2[align1[ise1]]; nSeqa++;
	}
      if(align1[ise1]==-1)
	{
	  *(*seq2a+nSeqa)=21; nSeqa++;
	}
    }

  if(endsMode==1)
    {
      for(ise1=seq1MatchTo+1; ise1<nse1; ise1++)
	{
	  *(*seq1a+nSeqaPre)=seq1[ise1]; 
	  *(*seq2a+nSeqaPre)=21; 
	  nSeqa++;
	}
      for(ise1=seq2MatchTo+1; ise1<nse2; ise1++)
	{
	  *(*seq1a+nSeqaPre)=21; 
	  *(*seq2a+nSeqaPre)=seq2[ise1]; 
	  nSeqa++;
	}
    }
}
///////////////////////////////////////////////////////////////////////////
void alignToRef(int *align1, int nse1, int nse2, 
		int **seAlign1, int **seAlign2, int **alignSE1, int **alignSE2,
		int& nSeqa, int endsMode) {
  int i;
  *seAlign1=new int [nse1];
  *seAlign2=new int [nse2];
  *alignSE1=new int [nse1+nse2];
  *alignSE2=new int [nse1+nse2];
  int seq1MatchFrom=-1, seq1MatchTo=-1, seq2MatchFrom, seq2MatchTo;
  int ise1, nSeqaPre=0;
  
  for(ise1=0; ise1<nse1; ise1++)
    {
      if(align1[ise1]>=0 &&  seq1MatchFrom==-1) 
	seq1MatchFrom=ise1;
      if(align1[ise1]>=0) seq1MatchTo=ise1;
    }
  seq2MatchFrom=align1[seq1MatchFrom]; seq2MatchTo=align1[seq1MatchTo];

  if(endsMode==1)
    {
      for(ise1=0; ise1<seq1MatchFrom; ise1++)
	{
	  *(*alignSE1+nSeqaPre)=ise1;
	  *(*seAlign1+ise1)=nSeqaPre;
	  *(*alignSE2+nSeqaPre)=-1; 
	  nSeqaPre++;
	}
      for(ise1=0; ise1<seq2MatchFrom; ise1++)
	{
	  *(*alignSE1+nSeqaPre)=-1; 
	  *(*alignSE2+nSeqaPre)=ise1; 
	  *(*seAlign2+ise1)=nSeqaPre;
	  nSeqaPre++;
	}
    }

  nSeqa=nSeqaPre;

  for(ise1=seq1MatchFrom; ise1<=seq1MatchTo; ise1++)
    {
      if(ise1>0)
	if(align1[ise1-1]!=-1 && align1[ise1]-align1[ise1-1]>1)
	  {
	    for(i=0; i<align1[ise1]-align1[ise1-1]-1; i++) 
	      {
		*(*alignSE1+nSeqa)=-1; nSeqa++;
	      }
	  }
      *(*alignSE1+nSeqa)=ise1; 
      *(*seAlign1+ise1)=nSeqa;
      nSeqa++;
    }
  nSeqa=nSeqaPre;
  for(ise1=seq1MatchFrom; ise1<=seq1MatchTo; ise1++)
    {
      if(ise1>0)
	if(align1[ise1-1]!=-1 &&align1[ise1]-align1[ise1-1]>1)
	  {
	    for(i=0; i<align1[ise1]-align1[ise1-1]-1; i++)
	      {
		*(*alignSE2+nSeqa)=align1[ise1-1]+i+1;
		*(*seAlign2+align1[ise1-1]+i+1)=nSeqa;
		nSeqa++;
	      }
	  }
      if(align1[ise1]>=0)
	{
	  *(*alignSE2+nSeqa)=align1[ise1]; 
	  *(*seAlign2+align1[ise1])=nSeqa;
	  nSeqa++;
	}
      if(align1[ise1]==-1)
	{
	  *(*alignSE2+nSeqa)=-1; nSeqa++;
	}
    }

  if(endsMode==1)
    {
      for(ise1=seq1MatchTo+1; ise1<nse1; ise1++)
	{
	  *(*alignSE1+nSeqa)=ise1;
	  *(*seAlign1+ise1)=nSeqa;
	  *(*alignSE2+nSeqa)=-1; 
	  nSeqa++;
	}
      for(ise1=seq2MatchTo+1; ise1<nse2; ise1++)
	{
	  *(*alignSE1+nSeqa)=-1; 
	  *(*alignSE2+nSeqa)=ise1; 
	  *(*seAlign2+ise1)=nSeqa;
	  nSeqa++;
	}
    }
}
///////////////////////////////////////////////////////////////////////////
void align_se_to_se_align(int* align_se, int lcmp, int nse, int *se_align) {
  for(int ie=0; ie<nse; ie++) se_align[ie]=-1;
  for(int ia=0; ia<lcmp; ia++) {
    if(align_se[ia]!=-1) {
      if(align_se[ia]>nse) {
	printf("conversion error align_se -> se_align\n"); exit(0);
      }
      se_align[align_se[ia]]=ia;
    }
  }
}
///////////////////////////////////////////////////////////////////////////
void align_se_to_align_rep(int* align_se, int lcmp, int nse, int *align_rep) {
  for(int ie = 0; ie < nse; ie++) align_rep[ie] = -2;
  for(int ia = 0; ia < lcmp; ia++) {
    if(align_se[ia] != -1) {
      if(align_se[ia] > nse) {
	printf("align_se_to_align_rep - conversion error\n"); exit(0);
      }
      align_rep[align_se[ia]] = align_se[ia+lcmp];
    }
  }
}
///////////////////////////////////////////////////////////////////////////
void sortAlign(int item, int *align_se, int *ent, int lcmp, int ncmp) {
  if(ent[0]!=item) {
    for(int iit=0; iit<ncmp; iit++) {
      if(ent[iit]==item) {
	int tmp;
	tmp=ent[0];
	ent[0]=ent[iit];
	ent[iit]=tmp;
	
	int ibase=lcmp*iit;
	for(int ia=0; ia<lcmp; ia++) {
	  tmp=align_se[ia];
	  align_se[ia]=align_se[ibase+ia];
	  align_se[ibase+ia]=tmp;
	}
	return;
      }
    }
    printf("sortAlign - item %d not found\n", item); exit(0);
  }
}
////////////////////////////////////////////////////////////////////
/*
main()
{
  DB db("/misc/x1/moose/db96");
  Property name_enp("name.enp"), iseq_enp("iseq_flt.enp");
  int nse1, nse2, nSeqa, seqMatchFrom, seqMatchTo, i;
  int iEnp1=-1, iEnp2=-1; 
  for(i=0; i<name_enp.getNumObjects() && (iEnp1==-1 || iEnp2==-1); i++)
    {
      if(!strcmp(name_enp.items(i), "4HHB:A")) iEnp1=i;
      if(!strcmp(name_enp.items(i), "4HHB:B")) iEnp2=i;
    }
  if(iEnp1==-1 || iEnp2==-1) exit(0);
  
  int1 *seq1=iseq_enp.items1n(iEnp1, nse1), *seq2=iseq_enp.items1n(iEnp2, nse2), 
    *seq1a, *seq2a; 
  int *align1=new int [nse1], *align2=new int [nse2];

  strCmp(iEnp1, iEnp2, align1, align2, 20);

  alignToSeq(align1, seq1, nse1, seq2, nse2, &seq1a, &seq2a, nSeqa, 
	     seqMatchFrom, seqMatchTo);


  printf("(1) %d-%d from %d (2) %d-%d from %d\n\n", seqMatchFrom+1,
	 seqMatchTo+1, nse1, align1[seqMatchFrom]+1, align1[seqMatchTo]+1, 
	 nse2);
  for(i=0; i<nSeqa; i++) printf("%c", aa[seq1a[i]]);
  printf("\n");
  for(i=0; i<nSeqa; i++) printf("%c", aa[seq2a[i]]);
  printf("\n");
}
*/
///////////////////////////////////////////////////////////////////////////
/*
main() {
  DB db("/misc/x1/moose/db96");
  Property name_enp("name.enp"), ca_enp("c_a.enp");
  int pos;
  int iEnp=name_enp.find("4HHB:A", 0, pos);
  if(iEnp==-1) exit(0);
  flt4 *ca=ca_enp.itemsf(iEnp);
  int winSize=40, winStart1=0, winStart2=50;
  XYZ *mol1=new XYZ [winSize], *mol2=new XYZ [winSize];
  
  double d[20];
  for(int i=0; i<winSize; i++) {
    mol1[i].X=ca[winStart1+i*3];
    mol1[i].Y=ca[winStart1+i*3+1];
    mol1[i].Z=ca[winStart1+i*3+2];
    mol2[i].X=ca[winStart2+i*3];
    mol2[i].Y=ca[winStart2+i*3+1];
    mol2[i].Z=ca[winStart2+i*3+2];
  }

  for(int n=0; n<100; n++) {
    sup_str(mol1, mol2, winSize, d);
  }
} 
*/
///////////////////////////////////////////////////////////////////////////
int sup_str(XYZ *molA, XYZ *molB, int nAtom, double d[20] ) {
  kabsch(molA, molB, nAtom, d); 
  return(1);
}
///////////////////////////////////////////////////////////////////////////
double det3(double r[3][3]) {
  return(r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1])-
	 r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0])+
	 r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]));
}
///////////////////////////////////////////////////////////////////////////
double calc_rmsd(XYZ *molA, XYZ *molB, int nAtom, double d_[20] ) {
  double rmsd=0.0, dx, dy, dz, d;
  for(int l=0; l<nAtom; l++) {
    dx=molB[l].X; 
    dy=molB[l].Y;
    dz=molB[l].Z;
    d=molA[l].X-(dx*d_[0]+dy*d_[1]+dz*d_[2]+d_[9]);
    rmsd+=d*d;
    d=molA[l].Y-(dx*d_[3]+dy*d_[4]+dz*d_[5]+d_[10]);
    rmsd+=d*d;
    d=molA[l].Z-(dx*d_[6]+dy*d_[7]+dz*d_[8]+d_[11]);
    rmsd+=d*d;
  }
  return(sqrt(rmsd/nAtom));
}
///////////////////////////////////////////////////////////////////////////
void rot_mol(XYZ *molA, XYZ *molB, int nAtom, double d_[20] ) {
  double dx, dy, dz;
  for(int l=0; l<nAtom; l++) {
    if(molA[l].X<1e10) {
      dx=molA[l].X; 
      dy=molA[l].Y;
      dz=molA[l].Z;
      molB[l].X=dx*d_[0]+dy*d_[1]+dz*d_[2]+d_[9];
      molB[l].Y=dx*d_[3]+dy*d_[4]+dz*d_[5]+d_[10];
      molB[l].Z=dx*d_[6]+dy*d_[7]+dz*d_[8]+d_[11];
    }  
    else {
      molB[l]=molA[l];
    }
  }
}
///////////////////////////////////////////////////////////////////////////
XYZ *arrayToXYZ(float *array, int nXYZ){
  XYZ *xyz=new XYZ [nXYZ];
  for(int i=0; i<nXYZ; i++) {
    xyz[i].X=array[i*3];
    xyz[i].Y=array[i*3+1];
    xyz[i].Z=array[i*3+2];
  }
  return(xyz);
}
///////////////////////////////////////////////////////////////////////////
extern int aa_decode[25];
///////////////////////////////////////////////////////////////////////////
void seqCharToInt(char *seq1, char *seq2, int nseq)
{
  char seq_, aa_code_=20;
  for(int i=0; i<nseq; i++)
    {
      aa_code_=20;
      seq_=seq1[i];
      if(seq_!='*')
	{
	  seq_-=65;
	  if(seq_<0 || seq_>24)
	    {
	      printf("seqCharToInt: error at %d %d\n", i, seq_); 
	      exit(0);
	    }
	  aa_code_=aa_decode[seq_]-1;
	  if(aa_code_<0 || aa_code_>19)
	    {
	      printf("seqCharToInt: error at %d %d\n", i, seq_); 
	      exit(0);
	    }
	}
      seq2[i]=aa_code_;
    }
}
///////////////////////////////////////////////////////////////////////////
double dAngle(XYZ &A, XYZ &B, XYZ &C, XYZ &D)
{
  double return_angle;
  int ichk=-1;
  XYZ zero;
  XYZ X1, X2, X3, X4;
  X1=A; X2=B; X3=C; X4=D;

//  printf("%.3f %.3f %.3f %.3f %.3f %.3f\n", A.X, A.Y, A.Z,
//	 X1.X, X1.Y, X1.Z);
  zero=0.0;

  Vecdif(A,B,X3);
//  DisChk(X3,&ichk);
  if(ichk > 0)
    return(999.0);
  Vecdif(C,B,X4);
//  DisChk(X4,&ichk);
  if(ichk > 0)
    return(999.0);
  Cross(X3,X4,X1);
  Vecdif(B,C,X3);
  Vecdif(D,C,X4);
//  DisChk(X4,&ichk);
  if(ichk > 0)
    return(999.0);
  Cross(X3,X4,X2);


  double ang;
  double Q;

  Vecdif(X1,zero,X1);
  Vecdif(X2,zero,X2);
  
  Q=Amag(X1)*Amag(X2);
  if(Q < .0000001)
    Q=.0000001;
  ang=Dot(X1,X2)/Q;
  if(ang > 1.0)
    ang=1.0;
  if(ang < -1.0)
    ang=-1.0;
 return_angle=57.29577951*acos(ang);

  Cross(X1,X2,X4);
  X1.X=Dot(X3,X4);
  if(X1.X > 0.0)
    return_angle=-return_angle;
  return(return_angle);
}
//////////////////////////////////////////////////////////////////
void Cross(XYZ &a, XYZ &b, XYZ &c)
{
  c.X=a.Y*b.Z-b.Y*a.Z;
  c.Y=-a.X*b.Z+b.X*a.Z;
  c.Z=a.X*b.Y-b.X*a.Y;
}
//////////////////////////////////////////////////////////////////
double Dot(XYZ &a, XYZ &b)
{
  double val=0.0;
  
  val+=a.X*b.X;
  val+=a.Y*b.Y;
  val+=a.Z*b.Z;
  
  return(val);
}
//////////////////////////////////////////////////////////////////
void Vecdif(XYZ &a, XYZ &b, XYZ &c)
{
  c.X=a.X-b.X;
  c.Y=a.Y-b.Y;
  c.Z=a.Z-b.Z;
}
//////////////////////////////////////////////////////////////////
void DisChk(XYZ &x, int *ichk)
{
  double DSSQ;

  *ichk=0;
  DSSQ=Dot(x,x);
  if((DSSQ < .50) || (DSSQ > 9.00))
    *ichk=1;
  return;
}
//////////////////////////////////////////////////////////////////
double Amag(XYZ &x)
{
  return(sqrt(Dot(x,x)));
}
//////////////////////////////////////////////////////////////////
int lali_of_trace(int *trace, int nTrace) {
  int lali = 0;
  for(int it = 0; it < nTrace; it++)
    lali += trace[it*3+2];
  return(lali);
}
//////////////////////////////////////////////////////////////////
void align_se_to_trace(int* align_se, int lcmp, int4 **trace, int& nTrace) {
  align_se_to_trace(align_se, align_se+lcmp, lcmp, trace, nTrace);
}
//////////////////////////////////////////////////////////////////
void align_se_to_trace(int* align_se1, int* align_se2, int lcmp, int4 **trace, 
		       int& nTrace) {
  if(lcmp==0) {
    nTrace=0; return;
  }
  int k=0, n=0;
  int4 *trace_=new int4 [lcmp*3];
  for(int ia=0; ia<lcmp; ia++) 
    if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
      if(k==0) {
	trace_[n*3]=align_se1[ia];
	trace_[n*3+1]=align_se2[ia];
      }
      k++;
    }
    else {
      if(k>0) {
	trace_[n*3+2]=k;
	n++;
	k=0;
      }
    }
  if(k>0) {
    trace_[n*3+2]=k;
    n++;
  }
  (*trace)=new int4 [n*3];
  for(int it=0; it<n*3; it++) *(*trace+it)=trace_[it];
  nTrace=n;
  delete [] trace_;
}
//////////////////////////////////////////////////////////////////
void trace_to_align_se(int *trace, int nTrace, int** align_se, int& lcmp) {
  if(nTrace==0) {
    lcmp=0; return;
  }

  int lcmp_=0;
  for(int it=0; it<nTrace; it++) {
    lcmp_+=trace[it*3+2];
    if(it<nTrace-1)
      lcmp_+=trace[(it+1)*3]-(trace[it*3]+trace[it*3+2])+
	trace[(it+1)*3+1]-(trace[it*3+1]+trace[it*3+2]);
  }

  int *align_se_=new int [lcmp_*2];
  int l=0;
  int is1=trace[0], is2=trace[1];
  
  for(int it=0; it<nTrace; it++) {
    for(int k=0; k<trace[it*3+2]; k++) {
      align_se_[l]=is1;
      align_se_[l+lcmp_]=is2;
      is1++; is2++; l++;
      if(l>lcmp_) goto err;
    }
    if(it<nTrace-1) {
      for(int is=is1; is<trace[(it+1)*3]; is++) {
	align_se_[l]=is1; 
	align_se_[l+lcmp_]=-1; 
	is1++; l++;
	if(l>lcmp_) goto err;
      }
      for(int is=is2; is<trace[(it+1)*3+1]; is++) {
	align_se_[l]=-1; 
	align_se_[l+lcmp_]=is2; 
	is2++; l++;
	if(l>lcmp_) goto err;
      }
    }
  }
  lcmp=l;
  (*align_se)=new int [lcmp*2];
  for(int ia=0; ia<lcmp; ia++) {
    *(*align_se+ia)=align_se_[ia];
    *(*align_se+ia+lcmp)=align_se_[ia+lcmp_];
  }
  delete [] align_se_;
  return;
 err:
  printf("trace_to_align_se error\n"); exit(0);
}
//////////////////////////////////////////////////////////////////


/* This code does a basic structure based alignment */
/* It is based on the method of Gerstein & Levitt but does not follow */
/* their method exactly */
/* This program is probably only useful if the sequence ID */
/* between the structures is quite high */
/* This code is available AS IS and at your own risk */
/* Mansoor Saqi, Mike Hartshorn */
/* Some code from Roger Sayle */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define debug 0 

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif

#define     EQN_EPS     1.e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)

#define cross(a,b,c) (a[0] = (b[1] * c[2]) - (b[2] * c[1]), \
		      a[1] = (b[2] * c[0]) - (b[0] * c[2]), \
                      a[2] = (b[0] * c[1]) - (b[1] * c[0]))


double kabsch(XYZ *strBuf1, XYZ *strBuf2, int lali_str, double d[20]) {

/*
 * Fit pairs of points using the method described in:
 *
 * `A discussion of the solution for the best rotation to relate
 *  two sets of vectors', W.Kabsch, Acta Cryst. (1978), A34, 827-828.
 *
 * The method generates the 4x4 matrix rot which will map the
 * coordinates y on to the coordinates x.
 *
 * The matrix is arranged such that rot[3][0], rot[3][1], rot[3][2]
 * are the translational components.
 *
 * The function returns the root mean square distance between
 * x and y.  The coordinates in x and y are not modified.
 */


  float rot[4][4];
  float xc[3], yc[3];
  float r[3][3], rtrans[3][3], rr[3][3];
  float mu[3], a[3][3], atrans[3][3], b[3][3], u[3][3];
  float yy[3];
  float rmsd;
  int i, j, k;
  
  int n;

  for(i = 0; i < 3; i++) {
    xc[i] = 0.0;
    yc[i] = 0.0;
  } 

  n = lali_str;
  

  /* Find center of each set of coordinates. */

  for (i = 0; i < n; i++)
    {
      xc[0] += strBuf1[i].X;
      yc[0] += strBuf2[i].X;
      xc[1] += strBuf1[i].Y;
      yc[1] += strBuf2[i].Y;
      xc[2] += strBuf1[i].Z;
      yc[2] += strBuf2[i].Z;
      
      //for (j = 0; j < 3; j++)
      //{
      //  xc[j] += MX[i][j];
      //  yc[j] += MY[i][j];
      //}
    }

  for (j = 0; j < 3; j++)
    {
      xc[j] /= (float) n;
      yc[j] /= (float) n;
    }

  /*
   * Initialise and then fill the r matrix.
   * Note that centre is subtracted at this stage.
   */

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      r[i][j] = 0.0;

  for (k = 0; k < n; k++)
    {

      r[0][0] += (strBuf2[k].X - yc[0]) * (strBuf1[k].X - xc[0]);
      r[0][1] += (strBuf2[k].X - yc[0]) * (strBuf1[k].Y - xc[1]);
      r[0][2] += (strBuf2[k].X - yc[0]) * (strBuf1[k].Z - xc[2]);
      r[1][0] += (strBuf2[k].Y - yc[1]) * (strBuf1[k].X - xc[0]);
      r[1][1] += (strBuf2[k].Y - yc[1]) * (strBuf1[k].Y - xc[1]);
      r[1][2] += (strBuf2[k].Y - yc[1]) * (strBuf1[k].Z - xc[2]);
      r[2][0] += (strBuf2[k].Z - yc[2]) * (strBuf1[k].X - xc[0]);
      r[2][1] += (strBuf2[k].Z - yc[2]) * (strBuf1[k].Y - xc[1]);
      r[2][2] += (strBuf2[k].Z - yc[2]) * (strBuf1[k].Z - xc[2]);

      //for (i = 0; i < 3; i++)
      //{
      //  for (j = 0; j < 3; j++)
      //    {
      //      r[i][j] += (MY[k][i] - yc[i]) * (MX[k][j] - xc[j]);
      //    }
      //}

    }

  /* Generate the transpose of r and form rtrans x r */

  matrix_transpose (rtrans, r);

  matrix_multiply (rr, rtrans, r);

  /*
   * Get the eigenvalues and vectors.
   * Reform a[2] as cross product of a[0] and a[1] to ensure
   * right handed system.
   */

  eigen_values (rr, mu, a);

  cross (a[2], a[0], a[1]);

  /* Transform first two eigenvectors and normalise them. */

  for (i = 0; i < 2; i++)
    {
      transformpoint (b[i], r, a[i]);

      normalise (b[i]);
    }

  /* Make right handed set. */

  cross (b[2], b[0], b[1]);

  /* Form the rotation matrix. */

  matrix_transpose (atrans, a);

  matrix_multiply (u, b, atrans);

  /* Make rot the identity matrix. */

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      rot[i][j] = (i == j) ? 1. : 0.;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      rot[i][j] = u[i][j];

  /* Transform offset of y coordinates by the rotation. */

  transformpoint (yy, u, yc);

  /* Build translational component of rot from offsets. */

  for (i = 0; i < 3; i++)
    rot[3][i] = -yy[i] + xc[i];

  /* Figure out the rms deviation of the fitted coordinates. */

  d[0] = rot[0][0]; d[1] = rot[1][0]; d[2] = rot[2][0]; 
  d[3] = rot[0][1]; d[4] = rot[1][1]; d[5] = rot[2][1]; 
  d[6] = rot[0][2]; d[7] = rot[1][2]; d[8] = rot[2][2]; 
  d[9] = rot[3][0]; d[10] = rot[3][1]; d[11] = rot[3][2]; 


  rmsd = 0.0;

  for (i = 0; i < n; i++)
    {
      float xt[3], yt[3];

      /* Translate the pairs of coordinates to origin. */

      xt[0] = strBuf1[i].X - xc[0];
      xt[1] = strBuf1[i].Y - xc[1];
      xt[2] = strBuf1[i].Z - xc[2];

      yt[0] = strBuf2[i].X - yc[0];
      yt[1] = strBuf2[i].Y - yc[1];
      yt[2] = strBuf2[i].Z - yc[2];

      //for (j = 0; j < 3; j++)
      //{
      //  xt[j] = MX[i][j] - xc[j];
      //  yt[j] = MY[i][j] - yc[j];
      //}

      /* Transform the y coordinate by u. */

      transformpoint (yy, u, yt);

      /* Accumulate the difference. */

      rmsd +=
	(yy[0] - xt[0]) * (yy[0] - xt[0]) +
	(yy[1] - xt[1]) * (yy[1] - xt[1]) +
	(yy[2] - xt[2]) * (yy[2] - xt[2]);
    }

  rmsd /= (float) n;

  rmsd = sqrt (rmsd);

  return rmsd;
}


int cubic_roots (double c[4], double s[3])
{
  int i, num;
  double sub;
  double A, B, C;
  double A2, p, q;
  double p3, D;

  /* normal form: x^3 + Ax^2 + Bx + C = 0 */

  A = c[2] / c[3];
  B = c[1] / c[3];
  C = c[0] / c[3];

  /*  substitute x = y - A/3 to eliminate quadric term:
     x^3 +px + q = 0 */

  A2 = A * A;
  p = 1.0 / 3 * (-1.0 / 3 * A2 + B);
  q = 1.0 / 2 * (2.0 / 27 * A * A2 - 1.0 / 3 * A * B + C);

  /* use Cardano's formula */

  p3 = p * p * p;
  D = q * q + p3;

  if (IsZero (D))
    {				/* one triple solution */
      if (IsZero (q))
	{
	  s[0] = 0;
	  num = 1;
	}
      else
	{			/* one single and one double solution */
	  double u = cbrt (-q);
	  s[0] = 2 * u;
	  s[1] = -u;
	  num = 2;
	}
    }
  else if (D < 0)
    {				/* Casus irreducibilis: three real solutions */

      double phi = 1.0 / 3 * acos (-q / sqrt (-p3));
      double t = 2 * sqrt (-p);

      s[0] = t * cos (phi);
      s[1] = -t * cos (phi + M_PI / 3);
      s[2] = -t * cos (phi - M_PI / 3);
      num = 3;
    }
  else
    {				/* one real solution */

      double sqrt_D = sqrt (D);
      double u = cbrt (sqrt_D - q);
      double v = -cbrt (sqrt_D + q);

      s[0] = u + v;
      num = 1;
    }

  /* resubstitute */

  sub = 1.0 / 3 * A;

  for (i = 0; i < num; ++i)
    s[i] -= sub;

  return num;
}


void eigen_values (float m[3][3], float values[3], float vectors[3][3])
/*
 * Find the eigen values and vectors for the matrix m.
 */
{
  double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2];
  double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2];
  double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2];
  double c[4];
  double x, y, z, norm;
  double roots[3], l;
  int nroots, iroot;
  int i;

  /*
   * Expanding the characteristic equation of a 3x3 matrix...
   * Maple tells us the following is the cubic in l
   */

  c[0] = a1 * (b2 * c3 - b3 * c2) + c1 * (a2 * b3 - a3 * b2) - b1 * (a2 * c3 - a3 * c2);
  c[1] = (a1 * (-b2 - c3) - b2 * c3 + b3 * c2 + c1 * a3 + b1 * a2);
  c[2] = (a1 + b2 + c3);
  c[3] = -1.;

  nroots = cubic_roots (c, roots);

  /* Degenerate roots are not returned individually. */

  for (i = 0; i < nroots; i++)
    {
      iroot = i > nroots ? nroots : i;

      l = roots[iroot];

      values[i] = l;

      /*
       * Find the eigen vectors by solving pairs of the
       * three simultaneous equations.  From `Mathematical Methods
       * in Science and Engineering', Heiding, pg.19
       *
       * Sometimes we get x = y = z = 0.0, so try the other two
       * pairs of equations and hope that one of them gives a solution.
       */

      x = b1 * c2 - (b2 - l) * c1;
      y = -((a1 - l) * c2 - a2 * c1);
      z = ((a1 - l) * (b2 - l) - a2 * b1);

      if (IsZero (x) && IsZero (y) && IsZero (z))
	{
	  x = b1 * (c3 - l) - b3 * c1;
	  y = -((a1 - l) * (c3 - l) - a3 * c1);
	  z = ((a1 - l) * b3 - a3 * b1);

	  if (IsZero (x) && IsZero (y) && IsZero (z))
	    {
	      x = (b2 - l) * (c3 - l) - b3 * c2;
	      y = -(a2 * (c3 - l) - a3 * c2);
	      z = (a2 * b3 - a3 * (b2 - l));
	      if (IsZero (x) && IsZero (y) && IsZero (z))
		{
		  printf ("eigen: no solution for eigen vector %d\n", i);
		  exit(0);
		}
	    }
	}

      norm = sqrt (x * x + y * y + z * z);

      if (!IsZero (norm))
	{
	  vectors[i][0] = x / norm;
	  vectors[i][1] = y / norm;
	  vectors[i][2] = z / norm;
	}
    }
}

void transformpoint (float p_new[3], float m[3][3], float p[3])
{
  int i, column;

  for (i = 0; i < 3; i++)
    {
      p_new[i] = 0.0;
      for (column = 0; column < 3; column++)
	p_new[i] += m[column][i] * p[column];
    }
}

void matrix_multiply (float m[3][3], float a[3][3], float b[3][3])
{
  int row, column, i;

  for (row = 0; row < 3; row++)
    for (column = 0; column < 3; column++)
      {
	m[row][column] = 0.0;
	for (i = 0; i < 3; i++)
	  m[row][column] += b[row][i] * a[i][column];
      }
}

void matrix_transpose (float result[3][3], float a[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      result[i][j] = a[j][i];
}

void normalise (float a[3])
{
  float norm = sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  if (norm > 1.e-9)
    {
      a[0] /= norm;
      a[1] /= norm;
      a[2] /= norm;
    }
}


////////////////////////////////////////////////////////////////////

