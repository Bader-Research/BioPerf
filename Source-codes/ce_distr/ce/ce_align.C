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
///////////////////////////////////////////////////////////////////////////
double ce_1(char *name1, char *name2,
		XYZ *ca1, XYZ *ca2, int nse1, int nse2, int *align_se1, 
		int *align_se2, int& lcmp, int winSize_,
		double rmsdThr, double rmsdThrJoin, double d_[20], 
	    int isPrint) {

  int distAll=0;

  if(isPrint == 2) printf("CE Algorithm, version 1.00, 1998.");
  else if (isPrint == 1) for(int ip = 0; ip < 85; ip++) printf("-");

  if(isPrint) printf("\n\nChain 1: %s (Size=%d) \nChain 2: %s (Size=%d)\n\n", name1, nse1, name2, nse2);

  tms time1, time2;
  times(&time1);
  long nTraces, tracesLimit=(long)5e7;
  
  double rmsdThrSq=rmsdThr*rmsdThr;
  double bestTraceScore, oldBestTraceScore=10000.0, bestTraceZScore;

  double d, dd, d1, d2, dx, dy, dz, score, score0, score1, score2, 
    rmsd, rmsdNew, z; 

  int ise1, ise2, jse1, jse2, kse1, kse2, mse1, mse2, l1, l2, i, j, nd;

  int traceMaxSize=nse1<nse2?nse1:nse2;

  int bestTracesMax=30;
  int **bestTraces1, **bestTraces2, *bestTracesN=new int [bestTracesMax];
  double *bestTracesScores=new double [bestTracesMax];
  NewArray(&bestTraces1, bestTracesMax, traceMaxSize);
  NewArray(&bestTraces2, bestTracesMax, traceMaxSize);
  int nBestTraces=0, newBestTrace=0;

  for(int it=0; it<bestTracesMax; it++) {
    bestTracesN[it]=0;
    bestTracesScores[it]=100;
  }

  int *f1=new int [nse1], *f2=new int [nse2];
  double **mat;
  double **dist1, **dist2;
  NewArray(&mat, nse1, nse2);
  NewArray(&dist1, nse1, nse1);
  NewArray(&dist2, nse2, nse2);
  int gapMax=30;
  int iterDepth=gapMax*2+1;
  int *trace1=new int [traceMaxSize], *trace2=new int [traceMaxSize];
  int *traceIterLevel=new int [traceMaxSize],
    *traceIndex=new int [traceMaxSize];
  double **traceScore;
  int **traceUsage;
  NewArray(&traceScore, traceMaxSize, iterDepth);
  NewArray(&traceUsage, traceMaxSize, iterDepth);
  
  int nTrace, nGaps, igap, idir, jgap, jdir, nBestTrace=0, nBestTrace0;
  XYZ *strBuf1=new XYZ [traceMaxSize], *strBuf2=new XYZ [traceMaxSize];
  int *bestTrace1=new int [traceMaxSize], *bestTrace2=new int [traceMaxSize];

  int winSize=winSize_;
  
  for(int is=0; is<nse1-winSize+1; is++) {
    f1[is]=0;
    for(int ir=0; ir<winSize; ir++) if(ca1[is+ir].X>=1e10) goto cont1;
    f1[is]=1;
  cont1: ;
  }

  for(int is=0; is<nse2-winSize+1; is++) {
    f2[is]=0;
    for(int ir=0; ir<winSize; ir++) if(ca2[is+ir].X>=1e10) goto cont2;
    f2[is]=1;
  cont2: ;
  }

  for(ise1=0; ise1<nse1; ise1++) 
    for(ise2=0; ise2<nse1; ise2++) 
      *(dist1[ise1]+ise2)=
	ca1[ise1].X<1e10&&ca1[ise2].X<1e10?ca1[ise1].dist(ca1[ise2]):2e10;
  
  for(ise1=0; ise1<nse2; ise1++) 
    for(ise2=0; ise2<nse2; ise2++) 
      *(dist2[ise1]+ise2)=
	ca2[ise1].X<1e10&&ca2[ise2].X<1e10?ca2[ise1].dist(ca2[ise2]):2e10;
  
  int isGood;
  
  int ise11, ise12, ise21, ise22;
  
  int winSizeComb1=(winSize-1)*(winSize-2)/2;
  int winSizeComb2=distAll?winSize*winSize:winSize;

  int *a=new int [traceMaxSize];

  if(distAll)
    for(i=0; i<traceMaxSize; i++)
      a[i]=(i+1)*i*winSize*winSize/2+(i+1)*winSizeComb1;
  else
    for(i=0; i<traceMaxSize; i++)
      a[i]=(i+1)*i*winSize/2+(i+1)*winSizeComb1;
  
  for(ise1=0; ise1<nse1; ise1++) {
    for(ise2=0; ise2<nse2; ise2++) {
      *(mat[ise1]+ise2)=-1.0;
      if(ise1>nse1-winSize || ise2>nse2-winSize) continue;
      
      isGood=1;
      
      for(l1=0; l1<winSize; l1++) 
	if(ca1[ise1+l1].X>1e10 || ca2[ise2+l1].X>1e10) {
	    isGood=0; break;
	}
      
      if(!isGood) continue;
      
      d=0.0; 
      for(int is1=0; is1<winSize-2; is1++)
	for(int is2=is1+2; is2<winSize; is2++) {
	  d+=fabs(*(dist1[ise1+is1]+ise1+is2)-
		      *(dist2[ise2+is1]+ise2+is2));
	}
      *(mat[ise1]+ise2)=d/winSizeComb1; 
      
    }
  }
    
  /*
    printf("  ");  
    for(i=0; i<nse2; i++) printf("%c", i%10==0?(i%100)/10+48:32);
    printf("\n");  
    for(i=0; i<nse1; i++) {
    printf("%c ", i%10==0?(i%100)/10+48:32);
    for(j=0; j<nse2; j++) 
    printf("%c", *(mat[i]+j)<rmsdThr?'+':'X');
    //printf("%c", ((int)*(mat[i]+j)/40)>9?'*':((int)*(mat[i]+j)/40)+48);
    printf("\n");
    }
    printf("\n");
    */
  
  // tracing trough fragment matrix
  
  int nIter=1;
  int iterLevel=1;

  nBestTrace=0;    
  bestTraceScore=100.0;
  bestTraceZScore=-1.0;
  nTraces=0;

  double z0, zThr=-0.01;

  for(int iter=0; iter<nIter; iter++) {

    if(iter>2) {
      if(oldBestTraceScore<=bestTraceScore) break;
    }
    oldBestTraceScore=bestTraceScore;

    if(iter==1) {
      z0=zStrAlign(winSize, nBestTrace, bestTraceScore, 
		   bestTrace1[nBestTrace]+winSize-bestTrace1[0]+
		   bestTrace2[nBestTrace]+winSize-bestTrace2[0]-
		   nBestTrace*2*winSize);
      if(z0<zThr) break;
      nBestTrace0=nBestTrace;
      nBestTrace=0; 
      bestTraceScore=100.0;
      bestTraceZScore=-1.0;
      nTraces=0;
    }

    int iseStep;
    if(iter==0) {
      ise11=0; ise12=nse1;
      ise21=0; ise22=nse2;
      iseStep=1;
    }
    else {
      if(iter==1) {
	ise11=bestTrace1[0]; ise12=bestTrace1[0]+1;
	ise21=bestTrace2[0]; ise22=bestTrace2[0]+1;
      }
      else {
	ise11=bestTrace1[0]-1; ise12=bestTrace1[0]+2;
	ise21=bestTrace2[0]-1; ise22=bestTrace2[0]+2;
      }
      if(ise11<0) ise11=0;
      if(ise12>nse1) ise12=nse1;
      if(ise21<0) ise21=0;
      if(ise22>nse2) ise22=nse2;
      iseStep=1;
      
    }


    for(int ise1_=ise11; ise1_<ise12; ise1_++) {
      for(int ise2_=ise21; ise2_<ise22; ise2_++) {

	ise1=ise1_;
	ise2=ise2_;
	if(iter>1 && ise1==ise11+1 && ise2==ise21+1) continue;

	//if(ise2==ise21) printf("(%d, %d)\n",ise1, nTraces);

	
	if(iter==0 && (ise1>nse1-winSize*(nBestTrace-1) || 
	ise2>nse2-winSize*(nBestTrace-1))) continue;

	if(*(mat[ise1]+ise2)<0.0) continue;
	if(*(mat[ise1]+ise2)>rmsdThr) continue;
	
	nTrace=0;
	trace1[nTrace]=ise1; 
	trace2[nTrace]=ise2;
	traceIndex[nTrace]=0;
	traceIterLevel[nTrace]=0;
	
	score0=*(mat[ise1]+ise2);
	
	/*
	if(nBestTrace==0 || 
	   (nBestTrace==1 && score0<bestTraceScore)) {
	  bestTraceScore=score0;
	  bestTrace1[0]=ise1;
	  bestTrace2[0]=ise2;
	  nBestTrace=1;
	}
	*/
	
	/*
	z=zStrAlign(winSize, 1, score0, 0);
	if(z>bestTraceZScore) {
	  bestTraceZScore=z;
	  bestTraceScore=score0;
	  bestTrace1[0]=ise1;
	  bestTrace2[0]=ise2;
	  nBestTrace=1;
	}
	*/
	
	
	nTrace++;
	int isTraceUp=1;
	int traceIndex_=0;
	
	while(nTrace>0) {
	  
	  iterLevel=1;
	  
	  //if(iter==0) {
	  //if(nBestTrace<4) iterLevel=3;
	  //if(nBestTrace-nTrace<4 && nTrace>4) iterLevel=3;
	  //}

	  
	  kse1=trace1[nTrace-1]+winSize;
	  kse2=trace2[nTrace-1]+winSize;
     
	  while(1) {
	    if(kse1>nse1-winSize-1) break;
	    if(kse2>nse2-winSize-1) break;
	    if(*(mat[kse1]+kse2)>=0.0) break;
	    if(f1[kse1]==0) kse1++;
	    if(f2[kse2]==0) kse2++;
	  }
	  

	  traceIndex_=-1; 

	  if(isTraceUp) {
	    int nBestExtTrace=nTrace;
	    double bestExtScore=100.0;

	    for(int it=0; it<iterDepth; it++) {
	      jgap=(it+1)/2;
	      jdir=(it+1)%2;
	      if(jdir==0) {
		mse1=kse1+jgap;
		mse2=kse2;
	      }
	      else {
		mse1=kse1;
		mse2=kse2+jgap;
	      }
	      
	      if(mse1>nse1-winSize-1) continue;
	      if(mse2>nse2-winSize-1) continue;
	      
	      if(*(mat[mse1]+mse2)<0.0) continue;
	      if(*(mat[mse1]+mse2)>rmsdThr) continue;
	      
	      nTraces++;
	      if(nTraces>tracesLimit) {
		if(isPrint) printf("exited ");
		goto exit_tracing;
	      }
	      
	      score=0.0;
	      if(!distAll) {
		// (winSize) "best" dist
		for(int itrace=0; itrace<nTrace; itrace++) {
		  score+=fabs(*(dist1[trace1[itrace]]+mse1)-
			      *(dist2[trace2[itrace]]+mse2));
		  score+=fabs(*(dist1[trace1[itrace]+winSize-1]+
				mse1+winSize-1)-
			      *(dist2[trace2[itrace]+winSize-1]+
				mse2+winSize-1));
		  
		  for(int id=1; id<winSize-1; id++) 
		    score+=fabs(*(dist1[trace1[itrace]+id]+mse1+winSize-1-id)-
				*(dist2[trace2[itrace]+id]+mse2+winSize-1-id));
		  
		}
		score1=score/(nTrace*winSize);
	      }
	      else {
		// all dist
		for(int itrace=0; itrace<nTrace; itrace++) {
		  for(int is1=0; is1<winSize; is1++)
		    for(int is2=0; is2<winSize; is2++)
		      score+=fabs(*(dist1[trace1[itrace]+is1]+mse1+is2)-
				  *(dist2[trace2[itrace]+is1]+mse2+is2));
		}
		score1=score/(nTrace*winSize*winSize);
	      }
	      
	      
	      if(score1>rmsdThrJoin) continue;
	      
	      score2=score1;
	      
	      /*
		printf("[%2d %3.1f %3.1f %d %d %d] ", nTrace, score1, score2, 
		extIter[0], extIter[1], traceIndex_);
		for(int k=0; k<nTrace; k++)
		printf("(%d,%d) ", trace1[k], trace2[k]);
		printf("(%d,%d)\n", mse1, mse2);
		*/
	      
	      if(score2>rmsdThrJoin) continue;
	      
	      if(nTrace>nBestExtTrace || (nTrace==nBestExtTrace &&
					   score2<bestExtScore)) {
		bestExtScore=score2;
		nBestExtTrace=nTrace;
		traceIndex_=it;
		*(traceScore[nTrace-1]+traceIndex_)=score1;
	      }

	    }
	  }

	  double rmsdScore, traceTotalScore, traceScoreMax;

	  if(traceIndex_!=-1) {
	    jgap=(traceIndex_+1)/2;
	    jdir=(traceIndex_+1)%2;
	    if(jdir==0) {
	      jse1=kse1+jgap;
	      jse2=kse2;
	    }
	    else {
	      jse1=kse1;
	      jse2=kse2+jgap;
	    }

	    if(iter==0){
	      score1=(*(traceScore[nTrace-1]+traceIndex_)*winSizeComb2*nTrace+
		      *(mat[jse1]+jse2)*winSizeComb1)/(winSizeComb2*nTrace+
						       winSizeComb1);
	      
	      score2=
		((nTrace>1?*(traceScore[nTrace-2]+traceIndex[nTrace-1]):score0)
		 *a[nTrace-1]+score1*(a[nTrace]-a[nTrace-1]))/a[nTrace];
	      
	      if(score2>rmsdThrJoin) traceIndex_=-1;
	      else traceTotalScore=*(traceScore[nTrace-1]+traceIndex_)=score2;
	    }
	    else {
	      if(traceScoreMax>rmsdThrJoin && nBestTrace>=nBestTrace0) 
		traceIndex_=-1;
	      traceTotalScore=traceScoreMax;
	    }
	  }


	  if(traceIndex_==-1) {
	    //if(iterLevel==1) break;
	    nTrace--;
	    isTraceUp=0;
	    continue;
	  }
	  else {
	    traceIterLevel[nTrace-1]++;
	    trace1[nTrace]=jse1;
	    trace2[nTrace]=jse2;
	    traceIndex[nTrace]=traceIndex_;
	    traceIterLevel[nTrace]=0;
	    nTrace++;
	    isTraceUp=1;
	    
	      if(nTrace>nBestTrace || 
		 (nTrace==nBestTrace  && 
		  bestTraceScore>traceTotalScore)) {
		
		for(int itrace=0; itrace<nTrace; itrace++) {
		  bestTrace1[itrace]=trace1[itrace];
		  bestTrace2[itrace]=trace2[itrace];
		}
		bestTraceScore=traceTotalScore;
		nBestTrace=nTrace;
              }

	    if(iter==0)
	      if(nTrace>bestTracesN[newBestTrace] ||
		 (nTrace==bestTracesN[newBestTrace] && 
		  bestTracesScores[newBestTrace]>traceTotalScore)) {
		for(int itrace=0; itrace<nTrace; itrace++) {
		  *(bestTraces1[newBestTrace]+itrace)=trace1[itrace];
		  *(bestTraces2[newBestTrace]+itrace)=trace2[itrace];
		  bestTracesN[newBestTrace]=nTrace;
		  bestTracesScores[newBestTrace]=traceTotalScore;
		}
		
		if(nTrace>nBestTrace) nBestTrace=nTrace;
		
		if(nBestTraces<bestTracesMax) {
		  nBestTraces++;
		  newBestTrace++;
		}
		
		if(nBestTraces==bestTracesMax) {
		  newBestTrace=0;
		  double scoreTmp=bestTracesScores[0];
		  int nTraceTmp=bestTracesN[0];
		  for(int ir=1; ir<nBestTraces; ir++) {
		    if(bestTracesN[ir]<nTraceTmp || 
		       (bestTracesN[ir]==nTraceTmp && 
			scoreTmp<bestTracesScores[ir])) {
		      nTraceTmp=bestTracesN[ir];
		      scoreTmp=bestTracesScores[ir];
		      newBestTrace=ir;
		    }
		  }
		}
	      }
	    
	    
	    /*
	    z=zStrAlign(winSize, nTrace, traceTotalScore, 
			trace1[nTrace-1]-trace1[0]+trace2[nTrace-1]-trace2[0]-
			2*(nTrace-1)*winSize);
	    if(z>bestTraceZScore) {
	      for(int itrace=0; itrace<nTrace; itrace++) {
		bestTrace1[itrace]=trace1[itrace];
		bestTrace2[itrace]=trace2[itrace];
	      }
	      bestTraceZScore=z;
	      bestTraceScore=*(traceScore[nTrace-2]+traceIndex_);
	      nBestTrace=nTrace;
	    }
	    */
	  }
	}
      }
    }
  }
  exit_tracing: ;
  lcmp=0;
  z=0.0;
  if(nBestTrace>0) {
    int is=0;
    for(int jt=0; jt<nBestTrace; jt++) {
      for(i=0; i<winSize; i++) {
	strBuf1[is+i]=ca1[bestTrace1[jt]+i];
	strBuf2[is+i]=ca2[bestTrace2[jt]+i];
      }
      is+=winSize;
    }
    
    sup_str(strBuf1, strBuf2, nBestTrace*winSize, d_);
    rmsdNew=calc_rmsd(strBuf1, strBuf2, nBestTrace*winSize, d_);
    /*
      printf("rmsd=%.2f ", rmsdNew);
      for(int k=0; k<nBestTrace; k++)
      printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
      8);
      printf("\n");
    */
    
    rmsd=100.0;
    int iBestTrace=0;
    for(int ir=0; ir<nBestTraces; ir++) {
      if(bestTracesN[ir]!=nBestTrace) continue;
      int is=0;
      for(int jt=0; jt<bestTracesN[ir]; jt++) {
	for(i=0; i<winSize; i++) {
	  strBuf1[is+i]=ca1[*(bestTraces1[ir]+jt)+i];
	  strBuf2[is+i]=ca2[*(bestTraces2[ir]+jt)+i];
	}
	is+=winSize;
      }
      sup_str(strBuf1, strBuf2, bestTracesN[ir]*winSize, d_);
      rmsdNew=calc_rmsd(strBuf1, strBuf2, bestTracesN[ir]*winSize, d_);
      //printf("%d %d %d %.2f\n", ir, bestTracesN[ir], nBestTrace, rmsdNew);
      if(rmsd>rmsdNew) {
	iBestTrace=ir;
	rmsd=rmsdNew;
      }
    }
    for(int it=0; it<bestTracesN[iBestTrace]; it++) {
      bestTrace1[it]=*(bestTraces1[iBestTrace]+it);
      bestTrace2[it]=*(bestTraces2[iBestTrace]+it);
    }
    nBestTrace=bestTracesN[iBestTrace];
    bestTraceScore=bestTracesScores[iBestTrace];
    
    //printf("\nOptimizing gaps...\n");
    
    int *traceLen=new int [traceMaxSize], *bestTraceLen=new int [traceMaxSize];
    
    
    int strLen=0;
    
    int jt;
    strLen=0;
    nTrace=nBestTrace;  
    nGaps=0;    
    
    for(jt=0; jt<nBestTrace; jt++) {
      trace1[jt]=bestTrace1[jt];
      trace2[jt]=bestTrace2[jt];
      traceLen[jt]=winSize;
      
      if(jt<nBestTrace-1) {
	nGaps+=bestTrace1[jt+1]-bestTrace1[jt]-winSize+
	  bestTrace2[jt+1]-bestTrace2[jt]-winSize;
      }
    }    
    nBestTrace=0;
    for(int it=0; it<nTrace; ) {
      int cSize=traceLen[it];
      for(jt=it+1; jt<nTrace; jt++) {
	if(trace1[jt]-trace1[jt-1]-traceLen[jt-1]!=0 ||
	   trace2[jt]-trace2[jt-1]-traceLen[jt-1]!=0) break;
	cSize+=traceLen[jt]; 
      }
      bestTrace1[nBestTrace]=trace1[it];
      bestTrace2[nBestTrace]=trace2[it];
      bestTraceLen[nBestTrace]=cSize;
      nBestTrace++;
      strLen+=cSize;
      it=jt;
    }
    
    is=0;
    for(jt=0; jt<nBestTrace; jt++) {
      for(i=0; i<bestTraceLen[jt]; i++) {
	strBuf1[is+i]=ca1[bestTrace1[jt]+i];
	strBuf2[is+i]=ca2[bestTrace2[jt]+i];
      }
      is+=bestTraceLen[jt];
    }
    sup_str(strBuf1, strBuf2, strLen, d_);
    rmsd=calc_rmsd(strBuf1, strBuf2, strLen, d_);
    
    int isCopied=0;
    
    for(int it=1; it<nBestTrace; it++) {
      int igap;
      if(bestTrace1[it]-bestTrace1[it-1]-bestTraceLen[it-1]>0) igap=0;
      if(bestTrace2[it]-bestTrace2[it-1]-bestTraceLen[it-1]>0) igap=1;
      
      int wasBest=0;
      for(idir=-1; idir<=1; idir+=2) {
	if(wasBest) break;
	for(int idep=1; idep<=winSize/2; idep++) {
	  
	  if(!isCopied)
	    for(jt=0; jt<nBestTrace; jt++) {
	      trace1[jt]=bestTrace1[jt];
	      trace2[jt]=bestTrace2[jt];
	      traceLen[jt]=bestTraceLen[jt];
	    }
	  isCopied=0;
	  
	  traceLen[it-1]+=idir;
	  traceLen[it]-=idir;
	  trace1[it]+=idir;
	  trace2[it]+=idir;
	  
	  is=0;
	  for(jt=0; jt<nBestTrace; jt++) {
	    for(i=0; i<traceLen[jt]; i++) {
	      if(ca1[trace1[jt]+i].X>1e10 || ca2[trace2[jt]+i].X>1e10) 
		goto bad_ca;
	      strBuf1[is+i]=ca1[trace1[jt]+i];
	      strBuf2[is+i]=ca2[trace2[jt]+i];
	    }
	    is+=traceLen[jt];
	  }
	  sup_str(strBuf1, strBuf2, strLen, d_);
	  rmsdNew=calc_rmsd(strBuf1, strBuf2, strLen, d_);
	  //printf("step %d %d %d %.2f\n", it, idir, idep, rmsdNew);
	  if(rmsdNew<rmsd) {
	    for(jt=0; jt<nBestTrace; jt++) {
	      bestTrace1[jt]=trace1[jt];
	      bestTrace2[jt]=trace2[jt];
	      bestTraceLen[jt]=traceLen[jt];
	    }
	    isCopied=1;
	    wasBest=1;
	    rmsd=rmsdNew;
	    continue;
	  }
	bad_ca: break;
	}
      }
    }
    
    z=zStrAlign(winSize, strLen/winSize, bestTraceScore, nGaps);
    
    int nAtom=strLen;
    
    if(isPrint) {
      /*
      printf("size=%d rmsd=%.2f z=%.1f gaps=%d(%.1f%%) comb=%d\n", 
	     nAtom, rmsd, z, nGaps, nGaps*100.0/nAtom,
	     nTraces); 
      
      for(int k=0; k<nBestTrace; k++)
	printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
	       bestTraceLen[k]);
      printf("\n");
      */
    }
    
    if(z>=zThr) {
      
      // optimization on superposition
      XYZ *ca3=new XYZ [nse2];
      
      double oRmsdThr=2.0, rmsdLen;
      double _d[20];
      int isRmsdLenAssigned=0, nAtomPrev=-1;
      
      nAtom=0;
      while(nAtom<strLen*0.95 || 
	    (isRmsdLenAssigned && rmsd<rmsdLen*1.1 && nAtomPrev!=nAtom)) {
	nAtomPrev=nAtom;
	oRmsdThr+=0.5;
	rot_mol(ca2, ca3, nse2, d_);
	
	for(ise1=0; ise1<nse1; ise1++) {
	  for(ise2=0; ise2<nse2; ise2++) {
	    *(mat[ise1]+ise2)=-0.001;
	    
	    if(ca1[ise1].X<1e10 && ca3[ise2].X<1e10)
	      *(mat[ise1]+ise2)=oRmsdThr-ca1[ise1].dist(ca3[ise2]);
	  }
	}
	dpAlign(mat, nse1, nse2, align_se1, align_se2, lcmp, 5.0, 0.5, 0, 0);

	nAtom=0; nGaps=0; 
	for(int ia=0; ia<lcmp; ia++)
	  if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	    if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	      strBuf1[nAtom]=ca1[align_se1[ia]];
	      strBuf2[nAtom]=ca2[align_se2[ia]];
	      nAtom++;
	    }
	  }
	  else {
	    nGaps++;
	  }
	if(nAtom<4) continue;
	
	sup_str(strBuf1, strBuf2, nAtom, _d);
	rmsd=calc_rmsd(strBuf1, strBuf2, nAtom, _d);
	if(!(nAtom<strLen*0.95) && isRmsdLenAssigned==0) { 
	  rmsdLen=rmsd;
	  isRmsdLenAssigned=1;
	  }	
	//printf("nAtom %d %d rmsd %.1f\n", nAtom, nAtomPrev, rmsd);
	
      }
      /*
	nAtom=0; nGaps=0; 
	for(int ia=0; ia<lcmp; ia++)
	if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	strBuf1[nAtom]=ca1[align_se1[ia]];
	strBuf2[nAtom]=ca2[align_se2[ia]];
	nAtom++;
	}
	}
	else {
	nGaps++;
	}
	
	sup_str(strBuf1, strBuf2, nAtom, _d);
	rmsd=calc_rmsd(strBuf1, strBuf2, nAtom, _d);
      */
      nBestTrace=0;
      int newBestTrace=1;
      for(int ia=0; ia<lcmp; ia++)
	if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	  if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	    if(newBestTrace) {
	      bestTrace1[nBestTrace]=align_se1[ia];
	      bestTrace2[nBestTrace]=align_se2[ia];
	      bestTraceLen[nBestTrace]=0;
	      newBestTrace=0;
	      nBestTrace++;
	    }
	    bestTraceLen[nBestTrace-1]++;
	  }
	}
	else {
	  newBestTrace=1;
	}
      
      delete [] ca3;
      
      // end of optimization on superposition
      
      if(isPrint) {
	/*
	FILE *f=fopen("homologies", "a");
	fprintf(f, "%s(%d) %s(%d) %3d %4.1f %4.1f %d(%d) ", 
		name1, nse1, name2, nse2, nAtom, rmsd, z, 
		nGaps, nGaps*100/nAtom);
	for(int k=0; k<nBestTrace; k++)
	  fprintf(f, "(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1, 
		  bestTraceLen[k]);
	fprintf(f, "\n");
	fclose(f);
	*/
      }
    }
    else {
      for(int k=0; k<nBestTrace; k++) {
	for(int l=0; l<bestTraceLen[k]; l++) {
	  align_se1[lcmp+l]=bestTrace1[k]+l;
	  align_se2[lcmp+l]=bestTrace2[k]+l;
	}
	lcmp+=bestTraceLen[k];
	if(k<nBestTrace-1) {
	  if(bestTrace1[k]+bestTraceLen[k]!=bestTrace1[k+1])
	    for(int l=bestTrace1[k]+bestTraceLen[k]; l<bestTrace1[k+1]; l++) {
	      align_se1[lcmp]=l;
	      align_se2[lcmp]=-1;
	      lcmp++;
	    }
	  if(bestTrace2[k]+bestTraceLen[k]!=bestTrace2[k+1])
	    for(int l=bestTrace2[k]+bestTraceLen[k]; l<bestTrace2[k+1]; l++) {
	      align_se1[lcmp]=-1;
	      align_se2[lcmp]=l;
	      lcmp++;
	    }
	}
      }
      nAtom=lcmp;
    }
    
    times(&time2);
    long clk_tck=sysconf(_SC_CLK_TCK);
    long time_q=(time2.tms_stime-time1.tms_stime+
		 time2.tms_utime-time1.tms_utime);
    
    if(isPrint) {
      printf("Alignment length = %d Rmsd = %.2fA Z-Score = %.1f Gaps = %d(%.1f%%) CPU = %ds ", 
	     nAtom, rmsd, z, nGaps, nGaps*100.0/nAtom, 
	     (int)(((double)time_q)/((double)clk_tck))); 
      /*      
      for(int k=0; k<nBestTrace; k++)
	printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
	       bestTraceLen[k]);
      printf("\n");
      */
    }
    
    delete [] traceLen;
    delete [] bestTraceLen;
  }
  else {
    if(isPrint) {
      times(&time2);
      long clk_tck=sysconf(_SC_CLK_TCK);
      long time_q=(time2.tms_stime-time1.tms_stime+
		   time2.tms_utime-time1.tms_utime);
      printf("size=0 time=%d comb=%d\n", 
	     (int)(((double)time_q)/((double)clk_tck)), nTraces);
    }
  }
  
  delete [] a;
  
  DeleteArray(&mat, nse1);
  DeleteArray(&dist1, nse1);
  DeleteArray(&dist2, nse2);
  delete [] f1;
  delete [] f2;
  delete [] trace1;
  delete [] trace2;
  delete [] bestTrace1;
  delete [] bestTrace2;
  delete [] traceIterLevel;
  delete [] traceIndex;
  DeleteArray(&traceScore, traceMaxSize);
  DeleteArray(&traceUsage, traceMaxSize);
  delete [] strBuf1;
  delete [] strBuf2;

  delete [] bestTracesN;
  delete [] bestTracesScores;
  DeleteArray(&bestTraces1, bestTracesMax);
  DeleteArray(&bestTraces2, bestTracesMax);
  return(z);
}
///////////////////////////////////////////////////////////////////////////
double trace_rmsd(int *trace1, int *trace2, int nTrace, int winSize,
		  XYZ *ca1, XYZ *ca2) {

  int nStr=nTrace*winSize;
  XYZ *strBuf1=new XYZ [nStr], *strBuf2=new XYZ [nStr];

  int is=0;
  for(int it=0; it<nTrace; it++) {
    for(int i=0; i<winSize; i++) {
      strBuf1[is+i]=ca1[trace1[it]+i];
      strBuf2[is+i]=ca2[trace2[it]+i];
    }
    is+=winSize;
  }
  double d_[20];
  sup_str(strBuf1, strBuf2, nStr, d_);
  double rmsd=calc_rmsd(strBuf1, strBuf2, nStr, d_);

  delete [] strBuf1;
  delete [] strBuf2;

  return(rmsd-2.0);
}
///////////////////////////////////////////////////////////////////////////
void origin3(XYZ *p, double *d) {
  XYZ t[3], c;
  double r1[3][3], r2[3][3];
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) r1[i][j]=r2[i][j]=i==j?1.0:0.0;
  
  for(int i=0; i<3; i++) t[i]=p[i];
  for(int i=0; i<3; i++) t[i]-=p[1];

  double p_r=sqrt(t[0].X*t[0].X+t[0].Y*t[0].Y);
  double q_r=sqrt(t[0].X*t[0].X+t[0].Y*t[0].Y+t[0].Z*t[0].Z);
  
  double p_cos=t[0].X/p_r, p_sin=-t[0].Y/p_r;

  r1[0][0]=p_cos; r1[0][1]=-p_sin;
  r1[1][0]=p_sin; r1[1][1]=p_cos;
  
  double q_cos=p_r/q_r, q_sin=-t[0].Z/q_r;
  
  r2[0][0]=q_cos; r2[0][2]=-q_sin;
  r2[2][0]=q_sin; r2[2][2]=q_cos;

  rot(r2, r1);

  t[0]=rotXYZ(t[0], r1);
  t[2]=rotXYZ(t[2], r1);

  double s_r=sqrt(t[2].Y*t[2].Y+t[2].Z*t[2].Z);
  
  double s_cos=t[2].Z/s_r, s_sin=t[2].Y/s_r;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) r2[i][j]=i==j?1.0:0.0;

  r2[1][1]=s_cos; r2[1][2]=-s_sin;
  r2[2][1]=s_sin; r2[2][2]=s_cos;
  
  rot(r2, r1);

  
  for(int i=0; i<3; i++) t[i]=p[i];
  for(int i=0; i<3; i++) t[i]-=p[1];

  t[0]=rotXYZ(t[0], r1);
  t[2]=rotXYZ(t[2], r1);


  d[0]=r1[0][0]; d[1]=r1[0][1]; d[2]=r1[0][2];
  d[3]=r1[1][0]; d[4]=r1[1][1]; d[5]=r1[1][2];
  d[6]=r1[2][0]; d[7]=r1[2][1]; d[8]=r1[2][2];

  c=0.0;
  c-=p[1];

  c=rotXYZ(c, r1);
  
  d[9]=c.X; d[10]=c.Y; d[11]=c.Z;

}
///////////////////////////////////////////////////////////////////////////
void rot(double m1[3][3], double m2[3][3]) {
  double m3[3][3], m; int i,j,k;
  for(i=0; i<3; i++)
    for(j=0; j<3; j++) {
      m=0.0;
      for(k=0; k<3; k++) m+=(m1[i][k]*m2[k][j]);
      m3[i][j]=m;
    }
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)  m2[i][j]=m3[i][j];
}
///////////////////////////////////////////////////////////////////////////
XYZ rotXYZ(XYZ p, double r[3][3]) {
  XYZ p_;
  p_.X=r[0][0]*p.X+r[0][1]*p.Y+r[0][2]*p.Z;
  p_.Y=r[1][0]*p.X+r[1][1]*p.Y+r[1][2]*p.Z;
  p_.Z=r[2][0]*p.X+r[2][1]*p.Y+r[2][2]*p.Z;
  return(p_);
}
///////////////////////////////////////////////////////////////////////////
XYZ transXYZ(XYZ p, double d[20]) {
  XYZ p_;
  p_.X=d[0]*p.X+d[1]*p.Y+d[2]*p.Z+d[9];
  p_.Y=d[3]*p.X+d[4]*p.Y+d[5]*p.Z+d[10];
  p_.Z=d[6]*p.X+d[7]*p.Y+d[8]*p.Z+d[11];
  return(p_);
}
///////////////////////////////////////////////////////////////////////////
double scoreAv8[]={2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74,
		  4.99, 5.22, 5.46, 5.70, 5.94, 6.13, 6.36, 6.52, 6.68, 6.91};
double scoreSd8[]={1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04,
		  1.08, 1.10, 1.15, 1.19, 1.23, 1.25, 1.32, 1.34, 1.36, 1.45};
double gapsAv8[]={0.00, 11.50, 23.32, 35.95, 49.02, 62.44, 76.28, 90.26, 
		 104.86, 119.97, 134.86, 150.54, 164.86, 179.57, 194.39, 
		 209.38, 224.74, 238.96, 253.72, 270.79};
double gapsSd8[]={0.00, 9.88, 14.34, 17.99, 21.10, 23.89, 26.55, 29.00, 31.11,
		 33.10, 35.02, 36.03, 37.19, 38.82, 41.04, 43.35, 45.45, 
		 48.41, 50.87, 52.27};
double scoreAv6[]={1.98, 1.97, 2.22, 2.54, 2.87, 3.18, 3.48, 3.77, 4.05, 4.31,
		   4.57, 4.82, 5.03, 5.24, 5.43, 5.64, 5.82, 6.02, 6.21, 6.42};
double scoreSd6[]={1.15, 0.73, 0.63, 0.64, 0.71, 0.80, 0.87, 0.95, 1.01, 1.07,
		   1.13, 1.19, 1.22, 1.25, 1.28, 1.32, 1.35, 1.39, 1.45, 1.50};
double gapsAv6[]={0.00, 10.12, 20.25, 31.29, 42.95, 55.20, 67.53, 80.15, 
		  93.30, 106.47, 120.52, 134.38, 148.59, 162.58, 176.64, 
		  191.23, 204.12, 218.64, 231.82, 243.43};
double gapsSd6[]={0.00, 9.80, 14.44, 18.14, 21.35, 24.37, 27.00, 29.68, 32.22,
		  34.37, 36.65, 38.63, 40.31, 42.16, 43.78, 44.98, 47.08, 
		  49.09, 50.78, 52.15};
///////////////////////////////////////////////////////////////////////////
double zScore(int winSize, int nTrace, double score) {
  
  if(winSize==8) {
    
    if(nTrace<1) return(0.0);
    
    double scoreAv_, scoreSd_;
    if(nTrace<21) {
      scoreAv_=scoreAv8[nTrace-1];
      scoreSd_=scoreSd8[nTrace-1];
    }
    else {
      scoreAv_=0.209874*nTrace+2.944714;
      scoreSd_=0.039487*nTrace+0.675735;
    }
    if(score>scoreAv_) return(0.0);
    return((scoreAv_-score)/scoreSd_);
  }

  if(winSize==6) {
    
    if(nTrace<1) return(0.0);
    
    double scoreAv_, scoreSd_;
    if(nTrace<21) {
      scoreAv_=scoreAv6[nTrace-1];
      scoreSd_=scoreSd6[nTrace-1];
    }
    else {
      scoreAv_=0.198534*nTrace+2.636477;
      scoreSd_=0.040922*nTrace+0.715636;
    }
    if(score>scoreAv_) return(0.0);
    return((scoreAv_-score)/scoreSd_);
  }

  return(0.0);

}
///////////////////////////////////////////////////////////////////////////
double zGaps(int winSize, int nTrace, int nGaps) {

  if(nTrace<2) return(0.0);
  double scoreAv_, scoreSd_;

  if(winSize==8) {
    if(nTrace<21) {
      scoreAv_=gapsAv8[nTrace-1];
      scoreSd_=gapsSd8[nTrace-1];
    }
    else {
      scoreAv_=14.949173*nTrace-14.581193;
      scoreSd_=2.045067*nTrace+13.191095;
    }
    if(nGaps>scoreAv_) return(0.0);
    return((scoreAv_-nGaps)/scoreSd_);
  }

  if(winSize==6) {
    if(nTrace<21) {
      scoreAv_=gapsAv6[nTrace-1];
      scoreSd_=gapsSd6[nTrace-1];
    }
    else {
      scoreAv_=13.574490*nTrace-13.977223;
      scoreSd_=1.719977*nTrace+19.615014;
    }
    if(nGaps>scoreAv_) return(0.0);
    return((scoreAv_-nGaps)/scoreSd_);
  }

  return(0.0);
}
///////////////////////////////////////////////////////////////////////////
double zStrAlign(int winSize, int nTrace, double score, int nGaps) {
  double z1=zScore(winSize, nTrace, score);
  double z2=zGaps(winSize, nTrace, nGaps);
  return(zByZ(z1, z2));
}
///////////////////////////////////////////////////////////////////////////
double tableZtoP[]={
1.0, 9.20e-01,8.41e-01,7.64e-01,6.89e-01,6.17e-01,5.49e-01,4.84e-01,4.24e-01,3.68e-01,
3.17e-01,2.71e-01,2.30e-01,1.94e-01,1.62e-01,1.34e-01,1.10e-01,8.91e-02,7.19e-02,5.74e-02,
4.55e-02,3.57e-02,2.78e-02,2.14e-02,1.64e-02,1.24e-02,9.32e-03,6.93e-03,5.11e-03,3.73e-03,
2.70e-03,1.94e-03,1.37e-03,9.67e-04,6.74e-04,4.65e-04,3.18e-04,2.16e-04,1.45e-04,9.62e-05,
6.33e-05,4.13e-05,2.67e-05,1.71e-05,1.08e-05,6.80e-06,4.22e-06,2.60e-06,1.59e-06,9.58e-07,
5.73e-07,3.40e-07,1.99e-07,1.16e-07,6.66e-08,3.80e-08,2.14e-08,1.20e-08,6.63e-09,3.64e-09,
1.97e-09,1.06e-09,5.65e-10,2.98e-10,1.55e-10,8.03e-11,4.11e-11,2.08e-11,1.05e-11,5.20e-12,
2.56e-12,1.25e-12,6.02e-13,2.88e-13,1.36e-13,6.38e-14,2.96e-14,1.36e-14,6.19e-15,2.79e-15,
1.24e-15,5.50e-16,2.40e-16,1.04e-16,4.46e-17,1.90e-17,7.97e-18,3.32e-18,1.37e-18,5.58e-19,
2.26e-19,9.03e-20,3.58e-20,1.40e-20,5.46e-21,2.10e-21,7.99e-22,3.02e-22,1.13e-22,4.16e-23,
1.52e-23,5.52e-24,1.98e-24,7.05e-25,2.48e-25,8.64e-26,2.98e-26,1.02e-26,3.44e-27,1.15e-27,
3.82e-28,1.25e-28,4.08e-29,1.31e-29,4.18e-30,1.32e-30,4.12e-31,1.27e-31,3.90e-32,1.18e-32,
3.55e-33,1.06e-33,3.11e-34,9.06e-35,2.61e-35,7.47e-36,2.11e-36,5.91e-37,1.64e-37,4.50e-38,
1.22e-38,3.29e-39,8.77e-40,2.31e-40,6.05e-41,1.56e-41,4.00e-42,1.02e-42,2.55e-43,6.33e-44,
1.56e-44,3.80e-45,9.16e-46,2.19e-46,5.17e-47,1.21e-47,2.81e-48,6.45e-49,1.46e-49,3.30e-50};
double tablePtoZ[]={
0.00,0.73,1.24,1.64,1.99,2.30,2.58,2.83,3.07,3.29,
3.50,3.70,3.89,4.07,4.25,4.42,4.58,4.74,4.89,5.04,
5.19,5.33,5.46,5.60,5.73,5.86,5.99,6.11,6.23,6.35,
6.47,6.58,6.70,6.81,6.92,7.02,7.13,7.24,7.34,7.44,
7.54,7.64,7.74,7.84,7.93,8.03,8.12,8.21,8.30,8.40,
8.49,8.57,8.66,8.75,8.84,8.92,9.01,9.09,9.17,9.25,
9.34,9.42,9.50,9.58,9.66,9.73,9.81,9.89,9.97,10.04,
10.12,10.19,10.27,10.34,10.41,10.49,10.56,10.63,10.70,10.77,
10.84,10.91,10.98,11.05,11.12,11.19,11.26,11.32,11.39,11.46,
11.52,11.59,11.66,11.72,11.79,11.85,11.91,11.98,12.04,12.10,
12.17,12.23,12.29,12.35,12.42,12.48,12.54,12.60,12.66,12.72,
12.78,12.84,12.90,12.96,13.02,13.07,13.13,13.19,13.25,13.31,
13.36,13.42,13.48,13.53,13.59,13.65,13.70,13.76,13.81,13.87,
13.92,13.98,14.03,14.09,14.14,14.19,14.25,14.30,14.35,14.41,
14.46,14.51,14.57,14.62,14.67,14.72,14.77,14.83,14.88,14.93};
///////////////////////////////////////////////////////////////////////////
double zToP(double z) {
  int ind=(int)(z/0.1);
  if(ind<0) ind=0;
  if(ind>149) ind=149;
  return(tableZtoP[ind]);
}
///////////////////////////////////////////////////////////////////////////
double pToZ(double p) {
  int ind=(int)(-log10(p)*3.0);
  if(ind<0) ind=0;
  if(ind>149) ind=149;
  return(tablePtoZ[ind]);
}
///////////////////////////////////////////////////////////////////////////
double zByZ(double z1, double z2) {
  double p1=zToP(z1), p2=zToP(z2);
  return(pToZ(p1*p2));
}
///////////////////////////////////////////////////////////////////////////
double ce_a1(char *name1, char *name2,
	      XYZ *ca1, XYZ *ca2, int nse1, int nse2, int *align_se1, 
	      int *align_se2, int& lcmp, double rmsdThrJoin, int isPrint) {
  /*
  int distAll=0;
  int winSize_=8;
  double rmsdThr=3.0;
  double zThr=3.5;
  int bestTracesMax=30;
  int gapMax=30;
  long tracesLimit=5e7;

  if(isPrint) printf("%s(%d) %s(%d) ", name1, nse1, name2, nse2);

  tms time1, time2;
  times(&time1);
  long nTraces;
  
  double rmsdThrSq=rmsdThr*rmsdThr;
  double bestTraceScore=10000.0;

  double d, dd, d1, d2, dx, dy, dz, score, score0, score1, score2, 
    d_[20], rmsd, rmsdNew, z; 

  int jse1, jse2, kse1, kse2, mse1, mse2, l1, l2, i, j, nd;

  int traceMaxSize=nse1<nse2?nse1:nse2;

  int **bestTraces1, **bestTraces2, *bestTracesN=new int [bestTracesMax];
  double *bestTracesScores=new double [bestTracesMax];
  NewArray(&bestTraces1, bestTracesMax, traceMaxSize);
  NewArray(&bestTraces2, bestTracesMax, traceMaxSize);
  int nBestTraces=0, newBestTrace=0;

  for(int it=0; it<bestTracesMax; it++) {
    bestTracesN[it]=0;
    bestTracesScores[it]=100;
  }

  int *f1=new int [nse1], *f2=new int [nse2];
  double **mat;
  double **dist1, **dist2;
  NewArray(&mat, nse1, nse2);
  NewArray(&dist1, nse1, nse1);
  NewArray(&dist2, nse2, nse2);
  int iterDepth=gapMax*2+1;
  int *trace1=new int [traceMaxSize], *trace2=new int [traceMaxSize];
  int *traceIterLevel=new int [traceMaxSize],
    *traceIndex=new int [traceMaxSize];
  double **traceScore;
  int **traceUsage;
  NewArray(&traceScore, traceMaxSize, iterDepth);
  NewArray(&traceUsage, traceMaxSize, iterDepth);
  
  int nTrace, nGaps, igap, idir, jgap, jdir, nBestTrace=0;
  XYZ *strBuf1=new XYZ [traceMaxSize], *strBuf2=new XYZ [traceMaxSize];
  int *bestTrace1=new int [traceMaxSize], *bestTrace2=new int [traceMaxSize];

  int winSize=winSize_;
  
  for(int is=0; is<nse1-winSize+1; is++) {
    f1[is]=0;
    for(int ir=0; ir<winSize; ir++) if(ca1[is+ir].X>=1e10) goto cont1;
    f1[is]=1;
  cont1: ;
  }

  for(int is=0; is<nse2-winSize+1; is++) {
    f2[is]=0;
    for(int ir=0; ir<winSize; ir++) if(ca2[is+ir].X>=1e10) goto cont2;
    f2[is]=1;
  cont2: ;
  }

  for(int ise1=0; ise1<nse1; ise1++) 
    for(int ise2=0; ise2<nse1; ise2++) 
      *(dist1[ise1]+ise2)=
	ca1[ise1].X<1e10&&ca1[ise2].X<1e10?ca1[ise1].dist(ca1[ise2]):2e10;
  
  for(int ise1=0; ise1<nse2; ise1++) 
    for(int ise2=0; ise2<nse2; ise2++) 
      *(dist2[ise1]+ise2)=
	ca2[ise1].X<1e10&&ca2[ise2].X<1e10?ca2[ise1].dist(ca2[ise2]):2e10;
  
  times(&time2);
  
  long clk_tck=sysconf(_SC_CLK_TCK);
  long time_q=(time2.tms_stime-time1.tms_stime+
	       time2.tms_utime-time1.tms_utime);
  
  int t=(int)(((double)time_q)/((double)clk_tck));

  int isGood;
  
  int winSizeComb1=(winSize-1)*(winSize-2)/2;
  int winSizeComb2=distAll?winSize*winSize:winSize;

  int *a=new int [traceMaxSize];

  if(distAll)
    for(i=0; i<traceMaxSize; i++)
      a[i]=(i+1)*i*winSize*winSize/2+(i+1)*winSizeComb1;
  else
    for(i=0; i<traceMaxSize; i++)
      a[i]=(i+1)*i*winSize/2+(i+1)*winSizeComb1;
  
  for(int ise1=0; ise1<nse1; ise1++) {
    for(int ise2=0; ise2<nse2; ise2++) {
      *(mat[ise1]+ise2)=-1.0;
      if(ise1>nse1-winSize || ise2>nse2-winSize) continue;
      
      isGood=1;
      
      for(l1=0; l1<winSize; l1++) 
	if(ca1[ise1+l1].X>1e10 || ca2[ise2+l1].X>1e10) {
	    isGood=0; break;
	}
      
      if(!isGood) continue;
      
      d=0.0; 
      for(int is1=0; is1<winSize-2; is1++)
	for(int is2=is1+2; is2<winSize; is2++) {
	  d+=fabs(*(dist1[ise1+is1]+ise1+is2)-
		      *(dist2[ise2+is1]+ise2+is2));
	}
      *(mat[ise1]+ise2)=d/winSizeComb1; 
      
    }
  }

  // processing afp matrix - begin
  /*
  int **rep;
  NewArray(&rep, nse1, nse2);
  for(int ise1=0; ise1<nse1; ise1++) 
    for(int ise2=0; ise2<nse2; ise2++) *(rep[ise1]+ise2)=0;
      
   for(int ise1=0; ise1<nse1-4; ise1++) 
    for(int ise2=0; ise2<nse2-4; ise2++) 
      if(*(mat[ise1]+ise2)>=0.0 && *(rep[ise1]+ise2)==0) 
	for(int im=1; im<=4; im++)
	  if(*(mat[ise1+im]+ise2+im)>0.0 &&
	     *(mat[ise1+im]+ise2+im)<rmsdThr) {
	    *(rep[ise1+im]+ise2+im)=1;
	    *(mat[ise1]+ise2)=100.0;
	    break;
	  }
   DeleteArray(&rep, nse1);
  */
  // processing afp matrix - end
    
  /*
    printf("  ");  
    for(i=0; i<nse2; i++) printf("%c", i%10==0?(i%100)/10+48:32);
    printf("\n");  
    for(i=0; i<nse1; i++) {
    printf("%c ", i%10==0?(i%100)/10+48:32);
    for(j=0; j<nse2; j++) 
    printf("%c", *(mat[i]+j)<rmsdThr?'+':'X');
    //printf("%c", ((int)*(mat[i]+j)/40)>9?'*':((int)*(mat[i]+j)/40)+48);
    printf("\n");
    }
    printf("\n");
    */
  
  // tracing trough fragment matrix
  /*
  nBestTrace=0;    
  bestTraceScore=100.0;
  nTraces=0;

  for(int ise1=0; ise1<nse1; ise1++) {
    for(int ise2=0; ise2<nse2; ise2++) {
      
      if((ise1>nse1-winSize*(nBestTrace-1) || 
	  ise2>nse2-winSize*(nBestTrace-1))) continue;
      
      if(*(mat[ise1]+ise2)<0.0) continue;
      if(*(mat[ise1]+ise2)>rmsdThr) continue;
      
      nTrace=0;
      trace1[nTrace]=ise1; 
      trace2[nTrace]=ise2;
      traceIndex[nTrace]=0;
      traceIterLevel[nTrace]=0;
      
      score0=*(mat[ise1]+ise2);
      
      nTrace++;
      int isTraceUp=1;
      int traceIndex_=0;
      
      while(nTrace>0) {
	
	kse1=trace1[nTrace-1]+winSize;
	kse2=trace2[nTrace-1]+winSize;
	
	while(1) {
	  if(kse1>nse1-winSize-1) break;
	  if(kse2>nse2-winSize-1) break;
	  if(*(mat[kse1]+kse2)>=0.0) break;
	  if(f1[kse1]==0) kse1++;
	  if(f2[kse2]==0) kse2++;
	}
	  
	traceIndex_=-1; 
	
	if(isTraceUp) {
	  int nBestExtTrace=nTrace;
	  double bestExtScore=100.0;
	  
	  for(int it=0; it<iterDepth; it++) {
	    jgap=(it+1)/2;
	    jdir=(it+1)%2;
	    if(jdir==0) {
	      mse1=kse1+jgap;
	      mse2=kse2;
	    }
	    else {
	      mse1=kse1;
	      mse2=kse2+jgap;
	    }
	    
	    if(mse1>nse1-winSize-1) continue;
	    if(mse2>nse2-winSize-1) continue;
	    
	    if(*(mat[mse1]+mse2)<0.0) continue;
	    if(*(mat[mse1]+mse2)>rmsdThr) continue;
	    
	    nTraces++;
	    if(nTraces>tracesLimit) {
	      if(isPrint) printf("exited ");
	      goto exit_tracing;
	    }
	    
	    score=0.0;
	    if(!distAll) {
	      // (winSize) "best" dist
	      for(int itrace=0; itrace<nTrace; itrace++) {
		score+=fabs(*(dist1[trace1[itrace]]+mse1)-
			    *(dist2[trace2[itrace]]+mse2));
		score+=fabs(*(dist1[trace1[itrace]+winSize-1]+
			      mse1+winSize-1)-
			    *(dist2[trace2[itrace]+winSize-1]+
			      mse2+winSize-1));
		
		for(int id=1; id<winSize-1; id++) 
		  score+=fabs(*(dist1[trace1[itrace]+id]+mse1+winSize-1-id)-
			      *(dist2[trace2[itrace]+id]+mse2+winSize-1-id));
		
	      }
	      score1=score/(nTrace*winSize);
	    }
	    else {
	      // all dist
	      for(int itrace=0; itrace<nTrace; itrace++) {
		for(int is1=0; is1<winSize; is1++)
		  for(int is2=0; is2<winSize; is2++)
		    score+=fabs(*(dist1[trace1[itrace]+is1]+mse1+is2)-
				*(dist2[trace2[itrace]+is1]+mse2+is2));
	      }
	      score1=score/(nTrace*winSize*winSize);
	    }
	      
	    if(score1>rmsdThrJoin) continue;
	    
	    if(nTrace>nBestExtTrace || (nTrace==nBestExtTrace &&
					score1<bestExtScore)) {
	      bestExtScore=score1;
	      nBestExtTrace=nTrace;
	      traceIndex_=it;
	      *(traceScore[nTrace-1]+traceIndex_)=score1;
	    }
	  }
	}
	
	double rmsdScore, traceTotalScore, traceScoreMax;
	
	if(traceIndex_!=-1) {
	  jgap=(traceIndex_+1)/2;
	  jdir=(traceIndex_+1)%2;
	  if(jdir==0) {
	    jse1=kse1+jgap;
	    jse2=kse2;
	  }
	  else {
	    jse1=kse1;
	    jse2=kse2+jgap;
	  }
	  
	  score1=(*(traceScore[nTrace-1]+traceIndex_)*winSizeComb2*nTrace+
		  *(mat[jse1]+jse2)*winSizeComb1)/(winSizeComb2*nTrace+
						   winSizeComb1);
	  
	  score2=
	    ((nTrace>1?*(traceScore[nTrace-2]+traceIndex[nTrace-1]):score0)
	     *a[nTrace-1]+score1*(a[nTrace]-a[nTrace-1]))/a[nTrace];
	  
	  if(score2>rmsdThrJoin) traceIndex_=-1;
	  else traceTotalScore=*(traceScore[nTrace-1]+traceIndex_)=score2;
	}

	if(traceIndex_==-1) {
	  nTrace--;
	  isTraceUp=0;
	  continue;
	}
	else {
	  traceIterLevel[nTrace-1]++;
	  trace1[nTrace]=jse1;
	  trace2[nTrace]=jse2;
	  traceIndex[nTrace]=traceIndex_;
	  traceIterLevel[nTrace]=0;
	  nTrace++;
	  isTraceUp=1;
	  
	  if(nTrace>nBestTrace || 
	     (nTrace==nBestTrace  && 
	      bestTraceScore>traceTotalScore)) {
	    
	    for(int itrace=0; itrace<nTrace; itrace++) {
	      bestTrace1[itrace]=trace1[itrace];
	      bestTrace2[itrace]=trace2[itrace];
	    }
	    bestTraceScore=traceTotalScore;
	    nBestTrace=nTrace;
	  }
	  
	  if(nTrace>bestTracesN[newBestTrace] ||
	     (nTrace==bestTracesN[newBestTrace] && 
	      bestTracesScores[newBestTrace]>traceTotalScore)) {
	    for(int itrace=0; itrace<nTrace; itrace++) {
	      *(bestTraces1[newBestTrace]+itrace)=trace1[itrace];
	      *(bestTraces2[newBestTrace]+itrace)=trace2[itrace];
	      bestTracesN[newBestTrace]=nTrace;
	      bestTracesScores[newBestTrace]=traceTotalScore;
	    }
	    
	    if(nTrace>nBestTrace) nBestTrace=nTrace;
	    
	    if(nBestTraces<bestTracesMax) {
	      nBestTraces++;
	      newBestTrace++;
	    }
	    
	    if(nBestTraces==bestTracesMax) {
	      newBestTrace=0;
	      double scoreTmp=bestTracesScores[0];
	      int nTraceTmp=bestTracesN[0];
	      for(int ir=1; ir<nBestTraces; ir++) {
		if(bestTracesN[ir]<nTraceTmp || 
		   (bestTracesN[ir]==nTraceTmp && 
		    scoreTmp<bestTracesScores[ir])) {
		  nTraceTmp=bestTracesN[ir];
		  scoreTmp=bestTracesScores[ir];
		  newBestTrace=ir;
		}
	      }
	    }
	  }
	  
	}
      }
    }
  }
  
 exit_tracing: ;
  int is=0;
  for(int jt_=0; jt_<nBestTrace; jt_++) {
    for(i=0; i<winSize; i++) {
      strBuf1[is+i]=ca1[bestTrace1[jt_]+i];
      strBuf2[is+i]=ca2[bestTrace2[jt_]+i];
    }
    is+=winSize;
  }
  
  sup_str(strBuf1, strBuf2, nBestTrace*winSize, d_);
  rmsdNew=calc_rmsd(strBuf1, strBuf2, nBestTrace*winSize, d_);
  /*
    printf("rmsd=%.2f ", rmsdNew);
    for(int k=0; k<nBestTrace; k++)
    printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
    8);
    printf("\n");
  */
  /**
  rmsd=100.0;
  int iBestTrace=0;
  for(int ir=0; ir<nBestTraces; ir++) {
    if(bestTracesN[ir]!=nBestTrace) continue;
    int is=0;
    for(int jt_=0; jt<bestTracesN[ir]; jt++) {
      for(i=0; i<winSize; i++) {
	strBuf1[is+i]=ca1[*(bestTraces1[ir]+jt_)+i];
	strBuf2[is+i]=ca2[*(bestTraces2[ir]+jt_)+i];
      }
      is+=winSize;
    }
    sup_str(strBuf1, strBuf2, bestTracesN[ir]*winSize, d_);
    rmsdNew=calc_rmsd(strBuf1, strBuf2, bestTracesN[ir]*winSize, d_);
    //printf("%d %d %d %.2f\n", ir, bestTracesN[ir], nBestTrace, rmsdNew);
    if(rmsd>rmsdNew) {
      iBestTrace=ir;
      rmsd=rmsdNew;
    }
  }
  for(int it=0; it<bestTracesN[iBestTrace]; it++) {
    bestTrace1[it]=*(bestTraces1[iBestTrace]+it);
    bestTrace2[it]=*(bestTraces2[iBestTrace]+it);
  }
  nBestTrace=bestTracesN[iBestTrace];
  bestTraceScore=bestTracesScores[iBestTrace];
  
  //printf("\nOptimizing gaps...\n");
  
  int *traceLen=new int [traceMaxSize], *bestTraceLen=new int [traceMaxSize];
  
  lcmp=0;
  
  int strLen=0;
  if(nBestTrace>0) {
    
    int jt;
    strLen=0;
    nTrace=nBestTrace;  
    nGaps=0;    
    
    for(jt=0; jt<nBestTrace; jt++) {
      trace1[jt]=bestTrace1[jt];
      trace2[jt]=bestTrace2[jt];
      traceLen[jt]=winSize;
      
      if(jt<nBestTrace-1) {
	nGaps+=bestTrace1[jt+1]-bestTrace1[jt]-winSize+
	  bestTrace2[jt+1]-bestTrace2[jt]-winSize;
      }
    }    
    nBestTrace=0;
    for(int it=0; it<nTrace; ) {
      int cSize=traceLen[it];
      for(jt=it+1; jt<nTrace; jt++) {
	if(trace1[jt]-trace1[jt-1]-traceLen[jt-1]!=0 ||
	   trace2[jt]-trace2[jt-1]-traceLen[jt-1]!=0) break;
	cSize+=traceLen[jt]; 
      }
      bestTrace1[nBestTrace]=trace1[it];
      bestTrace2[nBestTrace]=trace2[it];
      bestTraceLen[nBestTrace]=cSize;
      nBestTrace++;
      strLen+=cSize;
      it=jt;
    }
    
    int is=0;
    for(jt=0; jt<nBestTrace; jt++) {
      for(i=0; i<bestTraceLen[jt]; i++) {
	strBuf1[is+i]=ca1[bestTrace1[jt]+i];
	strBuf2[is+i]=ca2[bestTrace2[jt]+i];
      }
      is+=bestTraceLen[jt];
    }
    sup_str(strBuf1, strBuf2, strLen, d_);
    rmsd=calc_rmsd(strBuf1, strBuf2, strLen, d_);
    
    int isCopied=0;
    
    for(int it=1; it<nBestTrace; it++) {
      int igap;
      if(bestTrace1[it]-bestTrace1[it-1]-bestTraceLen[it-1]>0) igap=0;
      if(bestTrace2[it]-bestTrace2[it-1]-bestTraceLen[it-1]>0) igap=1;
      
      int wasBest=0;
      for(idir=-1; idir<=1; idir+=2) {
	if(wasBest) break;
	for(int idep=1; idep<=winSize/2; idep++) {
	  
	  if(!isCopied)
	    for(jt=0; jt<nBestTrace; jt++) {
	      trace1[jt]=bestTrace1[jt];
	      trace2[jt]=bestTrace2[jt];
	      traceLen[jt]=bestTraceLen[jt];
	    }
	  isCopied=0;
	  
	  traceLen[it-1]+=idir;
	  traceLen[it]-=idir;
	  trace1[it]+=idir;
	  trace2[it]+=idir;
	  
	  is=0;
	  for(jt=0; jt<nBestTrace; jt++) {
	    for(i=0; i<traceLen[jt]; i++) {
	      if(ca1[trace1[jt]+i].X>1e10 || ca2[trace2[jt]+i].X>1e10) 
		goto bad_ca;
	      strBuf1[is+i]=ca1[trace1[jt]+i];
	      strBuf2[is+i]=ca2[trace2[jt]+i];
	    }
	    is+=traceLen[jt];
	  }
	  sup_str(strBuf1, strBuf2, strLen, d_);
	  rmsdNew=calc_rmsd(strBuf1, strBuf2, strLen, d_);
	  //printf("step %d %d %d %.2f\n", it, idir, idep, rmsdNew);
	  if(rmsdNew<rmsd) {
	    for(jt=0; jt<nBestTrace; jt++) {
	      bestTrace1[jt]=trace1[jt];
	      bestTrace2[jt]=trace2[jt];
	      bestTraceLen[jt]=traceLen[jt];
	    }
	    isCopied=1;
	    wasBest=1;
	    rmsd=rmsdNew;
	    continue;
	  }
	bad_ca: break;
	}
      }
    }
    
    z=zStrAlign(winSize, strLen/winSize, bestTraceScore, nGaps);

    int nAtom=strLen;
    
    if(isPrint) {
      printf("size=%d rmsd=%.2f z=%.1f gaps=%d(%.1f%%) comb=%d\n", 
	     nAtom, rmsd, z, nGaps, nGaps*100.0/nAtom,
	     nTraces); 
      
      for(int k=0; k<nBestTrace; k++)
	printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
	       bestTraceLen[k]);
      printf("\n");
    }
   
    if(z>=zThr) {

      // optimization on superposition
      XYZ *ca3=new XYZ [nse2];
      
      double oRmsdThr=2.0, rmsdLen;
      double _d[20];
      int isRmsdLenAssigned=0, nAtomPrev=-1;

      nAtom=0;
      while(nAtom<strLen*0.95 || 
	    (isRmsdLenAssigned && rmsd<rmsdLen*1.1 && nAtomPrev!=nAtom)) {
	nAtomPrev=nAtom;
	oRmsdThr+=0.5;
	rot_mol(ca2, ca3, nse2, d_);
	
	for(int ise1=0; ise1<nse1; ise1++) {
	  for(int ise2=0; ise2<nse2; ise2++) {
	    *(mat[ise1]+ise2)=-0.001;
	    
	    if(ca1[ise1].X<1e10 && ca3[ise2].X<1e10)
	      *(mat[ise1]+ise2)=oRmsdThr-ca1[ise1].dist(ca3[ise2]);
	  }
	}
	dpAlign(mat, nse1, nse2, align_se1, align_se2, lcmp, 5.0, 0.5, 0, 0);

	nAtom=0; 
	for(int ia=0; ia<lcmp; ia++)
	  if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	    if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10)
	      nAtom++;
	  }

	nAtom=0; nGaps=0; 
	for(int ia=0; ia<lcmp; ia++)
	  if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	    if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	      strBuf1[nAtom]=ca1[align_se1[ia]];
	      strBuf2[nAtom]=ca2[align_se2[ia]];
	      nAtom++;
	    }
	  }
	  else {
	    nGaps++;
	  }
      
	sup_str(strBuf1, strBuf2, nAtom, _d);
	rmsd=calc_rmsd(strBuf1, strBuf2, nAtom, _d);
	if(!(nAtom<strLen*0.95) && isRmsdLenAssigned==0) { 
	  rmsdLen=rmsd;
	  isRmsdLenAssigned=1;
	}	
	//printf("nAtom %d %d rmsd %.1f\n", nAtom, nAtomPrev, rmsd);
	
      }
      /*
      nAtom=0; nGaps=0; 
      for(int ia=0; ia<lcmp; ia++)
	if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	  if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	    strBuf1[nAtom]=ca1[align_se1[ia]];
	    strBuf2[nAtom]=ca2[align_se2[ia]];
	    nAtom++;
	  }
	}
	else {
	  nGaps++;
	}
      
      sup_str(strBuf1, strBuf2, nAtom, _d);
      rmsd=calc_rmsd(strBuf1, strBuf2, nAtom, _d);
      */
  /*
      nBestTrace=0;
      int newBestTrace=1;
      for(int ia=0; ia<lcmp; ia++)
	if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
	  if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
	    if(newBestTrace) {
	      bestTrace1[nBestTrace]=align_se1[ia];
	      bestTrace2[nBestTrace]=align_se2[ia];
	      bestTraceLen[nBestTrace]=0;
	      newBestTrace=0;
	      nBestTrace++;
	    }
	    bestTraceLen[nBestTrace-1]++;
	  }
	}
	else {
	  newBestTrace=1;
	}
      
      delete [] ca3;
      
      // end of optimization on superposition
      
      if(isPrint) {
	FILE *f=fopen("homologies", "a");
	fprintf(f, "%s(%d) %s(%d) %3d %4.1f %4.1f %d(%d) ", 
		name1, nse1, name2, nse2, nAtom, rmsd, z, 
		nGaps, nGaps*100/nAtom);
	for(int k=0; k<nBestTrace; k++)
	  fprintf(f, "(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1, 
		  bestTraceLen[k]);
	fprintf(f, "\n");
	fclose(f);
      }
    }
    else {
      for(int k=0; k<nBestTrace; k++) {
	for(int l=0; l<bestTraceLen[k]; l++) {
	  align_se1[lcmp+l]=bestTrace1[k]+l;
	  align_se2[lcmp+l]=bestTrace2[k]+l;
	}
	lcmp+=bestTraceLen[k];
	if(k<nBestTrace-1) {
	  if(bestTrace1[k]+bestTraceLen[k]!=bestTrace1[k+1])
	    for(int l=bestTrace1[k]+bestTraceLen[k]; l<bestTrace1[k+1]; l++) {
	      align_se1[lcmp]=l;
	      align_se2[lcmp]=-1;
	      lcmp++;
	    }
	  if(bestTrace2[k]+bestTraceLen[k]!=bestTrace2[k+1])
	    for(int l=bestTrace2[k]+bestTraceLen[k]; l<bestTrace2[k+1]; l++) {
	      align_se1[lcmp]=-1;
	      align_se2[lcmp]=l;
	      lcmp++;
	    }
	}
      }
      nAtom=lcmp;
    }

    times(&time2);
    long clk_tck=sysconf(_SC_CLK_TCK);
    long time_q=(time2.tms_stime-time1.tms_stime+
		 time2.tms_utime-time1.tms_utime);

    if(isPrint) {
      printf("size=%d rmsd=%.2f z=%.1f gaps=%d(%.1f%%) time=%d comb=%d\n", 
	     nAtom, rmsd, z, nGaps, nGaps*100.0/nAtom,
	     (int)(((double)time_q)/((double)clk_tck)), nTraces); 
      
      for(int k=0; k<nBestTrace; k++)
	printf("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
	       bestTraceLen[k]);
      printf("\n");
    }
    
  }
  else {
    if(isPrint) {
      times(&time2);
      long clk_tck=sysconf(_SC_CLK_TCK);
      long time_q=(time2.tms_stime-time1.tms_stime+
		   time2.tms_utime-time1.tms_utime);
      printf("size=0 time=%d comb=%d\n", 
	     (int)(((double)time_q)/((double)clk_tck)), nTraces);
    }
  }

  
  delete [] traceLen;
  delete [] bestTraceLen;
  delete [] a;

  DeleteArray(&mat, nse1);
  DeleteArray(&dist1, nse1);
  DeleteArray(&dist2, nse2);
  delete [] f1;
  delete [] f2;
  delete [] trace1;
  delete [] trace2;
  delete [] bestTrace1;
  delete [] bestTrace2;
  delete [] traceIterLevel;
  delete [] traceIndex;
  DeleteArray(&traceScore, traceMaxSize);
  DeleteArray(&traceUsage, traceMaxSize);
  delete [] strBuf1;
  delete [] strBuf2;

  delete [] bestTracesN;
  delete [] bestTracesScores;
  DeleteArray(&bestTraces1, bestTracesMax);
  DeleteArray(&bestTraces2, bestTracesMax);
  return(z);
  */
}
///////////////////////////////////////////////////////////////////////////
double zStrAlign(int4 *align_se, int lcmp, XYZ *ca1, XYZ *ca2) {
		
  double d=0.0; 
  int nd=0, lali=0, lgap=0;
  for(int l1=0; l1<lcmp; l1++) 
    if(align_se[l1]!=-1 && align_se[l1+lcmp]!=-1) {
      if(ca1[align_se[l1]].X<1e10 && ca2[align_se[l1+lcmp]].X<1e10) {
	lali++;
      }
    }
    else {
      lgap++;
    }

  for(int l1=0; l1<lcmp-2; l1++) 
    if(align_se[l1]!=-1 && align_se[l1+lcmp]!=-1) {

      if(ca1[align_se[l1]].X<1e10 && ca2[align_se[l1+lcmp]].X<1e10) 
	
	for(int l2=l1+1; l2<lcmp; l2++) 
	  if(align_se[l2]!=-1 && align_se[l2+lcmp]!=-1) 
	    if(abs(align_se[l1]-align_se[l2])>1) 
	      if(ca1[align_se[l2]].X<1e10 && ca2[align_se[l2+lcmp]].X<1e10) {
		d+=fabs(ca1[align_se[l1]].dist(ca1[align_se[l2]])-
			ca2[align_se[l1+lcmp]].dist(ca2[align_se[l2+lcmp]]));
		nd++;
	      }
    }
  if(lali>4) {
    return(zStrAlign(8, lali/8+1, d/nd, lgap));
  }
  return(0.0);
}
////////////////////////////////////////////////////////////////////
