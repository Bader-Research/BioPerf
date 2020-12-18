/*	dispn.c	associated subroutines for matching sequences */

#include <stdio.h>

#define XTERNAL
#include "upam.gbl"

#define TRUE 1
#define FALSE 0
extern FILE *outfd;
extern dnaseq;
extern int smark[4];
extern int min0,min1,max0,max1;
extern int smin0, smin1, smins;
extern long loffset;
extern long sq0off, sq1off;
extern char name0[],name1[];
extern int llen, markx, showall;

#define MAXOUT 201

discons(seqc0, seqc1, nc)
     char *seqc0, *seqc1;
     int nc;
{
  char line[3][MAXOUT], cline[2][MAXOUT+10];
  int il, i, lend, loff, il1, il2;
  int nm;
  int del0, del1, ic, ll0, ll1, ll01, cl0, cl1, rl0, rl1;
  int i00, i0n, i10, i1n;
  int ic0, ic1;
  int ioff0, ioff1;
  long qqoff, lloff;
  int tmp;
  char *name01;

  if (markx==2) name01=name1; else name01 = "\0";

  i00 = smark[0];
  i0n = smark[1];
  i10 = smark[2];
  i1n = smark[3];
  
  ioff0=smin0-smins;
  ioff1=smin1-smins;
  
  nm=0;
  
  
  fprintf ( stderr, "**%d %d %d %d", min0,max0, min1, max1);
  if (markx==4) return;

  if (markx==3) {
    fprintf(outfd,">%s ..\n",name0);
    for (i=0; i<nc; i++) {
      if (seqc0[i]=='-') fputc('.',outfd);
      else fputc(seqc0[i],outfd);
      if (i%50 == 49) fputc('\n',outfd);
    }
    fputc('\n',outfd);
    fprintf(outfd,">%s ..\n",name1);
    for (i=0; i<nc; i++) {
      if (seqc1[i]=='-') fputc('.',outfd);
      else fputc(seqc1[i],outfd);
      if (i%50 == 49) fputc('\n',outfd);
    }
    fputc('\n',outfd);
    return;
  }

  for (i=0; i<3; i++) memset(line[i],' ',MAXOUT);

  ic = 0; del0=del1=0;
  for (il=0; il<(nc+llen-1)/llen; il++) {
    loff=il*llen;
    lend=min(llen,nc-loff);

    ll0 = FALSE; ll1 = FALSE;
    ic0 = ic-del0;	/* save for labeling line */
    ic1 = ic-del1;


    
    for (i=0; i<2; i++) memset(cline[i],' ',MAXOUT);
    for (i=0; i<lend; i++, ic++,ioff0++,ioff1++) {
      cl0 =  cl1 = rl0 = rl1 = TRUE;
      if ((line[0][i]=seqc0[ic])=='-') {
	del0++; cl0=rl0=FALSE;
      }
      if ((line[2][i]=seqc1[ic])=='-') {
	del1++; cl1=rl1=FALSE;
      }
		
      if (seqc0[ic]==' ') {del0++; cl0=rl0=FALSE;}
      else ll0 = TRUE;
      if (seqc1[ic]==' ') {del1++; cl1=rl1=FALSE;}
      else ll1 = TRUE;
      qqoff = sq0off - 1 + (long)(ioff0-del0);
      if (cl0 && qqoff%10 == 9)  {
	sprintf(&cline[0][i],"%8ld",qqoff+1l);
	cline[0][i+8]=' ';
	rl0 = FALSE;
      }
      else if (cl0 && qqoff== -1) {
	sprintf(&cline[0][i],"%8ld",0l);
	cline[0][i+8]=' ';
	rl0 = FALSE;
      }
      else if (rl0 && (qqoff+1)%10 == 0) {
	sprintf(&cline[0][i],"%8ld",qqoff+1);
	cline[0][i+8]=' ';
      }
      
      lloff = sq1off-1 + loffset + (long)(ioff1-del1);
      if (cl1 && lloff%10 == 9)  {
	sprintf(&cline[1][i],"%8ld",lloff+1l);
	cline[1][i+8]=' ';
	rl1 = FALSE;
      }
      else if (cl1 && lloff== -1) {
	sprintf(&cline[1][i],"%8ld",0l);
	cline[1][i+8]=' ';
	rl1 = FALSE;
      }
      else if (rl1 && (lloff+1)%10 == 0) {
	sprintf(&cline[1][i],"%8ld",lloff+1);
	cline[1][i+8]=' ';
      }
      

      line[1][i] = ' ';
      if (ioff0-del0 >= min0 && ioff0-del0 <= max0) {
	if (toupper(line[0][i])==toupper(line[2][i]) || (dnaseq && (
	    (toupper(line[0][i])=='T' && toupper(line[2][i])=='U') ||
	    (toupper(line[0][i])=='U' && toupper(line[2][i])=='T'))))
	  switch (markx) {
	  case 0: line[1][i]= '|';
	    break;
	  case 1: line[1][i]= ' ';
	    break;
	   case 2: line[1][i]= '.';
	    break;
	  }
	else if (markx==2) line[1][i]=line[2][i];
	else if ((il1=iidex(sq,line[0][i]))>=0 &&
		 (il2=iidex(sq,line[2][i]))>=0 &&
		 pam2[il1][il2]>= 0)
	    line[1][i]= (markx) ? 'x':'.';
	else if ((il1=iidex(sq,line[0][i]))>=0 &&
		 (il2=iidex(sq,line[2][i]))>=0)
	    line[1][i]= (markx) ? 'X':' ';
      }
      else if (markx==2) line[1][i]=line[2][i];

      if (markx==0) {
	if (ioff0-del0 == i00 && ioff1-del1 == i10) {
	  line[1][i]='X';
	  i00 = i10 = -1;
	}
	if (ioff0-del0 == i0n && ioff1-del1 == i1n) {
	  line[1][i]='X';
	  i0n = i1n = -1;
	}
	if ((ioff0-del0 == i00) || (ioff0-del0 == i0n)) {
	  line[1][i]='^';
	  if(ioff0-del0 == i00) i00= -1;
	  else i0n = -1;
	}
	if (ioff1-del1 == i10 || ioff1-del1 == i1n) {
	  line[1][i]='v';
	  if(ioff1-del1 == i10) i10= -1;
	  else i1n = -1;
	}
      }
    }
    
    for (i=0; i<3; i++) {line[i][lend]=0;}
    for (i=0; i<2; i++) {cline[i][lend+7]=0;}
    
    ll01 = ll0&&ll1;
    if (markx==2 && (!showall || ll0)) ll1=0;
    fprintf(outfd,"\n");
    if (ll0) fprintf(outfd,"%s\n",cline[0]);
    if (ll0) fprintf(outfd,"%-6s %s\n",name0,line[0]);
    if (ll01) fprintf(outfd,"%-6s %s\n",name01,line[1]);
    if (ll1) fprintf(outfd,"%-6s %s\n",name1,line[2]);
    if (ll1) fprintf(outfd,"%s\n",cline[1]);
  }
}

static float gscale= -1.0;

disgraph(n0,n1,percent,score,min0,min1,max0,max1)
     int n0, n1, min0, min1, max0, max1;
     float percent;
     int score;
{
  int i, gstart, gstop, gend;
  char line[MAXOUT+1];

  memset(line,' ',llen);

  line[llen-1]='\0';
  if (gscale < 0.0) {
    gscale = (float)llen/(float)n0;
    fprintf(outfd,"%-6s %4d-%4d:     %5.1f%%:%s:\n",name0,1,n0,100.0,line);
  }

  gstart = (int)(gscale*(float)min0+0.5);
  gstop = (int)(gscale*(float)max0+0.5);
  gend = gstop+(int)(gscale*(float)(n1-max1));

  if (gstop >= llen) gstop = llen-1;
  if (gend >= llen) gend = llen-1;
  for (i=0; i<gstart; i++) line[i]=' ';
  for (; i<gstop; i++) line[i]='-';
  for (; i<llen; i++) line[i]=' ';

  line[gend]=':';
  line[llen]='\0';

  fprintf(outfd,"%-6s %4d-%4d:%4d %5.1f%%:%s\n",
	  name1,min0+1,max0,score,percent,line);
}

aancpy(to,from,count)
     char *to, *from;
     int count;
{
  char *tp;
  tp=to;
  while (count--&& *from>=0) {
    if (*from<nsq) *tp++ = sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp=0;
}

r_memcpy(dest,src,cnt)
     char *dest, *src;
     int cnt;
{
  while (cnt--) *dest++ = *src++;
}

l_memcpy(dest,src,cnt)
     char *dest, *src;
     int cnt;
{
  dest = dest+cnt;
  src = src+cnt;
  while (cnt--) *--dest = *--src;
}

iidex(str, chr)
	char *str, chr;
{
	int i;
	for (i=0; str[i]; i++) if (str[i]==chr) return i;
	return (-1);
	}

opnline() {}

clsline() {}



min(arg1, arg2)
	int arg1, arg2;
{
	return (arg1<=arg2) ? arg1 : arg2;
	}

#ifdef VMS
memset(str, c, cnt)
	char *str; int cnt; char c;
{
	while (cnt--) *str++ = c;
	}
#endif
