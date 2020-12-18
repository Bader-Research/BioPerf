/*	pam.c	19-June-86

	copyright (c) 1987 William R. Pearson

	read in the alphabet and pam matrix data

	designed for universal matcher
*/

#include <stdio.h>

#define XTERNAL
#include "uascii.gbl"
#include "upam.gbl"

extern int bestoff, bestscale, bkfact, scfact, bktup, bestmax, histint;

initpam(mfname)
	char *mfname;
{
	char line[512], *lp;
	int ggaptmp, gdeltmp;
	int i, iaa, ipam;
	FILE *fmat;

	if (strcmp(mfname,"250")==0) {
	  pam = apam250;
	  strcpy(mfname,"PAM250");
	  bestscale = 200;
	  bestoff = 27;
	  bkfact = 5;
	  scfact = 4;
	  return 1;
	}

	if (strcmp(mfname,"BL50")==0) {
	  pam = abl50;
	  strcpy(mfname,"BLOSUM50");
	  return 1;
	}

	if ((fmat=fopen(mfname,"r"))==NULL) {
		printf(" cannot open scoring matrix file %s\n",mfname);
		return 0;
		}

l1:	if (fgets(line,512,fmat)==NULL) {
		printf(" pam - cannot read first line of SMATRIX file\n");
		return(0);
		}

	if (line[0]==';') {
		if (line[1]=='P') {
			strcpy(sqnam,"aa"); strcpy(sqtype,"protein");}
		else if (line[1]=='D') {
			strcpy(sqnam,"nt"); strcpy(sqtype,"DNA");}
		goto l1;
		}

	if (sscanf(line," %d %d %d %d %d %d %d",
	    &scfact,&bestoff,&bestscale,&bkfact,&bktup,&bestmax,&histint)!=7) {
	    	printf("  bestcut parameters - bad format\n");
		exit(1);
		}

	if (fgets(line,512,fmat)==NULL) {
		printf(" pam - cannot read DELVAL line\n");
		return 0;
		}

	else if (sscanf(line," %d %d",&gdeltmp,&ggaptmp)!=2) {
		printf(" DELVAL parameters - bad format\n");
		exit(1);
		}
		
	if (!del_set) gdelval = gdeltmp;
	if (!gap_set) ggapval = ggaptmp;

	if (fgets(line,512,fmat)==NULL) {
		printf(" pam - cannot read EOS line\n");
		return 0;
		}

/*	clear out sascii	*/
	for (i=0; i<=AAMASK; i++) sascii[i]= NA;

/*	set end of line stop	*/
	sascii[0]=sascii['\r']=sascii['\n']= EL;

/*	set end of sequence stop */
	for (i=0; line[i]; i++) if (line[i]>' ') sascii[line[i]]= ES;
	
	if (fgets(line,512,fmat)==NULL) {
		printf(" pam - cannot read aa line\n");
		exit(1);
		}

/* read the alphabet */
	for (i=0,nsq=0; line[i]; i++) if (line[i]>' ')
		sq[nsq++]=toupper(line[i]);

/* initialize sascii */
	for (iaa=0; iaa<nsq; iaa++) {
		sascii[sq[iaa]]=iaa;
		if (sascii[aa[iaa]]<NA && sq[iaa]>='A' && sq[iaa]<='Z')
			sascii[aa[iaa]-'A'+'a']=sascii[aa[iaa]];
		}

/* read in hnt values */
	for (iaa=0; iaa<nsq; iaa++)
		if (fscanf(fmat,"%d",&hsq[iaa])!=1) {
			printf(" error reading hsq values\n");
			exit(1);
			}

	for (iaa=ipam=0; iaa<nsq; iaa++)
		for (i=0; i<=iaa; i++)
			if (fscanf(fmat,"%d",&pam[ipam++])!=1) {
				printf(" error reading pam matrix\n");
				exit(1);
				}

	fclose(fmat);
	return 1;
	}
