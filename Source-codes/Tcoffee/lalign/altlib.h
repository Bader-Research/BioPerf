
#ifdef BIGMEM
#define LASTLIB 12
#define NCBIBL13 11
#else
#define LASTLIB 11
#endif

#define NCBISRC 10
#define BINARYGB 9
#define DEFAULT 0
#define FULLGB 1
#define UNIXPIR 2
#define EMBLSWISS 3
#define INTELLIG 4
#define VMSPIR 5
#define LASTTXT 5

extern int (*getlib)(), (*ranlib)();

int agetlib(),aranlib();	/* pearson fasta format */
int ggetlib(),granlib();	/* compressed genbank format BINARYGB */
int lgetlib(),lranlib();	/* full uncompressed GB FULLGB*/
int pgetlib(),pranlib();	/* PIR UNIX protein UNIXPIR */
int egetlib(),eranlib();	/* EMBL/SWISS-PROT EMBLSWISS */
int igetlib(),iranlib();	/* Intelligenetics INTELLIG */
int vgetlib(),vranlib();	/* PIR VMS format */
extern int ncbi_getliba(), ncbi_ranlib(); /* ncbi search-cd format */

#ifdef BIGMEM
extern int ncbl_getliba(), ncbl_ranlib(); /* ncbi blast 1.3 format */
#endif

int (*getliba[LASTLIB])()={
	agetlib,lgetlib,pgetlib,egetlib,
	igetlib,vgetlib,agetlib,agetlib,
	agetlib,ggetlib,ncbi_getliba
#ifdef BIGMEM
	,ncbl_getliba
#endif
        };

int (*ranliba[LASTLIB])()={
	aranlib,lranlib,pranlib,eranlib,
	iranlib,vranlib,aranlib,aranlib,
	aranlib,granlib,ncbi_ranlib
#ifdef BIGMEM
	,ncbl_ranlib
#endif
        };

