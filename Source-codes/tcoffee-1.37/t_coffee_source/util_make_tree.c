#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"


static int nseqs;
static char **names;       /* the seq. names */

static double   **tmat;
static double 	*av;
static double 	*left_branch, *right_branch;
static int 	*tkill;






NT_node ** make_nj_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode)
        {
	  int a,b;
	  double **d;
	  NT_node **T;
	  
	  d=declare_double( out_nseq, out_nseq);

	  for ( a=0; a< out_nseq; a++)
	    for ( b=0; b< out_nseq; b++)
	      d[a][b]=distances[a][b];
	  T=dist2nj_tree ( d, out_seq_name, out_nseq, tree_file);
	  free_double (d, -1);
	  return T;
	}

NT_node ** int_dist2nj_tree (int **distances, char **out_seq_name, int out_nseq,  char *tree_file) 
        {
	  int a, b;
	  double **d;
	  NT_node **T;

	  d=declare_double( out_nseq, out_nseq);
	  for ( a=0; a< out_nseq; a++)
	    for ( b=0; b< out_nseq; b++)
	      d[a][b]=distances[a][b];
	  T=dist2nj_tree ( d, out_seq_name, out_nseq, tree_file);
	  free_double (d, -1);
	  return T;
	}
NT_node ** float_dist2nj_tree (float **distances, char **out_seq_name, int out_nseq,  char *tree_file)
        {
	  int a, b;
	  double **d;
	  NT_node **T;
	  
	  d=declare_double( out_nseq, out_nseq);
	  for ( a=0; a< out_nseq; a++)
	    for ( b=0; b< out_nseq; b++)
	      d[a][b]=distances[a][b];
	  T=dist2nj_tree ( d, out_seq_name, out_nseq, tree_file);
	  free_double (d, -1);
	  return T;
	}
NT_node ** dist2nj_tree (double **distances, char **out_seq_name, int out_nseq,  char *tree_file)
	{
	int max, min, a, b;
	double **d_dist;
	int tot_node=0;
	
	if ( !tree_file)tree_file=vtmpnam(NULL);

	d_dist=declare_double( out_nseq+1, out_nseq+1);
	
	min=max=distances[0][1];
	for ( a=0; a<out_nseq-1; a++)
		for ( b=a+1; b< out_nseq; b++)
			{
			max=(distances[a][b]>max)?distances[a][b]:max;
			min=(distances[a][b]<min)?distances[a][b]:min;
			}
	
	for ( a=0; a<out_nseq; a++)
		for ( b=0; b< out_nseq; b++)
			{
			if ( a!=b)
				d_dist[a+1][b+1]=(double)(max-distances[a][b]+min)/(double)max;
			else d_dist[a+1][b+1]=0;
			}
	guide_tree ( tree_file, d_dist, out_seq_name, out_nseq);
	free_double (d_dist,-1);
	return read_tree (tree_file,&tot_node,out_nseq,  out_seq_name);	
	}

void nj_tree(char **tree_description)
{
	register int i;
	int l[4],nude,k;
	int nc,mini,minj,j,ii,jj;
	double fnseqs,fnseqs2,sumd;
	double diq,djq,dij,d2r,dr,dio,djo,da;
	double tmin,total,dmin;
	double bi,bj,b1,b2,b3,branch[4];
	int typei,typej;             /* 0 = node; 1 = OTU */
	
	fnseqs = (double)nseqs;

/*********************** First initialisation ***************************/
	
	

	mini = minj = 0;

	left_branch 	= (double *) malloc( (nseqs+2) * sizeof (double)   );
	right_branch    = (double *) malloc( (nseqs+2) * sizeof (double)   );
	tkill 		= (int *) malloc( (nseqs+1) * sizeof (int) );
	av   		= (double *) malloc( (nseqs+1) * sizeof (double)   );

	for(i=1;i<=nseqs;++i) 
		{
		tmat[i][i] = av[i] = 0.0;
		tkill[i] = 0;
		}

/*********************** Enter The Main Cycle ***************************/

       	/**start main cycle**/
	for(nc=1; nc<=(nseqs-3); ++nc) {
		sumd = 0.0;
		for(j=2; j<=nseqs; ++j)
			for(i=1; i<j; ++i) {
				tmat[j][i] = tmat[i][j];
				sumd = sumd + tmat[i][j];
			}

		tmin = 99999.0;

/*.................compute SMATij values and find the smallest one ........*/

		for(jj=2; jj<=nseqs; ++jj) 
			if(tkill[jj] != 1) 
				for(ii=1; ii<jj; ++ii)
					if(tkill[ii] != 1) {
						diq = djq = 0.0;

						for(i=1; i<=nseqs; ++i) {
							diq = diq + tmat[i][ii];
							djq = djq + tmat[i][jj];
						}

						dij = tmat[ii][jj];
						d2r = diq + djq - (2.0*dij);
						dr  = sumd - dij -d2r;
						fnseqs2 = fnseqs - 2.0;
					        total= d2r+ fnseqs2*dij +dr*2.0;
						total= total / (2.0*fnseqs2);

						if(total < tmin) {
							tmin = total;
							mini = ii;
							minj = jj;
						}
					}
		

/*.................compute branch lengths and print the results ........*/


		dio = djo = 0.0;
		for(i=1; i<=nseqs; ++i) {
			dio = dio + tmat[i][mini];
			djo = djo + tmat[i][minj];
		}

		dmin = tmat[mini][minj];
		dio = (dio - dmin) / fnseqs2;
		djo = (djo - dmin) / fnseqs2;
		bi = (dmin + dio - djo) * 0.5;
		bj = dmin - bi;
		bi = bi - av[mini];
		bj = bj - av[minj];

		if( av[mini] > 0.0 )
			typei = 0;
		else
			typei = 1;
		if( av[minj] > 0.0 )
			typej = 0;
		else
			typej = 1;

		

/* 
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/		if( fabs(bi) < 0.0001) bi = 0.0;
		if( fabs(bj) < 0.0001) bj = 0.0;

	    	


	    	left_branch[nc] = bi;
	    	right_branch[nc] = bj;

		for(i=1; i<=nseqs; i++)
			tree_description[nc][i] = 0;

	     	if(typei == 0) { 
			for(i=nc-1; i>=1; i--)
				if(tree_description[i][mini] == 1) {
					for(j=1; j<=nseqs; j++)  
					     if(tree_description[i][j] == 1)
						    tree_description[nc][j] = 1;
					break;
				}
		}
		else
			tree_description[nc][mini] = 1;

		if(typej == 0) {
			for(i=nc-1; i>=1; i--) 
				if(tree_description[i][minj] == 1) {
					for(j=1; j<=nseqs; j++)  
					     if(tree_description[i][j] == 1)
						    tree_description[nc][j] = 1;
					break;
				}
		}
		else
			tree_description[nc][minj] = 1;
			

/* 
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
                if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;

/*........................Re-initialisation................................*/

		fnseqs = fnseqs - 1.0;
		tkill[minj] = 1;

		for(j=1; j<=nseqs; ++j) 
			if( tkill[j] != 1 ) {
				da = ( tmat[mini][j] + tmat[minj][j] ) * 0.5;
				if( (mini - j) < 0 ) 
					tmat[mini][j] = da;
				if( (mini - j) > 0)
					tmat[j][mini] = da;
			}

		for(j=1; j<=nseqs; ++j)
			tmat[minj][j] = tmat[j][minj] = 0.0;


/****/	}						/**end main cycle**/

/******************************Last Cycle (3 Seqs. left)********************/

	nude = 1;

	for(i=1; i<=nseqs; ++i)
		if( tkill[i] != 1 ) {
			l[nude] = i;
			nude = nude + 1;
		}

	b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
	b2 =  tmat[l[1]][l[2]] - b1;
	b3 =  tmat[l[1]][l[3]] - b1;
 
	branch[1] = b1 - av[l[1]];
	branch[2] = b2 - av[l[2]];
	branch[3] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
	if( fabs(branch[3]) < 0.0001) branch[3] = 0.0;

	left_branch[nseqs-2] = branch[1];
	left_branch[nseqs-1] = branch[2];
	left_branch[nseqs]   = branch[3];

	for(i=1; i<=nseqs; i++)
		tree_description[nseqs-2][i] = 0;

	

	for(i=1; i<=3; ++i) {
	   if( av[l[i]] > 0.0) {
	
		for(k=nseqs-3; k>=1; k--)
			if(tree_description[k][l[i]] == 1) {
				for(j=1; j<=nseqs; j++)
				 	if(tree_description[k][j] == 1)
					    tree_description[nseqs-2][j] = i;
				break;
			}
	   }
	   else  {
	      
		tree_description[nseqs-2][l[i]] = i;
	   }
	   if(i < 3) {
	   }
	}
}

void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap)
{

	fprintf(tree,"(\n");
 
	two_way_split(tree_description, tree, nseqs-2,1,bootstrap);
	fprintf(tree,":%7.5f,\n",left_branch[nseqs-2]);
	two_way_split(tree_description, tree, nseqs-2,2,bootstrap);
	fprintf(tree,":%7.5f,\n",left_branch[nseqs-1]);
	two_way_split(tree_description, tree, nseqs-2,3,bootstrap);

	fprintf(tree,":%7.5f)",left_branch[nseqs]);
	if (bootstrap) fprintf(tree,"TRICHOTOMY");
	fprintf(tree,";\n");
}


void two_way_split
(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap)
{
	int row, new_row, col, test_col;
	int single_seq;

	if(start_row != nseqs-2) fprintf(tree,"(\n"); 

	for(col=1; col<=nseqs; col++) {
		if(tree_description[start_row][col] == flag) {
			test_col = col;
			break;
		}
	}

	single_seq = TRUE;
	for(row=start_row-1; row>=1; row--) 
		if(tree_description[row][test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		tree_description[start_row][test_col] = 0;
		fprintf(tree,"%s",names[test_col+0-1]);
	}
	else {
		for(col=1; col<=nseqs; col++) {
		    if((tree_description[start_row][col]==1)&&
		       (tree_description[new_row][col]==1))
				tree_description[start_row][col] = 0;
		}
		two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
	}

	if(start_row == nseqs-2) {
/*		if (bootstrap && (flag==1)) fprintf(tree,"[TRICHOTOMY]");
*/
		return;
	}

	fprintf(tree,":%7.5f,\n",left_branch[start_row]);

	for(col=1; col<=nseqs; col++) 
		if(tree_description[start_row][col] == flag) {
			test_col = col;
			break;
		}
	
	single_seq = TRUE;
	for(row=start_row-1; row>=1; row--) 
		if(tree_description[row][test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		tree_description[start_row][test_col] = 0;
		fprintf(tree,"%s",names[test_col+0-1]);
	}
	else {
		for(col=1; col<=nseqs; col++) {
		    if((tree_description[start_row][col]==1)&&
		       (tree_description[new_row][col]==1))
				tree_description[start_row][col] = 0;
		}
		two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
	}

	fprintf(tree,":%7.5f)\n",right_branch[start_row]);
	
	
}


void guide_tree(char *fname, double **saga_tmat, char **saga_seq_name, int saga_nseq)
/* 
   Routine for producing unrooted NJ trees from seperately aligned
   pairwise distances.  This produces the GUIDE DENDROGRAMS in
   PHYLIP format.
*/
{    

        int i;
        FILE *fp;
        static char **standard_tree;



  	nseqs=saga_nseq;
  	names=saga_seq_name;
  	tmat=saga_tmat;
        standard_tree   = (char **) malloc( (nseqs+1) * sizeof (char *) );
        for(i=0; i<nseqs+1; i++)
                standard_tree[i]  = (char *) malloc( (nseqs+1) * sizeof(char));
         	
        	
        nj_tree(standard_tree);
        
        fp=fopen ( fname, "w");
        print_phylip_tree(standard_tree,fp,FALSE);

        free(left_branch);
        free(right_branch);
        free(tkill);
        free(av);
        for (i=1;i<nseqs+1;i++)
                free(standard_tree[i]);
        free(standard_tree);

        fclose(fp);

        
}

