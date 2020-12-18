#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
static float dm[]={
0.178091,
1.18065, 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 0.015315, 0.048298, 0.053803, 0.020662, 0.023612, 0.216147, 0.147226, 0.065438, 0.003758, 0.009621,
0.056591,
1.35583, 0.021465, 0.0103, 0.011741, 0.010883, 0.385651, 0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 0.029156, 0.018153, 0.0361, 0.07177, 0.419641,
0.0960191,
6.66436 ,0.561459, 0.045448, 0.438366, 0.764167, 0.087364, 0.259114, 0.21494, 0.145928, 0.762204, 0.24732, 0.118662, 0.441564, 0.174822, 0.53084, 0.465529, 0.583402, 0.445586, 0.22705, 0.02951, 0.12109,
0.078123,
2.08141, 0.070143, 0.01114, 0.019479, 0.094657, 0.013162, 0.048038, 0.077, 0.032939, 0.576639, 0.072293, 0.02824, 0.080372, 0.037661, 0.185037, 0.506783, 0.073732, 0.071587, 0.042532, 0.011254, 0.028723,
0.0834977,
2.08101, 0.041103, 0.014794, 0.00561, 0.010216, 0.153602, 0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 0.012049, 0.035799, 0.180085, 0.012744, 0.026466,
0.0904123,
2.56819, 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 0.027535, 0.088333, 0.94434, 0.004373, 0.016741,
0.114468,
1.76606, 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 0.01004, 0.187869, 0.050018, 0.110039, 0.038668, 0.119471, 0.065802, 0.02543, 0.003215, 0.018742,
0.0682132,
4.98768, 0.452171, 0.114613, 0.06246, 0.115702, 0.284246, 0.140204, 0.100358, 0.55023, 0.143995, 0.700649, 0.27658, 0.118569, 0.09747, 0.126673, 0.143634, 0.278983, 0.358482, 0.66175, 0.061533, 0.199373,
0.234585,
0.0995, 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 0.016931, 0.003647, 0.002184, 0.005019, 0.00599, 0.001473, 0.004158, 0.009055, 0.00363, 0.006583, 0.003172, 0.00369, 0.002967, 0.002772,0.002686};

double int_logB (int *i, int n)
	{
	static double *array;
	int a;
	
	if ( array==NULL)array=vcalloc ( 1000, sizeof (double));
	
	for ( a=0; a< n; a++)
		array[a]=(double)i[a];
	return double_logB(array, n);
	}
double float_logB (float *i, int n)
	{
	static double *array;
	int a;
	 
	if ( array==NULL)array=vcalloc ( 1000, sizeof (double));
	for ( a=0; a< n; a++)
		array[a]=(double)i[a];
	return double_logB(array, n);
	}
	  
double double_logB (double *x, int n)
	{
	double vx=0;
	double result=0;
	int i;
	
	
	for ( i=0; i<n; i++)vx+=x[i];
	for ( i=0; i<n; i++)result+=lgamma(x[i]);
	return 	result-lgamma(vx);
	} 	
double *** make_lup_table ( Mixture *D)
	{
	int a, b, c;
	double ***lup;
	
	lup=vcalloc ( 9, sizeof (double**));
	for ( a=0; a<9; a++)
		{
		lup[a]=vcalloc ( 20, sizeof (double*));
		for ( b=0; b< 20; b++)
			lup[a][b]=vcalloc ( 200, sizeof (double));
		}
	
	for ( a=0; a< 9; a++)
		for ( b=0; b< 20; b++)
			for ( c=0; c< 100; c++)
				lup[a][b][c]=lgamma(D->ALPHA[a][b]+c);
	
	return lup;
	}
	
double  double_logB2(int j, double *n,Mixture *D)
	{
	double vx=0;
	double result=0;
	int i;
	
	static double ***lup;
	
	
	
	if ( lup==NULL)lup=make_lup_table (D);
	
       

	for ( i=0; i<D->n_aa; i++)vx+=(double)n[i]+D->ALPHA[j][i];
	for ( i=0; i<D->n_aa; i++)
	  {
	    
	    
	    result+=lup[j][i][(int)n[i]];
	  }
	return 	result-lgamma(vx);
	} 	
			
double compute_exponant ( double *n, int j, Mixture *D)
	{
	
	if ( j>=9)fprintf ( stderr, "\nPB: j=%d", j);
	
	return double_logB2(j, n,D)-D->double_logB_alpha[j];
	}


double *compute_matrix_p ( double *n,int Nseq)
        {
	  
	  /*
	    reads in a frquency list of various amino acids:
	    
	    sum freq(aa)=1 (gaps are ignored)
	    
	    aa[1]=x1
	    aa[2]=x2
	    ....

	    Outputs a similar list with frequencies 'Blurred' using a pam250 mt
	  */

	    
	  
	  static int **matrix;
	  static double *R;
	  int a, b;
	  double v,min, tot;
	  
	  
	  if ( !matrix) 
	    {
	      matrix=read_matrice ( "pam250mt");	  
	    }
	  
	  if (R){free(R); R=vcalloc ( 26, sizeof (double));}
	  else R=vcalloc ( 26, sizeof (double));

	  for ( a=0; a<26; a++)
	    {
	      if (!is_aa(a+'a'))continue;
	      if ( n[a]==0)continue;
	      
	      for ( b=0; b< 26; b++)
		{
		  if (!is_aa(b+'a'))continue;
		  v=n[a]*(matrix[a][b]);
		  if ( v>0)
		    {
		      R[b]+=v+(10*n[a]);
		    }
		}
	    }
	  
	  min=R[0];
	  for ( min=R[0],a=0; a< 26; a++)min=MIN(min,R[a]);
	  for ( tot=0,   a=0; a< 26; a++)         {R[a]-=min;tot+=R[a];}
	  for ( a=0; a< 26; a++)if ( is_aa(a+'a')){R[a]=R[a]*((float)(100)/(float)tot);}
	  return R;
	}
	      

double *compute_dirichlet_p ( double *n,int Nseq)
	{
	  /*
	    Given a list of frequenceies measured for the residues, this function returns 
	    the p_values associated with each residue in the column
	  */
	  
	int a, b;
	double X_LIST[100];
	double sum, log_sum, max;
	static Mixture *D;
	static double *R;



	if (!D)
		{
		D=read_dirichlet (NULL);
		
		D->n_aa=20;
		R=vcalloc ( D->n_aa, sizeof (double));
		D->double_logB_alpha=vcalloc (D->N_COMPONENT , sizeof (double));
		
		D->exponant_list=vcalloc (D->N_COMPONENT , sizeof (double));
		precompute_log_B ( D->double_logB_alpha,D);
		D->alpha_tot=vcalloc (D->N_COMPONENT , sizeof (double));
		for ( a=0; a<D->N_COMPONENT; a++)
			for ( b=0; b< D->n_aa; b++)
				D->alpha_tot[a]+=D->ALPHA[a][b];
		}
	
	

	for ( D->tot_n=0,a=0; a< D->n_aa; a++)D->tot_n+=(double)n[a];
	max=D->exponant_list[0]=compute_exponant ( n, 0, D);	
	for ( a=1; a<D->N_COMPONENT; a++)
		{
		D->exponant_list[a]=compute_exponant ( n, a,D);
		max= ( max< D->exponant_list[a])?D->exponant_list[a]:max;
		}
	for ( a=1; a<D->N_COMPONENT; a++)D->exponant_list[a]=D->exponant_list[a]-max;
	
	
	for ( sum=0,log_sum=0,a=0; a< D->n_aa; a++)
		{
		sum+=X_LIST[a]=compute_X (n, a,D);
		}
	log_sum=log(sum);

		
	for (a=0; a<D->n_aa; a++)
		{
		R[a]=(log(X_LIST[a])-log_sum);
		}
	
	
	/*
	printf ( "\n[");
	for ( a=0;a< n_aa; a++)printf ("%d ", n[a]);
	printf ("] score=%f",(float) result );
	
	fprintf ( stderr, "\nRESULT=%f", (float)result);
	exit(0);
	*/
	return R;

	}
	
void precompute_log_B ( double *table,Mixture *D)
	{
	int a;
	for ( a=0; a< D->N_COMPONENT; a++)
		{
		table[a]=double_logB ( D->ALPHA[a], D->n_aa);
		}
	}			
double compute_X (double *n,int i,Mixture *D)
	{
	int  j;
	double term1, term2,result;
	double **alpha;
	double *q;

	
	
	alpha=D->ALPHA;
	q=D->DM_Q;

	for (result=0, j=0; j<D->N_COMPONENT; j++)
		{
		term1=exp (D->exponant_list[j])*q[j];
		term2=(alpha[j][i]+(double)n[i])/(D->alpha_tot[j]+D->tot_n);
		result+=term1*term2;
		}
	return result;
	}
Mixture * read_dirichlet ( char *name)
	{
	FILE *fp;
	int a,b, c;
	float f;	
	Mixture *D;


	D=vcalloc ( 1, sizeof (Mixture));
	
	
	D->N_COMPONENT=9;
	D->ALPHA=vcalloc (9, sizeof (double*));
	for ( a=0; a< 9; a++)
		D->ALPHA[a]=vcalloc (20, sizeof (double));
	D->DM_Q=vcalloc (9, sizeof (double));
	
	if (name!=NULL)
	  {
	    fp=vfopen ( name, "r");
	    for ( a=0; a< 9; a++)
		{
		fscanf(fp, "%f\n", &f);
		D->DM_Q[a]=(double)f;
		fscanf(fp, "%f", &f);
		
		for ( b=0; b<20; b++)
			{
			fscanf(fp, "%f", &f);
			D->ALPHA[a][b]=(double)f;
			}
		fscanf(fp, "\n");
		}
	    for ( a=0; a< 9; a++)
	      {
		fprintf(stderr, "\n%f\n",(float)D->DM_Q[a] );
		
		for ( b=0; b<20; b++)
		  {
		    fprintf(stderr, "%f ", (float)D->ALPHA[a][b]);
		  }
		fprintf(stderr, "\n");
	      }
	    fprintf ( stderr, "\nN_C=%d",D->N_COMPONENT );	
	    vfclose ( fp);
	  }
	else
	  {
	    for (c=0, a=0; a< 9;a++)
	      {
		D->DM_Q[a]=dm[c++];	
		for (b=0; b<20; b++)
		  D->ALPHA[a][b]=dm[c++];
	      }
	  }
	
	return D;
	}
int dirichlet_code( char aa)
	{
	
	char x;
	
	x=tolower (aa);
	
	if ( (x<'a') || (x>'z'))
		crash ( "CODE UNDEFINED");
	else if ( x<='a')
	    return x-'a';
	else if ( x<='i')
	    return x-('a'+1);
	else if ( x<= 'n')
	    return x-('a'+2);
	else if ( x<='t')
	    return x-('a'+3);
	else if ( x<='w')
	    return x-('a'+4);
	else if ( x=='y')
	    return x-('a'+5);
	else 
	  {
	    crash ("ERROR in dirichlet_code");
	    return 0;
	  }
	return 0;
	
	}
