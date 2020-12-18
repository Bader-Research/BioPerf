#define FILE_CHECK 1
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include <math.h>
#include <stdarg.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

/*********************************************************************/
/*                                                                   */
/*                                   HEAPSORT                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE * hsort_file ( FILE *fp,int n,int len, size_t size,int first_comp_field, int n_comp_fields,int (*compare)(const void *, const void*,int, int, size_t),void * (*copy)(void *,void*, size_t))
     {
     unsigned long i, ir, j, l;
     void *rra, *rrb, *ra_j, *ra_j_1;
     void *tp;
     long  start;
     int FE=1;
     
    
    
     start=ftell(fp);
     rra   =vcalloc ( len, size);
     rrb   =vcalloc ( len, size);
     ra_j  =vcalloc ( len, size);
     ra_j_1=vcalloc ( len, size);

     if ( n<2)return fp;
     l=(n >>1)+1;
     ir=n;

     for (;;)
         {
	 if ( l>FE)
	    {
	    l--;
	    fseek( fp, start+(((l-1)*len)*size), SEEK_SET);
	    fread( rra, size, len,fp); /*rra=ra[--l]*/
	    }
	 else
	    {
	    fseek( fp, start+((ir-1)*len*size), SEEK_SET);
	    fread( rra, size, len,fp); /*rra=ra[ir]*/
	    
	    fseek( fp, start, SEEK_SET);
	    fread( rrb, size, len,fp); /*rrb=ra[0]*/
	    
	    fseek( fp, start+((ir-1)*len*size), SEEK_SET);
	    fwrite(rrb,size, len, fp); /*ra[ir]=rrb=ra[0]*/ 

	    if (--ir ==FE)
	       {
	       fseek ( fp,start, SEEK_SET);
	       fwrite(rra,size, len, fp); /*ra[0]=rra*/ 
	       break;
	       }
	    }	      
	 i=l;
	 j=l+l;
	 while ( j<=ir)
	       {
	       fseek ( fp, start+((j-1)*len*size), SEEK_SET);
	       fread (ra_j, size, len, fp);
	       
	       if ( j<ir)
	          {
		  fseek ( fp, start+(((j-1)+1)*len*size), SEEK_SET);
		  fread (ra_j_1, size, len, fp);
		  }

	       if ( j<ir && compare( ra_j, ra_j_1,first_comp_field,n_comp_fields, size )<0)
		   {
		   SWAPP(ra_j, ra_j_1, tp);
		   j++;
		   }

	       if (compare(rra, ra_j, first_comp_field,n_comp_fields, size)<0)
	               {
		       fseek ( fp, start+((i-1)*len*size), SEEK_SET);
		       fwrite(ra_j,size, len, fp);
		       i=j;
		       j<<= 1;
		       }
	       else
		       j=ir+1;
	       
	       }
	 fseek ( fp, start+((i-1)*len*size), SEEK_SET);
	 fwrite(rra,size, len, fp);	
	 }
     vfree (  rra);
     vfree ( rrb);

     vfree (  ra_j);
     vfree (  ra_j_1);
     return fp;
     }
void ** hsort_array ( void **ra,int n,int len, size_t size,int first_comp_field, int n_comp_fields,int (*compare)(const void *, const void*,int, int, size_t),void * (*copy)(void *,void*, size_t))
     {
     unsigned long i, ir, j, l;
     void *rra;
     int FE=1;
     
     if ( FE==1){ra--;}
     else
     {
     n--;
     }


     rra   =vcalloc ( len, size);
    
     
     if ( n<2)return ra;
     l=(n >>1)+1;
     ir=n;

     for (;;)
         {
	 if ( l>FE)
	    {
	    copy ( rra, ra[--l],len);
	    }
	 else
	    {	   
	    copy ( rra, ra[ir],len);
	    copy ( ra[ir], ra[FE], len);
	    if (--ir ==FE)
	       {
	       copy ( ra[FE],rra,len);
	       break;
	       }
	    }
	 i=l;
	 j=l+l;
	 while ( j<=ir)
	       {	       
	       if ( j<ir && compare( ra[j], ra[j+1],first_comp_field,n_comp_fields, size )<0)j++;
	       if (compare(rra, ra[j], first_comp_field,n_comp_fields, size)<0)
	               {copy(ra[i], ra[j],len);i=j;j<<= 1;}
	       else
		       j=ir+1;	       
	       }
	 copy( ra[i], rra,len);	 
	 }
     vfree (rra);
     ra+=FE;
     
     return ra;
     }
/*********************************************************************/
/*                                                                   */
/*                         CEDRIC BSEARCH                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void * bsearch_file ( const void *key,int *p,int comp_first,int comp_len, FILE *fp,int len, int entry_len,size_t el_size, int (*compare)(const void *, const void*,int, int, size_t))
       {
       int upper, lower, c, i;
       static void *key2;
       long start;
       static long key2_size;


       start=ftell(fp);

       upper=-1;
       lower=len;
       if ( key2==NULL){key2=vcalloc (entry_len, el_size);key2_size=entry_len* el_size;}
       else if (key2_size<  (entry_len* el_size)){vfree(key2);key2=vcalloc (entry_len, el_size);key2_size=entry_len* el_size;} 

       while ((lower-upper)>1)
             {
	     i=(lower+upper) >> 1;

	     fseek ( fp,start+(i*el_size*entry_len), SEEK_SET); 
	     fread ( key2, el_size, entry_len,fp);	     
	     c=compare(key2,key, comp_first, comp_len,el_size);
	     
	     if      ( c==0){p[0]=i;return key2;}
	     else if ( c< 0)upper=i;
	     else if ( c> 0)lower=i;
             }
       return NULL;
       }

void * bsearch_array ( const void *key,int *p,int comp_first, int comp_len,void**list,int len, int entry_len,size_t el_size, int (*compare)(const void *, const void*,int, int, size_t))
       {
       int upper, lower, c, i;
       void *key2;
       
       upper=-1;
       lower=len;
       while ((lower-upper)>1)
             {
	     i=(lower+upper) >>1;
	     key2=list[i];	     
	     c=compare(key2,key, comp_first,comp_len,el_size);
	     
	     if      ( c==0){p[0]=i;return key2;}
	     else if ( c< 0)upper=i;
	     else if ( c> 0)lower=i;
             }
       return NULL;
       }	             
		     
/**********************************************************************/
/*                                                                    */
/*                         HSORT/BSEARCH WRAPPERS                             */
/*                                                                    */
/*                                                                    */
/**********************************************************************/
void **search_in_list_file ( void *key, int *p,int comp_len,FILE *fp, int len, size_t size, int entry_len)  
      {
      static void **l;	
      

      if ( l==NULL)l=vcalloc ( 1, sizeof (int*));
	
      l[0]=bsearch_file (key,p,0,comp_len,fp,len,entry_len,size,hsort_cmp);  
      if (l[0]==NULL)return NULL;
      else return l;
      }
void **search_in_list_array ( void *key,int *p, int comp_len,void **L, int len, size_t size, int entry_len)  
      {
      static void **l;
 
      if ( l==NULL)l=vcalloc ( 1, sizeof (int*));	
      
      l[0]=bsearch_array (key,p,0,comp_len,L,len,entry_len,size,hsort_cmp);  
      if (l[0]==NULL)return NULL;
      else return l;
      }
void **hsort_list_array ( void **L, int len, size_t size, int entry_len, int first_comp_field, int n_comp_fields)  
       {
       return hsort_array (L, len,entry_len, size,first_comp_field, n_comp_fields,hsort_cmp , hsort_cpy);
       }
FILE  *hsort_list_file ( FILE*fp  , int len, size_t size, int entry_len, int first_comp_field, int n_comp_fields)  
       {
       return hsort_file (fp, len,entry_len, size,first_comp_field, n_comp_fields,hsort_cmp , hsort_cpy);
       }

int hsort_cmp ( const void *a, const void *b, int first, int clen, size_t size)
       {
       int*ax;
       int*bx;
       int p;
       
       ax=(int*)a;
       bx=(int*)b;
       for ( p=first; p<clen; p++)
	   {
	   if ( ax[p]<bx[p])return -1;
	   else if ( ax[p]==bx[p]);
	   else return 1;
	   }
       return 0;
       }
void *hsort_cpy(void*to, void *from, size_t size)
       {
       
       int *ax;
       int *bx;
       int p;
       ax=(int*)to;
       bx=(int*)from;
       for (p=0; p<(int)size; p++)
	   ax[p]=bx[p];
       
       
       return to;
      
       }
       
       
void test_hsort_list_array()
      {
      int **array;
      int a;
      int n=100;

      array=declare_int(n, 3);
      for ( a=0; a<n; a++)array[a][0]=a;
      
      hsort_list_array( (void**)array,n, sizeof (int), 3, 0, 1);
      for ( a=0; a<n; a++)fprintf ( stderr, "\n%d %d", array[a][0],a);
      exit(1);
      }
      

/*********************************************************************/
/*                                                                   */
/*                         B_SEARCH_FILE FUNCTIONS                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/		 

/*********************************************************************/
/*                                                                   */
/*                         SORT/COMPARE/SEARCH FUNCTIONS             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
static int sort_field;
int **search_in_list_int ( int *key, int k_len, int **list, int ne)
	{
	int **l;		
	sort_field=k_len;	
	l=bsearch (&key,list, ne, sizeof(int**),(int(*)(const void*,const void*))(cmp_list_int));  
	return l;
	}
void sort_float ( float **V,int N_F, int F, int left, int right)
	{
	sort_field=F;
	qsort ( V, right+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	}
int cmp_float ( const float **a, const float **b)
	{
	if ( a[0][sort_field]< b[0][sort_field])return-1;
	else if ( a[0][sort_field]==b[0][sort_field])return 0;
	else return 1;
	}

void sort_int_1D ( int *L, int n)
	{
	int **array;
	int a;
	
	array=declare_int ( n, 1);
	for ( a=0; a< n; a++)
		array[a][0]=L[a];
	sort_int ( array, 1, 0, 0, n-1);
	for ( a=0; a< n; a++)
		L[a]=array[a][0];
	free_int ( array, n);
	}

void sort_int ( int **V,int N_F, int F, int left, int right)
	{
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	}
void sort_list_int ( int **V,int N_F, int F, int left, int right)
	{
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_list_int));
	}

void sort_int_inv ( int **V,int N_F, int F, int left, int right)
	{
	int a,b;
	int **list;
	
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_int));
	
	list=declare_int ((right-left)+1, N_F);
	for ( a=left; a< (right-left)+1; a++)
		{
		for ( b=0; b< N_F; b++)
			{
			list[a-left][b]=V[a][b];
			}
		}
	for ( a=left; a< (right-left)+1; a++)
		{
		for ( b=0; b< N_F; b++)
			V[a][b]=list[(right-left)-a][b];
		}
	free_int (list, -1);
	}
	
void sort_list_int_inv ( int **V,int N_F, int F, int left, int right)
	{
	int a,b;
	int **list;
	
	sort_field=F;
	qsort ( V, (right-left)+1, sizeof(int**),(int(*)(const void*,const void*))(cmp_list_int));
	
	list=declare_int ((right-left)+1, N_F);	
	for ( a=left; a< (right-left)+1; a++)
		{
		for ( b=0; b< N_F; b++)
			{
			list[a-left][b]=V[a][b];
			}
		}
	for ( a=left; a< (right-left)+1; a++)
		{
		for ( b=0; b< N_F; b++)
			V[a][b]=list[(right-left)-a][b];
		}
	free_int (list, -1);
	}	
	



int cmp_int ( const int**a, const int**b)
	{
	if ( a[0][sort_field]< b[0][sort_field])return-1;
	else if ( a[0][sort_field]==b[0][sort_field])return 0;
	else return 1;
	}
int cmp_list_int (const int**a, const int**b)
	{
	int c;
	int undef=0;
	
	for ( c=0; c<=sort_field; c++)
		{
		if ( a[0][c]==UNDEFINED|| b[0][c]==UNDEFINED)return 0;
		else if ( a[0][c]>b[0][c])return 1;
		else if ( a[0][c]<b[0][c])return -1;	
		}
	if (undef==sort_field)
		{
		return 0;
		}
	return 0;
	}

	



int name_is_in_list ( char *name, char **name_list, int n_name, int len)
	{
	int a;
	int pos=-1;
	/*Note: RETURNS THE Offset of the LAST Occurence of name in name_list*/

	
	for ( a=0; a< n_name; a++)
		if ( len!=-1)
			if (strncmp ( name, name_list[a], len)==0)pos=a;			
		else if ( strm ( name, name_list[a]))pos=a;

	return pos;
	}
char * check_list_for_dup ( char **list, int ne)
        {
	int a, b;
	
	for ( a=0; a< ne-1; a++)
	    for ( b=a+1; b< ne; b++)if (strm ( list[a], list[b]))return list[a];
	return NULL;
	}
FILE *get_number_list_in_file ( FILE *fp, int *list, int *n, int *max_len)
	{
	
	int c;
	
	while ( isspace((c=fgetc (fp))));
	ungetc(c, fp);
	while ( c!='\n')
		{
		while ( isspace((c=fgetc (fp))) && c!='\n');
		
		if ( c!='\n')
			{
			ungetc(c, fp);
			if ( n[0]>=max_len[0])
				list=realloc ( list, (n[0]+100)*sizeof (int));
				max_len[0]=(n[0]+100);
			
			fscanf ( fp, "%d",&list[n[0]++]);
			}
		}
	return fp;
	}

/*********************************************************************/
/*                                                                   */
/*                         DUPLICATION                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
short* ga_memcpy_short ( short *array1, short *array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	
	return array2;
	}
int * ga_memcpy_int ( int *array1, int *array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	
	return array2;
	}			
			
float* ga_memcpy_float ( float *array1, float *array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	
	return array2;
	}
double*  ga_memcpy_double (double *array1, double*array2, int n)
	{
	int a;
	
	for ( a=0; a< n; a++)
		array2[a]=array1[a];
	
	return array2;
	}

/*********************************************************************/
/*                                                                   */
/*                          SIZES                                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
#define WRITE_SIZE(type,function)\
void function ( int x, type *array, int os)\
     {\
     int a,l;\
     char buf[SIZE_OF_INT+1];\
     array+=os*SIZE_OF_INT;\
     for ( a=0;a<SIZE_OF_INT; a++)array[a]=0;\
     sprintf ( buf, "%d", x);\
     l=strlen (buf);\
     array+=SIZE_OF_INT-l;\
     for (a=0; a<l; a++)array[a]=(type)buf[a];\
     }
WRITE_SIZE(short,write_size_short)
WRITE_SIZE(char,write_size_char)
WRITE_SIZE(int,write_size_int)
WRITE_SIZE(float,write_size_float)
WRITE_SIZE(double,write_size_double)

#define READ_ARRAY_SIZE(type, function)\
int function (type *array, int os)\
    {\
    int a, b;\
    char buf[SIZE_OF_INT+1];\
    a=b=0;\
    array+=os*SIZE_OF_INT;\
    while ( a!=SIZE_OF_INT && array[a]==0)a++;\
    while ( a!=SIZE_OF_INT)buf[b++]=(char)array[a++];\
    buf[b]='\0';\
    return atoi(buf);\
    }
READ_ARRAY_SIZE(short,read_size_short)
READ_ARRAY_SIZE(char,read_size_char)
READ_ARRAY_SIZE(int,read_size_int)
READ_ARRAY_SIZE(float,read_size_float)
READ_ARRAY_SIZE(double,read_size_double)
/*********************************************************************/
/*                                                                   */
/*                          DUPLICATION                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

#define SET_NUMBERS(type,function)\
type * function(type *list,int n, ...)\
     {\
     type *buf;\
     int  *index,i;\
     va_list ap;\
     int max;\
\
     va_start(ap, n);\
     buf=vcalloc(n, sizeof(type));\
     index=vcalloc(n, sizeof(int));\
     max=0;\
     for ( i=0; i< n; i++)\
         {\
	 buf[i]  =va_arg(ap,type);\
	 index[i]=va_arg(ap, int);\
	 if (index[i]>max)max=index[i];\
	 }\
     va_end(ap);\
     if (list==NULL)list=vcalloc ( max+1, sizeof (type));\
     for ( i=0; i<n; i++)list[index[i]]=buf[i];\
     vfree(buf);\
     vfree(index);\
     return list;\
     }
     /*SET_NUMBERS(short ,set_short)*/
     /*SET_NUMBERS(char  ,set_char)*/
SET_NUMBERS(int   ,set_int)
     /*SET_NUMBERS(float ,set_float)*/
SET_NUMBERS(double,set_double)

short ** duplicate_short ( short **array , int len, int field)
    {
    return copy_short (array ,declare_short ( len, field),  len, field);
    }
int ** duplicate_int ( int **array , int len, int field)
    {
    return copy_int (array ,declare_int ( len, field),  len, field);
    }
char ** duplicate_char ( char **array , int len, int field)
    {
    return copy_char (array ,declare_char ( len, field),  len, field);
    }
float ** duplicate_float ( float **array , int len, int field)
    {
    return copy_float (array ,declare_float ( len, field),  len, field);
    }
double ** duplicate_double ( double **array , int len, int field)
    {
    return copy_double (array ,declare_double ( len, field),  len, field);
    }



/*********************************************************************/
/*                                                                   */
/*                           COPY OF 2D ARRAY                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
short ** copy_short( short **array1, short  **array2, int len, int number_field)
    {
    int a;
    
    if ( len==-1)len=read_size_short (array1[-1], 0);
    if ( number_field==-1)number_field=read_size_short (array1[-1],1);
    if ( array2)free_short ( array2, -1);
    array2=declare_short ( len, number_field);

    for ( a=0; a< len; a++)
	ga_memcpy_short( array1[a],array2[a],number_field);
    
    return array2;
    }
char ** copy_char ( char **array1, char **array2, int len, int number_field)
    {
    int a;
    
    if ( len==-1)len=read_size_char (array1[-1], 0);
    if ( number_field==-1)
	{
	number_field=read_size_char (array1[-1],1);
	for ( a=0; a< len; a++)
	    number_field=MAX(number_field, strlen ( array1[a]))+1;
	}
    
    if ( array2)free_char (array2, -1);
    array2=declare_char(len, number_field);
    
    for ( a=0; a< len; a++)
      sprintf ( array2[a], "%s", array1[a]);

    return array2;
    } 
int ** copy_int ( int **array1, int **array2, int len, int number_field)
    {
    int a;

    if ( len==-1)len=read_size_int (array1[-1], 0);
    if ( number_field==-1)number_field=read_size_int (array1[-1],1);
    
    
    
    if (array2)free_int (array2, -1);
    array2=declare_int ( len, number_field);

    for ( a=0; a< len; a++)
	ga_memcpy_int( array1[a],array2[a],number_field);
    
    return array2;
    }

float ** copy_float ( float **array1, float **array2, int len, int number_field)
    {
    int a;
    if ( len==-1)len=read_size_float (array1[-1], 0);
    if ( number_field==-1)number_field=read_size_float (array1[-1],1);

    if ( array2)free_float (array2, -1);
    array2=declare_float ( len, number_field);
    
    for ( a=0; a< len; a++)
	ga_memcpy_float( array1[a],array2[a],number_field);
    return array2;
    }
double ** copy_double (double **array1, double **array2, int len, int number_field)
    {
    int a;
    if ( len==-1)len=read_size_double (array1[-1], 0);
    if ( number_field==-1)number_field=read_size_double (array1[-1],1);
    if ( array2)free_double (array2, -1);
    array2=declare_double ( len, number_field);

    for ( a=0; a< len; a++)
	ga_memcpy_double( array1[a],array2[a],number_field);
    return array2;
    }
/*********************************************************************/
/*                                                                   */
/*                        CONCATENATION                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/





Alignment ** cat_aln_list ( Alignment **list_to_cat,int first, int end, Alignment **rec_list)
    {
    int rec_list_start;
    int a, b;

    if ( list_to_cat==NULL)return rec_list;
    else
       {
       rec_list_start=(rec_list[-1])->nseq;       
       rec_list=realloc_aln_array ( rec_list, end-first);      
       for ( a=first, b=rec_list_start; a<end; a++, b++)copy_aln (list_to_cat[a], rec_list[b]);      
       free_aln_array ( list_to_cat);      
       return rec_list;
       }
    return rec_list;
    }
	   
/*********************************************************************/
/*                                                                   */
/*                         NUMBER ARRAY ANALYSE                      */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

FILE * output_array_int (int  **array, int len, int nf ,FILE *fp)
     {
     int a, b;
     for ( a=0; a<len; a++)
         {fprintf ( fp, "\n");
	    for (b=0; b<nf; b++)fprintf (fp, "%5d ", array[a][b]);
     
	 }
     fprintf ( fp, "\n");
     return fp;
     }
#define RETURN_MAX_COOR(type,wf,rf,function,comparison,undef)\
type function ( type ** array, int len_array, int field, int *coor)\
    {\
    type max;\
    int a;\
\
    if (array==NULL || len_array==0)return 0;\
    else\
        {\
        if (len_array==-1)len_array=rf(array[-1],0);\
        max=array[0][field];\
        coor[0]=0;\
        for ( a=1; a< len_array; a++)\
	      if ( max==undef)max=array[a][field];\
	      else if ( array[a][field]!=undef)\
		   if (array[a][field] comparison max)\
		       {max=array[a][field];\
                        coor[0]=a;\
	                }\
	 }\
    return max;\
    }

RETURN_MAX_COOR(short,write_size_short,read_size_short,return_max_coor_short,>, UNDEFINED_SHORT)
RETURN_MAX_COOR(char,write_size_char,read_size_char,return_max_coor_char,>, UNDEFINED_CHAR)
RETURN_MAX_COOR(int,write_size_int,read_size_int,return_max_coor_int,>, UNDEFINED_INT)
RETURN_MAX_COOR(float,write_size_float,read_size_float,return_max_coor_float,>, UNDEFINED_FLOAT)
RETURN_MAX_COOR(double,write_size_double,read_size_double,return_max_coor_double,>, UNDEFINED_DOUBLE)
RETURN_MAX_COOR(short,write_size_short,read_size_short,return_min_coor_short,<, UNDEFINED_SHORT)
RETURN_MAX_COOR(char,write_size_char,read_size_char,return_min_coor_char,<, UNDEFINED_CHAR)
RETURN_MAX_COOR(int,write_size_int,read_size_int,return_min_coor_int,<, UNDEFINED_INT)
RETURN_MAX_COOR(float,write_size_float,read_size_float,return_min_coor_float,<, UNDEFINED_FLOAT)
RETURN_MAX_COOR(double,write_size_double,read_size_double,return_min_coor_double,<, UNDEFINED_DOUBLE)
#define RETURN_MAX(type,wf,rf,function,comparison,undef)\
type function ( type ** array, int len_array, int field)\
    {\
    type max;\
    int a;\
\
    if (array==NULL || len_array==0)return 0;\
    else\
        {\
        if (len_array==-1)len_array=rf(array[-1],0);\
        max=array[0][field];\
        for ( a=1; a< len_array; a++)\
            if ( max==undef)max=array[a][field];\
	    else if ( array[a][field]!=undef)max=( array[a][field] comparison max)?array[a][field]:max;\
        }\
    return (max==undef)?0:max;\
    }

RETURN_MAX(short,write_size_short,read_size_short,return_max_short,>,UNDEFINED_SHORT)
RETURN_MAX(char,write_size_char,read_size_char,return_max_char,>,UNDEFINED_CHAR)
RETURN_MAX(int,write_size_int,read_size_int,return_max_int,>,UNDEFINED_INT)
RETURN_MAX(float,write_size_float,read_size_float,return_max_float,>,UNDEFINED_FLOAT)
RETURN_MAX(double,write_size_double,read_size_double,return_max_double,>,UNDEFINED_DOUBLE)
RETURN_MAX(short,write_size_short,read_size_short,return_min_short,<,UNDEFINED_SHORT)
RETURN_MAX(char,write_size_char,read_size_char,return_min_char,<,UNDEFINED_CHAR)
RETURN_MAX(int,write_size_int,read_size_int,return_min_int,<,UNDEFINED_INT)
RETURN_MAX(float,write_size_float,read_size_float,return_min_float,<,UNDEFINED_FLOAT)
RETURN_MAX(double,write_size_double,read_size_double,return_min_double,<,UNDEFINED_DOUBLE)



#define RETURN_2DMAX(type,wf,rf,function,comparison,undef)\
type function ( type ** array, int start, int len_array, int first_field, int number_field)\
    {\
    type max;\
    int a,b;\
    if (array==NULL || len_array==0 || first_field<0 || number_field==0)return 0;\
    else\
         {max=array[start][first_field];\
          for ( a=start; a< start+len_array; a++)\
	      for (b=first_field; b< first_field+number_field; b++)\
	             if (array[a][b]!=undef)max=( array[a][b] comparison max)?array[a][b]:max;\
         }\
    return max;\
    }
RETURN_2DMAX(short,write_size_short,read_size_short,return_2Dmax_short,>, UNDEFINED_SHORT)
RETURN_2DMAX(char,write_size_char,read_size_char,return_2Dmax_char,>,UNDEFINED_CHAR)
RETURN_2DMAX(int,write_size_int,read_size_int,return_2Dmax_int,>,UNDEFINED_INT)
RETURN_2DMAX(float,write_size_float,read_size_float,return_2Dmax_float,>,UNDEFINED_FLOAT)
RETURN_2DMAX(double,write_size_double,read_size_double,return_2Dmax_double,>,UNDEFINED_DOUBLE)
RETURN_2DMAX(short,write_size_short,read_size_short,return_2Dmin_short,<,UNDEFINED_SHORT)
RETURN_2DMAX(char,write_size_char,read_size_char,return_2Dmin_char,<,UNDEFINED_CHAR)
RETURN_2DMAX(int,write_size_int,read_size_int,return_2Dmin_int,<,UNDEFINED_INT)
RETURN_2DMAX(float,write_size_float,read_size_float,return_2Dmin_float,<,UNDEFINED_FLOAT)
RETURN_2DMAX(double,write_size_double,read_size_double,return_2Dmin_double,<,UNDEFINED_DOUBLE)

#define RETURN_2DMAX_COOR(type,wf,rf,function,compare,undef)\
type function ( type **array, int start1 , int end1, int start2, int end2,int *i, int *j)\
    {\
    int a, b;\
    double max=undef;\
    if ( start1==-1)start1=0;\
    if ( start2==-1)start2=0;\
    if ( end1==-1)end1=rf(array[-1],0);\
    if ( end2==-1)end2=rf(array[-1],1);\
    if ( array==NULL || (end1-start1)==0 || (end1-start1)>rf ( array[-1],0) || (end2-start2)==0)\
        {\
	return 0;\
        i[0]=0;\
        j[0]=0;\
        }\
    i[0]=0;\
    j[0]=0;\
    for ( a=start1; a<end1; a++)\
	for ( b=start2; b<end2; b++)\
	    {\
            if ( max==undef && array[a][b]!=undef)max=array[a][b];\
	    else if ( array[a][b]!=undef && (array[a][b] compare max))\
	       {\
	       max=array[a][b];\
	       i[0]=a;\
	       j[0]=b;\
	       }\
	    }\
    return (type)max;\
    }
RETURN_2DMAX_COOR(short,write_size_short,read_size_short,return_2Dmax_coor_short,>,UNDEFINED_SHORT)
RETURN_2DMAX_COOR(char,write_size_char,read_size_char,return_2Dmax_coor_char,>,UNDEFINED_CHAR)
RETURN_2DMAX_COOR(int,write_size_int,read_size_int,return_2Dmax_coor_int,>,UNDEFINED_INT)
RETURN_2DMAX_COOR(float,write_size_float,read_size_float,return_2Dmax_coor_float,>,UNDEFINED_FLOAT)
RETURN_2DMAX_COOR(double,write_size_double,read_size_double,return_2Dmax_coor_double,>,UNDEFINED_DOUBLE)
RETURN_2DMAX_COOR(short,write_size_short,read_size_short,return_2Dmin_coor_short,<,UNDEFINED_SHORT)
RETURN_2DMAX_COOR(char,write_size_char,read_size_char,return_2Dmin_coor_char,<,UNDEFINED_CHAR)
RETURN_2DMAX_COOR(int,write_size_int,read_size_int,return_2Dmin_coor_int,<,UNDEFINED_INT)
RETURN_2DMAX_COOR(float,write_size_float,read_size_float,return_2Dmin_coor_float,<,UNDEFINED_FLOAT)
RETURN_2DMAX_COOR(double,write_size_double,read_size_double,return_2Dmin_coor_double,<,UNDEFINED_DOUBLE)

#define RETURN_WMEAN(type,wf,rf,function,sum_function,undef)\
double function ( type **array, int len, int wfield,int sfield)\
    {\
    double b;\
    int a, c;\
    if ( len==0 ||array==NULL || len>rf ( array[-1],0))return 0;\
    else\
         {\
         if ( len==-1)len=rf(array[-1],0);\
         for ( b=0, c=0,a=0; a< len; a++)\
             {\
	     if (array[a][sfield]!=undef && array[a][wfield]!=undef )\
	        {\
		b+=array[a][sfield];\
		c+=array[a][wfield];\
		}\
             }\
         }\
    return (c==0)?0:(b/c);\
    }
RETURN_WMEAN(short,write_size_short,read_size_short,return_wmean_short, return_sum_short,UNDEFINED_SHORT)
RETURN_WMEAN(char,write_size_char,read_size_char, return_wmean_char,return_sum_char,UNDEFINED_CHAR)
RETURN_WMEAN(int,write_size_int,read_size_int,return_wmean_int,return_sum_int,UNDEFINED_INT)
RETURN_WMEAN(float,write_size_float,read_size_float,return_wmean_float,return_sum_float,UNDEFINED_FLOAT)
RETURN_WMEAN(double,write_size_double,read_size_double,return_wmean_double,return_sum_double,UNDEFINED_DOUBLE)

		     
#define RETURN_MEAN(type,wf,rf,function,sum_function,undef)\
double function ( type **array, int len, int field)\
    {\
    double b;\
    int a, c;\
    if ( len==0 ||array==NULL || len>rf ( array[-1],0))return 0;\
    else\
         {\
         for ( b=0, c=0,a=0; a< len; a++)\
             {\
	     if (array[a][field]!=undef)\
	        {\
		b+=array[a][field];\
		c++;\
		}\
             }\
         }\
    return (c==0)?0:(b/c);\
    }
RETURN_MEAN(short,write_size_short,read_size_short,return_mean_short, return_sum_short,UNDEFINED_SHORT)
RETURN_MEAN(char,write_size_char,read_size_char, return_mean_char,return_sum_char,UNDEFINED_CHAR)
RETURN_MEAN(int,write_size_int,read_size_int,return_mean_int,return_sum_int,UNDEFINED_INT)
RETURN_MEAN(float,write_size_float,read_size_float,return_mean_float,return_sum_float,UNDEFINED_FLOAT)
RETURN_MEAN(double,write_size_double,read_size_double,return_mean_double,return_sum_double,UNDEFINED_DOUBLE)

#define RETURN_SUM(type,wf,rf,function,undef)\
type function(type **array, int len, int field)\
{\
 int a;\
 type b=0;\
 if ( len==0 ||array==NULL)return 0;\
 else\
     {\
     if ( len==-1)len=rf ( array[-1],0);\
     for ( a=0; a< len; a++)\
          if ( array[a][field]!=undef)b+=array[a][field];\
     }\
  return b;\
  }
RETURN_SUM(short,write_size_short,read_size_short, return_sum_short,UNDEFINED_SHORT)
RETURN_SUM(char,write_size_char,read_size_char,return_sum_char,UNDEFINED_CHAR)
RETURN_SUM(int,write_size_int,read_size_int,return_sum_int,UNDEFINED_INT)
RETURN_SUM(float,write_size_float,read_size_float,return_sum_float,UNDEFINED_FLOAT)
RETURN_SUM(double,write_size_double,read_size_double,return_sum_double,UNDEFINED_DOUBLE)    

#define RETURN_SD(type,wf,rf,function,undef)\
type function ( type **array, int len, int field,type mean)\
    {\
    int a;\
    double c=0;\
    if ( len==0 ||array==NULL || len>rf ( array[-1],0))return 0;\
    else\
        {\
        for ( a=0; a< len; a++)\
    	    {\
	    if ((array[a][field]!=undef) && (mean-array[a][field])!=0)\
             c+=((double)mean-array[a][field])*((double)mean-array[a][field]);\
            }\
        c=sqrt(c)/(double)len;\
	return (type)MAX(c,1);\
	}\
    }
RETURN_SD(short,write_size_short,read_size_short, return_sd_short,UNDEFINED_SHORT)
RETURN_SD(char,write_size_char,read_size_char,return_sd_char,UNDEFINED_CHAR)
RETURN_SD(int,write_size_int,read_size_int,return_sd_int,UNDEFINED_INT)
RETURN_SD(float,write_size_float,read_size_float,return_sd_float,UNDEFINED_FLOAT)
RETURN_SD(double,write_size_double,read_size_double,return_sd_double,UNDEFINED_DOUBLE) 


double return_z_score( double x,double sum, double sum2, double n)
    {
    double sd;
    double avg;
    double z;
    
   
    sd=(n==0)?0:sqrt(sum2*n -sum*sum)/n;
    avg=(n==0)?0:(sum/n);
    z=(sd==0)?0:(x-avg)/sd;
    return z;
    }
#define INVERT_LIST(type,wf,rf,function,swap_function)\
type* function (type *list, int len)\
    {\
    int a, b;\
    for ( a=0, b=len-1; a<b; a++, b--)swap_function ( &list[a], &list[b], 1);\
    return list;\
    }
INVERT_LIST(short,write_size_short,read_size_short, invert_list_short,swap_short)
INVERT_LIST(char,write_size_char,read_size_char,invert_list_char,swap_char)
INVERT_LIST(int,write_size_int,read_size_int,invert_list_int,swap_int)
INVERT_LIST(float,write_size_float,read_size_float,invert_list_float,swap_float)
INVERT_LIST(double,write_size_double,read_size_double,invert_list_double,swap_double) 

#define SWAP_FUNCTION(type,wf,rf,function)\
void function(type *a, type *b, int n)\
    {\
    type t;\
    int c;\
    for ( c=0;c<n;c++)\
	{t=a[c];\
	 a[c]=b[c];\
	 b[c]=t;\
	}\
    }
SWAP_FUNCTION(short,write_size_short,read_size_short,swap_short)
SWAP_FUNCTION(char,write_size_char,read_size_char,swap_char)
SWAP_FUNCTION(int,write_size_int,read_size_int,swap_int)
SWAP_FUNCTION(float,write_size_float,read_size_float,swap_float)
SWAP_FUNCTION(double,write_size_double,read_size_double,swap_double) 

#define RETURN_MAX_HORIZ(type,wf,rf,function,comparison,undef)\
type function  (type ** array, int len_array, int field)\
    {\
    type max;\
    int a;\
    if ( len_array==0)return 0;\
    else\
        {\
	max=array[field][0];\
        for ( a=1; a< len_array; a++)\
	    if ( array[field][a]!=undef) max=( array[field][a] comparison max)?array[field][a]:max;\
        return (int)max;\
        }\
    }
RETURN_MAX_HORIZ(short,write_size_short,read_size_short,return_max_short_hor,>,UNDEFINED_SHORT)
RETURN_MAX_HORIZ(char,write_size_char,read_size_char,return_max_char_hor,>,UNDEFINED_CHAR)
RETURN_MAX_HORIZ(int,write_size_int,read_size_int,return_max_int_hor,>,UNDEFINED_INT)
RETURN_MAX_HORIZ(float,write_size_float,read_size_float,return_max_float_hor,>,UNDEFINED_FLOAT)
RETURN_MAX_HORIZ(double,write_size_double,read_size_double,return_max_double_hor,>,UNDEFINED_DOUBLE) 

RETURN_MAX_HORIZ(short,write_size_short,read_size_short,return_min_short_hor,<,UNDEFINED_SHORT)
RETURN_MAX_HORIZ(char,write_size_char,read_size_char,return_min_char_hor,<,UNDEFINED_CHAR)
RETURN_MAX_HORIZ(int,write_size_int,read_size_int,return_min_int_hor,<,UNDEFINED_INT)
RETURN_MAX_HORIZ(float,write_size_float,read_size_float,return_min_float_hor,<,UNDEFINED_FLOAT)
RETURN_MAX_HORIZ(double,write_size_double,read_size_double,return_min_double_hor,<,UNDEFINED_DOUBLE) 

#define BEST_OF_MANY(type,wf,rf,function,undef)\
type function (int n, ...)\
	{\
	va_list ap;\
	int *fop,a;\
	type v, best;\
	int maximise;\
	/*first Arg: number of values\
	  2nd   Arg: maximise(1)/minimise(0)\
	  3rd   Arg: *int contains the indice of the best value\
	  ...   Arg: n type values\
	*/\
	va_start (ap, n);\
	maximise=va_arg (ap, int);\
	fop=va_arg (ap, int*);\
	best=va_arg (ap, type);\
	fop[0]=0;\
	for ( a=1; a<n; a++)\
		{\
		v=va_arg (ap, type);\
		if (best==undef)\
			{\
			best=v;\
			fop[0]=a;\
			}\
		if ( best==undef || v==undef);\
		else if ( maximise==1 && v>best)\
			{\
			fop[0]=a;\
			best=v;\
			}\
		else if ( maximise==0 && v<best)\
			{\
			fop[0]=a;\
			best=v;\
			}\
		}\
	va_end (ap);\
	return best;\
	}
     /*BEST_OF_MANY(short,write_size_short,read_size_short, best_short,UNDEFINED_SHORT)*/
     /*BEST_OF_MANY(char,write_size_char,read_size_char,best_char,UNDEFINED_CHAR)*/
BEST_OF_MANY(int,write_size_int,read_size_int,best_int,UNDEFINED_INT)
     /*BEST_OF_MANY(float,write_size_float,read_size_float,best_float,UNDEFINED_FLOAT)*/
BEST_OF_MANY(double,write_size_double,read_size_double,best_double,UNDEFINED_DOUBLE)
#define IS_DEFINED(type,function,undef)\
int function(int n, ...)\
     {\
     int i;\
     va_list ap;\
\
     va_start(ap, n);\
     for ( i=0; i< n; i++)\
         {\
	 if(va_arg(ap,type)==undef)\
	      {\
		va_end(ap);\
		return 0;\
	      }\
	 }\
     va_end(ap);\
     return 1;\
     }	
     /*IS_DEFINED(short,is_defined_short,UNDEFINED_SHORT)*/
     /*IS_DEFINED(char,is_defined_char,  UNDEFINED_CHAR)*/
IS_DEFINED(int,is_defined_int,   UNDEFINED_INT)
     /*IS_DEFINED(float,is_defined_float, UNDEFINED_FLOAT)*/
IS_DEFINED(double,is_defined_double,UNDEFINED_DOUBLE)
int return_maxlen ( char ** array, int number)
    {
    int a;
    int max=0;
    for ( a=0; a< number; a++)
	max=( strlen ( array[a])>max)?strlen ( array[a]):max;

    return max;
    }

int return_minlen ( char ** array, int number)
    {
    int a;
    int min;

    min=strlen( array[0]);
    for ( a=1; a< number; a++)
	min=( strlen ( array[a])>min)?strlen ( array[a]):min;

    return min;
    }
 


float return_mean_diff_float ( float **array, int len, int field,float mean)
    {
    int a;
    float b=0;
    
    for ( a=0; a< len; a++)
    	{
    	if ( (mean-array[a][field])!=0) 
	 	b+=sqrt((double)((float) ( mean-array[a][field])*(float)(mean-array[a][field])));
	}
	
    return ((float)b/(float)len);
    }



void inverse_int ( int**array, int len, int field, int max, int min)
    {
    int a;
    for ( a=0; a< len; a++)
	array[a][field]=max-array[a][field]+min;
    }
void inverse_float ( float**array, int len, int field, int max, int min)
    {
    int a;
    for ( a=0; a< len; a++)
	array[a][field]=max-array[a][field]+min;
    }
void inverse_2D_float ( float **array, int start, int len, int start_field, int number_field, float max,float min)
    {
    int a, b;
    for ( a=0; a< start+len; a++)
	for ( b=start_field; b< start_field+ number_field; b++)
	    array[a][b]=max-array[a][b]+min;
    }
/*********************************************************************/
/*                                                                   */
/*                         SHELL INTERFACES                          */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char* get_env_variable ( const char *var, int mode)
        {
	    /*mode 0: return NULL if variable not set*/
	    /*mode 1: crash if variable not set*/
	    if ( !getenv (var))
	       {
		   if ( !mode)return NULL;
		   else
		      {
			  fprintf ( stderr, "\nYou must set the variable %s [FATAL]\n", var);
	 		  exit (1);
			  return NULL;
		      }
	       }
	    else return getenv (var);
	}
		

	    
void setenv_func ( char *string_name, char *string_value)
	{
	char command[1000];
	sprintf (command, "setenv %s %s", string_name,string_value);
	
	system ( command);
	}
void get_pwd ( char *name)
	{
	char *string;
	char command[1000];
	FILE *fp;
	 
	string=vcalloc ( L_tmpnam, sizeof (char));
	string=vtmpnam(string);
	sprintf ( command, "pwd > %s", string);
	system (command);
	fp=fopen ( string, "r");
	fscanf ( fp, "%s",name);
	fclose (fp);
	sprintf ( command, "rm %s", string);
	system ( command);
	vfree (string);
	
	}
int pg_is_installed ( char *pg)
        {
	char *fname;
	char command[1000];
	FILE *fp;
	int r=0;

	return 1;

	fname=vcalloc ( L_tmpnam+1, sizeof (char));
	fname= vtmpnam(fname);
	
	sprintf ( command, "which %s > %s", pg, fname);
	system ( command);
	
	
        if (fp=find_token_in_file ( fname, NULL, "Command")){r=1;vfclose(fp);}
	 
	free (fname);
	return r;
	
	}

	    
/*********************************************************************/
/*                                                                   */
/*                           MISC                                    */  
/*                                                                   */
/*********************************************************************/
char *num2plot (int value, int max, int line_len)
        {
	       int   len;
	       int   value_len;
	       char *buf;
	static char *string;

	if ( string==NULL)string=vcalloc (1000, sizeof(char));
 
	if ( line_len==-1)len=30;
	else len=line_len;
	
	value_len=((float)value/(float)max)*(float)len;
	if ( value==0)
	    sprintf ( string, "|");
	else
	    {
	    buf=generate_string(value_len, '*');
	    sprintf ( string,"%s", buf);
	    vfree(buf);
	    }
	return string;
	}
	
	   
float grep_function ( char *pattern, char *file)
	{
	char command [1000];
	int a, b, l;
	char buf1[100];
	char buf2[100];
	FILE*fp;
	char *s;
	float f;
	
	
	s=calloc ( L_tmpnam, sizeof (char));
	s=vtmpnam(s);
	
	sprintf ( command, "GREP_COMMAND %s %s > %s",pattern,file, s);
	system ( command);
	if ((fp=fopen (s, "r"))==NULL )return 0;
	else
		{
		fgets ( command, 900, fp);
		l=strlen ( command);
		while ( !isdigit (command[l]))l--;
		a=0;
		while ( isdigit (command[l]) || command[l]=='.')
			{
			buf1[a++]=command[l];
			l--;
			}
		buf1[a]='\0';
		l=strlen (buf1);
		for ( a=0, b=l-1;a< l; a++, b--)
			buf2[b]=buf1[a];
		buf2[l]='\0';
		
		sscanf ( buf2, "%f", &f);
		sprintf ( command,"rm %s", s);
		system ( command); 
		return f;
		}
	}    

void crash_if ( int val, char *s)
    {
    if ( val==0)crash(s);
    }
void crash ( char *s)
	{
	int *a;

	
	fprintf ( stderr, "\nCRASH\n");
	fprintf ( stderr, "%s\n",s);
	
	a=vcalloc ( 10, sizeof (int));
	a[20]=1;
	
	exit (1);
	}

static int *local_table;
int ** make_recursive_combination_table ( int tot_n_param, int *n_param, int *nc, int**table, int field)
    {
    int a, b, c;
    
    /* makes a table of all possible combinations*/

    if ( tot_n_param==0)
	{
	    nc[0]=1;
	    fprintf ( stderr, "\nNULL RETURNED");
	    return NULL;
	}
    if (table==NULL)
        {
        if ( local_table!=NULL)vfree (local_table);
	local_table=vcalloc ( tot_n_param, sizeof (int));
	field=0;
	for ( a=0; a< tot_n_param; a++)local_table[a]=-1;
	for ( a=0; a< tot_n_param; a++)nc[0]=nc[0]*n_param[a];
	

	table=declare_int ( nc[0],tot_n_param);
	nc[0]=0;
	}
    
    for ( b=0; b<n_param[field]; b++)
	       {
               
               local_table[field]=b;
	       if ( field<tot_n_param-1)
	          {
                  table=make_recursive_combination_table ( tot_n_param, n_param, nc, table, field+1);
		  }
	       else
	          {
                  for ( c=0; c< tot_n_param; c++)table[nc[0]][c]=local_table[c];
		  nc[0]++;
		  }
	       }
    return table;
    }
	          
/*********************************************************************/
/*                                                                   */
/*                         STRING PROCESSING                         */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
/*Replace by a gap the parts of the two strings that do not OVERLAP*/
/*returns the length of the overlap*/
FILE *print_array_char (FILE *out, char **array, int n, char *sep)
        {
	int a;
	if ( array==NULL || read_size_char (array[-1],0)<n)
	   {
	   fprintf ( stderr, "\nORB in print_array_char [FATAL]\n");
	   crash("");
	   }
	for ( a=0; a< n; a++)fprintf ( out, "%s%s", array[a],sep);
	return out;
	}

Fname* parse_fname ( char *array)
	 {
	 int l;
	 Fname *F;



	F=declare_fname ();
	
	sprintf ( F->path, "%s", array);	
	l=strlen (array);
	while (l!=-1 && (F->path)[l]!='/')(F->path)[l--]='\0';
	
	sprintf ( F->name, "%s", array+l+1);
	l=strlen (F->name);
	while (l!=-1)
	    {
	    if((F->name)[l]=='.')
	      {
	      F->name[l]='\0';
	      sprintf ( F->suffix, "%s", F->name+l+1);
	      break;
	      }
	    else l--;
	    }
	return F;
        }
char *extract_suffixe ( char *array)
        {
	int l;
	char *new_string;
	char *x;
	l=strlen (array);
	new_string=vcalloc ( l+1, sizeof (char));
	sprintf (new_string, "%s",array);
	
	x=new_string+l;
	while (x!=new_string && x[0]!='.' && x[0]!='/' )x--;
	if ( x[0]=='.')x[0]='\0';
	else if (x[0]=='/')return x+1;
 
	while ( x!=new_string && x[0]!='/')x--;
	
	return (x[0]=='/')?x+1:x;
	}
void string_array_upper ( char **string, int n)
     {
     int a;
     for ( a=0; a< n; a++)upper_string (string[a]);
     }
void string_array_lower ( char **string, int n)
     {
     int a;
     for ( a=0; a< n; a++)lower_string (string[a]);
     }

char *upper_string ( char *string)
	{
	int len, a;
	
	len=strlen ( string);
	for ( a=0; a< len; a++)string[a]=toupper ( string[a]);
	return string;
	}	
char *lower_string ( char *string)
	{
	int len, a;
	
	len=strlen ( string);
	for ( a=0; a< len; a++)string[a]=tolower ( string[a]);
	return string;
	}	
void string_array_convert ( char **array, int n_strings, int ns, char **sl)
        {
	int a;

	for ( a=0; a< n_strings; a++)string_convert ( array[a], ns, sl);
	}
void string_convert( char *string, int ns, char **sl)
        {
	int a, l;
	l=strlen ( string);
	for ( a=0; a< l; a++)
	    string[a]=convert(string[a], ns, sl);
	}
int convert ( char c, int ns, char **sl)
        {
	int a;
	int return_char;

	for ( a=0; a< ns; a++)
	    {
	    if ((return_char=convert2 ( c, sl[a]))!=-1)
		return return_char;
	    }
	return c;
	

	}
int convert2 ( char c, char *list)
    {
    int a;
    int l1;
    int return_char;

    l1=strlen ( list);
    
    return_char=(list[l1-1]=='#')?c:list[l1-1];

    for ( a=0; a< l1; a++)
	   if (list[a]=='#')return return_char;
           else if ( list[a]==c)return return_char;
    
    return -1;
    }

	


int str_overlap ( char *string1, char *string2, char x)
        {	
	int a, b;
	int l1, l2;
	int **array;
	int max1, max2;
	int score;
	int max_val=-1;
	int end1, end2;
	int start1, start2;

	if ( strm ( string1, string2))return 0;
	else
	    {
	    l1=strlen ( string1);
	    l2=strlen ( string2);
	    array=declare_int ( strlen ( string1), strlen ( string2));
	    for ( a=0; a< l1; a++)
		for ( b=0; b< l2; b++)
		    {
			if ( a==0 || b==0)array[a][b]=(string1[a]==string2[b])?1:0;
			else
			    if (string1[a]==string2[b])
			       {
			       score=array[a][b]=array[a-1][b-1]+1;
			       if ( max_val<score)
			          {
				  max_val=score;
				  max1=a;
				  max2=b;
				  }
			       }
		    }
	    start1=(max1+1)-max_val;
	    end1=  max1;
	    
	    start2=(max2+1)-max_val;
	    end2=  max2;
	    
	    for ( a=0; a< l1; a++)if ( a<start1 || a> end1)string1[a]=x;
	    for ( a=0; a< l2; a++)if ( a<start2 || a> end2)string2[a]=x;
	    
	    free_int ( array, l1);
	    
	    return max_val;
	    }
	}

int get_string_line ( int start, int n_lines, char *in, char *out)
	{
	int nl=0;
	int a=0;
	int c;
	
	while ( nl<n_lines)
		{
		while ( (c=in[start++])!='\n' && c!='\0')
			{
			out[a++]=c;
			}
		out[a++]='\n';
		nl++;
		}
	out[a]='\0';
	return (c=='\0')?-1:start;
	}
	

FILE * output_string_wrap ( int wrap,char *string, FILE *fp)
	{

	 int a, b,l;
	 
	 l=strlen ( string);
	 
	 for ( a=0, b=1; a< l; a++, b++)
	 	{
	 	fprintf ( fp, "%c", string[a]);
	 	if ( b==wrap)
	 		{
	 		fprintf ( fp, "\n");
	 		b=0;
	 		}
	 	}
	 return fp;
	 }

char * extract_char ( char * array, int first, int len)
    {
    char *array2;
    int a;

    len= ( len<0)?0:len;    

    array2=vmalloc ( sizeof ( char)*(len+1));
    
    for ( a=0; a<len; a++)
	array2[a]=array[first++];

    array2[a]='\0';

    return array2;
    }    


char** break_list ( char **argv, int *argc, char *separators)
       {
       int a, b;
       int n_in;
       char **out;
       char **ar;
       int n_ar;
       
       out=vcalloc ( 1000, sizeof (char*));
       n_in=argc[0];
       argc[0]=0;
       
       for ( a=0; a< n_in; a++)
           {
	   ar=get_list_of_tokens( argv[a], separators,&n_ar);
	   for ( b=0; b< n_ar; b++){out[argc[0]++]=ar[b];}	   
	   }
      
       return out;
       }
    
char ** get_list_of_tokens ( char *in_string, char *separators, int *n_tokens)
{
    char **list=NULL;
    char *p=NULL;
    char *string;
    
    n_tokens[0]=0;
    if ( in_string==NULL || strm(in_string, ""));
    else
        {
	list=declare_char (strlen ( in_string)+1, 1);	
	string=vcalloc ( strlen(in_string)+1, sizeof (char));
	sprintf ( string, "%s", in_string);	
	
	while ( (p=strtok ((p==NULL)?string:NULL, ((separators==NULL)?SEPARATORS:separators)))!=NULL)
           {
	   list[n_tokens[0]]=realloc ( list[n_tokens[0]], sizeof (char) *strlen (p)+1);
	   sprintf ( list[n_tokens[0]], "%s", p);
	   n_tokens[0]++;
	   }

	vfree (string);
	}
   return list;
   }

char **ungap_array ( char **array, int n)
        {
	int a;
	for ( a=0; a< n; a++)ungap(array[a]);
	return array;
	}

void ungap ( char *seq)
	{
	int a, b, l;
	
	l=strlen ( seq);
	for (b=0, a=0; a<=l; a++)
		if ( is_gap(seq[a]));
		else seq[b++]=seq[a];	
	seq[b]='\0';
	}
char **char_array2number ( char ** array, int n)
       {
       int a;
       for ( a=0; a< n; a++)array[a]=char2number(array[a]);
       return array;
       }
char *char2number ( char * array)
       {
       int a, l;
       
       
       l=strlen ( array);
       for ( a=0; a< l; a++)
	   {
	     if ( isdigit(array[a]) && array[a]!=NO_COLOR_RESIDUE && array[a]!=NO_COLOR_GAP )array[a]-='0';
	     else if ( array[a]==NO_COLOR_RESIDUE || array[a]==NO_COLOR_GAP)array[a]=NO_COLOR_RESIDUE;
	   }
           
       return array;
       }
char *mark_internal_gaps(char *seq, char symbol)
{
  int l, a, gap;
  int in_seq;
  char *cache_seq;
  
  l=strlen(seq);
  cache_seq=vcalloc ( l+1, sizeof (char));
  sprintf ( cache_seq, "%s", seq);
  
  for ( gap=0, in_seq=0,a=0; a< l; a++)
    {
      gap=is_gap(seq[a]);
      if ( !gap && !in_seq)in_seq=1;
      if (gap && in_seq)seq[a]=symbol;      
    }
  
  for (gap=0, in_seq=0,a=l-1; a>=0; a--)
    {
      gap=is_gap(seq[a]);
      if ( !gap && !in_seq)break;
      if (gap && !in_seq)seq[a]=cache_seq[a];
    }
  free(cache_seq);
  return seq;
}

void splice_out ( char *seq, char x)
    
        {
	int a, b, l;

	l=strlen ( seq);
	for (b=0, a=0; a<=l; a++)
		if ( seq[a]==x);
		else seq[b++]=seq[a];	
	seq[b]='\0';
	}
int isblanc ( char *buf)
   {
   int a, l;
   
   if ( buf==NULL)return 0;
   l=strlen (buf);
   for ( a=0; a< l; a++)
   	if (isalnum (buf[a]))return 0;
   return 1;
   }

int is_number ( char *num)
   {
   return array_is_in_set ( num, "0123456789.-+");
   }

int is_alnum_line ( char *buf)
   {
   int a, l;
   l=strlen (buf);
   for ( a=0; a< l; a++)
   	if (isalnum (buf[a]))return 1;
   return 0;
   } 
int is_alpha_line ( char *buf)
   {
   int a, l;
   l=strlen (buf);
   for ( a=0; a< l; a++)
   	if (isalpha (buf[a]))return 1;
   return 0;
   }
	
int get_string_sim ( char *string1, char *string2, char *ignore)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1, r2;
	
	
	len1=strlen (string1);
	len2=strlen (string2);
	
	if ( len1!=len2)return 0;
	
	for ( a=0; a< len1; a++)
		{
		r1=string1[a];
		r2=string2[a];
		if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if (r1==r2)
				{
				sim++;
				}
			}
		}
	
	if (pos==0)
		 return 0;
	else	
		return (int) (sim*100)/pos;
	
	}
int is_aa  ( char x)
         {
	 return (is_in_set (x, AA_ALPHABET) && !is_in_set (x, GAP_LIST));
	 }

int is_dna ( char x)
         {
	 return (is_in_set (x, DNA_ALPHABET)&& !is_in_set (x, GAP_LIST));
	 }
int is_gap ( char x)
         {
	 return ( is_in_set ( x, GAP_LIST));
	 }


char* get_alphabet   ( char *seq, char *alphabet)
/*This function builds an alphabet from the string seq*/
/*Alphabet does not need to be emppty. The total number of unique symbols is returned*/
/*Gaps, as defined in the GAP_LIST are ignored*/
    {
	int a;
	int len;
	
	if ( !alphabet) alphabet=vcalloc ( 200, sizeof (char));

	len=strlen (seq);
	
	for ( a=0; a<len; a++)
	    if ( !is_gap ( seq[a]) && !is_in_set ( seq[a], alphabet+1)){alphabet[++alphabet[0]]=seq[a];}
	return alphabet;
    }

int array_is_in_set ( char *x, char *set )
         {
	 int len;
	 int a;
	 
	 len=strlen (x);
	 for ( a=0; a< len; a++)
	     if ( !is_in_set ( x[a], set))return 0;
	 return 1;
	 }


int is_in_set ( char r, char *list)
	{
	static char s[2];
	
	s[0]=r;
	
	
	if ( strstr ( list,s)!=NULL)return 1;
	
	return 0;
	}
char * generate_void ( int x)
    {
    return generate_string (x, ' ');
    }
   
char * generate_null ( int x)
    {
    return generate_string ( x, '-');
    }
    	 	
char * generate_string ( int x, char y)
    {
	 int a;
    char *string;
    
    string = vcalloc ( x+1, sizeof ( char));

    for ( a=0; a< x; a++)
	string[a]=y;

    string[a]='\0';
    return string;
    }

void translate_name ( char *name)
	{
	int len;
	int a;
	
	len=strlen (name);
	
	for ( a=0; a<len; a++)
		{
		if ( isspace(name[a]))name[a]='\0';
		else if ( name[a]==';' ||name[a]==':' ||name[a]=='(' || name[a]==')')
			name[a]='_';
		}
	} 
int get_longest_string (char **array,int n, int *len, int *index)
     {
     int a, l;
     int max_len, max_index;

     if ( n==0 || array==NULL )return 0;     
     if ( n==-1)n=read_size_char(array[-1],0);
     

     if (read_size_char ( array[-1],0)<n)
         {
	 fprintf ( stderr, "\nBAD REQUEST: Number of strings=%d, expected number=%d", read_size_char ( array[-1],0),n);
	 crash ("[FATAL/util.c]");
	 }
     else
         {
	 max_len=strlen(array[0]);
	 max_index=0;
	 for ( a=1; a< n; a++)
	     {
	     l=strlen ( array[a]);
	     if ( l>max_len)
	        {
		max_len=l;
		max_index=a;
		}
	     }
	 if (index!=NULL)index[0]=max_index;
	 if (len!=NULL)len[0]=max_len;
	 }
     
     return max_len;
     }
int get_shortest_string (char **array,int n, int *len, int *index)
     {
     int a, l;
     int min_len;


     if ( n==0|| array==NULL || read_size_char ( array[-1],0)<n)return 0;
     else
         {
	 min_len=strlen(array[0]);

	 for ( a=1; a< n; a++)
	     {
	     l=strlen ( array[a]);
	     if ( l<min_len)
	        {
		min_len=l;

		}
	     }
	 if (index!=NULL)index[0]=a;
	 if (len!=NULL)len[0]=min_len;
	 }
     return min_len;
     }

char **pad_string_array ( char **array, int n, int len, char pad)
     {
     int a, b, l;
     for ( a=0; a< n; a++)
         {
	 l= strlen ( array[a]);
	 for (b=l; b<len; b++)array[a][b]=pad;
	 array[a][len]='\0';
	 }
     return array;
     }
char **right_pad_string_array ( char **array, int n, int len, char pad)
     {
     return pad_string_array(array, n, len, pad);
     }
char **left_pad_string_array ( char **array, int n, int len, char pad)
     {
     int a, b, c, l;
     char *buf;

     buf=vcalloc ( len+1, sizeof (char));
     for ( a=0; a< n; a++)
         {
	 l= strlen ( array[a]);
	 for (b=0; b<(len-l); b++)buf[b]=pad;
	 for (c=0,b=(len-l); b<=len; c++,b++)buf[b]=array[a][c];	 
	 sprintf (array[a], "%s", buf);

	 }
     vfree(buf);
     return array;
     }
char * crop_string (char *string, int start, int end)
     {
     static char *buf;
     
     if ( strlen (string)<end)
             {
	     fprintf ( stderr, "\nString=%s, start=%d end=%d", string, start, end);
	     crash ( "Wrong End in crop String [FATAL]");
	     }
     else
         {
	 buf=vcalloc (strlen (string)+1, sizeof (char));
	 string[end]='\0';
	 sprintf ( buf, "%s", string+start);
	 sprintf ( string,"%s", buf);
	 vfree (buf);
	 }

     return string;
     }
		  
int get_distance2char ( char *x, char*list)
    {

    char *y;
    y=x;
    if (x==NULL) return 0;
    while (!is_in_set (y[0], list) && y[0]!='\0')y++;
    return y-x;
    }
/*********************************************************************/
/*                                                                   */
/*                         TIME        FUNCTIONS                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
static long ref;
static int child;
static long ticks;
static long milli_sec_conv=1000;
int get_time ()
	{
	static long time;
        struct tms time_buf[1];
	long tms_stime, tms_utime;

	if ( child==1)return get_ctime();
	  
	if ( ticks==0)ticks = sysconf(_SC_CLK_TCK);
	times ( time_buf);

	tms_stime=(long)time_buf->tms_stime*milli_sec_conv;
	tms_utime=(long)time_buf->tms_utime*milli_sec_conv;
	


	
	if ( ref==0)
		{
		ref=(tms_stime+tms_utime);
		return 0;
		}
	else
		{
		time=(tms_utime+tms_stime)-ref;
		return (int) ((time)/ticks);
		}
	}
int get_ctime ()
	{
	static long time;
        struct tms time_buf[1];
	long   tms_cutime, tms_cstime; 
	
	if ( ticks==0)ticks = sysconf(_SC_CLK_TCK);
	times ( time_buf);



	tms_cstime=(long)time_buf->tms_cstime*milli_sec_conv;
	tms_cutime=(long)time_buf->tms_cutime*milli_sec_conv;

	if ( ref==0)
	        {
		child=1;
		ref=tms_cstime+tms_cutime;
		return 0;
		}
	else
		{
		time=(tms_cutime+tms_cstime)-ref;
		return (int)((time)/ticks);
		}
	}
int reset_time()
        {
	ref=0;
	return (int)get_time();
	}
int increase_ref_time(int increase)
        {
	if ( ref==0)get_time();
	
	ref-=(long)ticks*(long)increase;
	if (ref==0)ref++;
	return (int)ref;
	}
/*********************************************************************/
/*                                                                   */
/*                         SYSTEM CALLS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
int evaluate_sys_call_io ( char *out_file, char *com, char *fonc)
     {
     if ( check_file_exists (out_file)==1)return 1;
     else
        {
	fprintf ( stderr, "\nCommand\n%s\nFailed to produce File\n%s\n", com, out_file);
	return 0;
	}
     }

		

/*********************************************************************/
/*                                                                   */
/*                         IO FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char *program_name;
char ** standard_initialisation ( char **in_argv, int *in_argc)
      {
	static int standard_initialisation_done;
	
	if ( standard_initialisation_done==1)return NULL;
	else standard_initialisation_done=1;
		
/*Standard exit*/
      atexit (clean_function);	
      signal (SIGSEGV ,error_function);
      signal (SIGFPE , error_function);
      signal (SIGILL , error_function);
#ifndef WIN32
      signal (SIGBUS , error_function);
#endif
 
/*Process the parameters*/
      program_name=vcalloc ( 1000, sizeof (char));
      if (in_argv)
	  {
	    sprintf ( program_name, "%s", in_argv[0]);
	    in_argv=break_list ( in_argv, in_argc, "=;, \n");
	  }
      else
	{
	 sprintf ( program_name, "%s",PROGRAM); 
	}
     
/*Initialise the time*/
      get_time();	

      return in_argv;
      }
int   count_n_res_in_array  (char *array, int len)
      {
	return count_n_symbol_in_array(array, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", len);
      }
int   count_n_gap_in_array  (char *array, int len)
      {
	int l;
	if ( len<=0 ||len>strlen(array) )l=strlen(array);
	else l=len;
	
	return l- count_n_res_in_array (array,len);
      }
int count_n_symbol_in_array ( char *array, char *array_list, int len)
      {
	int a=0, t=0;
	int l;
	
	if ( len<=0 ||len>strlen(array) )l=strlen(array);
	else l=len;

	for ( a=0; a< l; a++)t+=is_in_set (array[a], array_list);
	return t;
      }
	
int   count_n_char_x_in_file(char *name, char x)
      {
      FILE *fp;
      int n, c;

      n=0;
      fp=vfopen(name, "r");
      while ( (c=fgetc(fp))!=EOF)n+=(c==x);
      vfclose (fp);
      return n;
      }
int   count_n_char_in_file(char *name)
      {
      int  c, n;
      FILE *fp;

      n=0;
      fp=vfopen(name, "r");
      while ( (c=fgetc(fp))!=EOF)n+=!isspace(c);
      vfclose (fp);
      return n;
      }
int count_n_line_in_file ( char *name )
     {
      int  c, n;
      FILE *fp;

      n=0;
      fp=vfopen(name, "r");
      while ( (c=fgetc(fp))!=EOF)n+=(c=='\n');
      vfclose (fp);
      return n;
      }
int measure_longest_line_in_file ( char *name )
     {
      int  c;
      FILE *fp;
      int longest=0;
      int current=0;


      fp=vfopen(name, "r");
      while ( (c=fgetc(fp))!=EOF)
	  {
	      if ( c=='\n'){longest=MAX(longest, current+1);current=0;}
	      else current++;
	  }
      
      vfclose (fp);
      return longest;
      }
char *input_name ()
	{
	char *string;
	int a;
	char ch;
	
	string= vcalloc ( 500, sizeof ( char));
	
	a=0;
    	while ( ( ch=getchar())!='\n')
    			string[a++]=ch;
    	string[a]='\0'; 
    	
    	if ( string[0]=='\0')
    		{
    		vfree (string);
    		return NULL;
    		}
    	else
    		return string;
    	}

FILE * vtmpfile()
    {
        return tmpfile();
    }

static char **tmpnam_list;
static int n_tmpnam;
static int max_tmpnam;
char *vtmpnam ( char *s)
     {
	 int name_size, a;
	 char **new_tmpnam_list;

	 standard_initialisation(NULL, NULL);

	 name_size=MAX( 2*L_tmpnam, MAXNAMES*2)+1;
	 
	 if ( s==NULL)
	     {
	     s=vcalloc (name_size, sizeof (char));
	     }
	 if ( !tmpnam_list)
	    {
		max_tmpnam=100;
		tmpnam_list=declare_char (max_tmpnam, name_size);
	    }

	 if ( max_tmpnam==n_tmpnam)
	   {
	       new_tmpnam_list=declare_char (max_tmpnam+100, name_size); 
	       for ( a=0; a< max_tmpnam; a++)sprintf ( new_tmpnam_list[a], "%s", tmpnam_list[a]);	       
	       max_tmpnam+=100;
	       
	       free_char ( tmpnam_list, -1);
	       tmpnam_list=new_tmpnam_list;
	       
	   }
	 
	 if ( n_tmpnam>TMP_MAX)
	    {
		crash ("The Function tmpnam cannot generate any more unique file [FATAL]\n");
	    }
	 s=tmpnam(s);
	 
	 sprintf ( tmpnam_list[n_tmpnam++], "%s", s);
	 return s;
     }

void error_function ( )
     {
	 
         fprintf( stderr, "\nAbnormal Program Termination:[%s, %s]",PROGRAM, VERSION);
	 fprintf( stderr, "\nPlease report the fault to: %s (Indicate the version number)\nThank you for your cooperation :-)\n", MAIL);	
	 clean_function();
	 exit (EXIT_FAILURE);
     }
void clean_function ( )
     {
	 int a;
	 
	 
	 
	 if ( getenv ( "CLEAN_TCOFFEE") && strm (getenv ( "CLEAN_TCOFFEE"), "NO"))return;
	 else if ( !tmpnam_list)return;
	 else 
	     {
		 for (a=0; a< n_tmpnam; a++)
		     {
			 
			 remove ( tmpnam_list[a]);
		     }
		 free_char ( tmpnam_list, -1);
		 tmpnam_list=NULL;		 
	     }
#ifdef __WIN32__
	 remove ( TO_NULL_DEVICE);
#endif
     return;
     }

FILE *NFP;/*Null file pointer: should only be open once*/
FILE * vfopen  ( char *name, char *mode)
    {
    FILE *fp;
    int get_new_name;
    int tolerate_mistake;    
    FILE *tmp_fp;
    char *tmp_file;
    int c;

    get_new_name=tolerate_mistake=0;    
    if ( mode[0]=='g'){get_new_name=1; mode++;}
    else if ( mode[0]=='t'){tolerate_mistake=1;mode++;}

  

    if (name==NULL ||strm5 ( name, "no","NO","No","NULL","/dev/null") || strm2 (name, "no_file", "NO_FILE"))
    		{       
 		if ( NFP==NULL)NFP=fopen (NULL_DEVICE, mode);
 		return NFP;
 		} 
    else if ( strm (name, "stderr"))return stderr; 
    else if ( strm ( name, "stdout"))return stdout;
    else if ( strm ( name, "stdin"))
	{
	    tmp_file=vtmpnam (NULL);
	    tmp_fp=vfopen ( tmp_file, "w");
	    while ( (c=fgetc(stdin))!=EOF)fprintf (tmp_fp, "%c", c);
	    vfclose ( tmp_fp);
	    sprintf ( name, "%s", tmp_file);
	    vfree ( tmp_file);
	    return vfopen (name, "r");
	}
	
    else if ( strm (name, "") && strm (mode, "w"))return stdout;
    else if ( strm (name, "") && strm (mode, "r"))return stdin;    
    else if ( (fp= fopen ( name, mode))==NULL)
    		{
    		if ( strcmp (mode, "r")==0)
    			{
    			fprintf (stderr, "\nCOULD NOT READ %s\n", name);
    			if ( get_new_name){fprintf ( stderr, "\nNew name: ");return vfopen (input_name(), mode-1);}
    			else if ( tolerate_mistake)return NULL;
			else 
			    {	
			    fprintf (stderr, "\nFORCED EXIT (NON INTERACTIVE MODE)\n");
			    exit (1);
			    }
			}
    		else
    			{
    			fprintf (stderr, "\nCANNOT WRITE %s\n", name);
    			if ( get_new_name==1){fprintf ( stderr, "\nNew name: ");return vfopen (input_name(), mode-1);}	
    			else if ( tolerate_mistake)return NULL;
			else 
			    {
			    fprintf (stderr, "\nFORCED EXIT (NON INTERACTIVE MODE): %s %s\n", (strcmp ( mode, "r")==0)?"READ":"WRITE", name);
			    exit(1);
			    }
			}		
    		} 
    else
	return fp;
    
    return NULL;
    }

FILE * vfclose ( FILE *fp)
       {
       if ( fp==NFP)return NULL;
       if ( fp==stdout)return stdout;
       if ( fp==stderr)return stderr;
       if ( fp==stdin) return stdin;
       if ( fp==NULL)return NULL;
       else fclose (fp);
       return NULL;
       }
	
int echo ( char *string, char *fname)
{
int a;
/*
description:
prints the content of string into file fname

in:
string= string to print
fname =name of the file to create
*/

FILE *fp;

    fp=vfopen ( fname, "w");
    fprintf (fp, "%s", string);
    a=fclose (fp);
    return a;

}

int get_cl_param (int argc, char **argv, FILE **fp,char *para_name, int *set_flag, char *type, int optional, int max_n_val,char *usage, ...)
        {
	/*
	usage:
	argc:       n_ arg
	argv        list *
	para_name   param
	set_flag    set to 1 if param set;
	para_type   F, I, S, R_FN (read_file, name), W_FN (written file, name), R_FP (pointer)
	max_n_val   maximum number of values;
	optional    1 for yes, 0 for no
	usage       usage list with optional value;
	val         pointer to the varaible holding the value(s)
	default1    default value (if value id not there)
	default2    default value if the flag is there but no value set ("")indicates an error
	range_left  min value ( "any" for any);
	range_right  max_value ( "any" for any);
	*/
	int pos;
	int a;
	va_list ap;
	
	int   *int_val;
	float *float_val;
	char  **string_val;

	
	char  *range_right;
	char  *range_left;
	
	
	char *default_value1;
	char *default_value2;
	int n_para=0;        
	double max, min;
	
	static char **parameter_list;
	static int    number_of_parameters;

	char **para_name_list;
	int    n_para_name;
	
	char **para_val;
	int    n_para_val;
	
	char **pv_l=NULL;
	int    n_pv_l;
	char **pv_r=NULL;
	int    n_pv_r;
	char   value[STRING];



/*CHECK THAT ALL THE PARAM IN ARG EXIST*/
	if ( para_name==NULL)
	   {
	   for ( a=1; a< argc; a++)
	       {
	       if ( is_parameter ( argv[a]))
		   if ( name_is_in_list ( argv[a], parameter_list, number_of_parameters, STRING)==-1)
		      {
		      fprintf ( stderr, "\n%s IS NOT A PARAMETER  OF %s [FATAL/%s]\n",argv[a], argv[0], argv[0]);
		      exit(1);
		      }
	       }
	   free_char (parameter_list,-1);
	   return 0;
	   }

	if ( parameter_list==NULL)parameter_list=declare_char(MAX_N_PARAM,STRING);
	para_name_list=get_list_of_tokens(para_name,NULL, &n_para_name);
	for ( a=0; a< n_para_name; a++)
	    {
	    sprintf ( parameter_list[number_of_parameters++],"%s", para_name_list[a]);
	    }
	free_char(para_name_list,-1);
        
	

	
	
	set_flag[0]=0;
	va_start (ap, usage);
	
	if (strm3 (type, "S","R_F","W_F"))
		string_val=va_arg(ap, char**);
	else if (strm2 (type, "D","FL"))
	        int_val=va_arg(ap, int*);
	else if (strm (type, "F"))
	        float_val=va_arg(ap, float*);
	else 
	    exit (1);
	
	
	
	default_value1=va_arg(ap, char*);
	default_value2=va_arg(ap, char*);
	range_left    =va_arg(ap, char*);
	range_right   =va_arg(ap, char*);
       	va_end(ap);


	para_name_list=get_list_of_tokens(para_name, NULL, &n_para_name);
	for ( a=0; a<n_para_name; a++) 
	    {
	    if ( (pos=name_is_in_list(para_name_list[a], argv,argc, STRING))!=-1)break;
	    }
	free_char (para_name_list,-1);

	
	if ( name_is_in_list("-h"      , argv,argc,  STRING)  !=-1  ||\
	     name_is_in_list("h"       , argv,argc,  STRING)  !=-1  ||\
	     name_is_in_list("-help"   , argv,argc  ,STRING)  !=-1  ||\
	     name_is_in_list("-options", argv,argc  ,STRING)  !=-1  ||\
	     name_is_in_list("help"    , argv,argc  ,STRING)  !=-1  ||\
	     name_is_in_list("Help"    , argv,argc  ,STRING)  !=-1  ||\
	     name_is_in_list("HELP"    , argv,argc  ,STRING)  !=-1  ||\
	     (name_is_in_list("-man"  , argv,argc  ,STRING)   !=-1 && name_is_in_list(para_name, argv,argc, STRING)!=-1)||\
	     (name_is_in_list("-info"  , argv,argc  ,STRING)  !=-1 && name_is_in_list(para_name, argv,argc, STRING)!=-1)||\
	     (name_is_in_list("-Info"  , argv,argc  ,STRING)  !=-1 && name_is_in_list(para_name, argv,argc, STRING)!=-1)||\
	     (name_is_in_list("-INFO"  , argv,argc  ,STRING)  !=-1 && name_is_in_list(para_name, argv,argc, STRING)!=-1) )
	
	   {
	       fprintf ( stderr, "PARAMETER: %s\n",  para_name);
	       fprintf ( stderr, "USAGE    : %s\n",       usage);
	       fprintf ( stderr, "DEFAULT  : %s OR %s (when flag set)\n", default_value1, default_value2);
	       fprintf ( stderr, "RANGE    : [%s]...[%s]\n", range_left,(strm(range_right,"any"))?"":range_right);
	       fprintf ( stderr, "TYPE     : %s\n\n", type);
	       return 0;
	   }
	else if( \
	     name_is_in_list("-man"  , argv,argc  ,STRING)   !=-1 ||\
	     name_is_in_list("-info"  , argv,argc  ,STRING)  !=-1 ||\
	     name_is_in_list("-Info"  , argv,argc  ,STRING)  !=-1 ||\
	     name_is_in_list("-INFO"  , argv,argc  ,STRING)  !=-1)
		 {
		     return 0;
		 }
	else if ( strm4 ( para_name,"-h","h","-man","man") || strm4( para_name, "-help", "help", "Help", "HELP" ) ||strm3( para_name, "-info", "INFO", "Info" ) )
	   {
	   exit (1);
	   }
	else if (para_name[0]!='-')
	   {
	   fprintf ( stderr, "\nWRONG PARAMETER DEFINITION %s Must Start with a dash", para_name);
	   exit (1);
	   }
     	else if (pos==-1)
	   {
	   if ( optional==OPTIONAL)
	      {
	      set_flag[0]=0;	      
	      para_val=get_list_of_tokens(default_value1, NULL, &n_para_val);
	      
	      for (n_para=0; n_para<n_para_val && !strm (default_value1, "NULL"); n_para++) 
	          {
		  if ( strm (para_val[n_para], ""))
		      {
		      set_flag[0]=0;
		      break;
		      }
		  else if ( strm (type, "FL"))
	              {
		      set_flag[0]=atoi(para_val[n_para]);
		      break;
		      }
		  else if (strm3 (type, "S", "R_F","W_F"))
		      {
		     
		      sprintf ( string_val[n_para], "%s",para_val[n_para]);
		      }		  
		  else if ( strm (type, "D"))
		      int_val[n_para]=atoi(para_val[n_para]);
		  else if ( strm (type, "F"))
		      float_val[n_para]=atof(para_val[n_para]);
		  }
	      free_char (para_val, -1);
               
	      if (n_para==0 && strm3(type, "S","W_F","R_F") && strm (default_value1, "NULL"))
		  {
		  vfree (string_val[0]);
		  string_val[0]=NULL;
		  
		  }
	      else if (n_para==0 && strm (type, "D") && strm (default_value1, "NULL"))int_val[0]=0;		   
	      else if (n_para==0 && strm (type, "F") && strm (default_value1, "NULL"))float_val[0]=0;
	 
	      }
	   else
	      {
	      fprintf ( stderr, "\nParameter %s is not optional",para_name);
	      exit (1);
	      }	   
	   }
	else if (pos!=-1)
	  {
	  set_flag[0]=1;
	  for (a=pos+1; a< argc; a++)
	      {
	      if ( is_parameter(argv[a]))break;
	      else
	          {
		  if ( n_para>=max_n_val)
		     {
		     fprintf ( stderr, "\nToo Many Values for %s (max=%d)[FATAL/%s]\n", para_name, max_n_val, argv[0]);
		     exit (1);
		     }
		  if ( !strm ( argv[a], "NULL"))		      	  
		      if ( strm3(type, "S", "R_F", "W_F"))
			  {
			  sprintf ( string_val[n_para],"%s", argv[a]);
			  }
		       else if (strm (type, "D"))
		          {
			  int_val[n_para]=atoi(argv[a]);
			  }
		       else if (strm ( type,"F"))
		          {
			  float_val[n_para]=atof(argv[a]);
			  }
		  n_para++;
		  }
	      }

	  if ( n_para==0 && !strm2(default_value2,"","NULL") && !strm(type, "FL"))
	      {
	      para_val=get_list_of_tokens(default_value2, NULL, &n_para_val);
	      for ( n_para=0; n_para<n_para_val; n_para++)
	          {
		  if ( strm3(type, "S", "R_F", "W_F"))sprintf ( string_val[n_para],"%s", para_val[n_para]);		  
		  else if (strm (type, "D"))int_val  [n_para]=atoi(para_val[n_para]);
		  else if (strm ( type,"F"))float_val[n_para]=atof(para_val[n_para]);		
		  }
	      free_char (para_val,-1);
	      }
	  else if (n_para==0 && strm (type, "FL"));
	  else if (n_para==0 && strm3(type, "S","W_F","R_F") && strm (default_value2, "NULL")){vfree (string_val[0]);string_val[0]=NULL;}
	  else if (n_para==0 && strm (type, "D") && strm (default_value2, "NULL"))int_val[0]=0;		   
	  else if (n_para==0 && strm (type, "F") && strm (default_value2, "NULL"))float_val[0]=0;
	  else if (n_para==0 && strm (default_value2, ""))
	         {
		 fprintf ( stderr, "\nParam %s needs a value [FATAL/%s]", para_name, argv[0]);
		 exit(1);
		 }
	  else;
	  }
	  
/*Check That The Parameters are in the Good Range*/
	
	pv_l=get_list_of_tokens( range_left , NULL, &n_pv_l);
	pv_r=get_list_of_tokens( range_right, NULL, &n_pv_r);

	for ( a=0; a< n_para; a++)
	    {	   
	    if ( strm (type, "R_F") && !check_file_exists(string_val[a]) && !check_file_exists(string_val[a]+1))
			{			
			fprintf ( stderr, "PARAM %s: File %s does not exist [FATAL/%s]\n",para_name,string_val[a], argv[0]);
			exit(1);
		        }	    
	    else if ( strm (pv_l[0], "any"));
	    else if ( strm (type, "D"))
	         {
		 if ( n_pv_l==1)
		    {
		    min=(double)atoi(pv_l[0]);
		    max=(double)atoi(pv_r[0]);
		    if ( int_val[a]<min || int_val[a]>max)
		       {
		       fprintf ( stderr, "\n%s out of range [%d %d] [FATAL/%s]\n", para_name, (int)min, (int)max,argv[0]);
		       exit (1);
		       }
		    }
		 else
		    {
		    sprintf ( value, "%d", int_val[a]);
		    if ( name_is_in_list(value, pv_l, n_pv_l, STRING)==-1)
			fprintf ( stderr, "\n%s out of range [%s: ", para_name, value);
		    print_array_char (stderr, pv_l, n_pv_l, " ");
		    fprintf ( stderr, "\n");
		    exit(1);
		    }
		 }
	    else if ( strm (type, "F"))
	         {
		  if ( n_pv_l==1)
		    {
		    min=(double)atof(range_left);
		    max=(double)atof(range_right);
		    if ( float_val[a]<min || float_val[a]>max)
		       {
		       fprintf ( stderr, "\n%s out of range [%f %f] [FATAL/%s]\n", para_name, (float)min, (float)max,argv[0]);
		       exit (1);
		       }
		     }
		  else
		     {
		     sprintf ( value, "%f", float_val[a]);
		     if ( name_is_in_list(value, pv_l, n_pv_l, STRING)==-1)
			fprintf ( stderr, "\n%s out of range [%s: ", para_name, value);
		     print_array_char (stderr, pv_l, n_pv_l, " ");
		     fprintf ( stderr, "\n");
		     exit(1);
		     }
		 }
	    }
         

	 if ( fp[0]!=NULL)
	      {
	      fprintf (fp[0], "%-15s\t%s\t[%d] ", para_name, type, set_flag[0]);
	      for (a=0; a<n_para; a++) 
		  {
		  if ( strm3 ( type, "S", "R_F", "W_F"))fprintf ( fp[0], "\t%s", string_val[a]);
		  else if ( strm  ( type, "D"))fprintf ( fp[0], "\t%d ", int_val[a]);
		  else if ( strm  ( type, "F"))fprintf ( fp[0], "\t%f ", float_val[a]);
	         }	    
	      if ( strm (type, "FL"))fprintf ( fp[0], "\t%d", int_val[0]);
	      fprintf ( fp[0], "\n");     
	      }

	free_char ( pv_l, -1);
	free_char ( pv_r, -1);
	return n_para;
	}

	
	
char ** get_parameter ( char *para_name, int *np, char *fname)
{
    /*
    In:
    para_name: the name of the parameter to look for
    fname: the name of the file containing the parameters
    np[0]: set to 0

    Out:
    char ** containing the np[0] values taken by para_name in fname.
    
    Special:
    if fname=NULL, para_name is searched using the last value taken by fp.
    
    Note: by default, the function keeps a file handle open until the first unsuccessful call.
    */

    static FILE *fp;
    static char *line;
    char ** return_value;
    
    if ( line==NULL)line=vcalloc ( VERY_LONG_STRING+1, sizeof (char));
    if ( fname!=NULL && fp!=NULL)fclose (fp);

    np[0]=0;

    if ((fp=find_token_in_file ( fname,(fname==NULL)?fp:NULL, para_name))==NULL) 
	{
	     return NULL;
	}
    else
        {
	fgets ( line, VERY_LONG_STRING,fp);
        return_value=get_list_of_tokens ( line, NULL, np);	
	return return_value;
	}
}

FILE * set_fp_id ( FILE *fp, char *id)
	{
/*Sets fp just after id, id needs to be at the begining of the line*/
	char string[10000];
	int cont=1;
	int c;
	
	while ( cont==1)	
		{
		c=fgetc(fp);
		if ( c!=EOF)
			{
			
			ungetc(c, fp);
			fscanf ( fp, "%s", string);
			
			if ( strcmp ( string, id)==0)
				return fp;
			else while ( c!='\n' && c!=EOF)
				c=fgetc(fp);			
			}
		else if ( c==EOF)
			{
			fclose ( fp);
			return NULL;
			}
		}	
	return fp;
	}
FILE * set_fp_after_char ( FILE *fp, char x)
	{
/*sets fp just after the first occurence of x*/
	


	int cont=1;
	int c;
	
	while ( cont==1)	
		{
		c=fgetc(fp);
		if ( c!=EOF)
			{
			if ( c==x)
				return fp;
			
			}
		else if ( c==EOF)
			{
			fclose ( fp);
			return NULL;
			}
		}	
	return fp;
	}
	
FILE * find_token_in_file_nlines ( char *fname, FILE * fp, char *token, int n_line)

        {
	  /*This function finds the string TOKEN (as a single word) in the n_line first lines of fname.
	    It returns NULL if not found or the position of the fp
	  */
	  
	  static char *tmp_name;
	  FILE *fp1;
	  FILE *fp2;
	  char buffer[10001];
	  int a;
	  if ( !fp && !check_file_exists(fname))return NULL;
	  if ( fp==NULL)
	    {
	      if ( tmp_name)
		{
		  remove(tmp_name);
		  free(tmp_name);
		}
	      tmp_name=vtmpnam ( NULL);
	
	      fp1=vfopen (fname, "r");
	      fp2=vfopen (tmp_name, "w");
	      
	      for ( a=0; a< n_line && fgets(buffer, 10000, fp1)!=NULL; a++)
		{
		  fprintf ( fp2, "%s", buffer);
		}
	      vfclose (fp1);
	      vfclose (fp2);
	    }
	  return find_token_in_file ( tmp_name,fp,token);
	}
	  

	  
FILE * find_token_in_file ( char *fname, FILE * fp, char *token)
	{
	int c;
	char *name;
	int token_len;

	int only_start;
	
	/*Note: Token: any string
	        If Token[0]=='\n' Then Token only from the beginning of the line
	*/

	if (!fp && !check_file_exists(fname))return NULL;

	if ( token[0]=='\n'){token++;only_start=1;}
	else only_start=0;

	token_len=strlen (token);


	name = vcalloc (((fname)?measure_longest_line_in_file (fname):10000)+1, sizeof (char));
	if (!fp)
	  {
	    fp=vfopen ( fname, "r");
	  }
	    
	while ( (fscanf ( fp, "%s", name))!=EOF)
		{
                
		if ( name[0]=='*')while ( ((c=fgetc (fp))!='\n')&& (c!=EOF));
		else if (strncmp ( name, token,token_len)==0){vfree (name);return fp;}
		else if (only_start) while ( ((c=fgetc (fp))!='\n')&& (c!=EOF));
		}
	vfree (name);
	fclose ( fp);
	return NULL;
	}
int **get_file_block_pattern (char *fname, int *n_blocks, int max_n_line)
        {
	int c;
	FILE *fp;
	char *line;
	int lline;
	int **l;
	int in_block;
	
	int max_block_size;
	int block_size;
	int x;
	int n_line;
	
	lline=measure_longest_line_in_file (fname)+1;			
	line=vcalloc ( sizeof (char),lline+1);
	
	fp=vfopen (fname, "r");
	max_block_size=block_size=0;
	in_block=1;
	n_blocks[0]=0;
	n_line=0;
	while ((c=fgetc(fp))!=EOF && (n_line<max_n_line || !max_n_line))
	    {
		  ungetc (c, fp);
		  fgets ( line, lline,fp);
		  n_line++;

		  if ( is_alnum_line (line) && !in_block){n_blocks[0]++;in_block=1;}
		  if ( is_alnum_line (line))
		     {
		     block_size++;
		     }
		  else
		      {
		      in_block=0;
		      max_block_size=MAX( max_block_size, block_size);
		      block_size=0;
		      }
	    }

	
	max_block_size=MAX( max_block_size, block_size);
	vfclose ( fp);
	
	l=declare_int (n_blocks[0]+1,max_block_size+1);  

	
	fp=vfopen (fname, "r");	
	in_block=1;
	n_blocks[0]=0;
	n_line=0;
	while ((c=fgetc(fp))!=EOF && (n_line<max_n_line || !(max_n_line)))
	      {
		  ungetc (c, fp);
		  fgets ( line, lline,fp);
		  n_line++;

		  if ( is_alnum_line (line) && !in_block){n_blocks[0]++;in_block=1;}
		  if ( is_alnum_line (line))
		      {
		      l[n_blocks[0]][0]++;
		      free_char (get_list_of_tokens (line, " \t\n*.:,", &x), -1); 
		      
		      if ( l[n_blocks[0]][0]> max_block_size)fprintf ( stderr, "\nERROR %d", l[n_blocks[0]][0]);

		      l[n_blocks[0]] [l[n_blocks[0]][0]]=x;
		      }
		  else
		      {
			  in_block=0;
		      }
	      }
	n_blocks[0]++;
	vfree(line);
	vfclose (fp);
	return l;
	}		  
		  
int check_file_exists ( char *fname)
	{
	FILE *fp;
	
	if ( strm5 (fname, "default", "stdin", "stdout","stderr", "/dev/null"))return 1;
	if ( strm5 (fname, "no", "NO", "No", "NO_FILE","no_file"))return 1;
	fp=fopen ( fname, "r");
	if ( fp==NULL) return 0;
	else
		vfclose (fp);
	return 1;
	}	


void create_file ( char *name)
	{
	FILE *fp;
	
	fp=fopen (name, "w");
	fclose (fp);
	}
void delete_file ( char *fname)
	{
	char command[1000];
	FILE * fp;
	
	fp=fopen ( fname, "w");
	fprintf ( fp, "x");
    	fclose ( fp);
	
	sprintf ( command, "rm %s", fname);
    	system ( command);
    		
    	}	

int util_rename ( char *from, char *to)
        {
	FILE *fp_from;
	FILE *fp_to;
	int c;

	
	if ( check_file_exists (from)==0)return 0;
        else if ( check_file_exists (to)==1 && remove (to)==0 && rename ( from, to)==0 );
        else
                {
		 
	        fp_from=vfopen ( from, "r");
		fp_to=vfopen ( to, "w");
		
		while ( (c=fgetc (fp_from))!=EOF)fprintf ( fp_to, "%c", c);

		fclose (fp_from);
		fclose ( fp_to);

		remove ( from);
		return 1;
		}
	return 0;
	}


int util_copy (  char *from, char *to)
        {
	FILE *fp_from;
	FILE *fp_to;
	int c;

	
	if ( check_file_exists (from)==0)return 0;
        else
                {
		 
	        fp_from=vfopen ( from, "r");
		fp_to=vfopen ( to, "w");
		
		while ( (c=fgetc (fp_from))!=EOF)fprintf ( fp_to, "%c", c);

		fclose (fp_from);
		fclose ( fp_to);
		return 1;
		}
	return 0;
	} 
FILE * output_completion ( FILE *fp,int n, int tot, int n_reports)
        {

	static int ref_val, flag;

        n++;
	if ( n==1)ref_val=flag=0;

	if ( !ref_val && !flag)
	   {
	   fprintf (fp, "\n\t\t[TOT=%5d][%3d %%]",tot,(tot==1)?100:0);
	   flag=1;
	   }
	else if ( n==tot)fprintf (fp, "\r\t\t[TOT=%5d][100 %%]", tot);
	else if ( ((n*100)/tot)>ref_val)
	      {
	      ref_val=((n*100)/tot);
	      fprintf (fp, "\r\t\t[TOT=%5d][%3d %%]", tot,ref_val);
	      flag=0;
	      }	   
	return fp;
	}
	   
	
void * null_function (int a,...)
{
  fprintf ( stderr, "\n[ERROR] Attempt to use the Null Function [FATAL:%s]", PROGRAM);
  crash ("");
  return NULL;
}

int  btoi ( int nc,...)
{
  va_list ap;
  int a, b;
  va_start (ap, nc);
  for ( a=0, b=0; a< nc; a++)
    {
      b+=pow(2,a)*va_arg (ap,int);
    }
  va_end(ap);
  return b;
}

