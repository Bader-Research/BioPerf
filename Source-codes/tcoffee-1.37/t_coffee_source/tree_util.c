#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"


int distance_tree;
int rooted_tree;
int tot_nseq;


int ** make_sub_tree_list ( NT_node **T, int nseq, int n_node)
        {


/*This function produces a list of all the sub trees*/
	 

/*  	  /A */
/*  	-* */
/*  	  \      /B */
/*  	   \    / */
/*  	    ---* */
/*                  \ */
/*                   *--C */
/*                    \ */
/*                     \D */
	
/*  	Contains 4 i_nodes */
/*  	       8   nodes (internal nodes +leaves) */
/*  	       8 sub trees: */
/*  	       ABCD */
/*                 1111 */
/*  	       0111 */
/*  	       1000 */
/*  	       0100 */
/*  	       0011 */
/*  	       0001 */
/*  	       0010 */
       
	int **sub_tree_list;
	int a, n=0;
	
	
	if (T)
	     {
		 sub_tree_list=declare_int ( (n_node), nseq);
		 make_all_sub_tree_list (T[3][0],sub_tree_list, &n);

	     }
	else
	     {
		 sub_tree_list=declare_int (nseq, nseq);
		 for ( a=0; a< nseq; a++)sub_tree_list[a][a]=1;
	     }
	
	return sub_tree_list;
	}

void make_all_sub_tree_list ( NT_node N, int **list, int *n)
        {
	make_one_sub_tree_list (N, list[n[0]++]);
	if (!N->leaf)
	    {
	    make_all_sub_tree_list (N->left , list, n);
	    make_all_sub_tree_list (N->right, list, n);
	    }
	return;
	}

void make_one_sub_tree_list ( NT_node T,int *list)
        {
	if (T->leaf)
	    {
	    
	    list[T->seq]=1;
	    }
	else
	    {
	    make_one_sub_tree_list(T->left , list);
	    make_one_sub_tree_list(T->right, list);
	    }
	return;
	}

NT_node** read_tree(char *treefile, int *tot_node,int nseq, char **seq_names)
	{

	    /*The Tree Root is in the TREE[3][0]...*/
	    /*TREE[0][ntot]--> pointer to each node and leave*/
  	char ch;
	int  a,b;

  	FILE *fp;
  	int nseq_read = 0;
  	int nnodes = 0;/*Number of Internal Nodes*/
  	int ntotal = 0;/*Number of Internal Nodes + Number of Leaves*/
  	int flag;
  	int c_seq;

  	



	NT_node **lu_ptr;
	NT_node seq_tree, root,p;
	
	tot_nseq=nseq;
	rooted_tree=distance_tree=TRUE;
	
  	fp = vfopen(treefile, "r");
	fp=skip_space(fp);
  	ch = (char)getc(fp);
  	if (ch != '(')
    		{
    	  	fprintf(stdout, "Error: Wrong format in tree file %s\n", treefile);
      		exit (1);
    		}
  	rewind(fp);
  
  	
  	lu_ptr=(NT_node **)vcalloc(4,sizeof(NT_node*));
  	     lu_ptr[0] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
  	     lu_ptr[1] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
  	     lu_ptr[2] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
  	lu_ptr[3] =(NT_node *) vcalloc(1,sizeof(NT_node));
  	
  	seq_tree =(NT_node) declare_tree_node();
  	
  	set_info(seq_tree, NULL, 0, "  ", 0.0);
	
	fp=create_tree(seq_tree,NULL,&nseq_read, &ntotal, &nnodes, lu_ptr, fp);
  	
  	fclose (fp);
  	if (nseq != tot_nseq)
     		{
        	fprintf(stderr," Error: tree not compatible with alignment (%d sequences in alignment and %d in tree\n", nseq,nseq_read);
         	exit (1);
     		}


  	if (distance_tree == FALSE)
     		{
  		if (rooted_tree == FALSE)
          		{
       	     		fprintf(stderr,"Error: input tree is unrooted and has no distances, cannot align sequences\n");
             		exit (1);
          		}
          	}
        
        if (rooted_tree == FALSE)
     		{
     		root = reroot(seq_tree, nseq,ntotal,nnodes, lu_ptr);
     		lu_ptr[1][nnodes++]=lu_ptr[0][ntotal++]=root;
     		}
  	else
     		{
        	root = seq_tree;
     		}    	
  	lu_ptr[3][0]=root;
  	tot_node[0]=nnodes;
	
	
	
  	for ( a=0; a< ntotal; a++)
  		{
  		(lu_ptr[0][a])->isseq=(lu_ptr[0][a])->leaf;
  		(lu_ptr[0][a])->dp=(lu_ptr[0][a])->dist;
  		}
  		
  	
  	for ( a=0; a< nseq; a++)
  		{
  		for ( flag=0,b=0; b<nseq; b++)
  			{
  			if ( strncmp ( (lu_ptr[2][a])->name, seq_names[b], MAXNAMES)==0)
  				{
  				flag=1;
  				
  				(lu_ptr[2][a])->order=(lu_ptr[2][a])->seq=b;
  				free ( (lu_ptr[2][a])->name);
  				(lu_ptr[2][a])->name= seq_names[b];
  				}
  			}
  		if ( flag==0  && (lu_ptr[0][a])->leaf==1)
  			{
  			fprintf ( stderr, "\n%s* not in tree",(lu_ptr[2][a])->name);
  			for ( a=0; a< ntotal; a++)
  				{
  			        fprintf ( stderr, "\n%d %s",(lu_ptr[2][a])->leaf, (lu_ptr[2][a])->name);
  				} 
  			} 
  		}
  	
  	for ( a=0; a< nseq; a++)
  		{
  		p=lu_ptr[2][a];
  		c_seq=p->seq;
  		
  		while ( p->parent!=NULL)
  			{
  			p->lseq[p->nseq]=c_seq;
  			p->nseq++;
  			p=p->parent;
  			}
  		}
  	
  	return lu_ptr;   	
 	}


FILE * create_tree(NT_node ptree, NT_node parent,int *nseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp)
	{

  	int i, type;
  	float dist=0;
  	char *name;
  	int ch;
	
	name=calloc ( MAXNAMES+1, sizeof (char));
	sprintf ( name, "   ");
  	fp=skip_space(fp);
  	ch = (char)getc(fp);
  	
  	if (ch == '(')
    		{
    		type = NODE;
      		name[0] = '\0';
      		lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
      		ntotal[0]++;
      		nnodes[0]++;
      		create_tree_node(ptree, parent);
      		fp=create_tree(ptree->left, ptree, nseq,ntotal,nnodes,lu,fp);
     		ch = (char)getc(fp);
                if ( ch == ',')
       			{
          		fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
			ch = (char)getc(fp);
          		if ( ch == ',')
            			{
            			
               			ptree = insert_tree_node(ptree);
               			lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
               			ntotal[0]++;
               			nnodes[0]++;
               			fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
               			rooted_tree = FALSE;
            			}
       			}

      		fp=skip_space(fp);
      		ch = (char)getc(fp);
    		}
 	else
    		{
     	 	type=LEAF;
      		lu[0][ntotal[0]] = lu[2][nseq[0]] = ptree;
      		ntotal[0]++;
      		nseq[0]++;
      		name[0] = ch;
      		i=1;
      		ch = (char)getc(fp);
      		while ((ch != ':') && (ch != ',') && (ch != ')'))
        		{
          		if (i < MAXNAMES) name[i++] = ch;
          		ch = (char)getc(fp);
        		}
      		name[i] = '\0';
      		if ( i>=(MAXNAMES+1))exit (1);
      		if (ch != ':')
         		{
           		distance_tree = FALSE;
         		}
         	}
	if (ch == ':')
     		{
       		fp=skip_space(fp);
       		fscanf(fp,"%f",&dist);
       		skip_space(fp);
       		}
       	else
       		{
       		ungetc ( ch, fp);
       		skip_space(fp);
       		}
        set_info(ptree, parent, type, name, dist);
  	
  	
  	free (name);
  	return fp;
  	}


NT_node declare_tree_node ()
	{
	NT_node p;
	
	
	p= (NT_node)vcalloc (1, sizeof ( Treenode)); 
	p->left = NULL;
   	p->right = NULL;
   	p->parent = NULL;
   	p->dist = 0.0;
   	p->leaf = 0;
   	p->order = 0;
   	p->name=(char*)vcalloc (MAXNAMES+1,sizeof (char));
   	p->name[0]='\0';
   	p->lseq=(int*)vcalloc ( tot_nseq, sizeof (int));
   	return p;
	
	}

void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist)
   	{
 	p->parent = parent;
   	p->leaf = pleaf;
   	p->dist = pdist;
   	p->order = 0;
   	
   	sprintf (p->name, "%s", pname);
   	if (pleaf ==1)
    	 	{
        	p->left = NULL;
        	p->right = NULL;
     		}

   }
NT_node insert_tree_node(NT_node pptr)
	{

   	NT_node newnode;

   	newnode = declare_tree_node();
   	create_tree_node(newnode, pptr->parent);

   	newnode->left = pptr;
   	pptr->parent = newnode;

   	set_info(newnode, pptr->parent, 0, "", 0.0);

   	return(newnode);
	}

void create_tree_node(NT_node pptr, NT_node parent)
	{
  	pptr->parent = parent;
  	pptr->left =declare_tree_node() ;
  	(pptr->left)->parent=pptr;
  	
  	pptr->right =declare_tree_node() ;    
	(pptr->right)->parent=pptr;
	}	

FILE * skip_space(FILE *fp)
	{
  	int   c;
  
  	do
     	c = getc(fp);
  	while(isspace(c));
	if ( c==EOF)
		{
		fprintf ( stderr, "\nEOF");
		exit (1);
		}
  	ungetc(c, fp);
  	return fp;
	}


NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu)
	{
	NT_node p, rootnode, rootptr;
 	float   diff, mindiff, mindepth = 1.0, maxdist;
	int   i;
	int first = TRUE;
	
	
	
	rootptr = ptree;
	
   	
   	for (i=0; i<ntotal; i++)
     		{
        	p = lu[0][i];
        	if (p->parent == NULL)
           		diff = calc_root_mean(p, &maxdist, nseq, lu);
        	else
           		diff = calc_mean(p, &maxdist, nseq, lu);

        	if ((diff == 0) || ((diff > 0) && (diff < 2 * p->dist)))
         		 {
              		if ((maxdist < mindepth) || (first == TRUE))
                 		{
                    		first = FALSE;
                    		rootptr = p;
                    		mindepth = maxdist;
                    		mindiff = diff;
                 		}
           		}

     		}
	if (rootptr == ptree)
     		{
        	mindiff = rootptr->left->dist + rootptr->right->dist;
        	rootptr = rootptr->right;
     		}
     		
   	rootnode = insert_root(rootptr, mindiff);
  	
   	diff = calc_root_mean(rootnode, &maxdist, nseq, lu);

   	return(rootnode);
	}


float calc_root_mean(NT_node root, float *maxdist, int nseq, NT_node **lu)
	{
   	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   	NT_node p;
   	int i;
   	int nl, nr;
   	int direction;
   	
   
   	dist = (*maxdist) = 0;
   	nl = nr = 0;
   	for (i=0; i< nseq; i++)
   	  {
          p = lu[2][i];
          dist = 0.0;
          while (p->parent != root)
         	{
               	dist += p->dist;
               	p = p->parent;
           	}
         if (p == root->left) direction = LEFT;
         else direction = RIGHT;
         dist += p->dist;

         if (direction == LEFT)
           	{
           	lsum += dist;
             	nl++;
           	}
         else
           	{
             	rsum += dist;
             	nr++;
           	}
        
        if (dist > (*maxdist)) *maxdist = dist;
     	}

      lmean = lsum / nl;
      rmean = rsum / nr;

      diff = lmean - rmean;
      return(diff);
      }

float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu)
	{
   	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   	NT_node p, *path2root;
   	float *dist2node;
   	int depth = 0, i,j , n;
   	int nl , nr;
   	int direction, found;
	
	
	path2root = (NT_node *)vcalloc(nseq,sizeof(Treenode));
	dist2node = (float *)vcalloc(nseq,sizeof(float));

   	depth = (*maxdist) = dist = 0;
   	nl = nr = 0;
   	p = nptr;
   	while (p != NULL)
     		{
         	path2root[depth] = p;
         	dist += p->dist;
         	dist2node[depth] = dist;
         	p = p->parent;
         	depth++;
     		}
 
/***************************************************************************
   *nl = *nr = 0;
   for each leaf, determine whether the leaf is left or right of the node.
   (RIGHT = descendant, LEFT = not descendant)
****************************************************************************/
   	for (i=0; i< nseq; i++)
     		{
       		p = lu[2][i];
       		if (p == nptr)
         		{
            		direction = RIGHT;
            		dist = 0.0;
         		}
       		else
         		{
         		direction = LEFT;
         		dist = 0.0;
			
        		found = FALSE;
         		n = 0;
         		while ((found == FALSE) && (p->parent != NULL))
           			{
               			for (j=0; j< depth; j++)
                 			if (p->parent == path2root[j])
                    				{ 
                      				found = TRUE;
                      				n = j;
                    				}
               			dist += p->dist;
               			p = p->parent;
           			}
         
         		if (p == nptr) direction = RIGHT;
        
			}
         	if (direction == LEFT)
           		{
             		lsum += dist;
             		lsum += dist2node[n-1];
             		nl++;
           		}
         	else
           		{
             		rsum += dist;
             		nr++;
           		}

        	if (dist > (*maxdist)) *maxdist = dist;
     		}

   free(dist2node);
   free(path2root);

   

   if ( nl==0 || nr==0)exit (1);	
   lmean = lsum / nl;
   rmean = rsum / nr;
   
   diff = lmean - rmean;
   return(diff);
}

NT_node insert_root(NT_node p, float diff)
{
   NT_node newp, prev, q, t;
   float dist, prevdist,td;

   
   newp = declare_tree_node();

   t = p->parent;
   prevdist = t->dist;

   p->parent = newp;

   dist = p->dist;

   p->dist = diff / 2;
   if (p->dist < 0.0) p->dist = 0.0;
   if (p->dist > dist) p->dist = dist;

   t->dist = dist - p->dist; 

   newp->left = t;
   newp->right = p;
   newp->parent = NULL;
   newp->dist = 0.0;
   newp->leaf = NODE;

   if (t->left == p) t->left = t->parent;
   else t->right = t->parent;

   prev = t;
   q = t->parent;

   t->parent = newp;

   while (q != NULL)
     {
        if (q->left == prev)
           {
              q->left = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->left;
           }
        else
           {
              q->right = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->right;
           }
    }

/*
   remove the old root node
*/
   q = prev;
   if (q->left == NULL)
      {
         dist = q->dist;
         q = q->right;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->right = NULL;
      }
   else
      {
         dist = q->dist;
         q = q->left;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->left = NULL;
      }

   return(newp);
}

