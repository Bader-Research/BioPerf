#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <vibrant.h>
#include <document.h>

#include "clustalw.h"
#include "xmenu.h"

static void build_profile(int prf_length,int first_seq,int last_seq,sint matrix[NUMRES][NUMRES],
					sint *weight,sint **profile);
static void calc_colscores(Boolean update_seqs,Boolean update_scores);
static void calc_panel_segment_exceptions(PaneL p);
static void calc_weights(PaneL p);
static void remove_short_segments(PaneL p);


extern Boolean aln_mode;
extern Boolean profile1_empty,profile2_empty;

extern int score_cutoff;    /* cutoff for residue exceptions */
extern int score_hwin;    /* half window for summing alignment column scores */
extern int score_scale;
extern int segment_dnascale;
extern int length_cutoff; /* cutoffs for segment exceptions */
extern Boolean residue_exceptions;
extern Boolean segment_exceptions;
extern int score_matnum;
extern char score_mtrxname[];
extern int segment_matnum;
extern char segment_mtrxname[];
extern int score_dnamatnum;
extern char score_dnamtrxname[];
extern int segment_dnamatnum;
extern char segment_dnamtrxname[];
extern IteM segment_item;

extern double   **tmat;

extern short score_matrix[];
extern short score_aa_xref[];
extern short segment_matrix[];
extern short segment_aa_xref[];
extern short score_dnamatrix[];
extern short score_dna_xref[];
extern short segment_dnamatrix[];
extern short segment_dna_xref[];

extern WindoW mainw;
extern FonT datafont;
extern short idmat[];
extern short   def_dna_xref[],def_aa_xref[];
extern short   swgapdnamt[],clustalvdnamt[];  /* used for alignment scores */
extern short   gon80mt[],gon120mt[],gon250mt[],gon350mt[];
extern Boolean  dnaflag;
extern sint     max_aa;
extern sint     *seqlen_array;
extern char     **seq_array;
extern sint     gap_pos1, gap_pos2;
extern sint	*output_index;
extern spanel  seq_panel;        /* data for multiple alignment area */
extern spanel  prf_panel[];       /* data for profile alignment areas */

extern PrompT residue_cutofftext;
extern PrompT length_cutofftext;
extern PrompT scorescaletext;
extern PrompT segmentdnascaletext;
extern PrompT scoremattext;
extern PrompT segmentmattext;
extern PrompT scorednamattext;
extern PrompT segmentdnamattext;
extern PopuP show_seg_toggle;

extern GrouP score_matrix_list,seg_matrix_list;
extern GrouP score_dnamatrix_list,seg_dnamatrix_list;

static Char filename[FILENAMELEN]; /* used in temporary file selection window */



void draw_colscores(PaneL p)
{ 
	RecT  block,r;
	int i, b, x, y;
	panel_data data;

	UseWindow(mainw);
	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	if(data.colscore == NULL) return;
	if(data.vcols<=0) return;

	ObjectRect (p, &r);
	InsetRect(&r,1,1);
	block.bottom=r.bottom;
	block.top=block.bottom-SCOREHEIGHT-1;
	block.left=r.left;
	block.right=r.right;
	data_colors();
	EraseRect(&block);

	Gray();
        r.left=block.left+data.charwidth;
	r.bottom=b=block.bottom;
	for(i=data.firstvcol;i<data.firstvcol+data.vcols && i<data.ncols;i++)
	{
		x=r.left;
		MoveTo(x,b);
		r.right=r.left+data.charwidth;
        	r.top=block.bottom-SCOREHEIGHT*(float)data.colscore[i]/100.0;
		PaintRect(&r);
		r.left+=data.charwidth;
	}
	black_on_white();
} 

 
void make_colscores(panel_data data)
{
/*	FILE *fd;*/
	int n,i,s,p,r,r1;
	short  *mat_xref, *matptr;
	float median,mean;
	float t,q1,q3,ul;
	float *seqdist,*sorteddist,diff;
	sint maxres;
	sint *seqvector;
	sint *freq,**profile;
        sint   matrix[NUMRES][NUMRES];
	Boolean include_gaps=FALSE;
	panel_data data1;

	if(dnaflag)
	{
		if (score_dnamatnum==1)
		{
			matptr = swgapdnamt;
			mat_xref = def_dna_xref;
		}
		else if (score_dnamatnum==2)
		{
			matptr = clustalvdnamt;
			mat_xref = def_dna_xref;
		}
		else
		{
			matptr = score_dnamatrix;
			mat_xref = score_dna_xref;
		}
	}
	else if (score_matnum==1)
	{
		matptr = idmat;
		mat_xref = def_aa_xref;
	}
	else if (score_matnum==2)
	{
		matptr = gon80mt;
		mat_xref = def_aa_xref;
	}
	else if (score_matnum==3)
	{
		matptr = gon120mt;
		mat_xref = def_aa_xref;
	}
	else if (score_matnum==4)
	{
		matptr = gon250mt;
		mat_xref = def_aa_xref;
	}
	else if (score_matnum==5)
	{
		matptr = gon350mt;
		mat_xref = def_aa_xref;
	}
	else
	{
		matptr = score_matrix;
		mat_xref = score_aa_xref;
	}
	maxres = get_matrix(matptr, mat_xref, matrix, FALSE, 100);
	if (maxres == 0)
	{
		error("matrix not found for aln score");
		return;
	}

	profile = (sint **) ckalloc( (data.ncols+2) * sizeof (sint *) );
	for(p=0; p<data.ncols; p++)
		profile[p] = (sint *) ckalloc( (max_aa+2) * sizeof(sint) );
	freq = (sint *) ckalloc( (max_aa+2) * sizeof (sint) );
 
	for(p=0;p<data.ncols;p++)
	{
		for(r=0;r<max_aa;r++)
			freq[r]=0;
		for(s=data.firstseq;s<data.firstseq+data.nseqs;s++)
			if(p<seqlen_array[s+1] && seq_array[s+1][p+1]>=0 && seq_array[s+1][p+1]<max_aa)
			{
				freq[seq_array[s+1][p+1]]++;
			}
		for(r=0;r<max_aa;r++)
		{
			profile[p][r]=0;
			for(r1=0;r1<max_aa;r1++)
				profile[p][r]+=freq[r1]*matrix[r1][r];
			profile[p][r]/=(float)data.nseqs;
		}
	}
/*
fprintf(fd,"Profile...\n");
for(r=0;r<max_aa;r++)
{
	for(p=0;p<data.ncols;p++)
		fprintf(fd,"%d\t",profile[p][r]);
	fprintf(fd,"\n");
}
*/
	seqvector = (sint *) ckalloc( (max_aa+2) * sizeof(sint) );
	seqdist=(float *)ckalloc((data.nseqs+1)*sizeof(float));
	sorteddist=(float *)ckalloc((data.nseqs+1)*sizeof(float));

    	for(p=0; p<data.ncols; p++)
	{
    		for(s=data.firstseq; s<data.firstseq+data.nseqs; s++)
		{
			if (p<seqlen_array[s+1])
				for (r=0;r<max_aa; r++)
					seqvector[r]=matrix[r][(int)seq_array[s+1][p+1]];
			else
				for (r=0;r<max_aa; r++)
					seqvector[r]=matrix[r][gap_pos1];
			seqdist[s-data.firstseq]=0.0;
			for(r=0;r<max_aa;r++)
			{
				diff=profile[p][r]-seqvector[r];
				diff/=1000.0;
				seqdist[s-data.firstseq]+=diff*diff;
			}
			seqdist[s-data.firstseq]=sqrt((double)seqdist[s-data.firstseq]);
		}
/*
fprintf(fd,"\n\nPosition %d:\n",p+1);
fprintf(fd,"Sequence Distances...\n");
for(s=0;s<data.nseqs;s++)
	fprintf(fd,"%.1f\t",seqdist[s]);
*/
/* calculate mean,median and rms of seq distances */
		mean=median=0.0;
		if(include_gaps)
		{
    			for(s=0; s<data.nseqs; s++)
				mean+=seqdist[s];
			mean/=data.nseqs;
			n=data.nseqs;
    			for(s=0; s<data.nseqs; s++)
					sorteddist[s]=seqdist[s];
		}
		else
		{
			n=0;
    			for(s=data.firstseq; s<data.firstseq+data.nseqs; s++)
				if(p<seqlen_array[s+1] && seq_array[s+1][p+1]>=0 && seq_array[s+1][p+1]<max_aa)
				{
					mean+=seqdist[s-data.firstseq];
					n++;
				}
			if(n>0) mean/=n;
    			for(s=data.firstseq,i=0; s<data.firstseq+data.nseqs; s++)
				if(p<seqlen_array[s+1] && seq_array[s+1][p+1]>=0 && seq_array[s+1][p+1]<max_aa)
					sorteddist[i++]=seqdist[s-data.firstseq];
		}
		sort_scores(sorteddist,0,n-1);
/*
fprintf(fd,"\nSorted:\n");
for(s=0;s<n;s++)
	fprintf(fd,"%.1f ",sorteddist[s]);
*/

		if(n == 0)
			median = 0;
		else if(n % 2 == 0)
			median=(sorteddist[n/2-1]+sorteddist[n/2])/2.0;
		else
			median=sorteddist[n/2];
		if(score_scale<=5)
			data.colscore[p]=exp((double)(-mean*(6-score_scale)/4.0))*100.0*n/data.nseqs;
		else
			data.colscore[p]=exp((double)(-mean/(4.0*(score_scale-4))))*100.0*n/data.nseqs;
/*
fprintf(fd,"\nMean %.1f Median %.1f Score %.1f\n",mean,median,data.colscore[p]);
*/
		if(n==0)
		{
			ul=0;
		}
		else
		{
			t = n/4.0 + 0.5;
			if(t - (int)t == 0.5)
			{
				q3=(sorteddist[(int)t]+sorteddist[(int)t+1])/2.0;
				q1=(sorteddist[n-(int)t]+sorteddist[n-(int)t-1])/2.0;
			}
			else if(t - (int)t > 0.5)
			{
				q3=sorteddist[(int)t+1];
				q1=sorteddist[n-(int)t-1];
			}
			else 
			{
				q3=sorteddist[(int)t];
				q1=sorteddist[n-(int)t];
			}
			if (n<4)ul=sorteddist[0];
			else ul=q3+(q3-q1)*((float)score_cutoff/2.0);
		}
/*
fprintf(fd,"\nMedian %.1f Q1 %.1f Q3 %.1f UL %.1f\n",median,q1,q3,ul);
fprintf(fd,"\nExceptions: ");
for(s=0;s<data.nseqs;s++)
	if(seqdist[s]>ul) fprintf(fd,"%d ",s+1);
*/
		for(s=data.firstseq;s<data.firstseq+data.nseqs;s++)
			if(seqdist[s-data.firstseq]>ul && p<seqlen_array[s+1] && seq_array[s+1][p+1]>=0 && seq_array[s+1][p+1]<max_aa)
				data.residue_exception[s-data.firstseq][p]=TRUE;
			else
				data.residue_exception[s-data.firstseq][p]=FALSE;
	}
/*
fclose(fd);
*/
	for(p=0;p<data.ncols;p++)
		ckfree(profile[p]);
	ckfree(profile);
	ckfree(freq);
	ckfree(seqvector);
	ckfree(seqdist);
	ckfree(sorteddist);


}

 
void sort_scores(float *scores,int f,int l)
{
	int i,last;

	if(f>=l) return;

	swap(scores,f,(f+l)/2);
	last=f;
	for(i=f+1;i<=l;i++)
	{
		if(scores[i]>scores[f])
			swap(scores,++last,i);
	}
	swap(scores,f,last);
	sort_scores(scores,f,last-1);
	sort_scores(scores,last+1,l);

}

void swap(float *scores,int s1, int s2)
{
	float temp;

	temp=scores[s1];
	scores[s1]=scores[s2];
	scores[s2]=temp;
}


void set_scorescale(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
        char str[FILENAMELEN];
	panel_data data;
 
        score_scale = newval+1;

	calc_colscores(FALSE,TRUE);

	sprintf(str,"Score Plot Scale:   %2d",score_scale);
	SetTitle(scorescaletext,str);
}

void set_scorecutoff(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
        char str[FILENAMELEN];
	int temp;
	panel_data data;
	temp=newval+1;

        score_cutoff = temp;

	calc_colscores(residue_exceptions,FALSE);
	sprintf(str,"Residue Exception Cutoff:   %2d",score_cutoff);
	SetTitle(residue_cutofftext,str);
}

 
 

void calc_segment_exceptions(IteM i)
{
	WatchCursor();
	segment_exceptions=TRUE;
	calc_seg_exceptions();
	show_segment_exceptions();
	SetValue(show_seg_toggle,1);
	SetStatus(segment_item,segment_exceptions);
	ArrowCursor();
}

void set_lengthcutoff(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
        char str[100];
	int temp;
 
	temp=newval+1;

        length_cutoff = temp;
	sprintf(str,"Minimum Length of Segments:   %2d",length_cutoff);

	if(aln_mode==MULTIPLEM)
	{
		remove_short_segments(seq_panel.seqs);
	}
	else
	{
		remove_short_segments(prf_panel[0].seqs);
		remove_short_segments(prf_panel[1].seqs);
	}
	if(segment_exceptions) show_segment_exceptions();
	SetTitle(length_cutofftext,str);

}
 
 
void set_score_user_matrix(ButtoN but)
{

	if(get_user_matrixname(score_mtrxname,score_matrix,score_aa_xref,6,&score_matnum,scoremattext))
	{
		calc_colscores(residue_exceptions,TRUE);
		SetValue(score_matrix_list,score_matnum);
	}
}
 
void set_score_matrix(GrouP g)
{
	int tmp;

        tmp = GetValue(g);
	if(tmp>0 && tmp<6)
	{
		score_matnum=tmp;
	}
	else
	{
                if (score_mtrxname[0]=='\0')
		{
			get_user_matrixname(score_mtrxname,score_matrix,score_aa_xref,6,&score_matnum,scoremattext);
		}
		else score_matnum=6;
	}
	calc_colscores(residue_exceptions,TRUE);

	SetValue(score_matrix_list,score_matnum);
}
 
void set_segment_user_matrix(ButtoN but)
{

	if(get_user_matrixname(segment_mtrxname,segment_matrix,segment_aa_xref,5,&segment_matnum,segmentmattext))
	{
		calc_seg_exceptions();
		if(segment_exceptions) show_segment_exceptions();
		SetValue(seg_matrix_list,segment_matnum);
	}
}
 
 
void set_segment_matrix(GrouP g)
{
	int tmp;

        tmp = GetValue(g);
	if(tmp>0 && tmp<5)
	{
		segment_matnum=tmp;
	}
	else
	{
		if (segment_mtrxname[0]=='\0')
		{
			get_user_matrixname(segment_mtrxname,segment_matrix,segment_aa_xref,5,&segment_matnum,segmentmattext);
		}
		else segment_matnum=5;
	}

	calc_seg_exceptions();
	if(segment_exceptions) show_segment_exceptions();

	SetValue(seg_matrix_list,segment_matnum);
}
 
 
void set_score_user_dnamatrix(ButtoN but)
{

	if(get_user_matrixname(score_dnamtrxname,score_dnamatrix,score_dna_xref,3,&score_dnamatnum,scorednamattext))
	{
		calc_colscores(residue_exceptions,TRUE);
		SetValue(score_dnamatrix_list,score_dnamatnum);
	}
}
 
void set_score_dnamatrix(GrouP g)
{
	int tmp;

        tmp = GetValue(g);
	if(tmp>0 && tmp<3)
	{
		score_dnamatnum=tmp;
	}
	else
	{
                if (score_dnamtrxname[0]=='\0')
		{
			get_user_matrixname(score_dnamtrxname,score_dnamatrix,score_dna_xref,3,&score_dnamatnum,scorednamattext);
		}
		else score_dnamatnum=3;
	}
	calc_colscores(residue_exceptions,TRUE);

	SetValue(score_dnamatrix_list,score_dnamatnum);
}
 
void set_segment_user_dnamatrix(ButtoN but)
{

	if(get_user_matrixname(segment_dnamtrxname,segment_dnamatrix,segment_dna_xref,3,&segment_dnamatnum,segmentdnamattext))
		calc_seg_exceptions();
	if(segment_exceptions) show_segment_exceptions();

	SetValue(seg_dnamatrix_list,segment_dnamatnum);
}

void set_segment_dnamatrix(GrouP g)
{
        int tmp;
 
        tmp = GetValue(g);
        if(tmp>0 && tmp<3)
        {
                segment_dnamatnum=tmp;
        }
        else
        {
                if (segment_dnamtrxname[0]=='\0')
                {
			get_user_matrixname(segment_dnamtrxname,segment_dnamatrix,segment_dna_xref,3,&segment_dnamatnum,segmentdnamattext);
                }
                else segment_dnamatnum=3;
        }
 
	calc_seg_exceptions();
        if(segment_exceptions) show_segment_exceptions();
 
        SetValue(seg_dnamatrix_list,segment_dnamatnum);
}

static void calc_colscores(Boolean update_seqs,Boolean update_scores)
{
	panel_data data;

	if(aln_mode==MULTIPLEM)
	{
		GetPanelExtra(seq_panel.seqs,&data);
		make_colscores(data);
		SetPanelExtra(seq_panel.seqs,&data);
		if(update_seqs) draw_seqs(seq_panel.seqs);
		if(update_scores) draw_colscores(seq_panel.seqs);
	}
	else
	{
		GetPanelExtra(prf_panel[0].seqs,&data);
		make_colscores(data);
		SetPanelExtra(prf_panel[0].seqs,&data);
		if(update_seqs) draw_seqs(prf_panel[0].seqs);
		if(update_scores) draw_colscores(prf_panel[0].seqs);
		GetPanelExtra(prf_panel[1].seqs,&data);
		make_colscores(data);
		SetPanelExtra(prf_panel[1].seqs,&data);
		if(update_seqs) draw_seqs(prf_panel[1].seqs);
		if(update_scores) draw_colscores(prf_panel[1].seqs);
	}
}
 
void calc_seg_exceptions(void)
{
	if(aln_mode==MULTIPLEM)
	{
		calc_panel_segment_exceptions(seq_panel.seqs);
	}
	else
	{
		calc_panel_segment_exceptions(prf_panel[0].seqs);
		calc_panel_segment_exceptions(prf_panel[1].seqs);
	}
}

void show_segment_exceptions(void)
{
	if(aln_mode==MULTIPLEM)
	{
		draw_seqs(seq_panel.seqs);
	}
	else
	{
		draw_seqs(prf_panel[0].seqs);
		draw_seqs(prf_panel[1].seqs);
	}
}

static void remove_short_segments(PaneL p)
{
	int i,j,k,start;
	panel_data data;

	GetPanelExtra(p,&data);
	if(data.nseqs<=0) return;

/* Reset all the exceptions - a value of 1 indicates an exception that
will be displayed. A value of -1 is used to remember exceptions that
are temporarily hidden in the display */
        for(i=0;i<data.nseqs;i++)
                for(j=0;j<data.ncols;j++)
			if(data.segment_exception[i][j] == -1)
				data.segment_exception[i][j] = 1;

        for(i=0;i<data.nseqs;i++)
        {
		start = -1;
                for(j=0;j<=data.ncols;j++)
                {
			if(start == -1)
			{
				if(data.segment_exception[i][j]==1)
					start=j;
			}
			else
			{
				if(j==data.ncols || data.segment_exception[i][j]==0)
				{
					if(j-start<length_cutoff)
						for(k=start;k<j;k++)
							data.segment_exception[i][k] = -1;
					start = -1;
				}
			}

		}
	}

	SetPanelExtra(p,&data);
}

static void calc_weights(PaneL p)
{
	int i,j;
	int status;
	sint *weight;
	float dscore;
	FILE *tree;
	panel_data data;

#ifdef UNIX
	char tree_name[FILENAMELEN]=".score.ph";
#else
	char tree_name[FILENAMELEN]="tmp.ph";
#endif

	GetPanelExtra(p,&data);
	if(data.nseqs<=0) return;

/* if sequence weights have been calculated before - don't bother
doing it again (it takes too long). data.seqweight is set to NULL when
 new sequences are loaded. */
	if(data.seqweight!=NULL) return;

	WatchCursor();
	info("Calculating sequence weights...");
/* count pairwise percent identities to make a phylogenetic tree */
	if(data.nseqs>=2)
	{
        	for (i=1;i<=data.nseqs;i++) {
                	for (j=i+1;j<=data.nseqs;j++) {
                        	dscore = countid(i+data.firstseq,j+data.firstseq);
                        	tmat[i][j] = (100.0 - dscore)/100.0;
                        	tmat[j][i] = tmat[i][j];
                	}
        	}

                if((tree = open_explicit_file(tree_name))==NULL) return;

		guide_tree(tree,data.firstseq+1,data.nseqs);

       		status = read_tree(tree_name, data.firstseq, data.firstseq+data.nseqs);
        	if (status == 0) return;
 
	}
 
	weight = (sint *) ckalloc( (data.firstseq+data.nseqs+1) * sizeof(sint) );
/* get the sequence weights */
 	calc_seq_weights(data.firstseq, data.firstseq+data.nseqs,weight);
	if(data.seqweight==NULL) data.seqweight=(sint *)ckalloc((data.nseqs+1) * sizeof(sint));
	for(i=data.firstseq;i<data.firstseq+data.nseqs;i++)
		data.seqweight[i-data.firstseq]=weight[i];

/* clear the memory for the phylogenetic tree */
   	if (data.nseqs >= 2)
	{
        	clear_tree(NULL);
		remove(tree_name);
	}
	ckfree(weight);
	SetPanelExtra(p,&data);
	info("Done.");
	ArrowCursor();
}

static void calc_panel_segment_exceptions(PaneL p)
{
        int i,j;
        float sum,prev_sum;
	float gscale;
	sint **profile;
	sint *weight,sweight;
	sint *gaps;
	sint maxres;
	int max=0,offset;
	short  *mat_xref, *matptr;
	sint matrix[NUMRES][NUMRES];
	float *fsum;
	float *bsum;
	float *pscore;
	panel_data data;
 
/* First, calculate sequence weights which will be used to build the
profile */
	calc_weights(p);

	GetPanelExtra(p,&data);
	if(data.nseqs<=0) return;
 
	WatchCursor();
	info("Calculating profile scores...");

        for(i=0;i<data.nseqs;i++)
                for(j=0;j<data.ncols;j++)
                     	data.segment_exception[i][j]=0;

/* get the comparison matrix for building the profile */
	if(dnaflag)
	{
		if (segment_dnamatnum==1)
		{
			matptr = swgapdnamt;
			mat_xref = def_dna_xref;
		}
		else if (segment_dnamatnum==2)
		{
			matptr = clustalvdnamt;
			mat_xref = def_dna_xref;
		}
		else
		{
			matptr = segment_dnamatrix;
			mat_xref = segment_dna_xref;
		}
/* get a positive matrix - then adjust it according to scale */
		maxres = get_matrix(matptr, mat_xref, matrix, FALSE, 100);
/* find the maximum value */
	for(i=0;i<=max_aa;i++)
		for(j=0;j<=max_aa;j++)
			if(matrix[i][j]>max) max=matrix[i][j];
/* subtract max*scale/2 from each matrix value */
	offset=(float)(max*segment_dnascale)/20.0;

	for(i=0;i<=max_aa;i++)
		for(j=0;j<=max_aa;j++)
			matrix[i][j]-=offset;
	}
	else
	{
		if (segment_matnum==1)
		{
			matptr = gon80mt;
			mat_xref = def_aa_xref;
		}
		else if (segment_matnum==2)
		{
			matptr = gon120mt;
			mat_xref = def_aa_xref;
		}
		else if (segment_matnum==3)
		{
			matptr = gon250mt;
			mat_xref = def_aa_xref;
		}
		else if (segment_matnum==4)
		{
			matptr = gon350mt;
			mat_xref = def_aa_xref;
		}
		else
		{
			matptr = segment_matrix;
			mat_xref = segment_aa_xref;
		}
/* get a negative matrix */
		maxres = get_matrix(matptr, mat_xref, matrix, TRUE, 100);
	}

	profile = (sint **) ckalloc( (data.ncols+2) * sizeof (sint *) );
	for(i=0; i<data.ncols+1; i++)
		profile[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

/* calculate the profile */
	gaps = (sint *) ckalloc( (data.ncols+1) * sizeof (sint) );
	for (j=1; j<=data.ncols; j++)
	{
		gaps[j-1] = 0;
		for(i=data.firstseq+1;i<data.firstseq+data.nseqs;i++)
			if (j<seqlen_array[i])
				if ((seq_array[i][j] < 0) || (seq_array[i][j] > max_aa))
					gaps[j-1]++;
	}
	weight = (sint *) ckalloc( (data.firstseq+data.nseqs+1) * sizeof(sint) );
	for(i=data.firstseq;i<data.firstseq+data.nseqs;i++)
		weight[i]=data.seqweight[i-data.firstseq];

	build_profile(data.ncols,data.firstseq,data.firstseq+data.nseqs,matrix,weight,profile);

	sweight=0;
        for(i=data.firstseq;i<data.firstseq+data.nseqs;i++)
		sweight+=weight[i];

/*Now, use the profile scores to mark segments of each sequence which score
badly. */

	fsum = (float *) ckalloc( (data.ncols+2) * sizeof (float) );
	bsum = (float *) ckalloc( (data.ncols+2) * sizeof (float) );
	pscore = (float *) ckalloc( (data.ncols+2) * sizeof (float) );
        for(i=data.firstseq+1;i<data.firstseq+data.nseqs+1;i++)
        {
/* In a forward phase, sum the profile scores. Mark negative sums as exceptions.
If the sum is positive, then it gets reset to 0. */
		sum=0.0;
                for(j=1;j<=seqlen_array[i];j++)
                {
           		gscale = (float)(data.nseqs-gaps[j-1]) / (float)data.nseqs;
			if(seq_array[i][j]<0 || seq_array[i][j]>=max_aa)
			{
				pscore[j-1]=0.0;
				sum=0.0;
			}
			else
                        	pscore[j-1]=(profile[j][seq_array[i][j]]-
			weight[i-1]*matrix[seq_array[i][j]][seq_array[i][j]])*gscale/sweight;
                        sum+=pscore[j-1];
			if(sum>0.0) sum=0.0;
			fsum[j-1]=sum;
                }
/* trim off any positive scoring residues from the end of the segments */
		prev_sum=0;
                for(j=seqlen_array[i]-1;j>=0;j--)
                {
			if(prev_sum>=0.0 && fsum[j]<0.0 && pscore[j]>=0.0)
				fsum[j]=0.0;
			prev_sum=fsum[j];
                }

/* Now, in a backward phase, do the same summing process. */
		sum=0.0;
                for(j=seqlen_array[i];j>=1;j--)
                {
			if(seq_array[i][j]<0 || seq_array[i][j]>=max_aa)
				sum=0;
			else
                        	sum+=pscore[j-1];
			if(sum>0.0) sum=0.0;
			bsum[j-1]=sum;
                }
/* trim off any positive scoring residues from the start of the segments */
		prev_sum=0;
                for(j=0;j<seqlen_array[i];j++)
                {
			if(prev_sum>=0.0 && bsum[j]<0.0 && pscore[j]>=0.0)
				bsum[j]=0.0;
			prev_sum=bsum[j];
                }
/*Mark residues as exceptions if they score negative in the forward AND backward directions. */
                for(j=1;j<=seqlen_array[i];j++)
			if(fsum[j-1]<0.0 && bsum[j-1]<0.0)
				if(seq_array[i][j]>=0 && seq_array[i][j]<max_aa)
				data.segment_exception[i-data.firstseq-1][j-1]=-1;
/*
if(i==5) {
fprintf(stderr,"%4d ",j);
fprintf(stderr,"\n");
for(j=0;j<seqlen_array[i];j++)
fprintf(stderr,"%4d ",(int)fsum[j]);
fprintf(stderr,"\n");
}
*/
        }
	for(i=0; i<data.ncols+1; i++)
		ckfree(profile[i]);
	ckfree(profile);
	ckfree(weight);
	ckfree(gaps);
	ckfree(pscore);
	ckfree(fsum);
	ckfree(bsum);

	SetPanelExtra(p,&data);

/* Finally, apply the length cutoff to the segments - removing segments shorter
than the cutoff */
	remove_short_segments(p);

	info("Done.");
	ArrowCursor();
}


void set_segment_dnascale(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
        char str[FILENAMELEN];
	panel_data data;
 
        segment_dnascale = newval+1;
	calc_seg_exceptions();
	if(segment_exceptions) show_segment_exceptions();
	sprintf(str,"DNA Marking Scale:   %2d",segment_dnascale);
	SetTitle(segmentdnascaletext,str);
}

 
static void build_profile(int prf_length,int first_seq,int last_seq,sint matrix[NUMRES][NUMRES],sint *weight,sint **profile)
{
  sint **weighting, d, i, res; 
  sint r, pos;
  int f;

  weighting = (sint **) ckalloc( (NUMRES+2) * sizeof (sint *) );
  for (i=0;i<NUMRES+2;i++)
    weighting[i] = (sint *) ckalloc( (prf_length+2) * sizeof (sint) );

  for (r=0; r<prf_length; r++)
   {
      for (d=0; d<=max_aa; d++)
        {
            weighting[d][r] = 0;
            for (i=first_seq; i<last_seq; i++)
		if (r+1<seqlen_array[i+1])
               		if (d == seq_array[i+1][r+1]) weighting[d][r] += weight[i];
        }
      weighting[gap_pos1][r] = 0;
      for (i=first_seq; i<last_seq; i++)
	if (r+1<seqlen_array[i+1])
         if (gap_pos1 == seq_array[i+1][r+1]) weighting[gap_pos1][r] += weight[i];
      weighting[gap_pos2][r] = 0;
      for (i=first_seq; i<last_seq; i++)
	if (r+1<seqlen_array[i+1])
         if (gap_pos2 == seq_array[i+1][r+1]) weighting[gap_pos2][r] += weight[i];
   }

  for (pos=0; pos< prf_length; pos++)
    {
           for (res=0; res<=max_aa; res++)
             {
                f = 0;
                for (d=0; d<=max_aa; d++)
                     f += (weighting[d][pos] * matrix[d][res]);
                f += (weighting[gap_pos1][pos] * matrix[gap_pos1][res]);
                f += (weighting[gap_pos2][pos] * matrix[gap_pos2][res]);
                profile[pos+1][res] = f;
             }
           f = 0;
           for (d=0; d<=max_aa; d++)
                f += (weighting[d][pos] * matrix[d][gap_pos1]);
           f += (weighting[gap_pos1][pos] * matrix[gap_pos1][gap_pos1]);
           f += (weighting[gap_pos2][pos] * matrix[gap_pos2][gap_pos1]);
           profile[pos+1][gap_pos1] = f;
           f = 0;
           for (d=0; d<=max_aa; d++)
                f += (weighting[d][pos] * matrix[d][gap_pos2]);
           f += (weighting[gap_pos1][pos] * matrix[gap_pos1][gap_pos2]);
           f += (weighting[gap_pos2][pos] * matrix[gap_pos2][gap_pos2]);
           profile[pos+1][gap_pos2] = f;
    }

  for (i=0;i<=max_aa;i++)
    weighting[i]=ckfree((void *)weighting[i]);
  weighting=ckfree((void *)weighting);

}


