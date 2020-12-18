#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <vibrant.h>
#include <document.h>

#include "clustalw.h"
#include "xmenu.h"

static void VscrollMulti(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollMultiN(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollMultiS(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void VscrollPrf1(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollPrf1N(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollPrf1S(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void VscrollPrf2(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollPrf2N(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);
static void HscrollPrf2S(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval);

static void NameClick(PaneL panel, PoinT pt);
static void NameDrag(PaneL panel, PoinT pt);
static void NameRelease(PaneL panel, PoinT pt);
static void SeqClick(PaneL panel, PoinT pt);
static void SeqDrag(PaneL panel, PoinT pt);
static void SeqRelease(PaneL panel, PoinT pt);

static void fit_seq_display(RecT wr,Boolean mv_message);
static void fit_prf_displays(RecT wr,int numseqs1,int numseqs2,Boolean mv_message);

static void vscrollnames(BaR bar, int newval, int oldval);
static void hscrollnames(BaR bar, int newval, int oldval);
static void vscrollseqs(BaR bar, int newval, int oldval);
static void hscrollseqs(BaR bar, int newval, int oldval);

static void correct_scrollbar(BaR b,int visible,int total,int value,Boolean reset);

static PaneL make_panel(int type,GrouP g,int width,int height,int firstline,int tlines);
static panel_data free_panel_data(panel_data data);
static panel_data alloc_name_data(panel_data data);
static panel_data alloc_seq_data(panel_data data);

extern sint max_names;

extern int    mheader; /* maximum header lines */
extern int    mfooter; /* maximum footer lines */
extern int max_plines;     /*   profile align display length */
extern int min_plines1;     /*   profile align display length */
extern int min_plines2;     /*   profile align display length */
extern int loffset,boffset,toffset;
extern int roffset;
extern int poffset;

extern Boolean aln_mode;
extern Boolean fixed_prf_scroll;
extern Boolean window_displayed;

extern PrompT   message;           /* used in temporary message window */

extern spanel  seq_panel;        /* data for multiple alignment area */
extern spanel  prf_panel[];       /* data for profile alignment areas */
extern spanel  active_panel;       /* 'in-use' panel -scrolling,clicking etc. */
extern FonT datafont;
extern WindoW mainw;
extern GrouP  seq_display,prf1_display,prf2_display;

extern int ncolors;
extern int    inverted;

extern Boolean  dnaflag;
extern sint     nseqs;
extern sint    profile1_nseqs;
extern sint     output_order;
extern sint     *output_index;
extern sint     *seqlen_array;
extern char     **seq_array;
extern char     **names, **titles;
extern char     *amino_acid_codes;
extern sint     gap_pos1, gap_pos2;
extern char *gap_penalty_mask1,*gap_penalty_mask2;
extern char *sec_struct_mask1,*sec_struct_mask2;
extern sint struct_penalties1,struct_penalties2;
extern sint output_struct_penalties;
extern Boolean use_ss1, use_ss2;

extern char *explicit_par_file;
extern char *par_file;
extern char def_protpar_file[];
extern char def_dnapar_file[];
extern sint     ncutseqs;
extern Boolean residue_exceptions;
extern Boolean segment_exceptions;
extern color color_lut[];
extern char *res_cat1[];
extern char *res_cat2[];

static range selected_seqs;           /* sequences selected by clicking on names */
static range selected_res;           /* residues selected by clicking on seqs */

 
static int fromvscroll,fromhscroll; /* set by scrolling functions,
                            used by DrawPanel, draw_names, draw_seqs */


void resize_multi_window(void)
{
	RecT r;

	ObjectRect(mainw,&r);
	fit_seq_display(r,FALSE);
}

void resize_prf_window(int numseqs1,int numseqs2)
{
	RecT r;

	SelectFont(datafont);
	stdCharWidth=CharWidth('A');
        stdLineHeight=LineHeight();

	if(numseqs1>max_plines)
		numseqs1=max_plines;
	else if(numseqs1<min_plines1)
		numseqs1=min_plines1;
	if(numseqs2>max_plines)
		numseqs2=max_plines;
	else if(numseqs2<min_plines2)
		numseqs2=min_plines2;
	ObjectRect(mainw,&r);

	fit_prf_displays(r,numseqs1,numseqs2,FALSE);
}

static void fit_seq_display(RecT wr,Boolean mv_message)
{
	int width,height,moffset;
	RecT nr,sr,mr;
	panel_data data;

	ObjectRect(seq_panel.names,&nr);
	ObjectRect(message,&mr);
	moffset=mr.top-nr.bottom;
	width=nr.right-nr.left;
	height=wr.bottom-wr.top-boffset-toffset;
	nr.top=toffset;
	nr.left=loffset;
	nr.bottom=nr.top+height;
	nr.right=nr.left+width;
	SetPosition(seq_panel.names,&nr);

	GetPanelExtra(seq_panel.names,&data);
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(seq_panel.names,&data);


	sr.top=nr.top;
	sr.left=nr.right;
	sr.bottom=sr.top+height;
	sr.right=wr.right-wr.left-roffset;
	width=sr.right-sr.left;
	SetPosition(seq_panel.seqs,&sr);

	GetPanelExtra(seq_panel.seqs,&data);
	data.vcols=width/data.charwidth - MARGIN*2;
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(seq_panel.seqs,&data);

	if(mv_message) {
		height=mr.bottom-mr.top;
		mr.top=nr.bottom+moffset;
		mr.bottom=mr.top+height;
		SetPosition(message,&mr);
	}


	position_scrollbars(seq_panel);
	correct_name_bars(TRUE);
	correct_seq_bars(TRUE);

}

static void fit_prf_displays(RecT wr,int numseqs1,int numseqs2,Boolean mv_message)
{
	int width,height,moffset;
	RecT mr,nr,sr;
	panel_data data;

	ObjectRect(prf_panel[1].names,&nr);
	ObjectRect(message,&mr);
	moffset=mr.top-nr.bottom;

	ObjectRect(prf_panel[0].names,&nr);
	width=nr.right-nr.left;

	nr.top=toffset;
	nr.left=loffset;
	height=(wr.bottom-wr.top-boffset-toffset-poffset)*numseqs1/(numseqs1+numseqs2);
	nr.bottom=nr.top+height;
	nr.right=nr.left+width;
	SetPosition(prf_panel[0].names,&nr);
	GetPanelExtra(prf_panel[0].names,&data);
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(prf_panel[0].names,&data);
	sr.top=nr.top;
	sr.left=nr.right;
	sr.bottom=sr.top+height;
	sr.right=wr.right-wr.left-roffset;
	width=sr.right-sr.left;
	SetPosition(prf_panel[0].seqs,&sr);
	GetPanelExtra(prf_panel[0].seqs,&data);
	data.vcols=width/data.charwidth - MARGIN*2;
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(prf_panel[0].seqs,&data);
	position_scrollbars(prf_panel[0]);



	nr.top=nr.bottom+poffset;
	height=(wr.bottom-wr.top-boffset-toffset-poffset)*numseqs2/(numseqs1+numseqs2);
	nr.bottom=nr.top+height;
	SetPosition(prf_panel[1].names,&nr);
	GetPanelExtra(prf_panel[1].names,&data);
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(prf_panel[1].names,&data);
	sr.top=nr.top;
	sr.bottom=sr.top+height;
	SetPosition(prf_panel[1].seqs,&sr);
	GetPanelExtra(prf_panel[1].seqs,&data);
	data.vcols=width/data.charwidth - MARGIN*2;
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	SetPanelExtra(prf_panel[1].seqs,&data);
	position_scrollbars(prf_panel[1]);

	if(mv_message) {
		height=mr.bottom-mr.top;
		mr.top=nr.bottom+moffset;
		mr.bottom=mr.top+height;
		SetPosition(message,&mr);
	}

	correct_name_bars(TRUE);
	correct_seq_bars(TRUE);
}

void ResizeWindowProc(WindoW w)
{
	int numseqs1,numseqs2;
	RecT wr;
	panel_data data;

	if(window_displayed==FALSE) return;

	ObjectRect(w,&wr);
	if (aln_mode==MULTIPLEM)
	{
/* if the window is too small, hide everything */
		if(wr.bottom-wr.top < toffset+boffset)
		{
			Hide(seq_display);
			Hide(message);
			return;
		}
		fit_seq_display(wr,TRUE);
		Show(seq_display);
		Show(message);
	}
	else
	{
/* if the window is too small, hide everything */
		if(wr.bottom-wr.top < toffset+boffset+2*poffset)
		{
			Hide(prf1_display);
			Hide(prf2_display);
			Hide(message);
			return;
		}
		GetPanelExtra(prf_panel[0].names,&data);
		numseqs1=data.nseqs;
		if(numseqs1<min_plines1)numseqs1=min_plines1;
		else if(numseqs1>max_plines)numseqs1=max_plines;
		GetPanelExtra(prf_panel[1].names,&data);
		numseqs2=data.nseqs;
		if(numseqs2<min_plines2)numseqs2=min_plines2;
		else if(numseqs2>max_plines)numseqs2=max_plines;

		fit_prf_displays(wr,numseqs1,numseqs2,TRUE);
		Show(prf1_display);
		Show(prf2_display);
		Show(message);
	}

}

void position_scrollbars(spanel p)
{
	int height;
	RecT hr,vr,nr,sr;
	panel_data data;

	ObjectRect(p.names,&nr);
	GetPanelExtra(p.names,&data);
	ObjectRect(data.hscrollbar,&hr);
	height=hr.bottom-hr.top;
	LoadRect(&hr,nr.left,nr.bottom,nr.right,nr.bottom+height);
	SetPosition(data.hscrollbar,&hr);
#ifdef WIN_MAC
	AdjustPrnt(data.hscrollbar,&hr,FALSE);
#endif
	ObjectRect(p.seqs,&sr);
	GetPanelExtra(p.seqs,&data);
	ObjectRect(data.hscrollbar,&hr);
	height=hr.bottom-hr.top;
	LoadRect(&hr,sr.left,sr.bottom,sr.right,sr.bottom+height);
	SetPosition(data.hscrollbar,&hr);
#ifdef WIN_MAC
	AdjustPrnt(data.hscrollbar,&hr,FALSE);
#endif
	ObjectRect(data.vscrollbar,&vr);
	LoadRect(&vr,vr.left,sr.top,vr.right,sr.bottom);
	SetPosition(data.vscrollbar,&vr);
#ifdef WIN_MAC
	AdjustPrnt(data.vscrollbar,&vr,FALSE);
#endif
}





void load_aln_data(spanel p,int fs,int ls,Boolean reset)
{
	int i,j,slength=0;
	int nhead;
	sint val;
	panel_data name_data,seq_data;

	WatchCursor();

	GetPanelExtra(p.names,&name_data);
	GetPanelExtra(p.seqs,&seq_data);
	name_data=free_panel_data(name_data);
	seq_data=free_panel_data(seq_data);
	SetPanelExtra(p.names,&name_data);
	SetPanelExtra(p.seqs,&seq_data);

	name_data.nseqs=ls-fs+1;
	seq_data.nseqs=name_data.nseqs;
	name_data.firstseq=fs;
	seq_data.firstseq=fs;

/* find the maximum length of sequence */
	for(i=fs;i<=ls;i++)
           	if (slength < seqlen_array[i+1]) slength = seqlen_array[i+1];
	name_data.ncols=max_names;
	seq_data.ncols=slength;

	if (name_data.nseqs>0)
	{
		name_data=alloc_name_data(name_data);
		seq_data=alloc_seq_data(seq_data);

	
		for(i=fs;i<=ls;i++)
        	{
                	strncpy(name_data.lines[i-fs],names[i+1],MAXNAMES);
                	name_data.lines[i-fs][MAXNAMES]='\0';
                	for(j=0;j<seqlen_array[i+1];j++)
                	{
                        	val = seq_array[i+1][j+1];
                        	if((val == -3) || (val == 253))
                                	break;
                        	else if((val == gap_pos1) || (val == gap_pos2))
                                	seq_data.lines[i-fs][j] = '-';
                        	else {
                                	seq_data.lines[i-fs][j] = amino_acid_codes[val];
                        	}
                	}
                	for(j=seqlen_array[i+1];j<slength;j++)
                        	seq_data.lines[i-fs][j] = ' ';
			seq_data.lines[i-fs][j]='\0';

			name_data.selected[i-fs]=FALSE;

        	}


		make_consensus(seq_data,name_data.header[0],seq_data.header[0]);
		nhead=make_struct_data(seq_data.prf_no,slength,name_data.header[1],seq_data.header[1]);
		if (nhead==0)
			nhead=make_gp_data(seq_data.prf_no,slength,name_data.header[1],seq_data.header[1]);
		seq_data.nhead=name_data.nhead=nhead+1;

		seq_data.nfoot=name_data.nfoot=1;
		seq_data.consensus=NULL;
		make_ruler(slength,name_data.footer[0],seq_data.footer[0]);
		make_colscores(seq_data);
	}
	else
	{
		seq_data.ncols=name_data.ncols=0;
	}

	if(reset==TRUE)
	{
		name_data.firstvline=0;
		name_data.firstvcol=0;
		seq_data.firstvline=0;
		seq_data.firstvcol=0;
	}
        name_data.vseqs=name_data.vlines-name_data.nhead-name_data.nfoot;
        seq_data.vseqs=seq_data.vlines-seq_data.nhead-seq_data.nfoot;

	if(seq_data.nseqs>0)
	{
/* try to find the user's color parameter file */
		if (explicit_par_file == NULL)
		{
			if (par_file != NULL)
				ckfree(par_file);
			if(dnaflag)
				par_file=find_file(def_dnapar_file);
			else
				par_file=find_file(def_protpar_file);
		}
        	init_color_parameters(par_file);
        	make_colormask(seq_data);
	}

	SetPanelExtra(p.names,&name_data);
	SetPanelExtra(p.seqs,&seq_data);

	ArrowCursor();
}

void load_aln(spanel p,int fs,int ls,Boolean reset)
{

	load_aln_data(p,fs,ls,reset);

	DrawPanel(p.names);
	DrawPanel(p.seqs);
	correct_name_bars(reset);
	correct_seq_bars(reset);

}

static panel_data alloc_name_data(panel_data data)
{
	int i;

	data.lines=(char **)ckalloc((data.nseqs+1)*sizeof(char *));
	data.colormask=NULL;
	data.selected=(int *)ckalloc((data.nseqs+1)*sizeof(int));

	for(i=0;i<data.nseqs;i++)
        {
		data.lines[i]=(char *)ckalloc((MAXNAMES+1)*sizeof(char));
		strncpy(data.lines[i],names[i+1],MAXNAMES);
		data.lines[i][MAXNAMES]='\0';
	}

	data.header=(char **)ckalloc((mheader+1)*sizeof(char *));
	for(i=0;i<mheader;i++)
		data.header[i]=(char *)ckalloc((MAXNAMES+1)*sizeof(char));
	data.footer=(char **)ckalloc((mfooter+1)*sizeof(char *));
	for(i=0;i<mfooter;i++)
		data.footer[i]=(char *)ckalloc((MAXNAMES+1)*sizeof(char));
	return(data);
}

static panel_data alloc_seq_data(panel_data data)
{
	int i;

	data.lines=(char **)ckalloc((data.nseqs+1)*sizeof(char *));
	data.colormask=(char **)ckalloc((data.nseqs+1)*sizeof(char *));
	data.firstsel=data.lastsel=-1;

	for(i=0;i<data.nseqs;i++)
        {
		data.lines[i]=(char *)ckalloc((data.ncols+1)*sizeof(char));
		data.colormask[i]=(char *)ckalloc((data.ncols+1)*sizeof(char));
	}

	data.selected=(int *)ckalloc((data.ncols+1)*sizeof(int));
	for(i=0;i<data.ncols;i++)
		data.selected[i]=FALSE;

	data.header=(char **)ckalloc((mheader+1)*sizeof(char *));
	for(i=0;i<mheader;i++)
		data.header[i]=(char *)ckalloc((data.ncols+1)*sizeof(char));

	data.colscore=(sint *)ckalloc((data.ncols+1)*sizeof(sint));
	data.residue_exception=(Boolean **)ckalloc((data.nseqs+1)*sizeof(Boolean *));
	for(i=0;i<data.nseqs;i++)
		data.residue_exception[i]=(Boolean *)ckalloc((data.ncols+1)*sizeof(Boolean));
	data.segment_exception=(short **)ckalloc((data.nseqs+1)*sizeof(short *));
	for(i=0;i<data.nseqs;i++)
		data.segment_exception[i]=(short *)ckalloc((data.ncols+1)*sizeof(short));

	data.footer=(char **)ckalloc((mfooter+1)*sizeof(char *));
	for(i=0;i<mfooter;i++)
		data.footer[i]=(char *)ckalloc((data.ncols+1)*sizeof(char));
	return(data);
}

void correct_name_bars(Boolean reset)
{
	panel_data data,data1;

	if(aln_mode==PROFILEM)
	{
		GetPanelExtra(prf_panel[0].names,&data);
		GetPanelExtra(prf_panel[1].names,&data1);
		if(reset==TRUE)
		{
			data.firstvcol=0;
			data1.firstvcol=0;
		}
		correct_scrollbar(data.hscrollbar,data.vcols,data.ncols,data.firstvcol,reset);
		correct_scrollbar(data1.hscrollbar,data1.vcols,data1.ncols,data.firstvcol,reset);
		if(reset==TRUE)
		{
			data.firstvline=0;
			data1.firstvline=0;
		}
		correct_scrollbar(data.vscrollbar,data.vseqs,data.nseqs,data.firstvline,reset);
		correct_scrollbar(data1.vscrollbar,data1.vseqs,data1.nseqs,data1.firstvline,reset);
		SetPanelExtra(prf_panel[0].names,&data);
		SetPanelExtra(prf_panel[1].names,&data1);
	}
	else
	{
		GetPanelExtra(seq_panel.names,&data);
		if(reset==TRUE)
		{
			data.firstvcol=0;
			data.firstvline=0;
		}
		correct_scrollbar(data.vscrollbar,data.vseqs,data.nseqs,data.firstvline,reset);
		correct_scrollbar(data.hscrollbar,data.vcols,data.ncols,data.firstvcol,reset);

		SetPanelExtra(seq_panel.names,&data);
	}
	
}

void correct_seq_bars(Boolean reset)
{
	int maxcols,m1,m2;
	panel_data data,data1;

	if(aln_mode==PROFILEM)
	{
		GetPanelExtra(prf_panel[0].seqs,&data);
		GetPanelExtra(prf_panel[1].seqs,&data1);
		if(fixed_prf_scroll==TRUE)
		{
			Hide(data.hscrollbar);
			m1=MAX(data.firstvcol,data1.firstvcol);
			m2=MAX(data.ncols-data.firstvcol,data1.ncols-data1.firstvcol);
			maxcols=m1+m2;
			if(reset==TRUE)
			{
				data.firstvcol=0;
				data1.firstvcol=0;
			}
			data.lockoffset= -MAX(data1.firstvcol-data.firstvcol,0);
			data1.lockoffset= -MAX(data.firstvcol-data1.firstvcol,0);
			correct_scrollbar(data1.hscrollbar,data1.vcols,maxcols,m1,TRUE);
		}
		else
		{
			Show(data.hscrollbar);
			if(reset==TRUE)
			{
				data.firstvcol=0;
				data1.firstvcol=0;
			}
			data.lockoffset=0;
			data1.lockoffset=0;
			correct_scrollbar(data.hscrollbar,data.vcols,data.ncols,data.firstvcol,reset);
			correct_scrollbar(data1.hscrollbar,data1.vcols,data1.ncols,data.firstvcol,reset);
		}
		if(reset==TRUE)
		{
			data.firstvline=0;
			data1.firstvline=0;
		}
		correct_scrollbar(data.vscrollbar,data.vseqs,data.nseqs,data.firstvline,reset);
		correct_scrollbar(data1.vscrollbar,data1.vseqs,data1.nseqs,data.firstvline,reset);
		SetPanelExtra(prf_panel[0].seqs,&data);
		SetPanelExtra(prf_panel[1].seqs,&data1);
	}
	else
	{
		GetPanelExtra(seq_panel.seqs,&data);
		if(reset==TRUE)
		{
			data.firstvcol=0;
			data.firstvline=0;
		}
		correct_scrollbar(data.vscrollbar,data.vseqs,data.nseqs,data.firstvline,reset);
		correct_scrollbar(data.hscrollbar,data.vcols,data.ncols,data.firstvcol,reset);

		SetPanelExtra(seq_panel.seqs,&data);
	}
	
}

static void correct_scrollbar(BaR b,int visible,int total,int value,Boolean reset)
{
	int max;

	if (b!=NULL)
	{
		if (visible > 0 && total > visible)
			max=total-visible;
		else
			max=0;
       		if(reset==TRUE) CorrectBarValue(b,0);
       		CorrectBarPage(b,visible,visible);
       		CorrectBarValue(b,value);
       		CorrectBarMax(b,max);
	}
}


void color_seqs(void)
{
	panel_data data;

	GetPanelExtra(seq_panel.seqs,&data);
	if (data.nseqs == 0) return;

	info("Coloring sequences...");
	make_colormask(data);
	DrawPanel(seq_panel.seqs);
	info("Done.");
}

void color_prf1(void)
{
	panel_data data;

	GetPanelExtra(prf_panel[0].seqs,&data);
	if (data.nseqs == 0) return;

	make_colormask(data);
	info("Coloring profile 1...");
	DrawPanel(prf_panel[0].seqs);
	info("Done.");
}

void color_prf2(void)
{
	panel_data data;

	GetPanelExtra(prf_panel[1].seqs,&data);
	if (data.nseqs == 0) return;

	make_colormask(data);
	info("Coloring profile 2...");
	DrawPanel(prf_panel[1].seqs);
	info("Done.");
}

void remove_gap_pos(int fseq, int lseq,int prf_no)
{
	int i,j,k,ngaps;


	if (fseq>=lseq) return;

	for (i=1;i<=seqlen_array[fseq];)
	{
		ngaps=0;
		for (j=fseq;j<=lseq;j++)
			if(seq_array[j][i]==gap_pos1 || seq_array[j][i]==gap_pos2) ngaps++;
		if (ngaps==lseq-fseq+1)
		{
			for (j=fseq;j<=lseq;j++)
			{
				for(k=i+1;k<=seqlen_array[j]+1;k++)
					seq_array[j][k-1]=seq_array[j][k];
				seqlen_array[j]--;
			}
			if(prf_no==1 && sec_struct_mask1 != NULL)
				for(k=i;k<=seqlen_array[fseq];k++)
					sec_struct_mask1[k-1]=sec_struct_mask1[k];
			if(prf_no==1 && gap_penalty_mask1 != NULL)
				for(k=i;k<=seqlen_array[fseq];k++)
					gap_penalty_mask1[k-1]=gap_penalty_mask1[k];
			if(prf_no==2 && sec_struct_mask2 != NULL)
				for(k=i;k<=seqlen_array[fseq];k++)
					sec_struct_mask2[k-1]=sec_struct_mask2[k];
			if(prf_no==2 && gap_penalty_mask2 != NULL)
				for(k=i;k<=seqlen_array[fseq];k++)
					gap_penalty_mask2[k-1]=gap_penalty_mask2[k];
			if(seqlen_array[fseq]<=0) break;
		}
		else i++;
	}
}

/* width and height passed here are in pixels */

static PaneL make_panel(int type,GrouP g,int width,int height,int firstseq,int nseqs)
{
	int i,l,length=0;
	PaneL p;
	panel_data data;

	data.type=type;
	SelectFont(datafont);
	data.lineheight=LineHeight();
	data.charwidth=CharWidth('A');
	if(type==NAMES)
	{
/* find the maximum length of sequence name */
        	for (i=firstseq;i<=firstseq+nseqs-1;i++)
		{
           		l = strlen(names[i]);
           		if (length < l) length = l;
		}
		data.vcols=width/data.charwidth - MARGIN*2 - DNUMBER;
	}
	else
	{
        	for (i=firstseq;i<=firstseq+nseqs-1;i++)
           		if (length < seqlen_array[i]) length = seqlen_array[i];
		data.vcols=width/data.charwidth - MARGIN*2;
	}
 
	data.lines=NULL;
	data.nhead=0;
	data.nfoot=0;
	data.header=NULL;
	data.footer=NULL;
	data.consensus=NULL;
	data.colormask=NULL;
	data.vlines=(height-SCOREHEIGHT)/data.lineheight - MARGIN;
	data.vseqs=data.vlines-data.nhead-data.nfoot;
	data.nseqs=nseqs;
	data.ncols=length;
	data.firstseq=firstseq-1;
	data.firstvline=0;
	data.firstvcol=0;
	data.lockoffset=0;
	data.ascent=Ascent();
	data.descent=Descent();
	data.selected=NULL;
	data.firstsel=-1;
	data.lastsel=-1;
	data.colscore=NULL;
	data.seqweight=NULL;
	data.subgroup=NULL;
	data.residue_exception=NULL;
	data.segment_exception=NULL;
	data.vscrollbar=NULL;
	data.hscrollbar=NULL;

	p=AutonomousPanel(g, width, height, DrawPanel, NULL,NULL,sizeof(panel_data), NULL, NULL);

	SetPanelExtra(p, &data);
	return p;

}
 
void DrawPanel(PaneL p)
{
	RecT r;
        panel_data data;
	int pixelwidth,pixelheight;

	UseWindow(mainw);
	Select(p);

	if (fromvscroll==0 && fromhscroll==0)
	{
		ObjectRect(p,&r);
        	pixelwidth=r.right-r.left;
        	pixelheight=r.bottom-r.top;

		SelectFont(datafont);
		GetPanelExtra(p, &data);
		data.lineheight=LineHeight();
		data.charwidth=CharWidth('A');
		if (data.type==NAMES)
			data.vcols=pixelwidth/data.charwidth-MARGIN*2-DNUMBER;
		else
			data.vcols=pixelwidth/data.charwidth-MARGIN*2;
		data.vlines=(pixelheight-SCOREHEIGHT)/data.lineheight - MARGIN;
		data.vseqs=data.vlines-data.nhead-data.nfoot;
		if(data.vseqs<0)data.vseqs=0;
		if(data.vcols<0)data.vcols=0;
		SetPanelExtra(p, &data);
/* draw the outside frame */
		ObjectRect (p, &r);
		Black();
		FrameRect(&r);
		InsetRect(&r,1,1);
		black_on_white();
		EraseRect(&r);
		if(data.nseqs == 0) return;
	}

/* draw the structure and gap penalty data */
/* draw the footer */
	if (fromvscroll==0)
	{
		draw_header(p);
		draw_footer(p);
		draw_colscores(p);
	}

/* draw the data lines */
	if (data.type==NAMES)
		draw_names(p);
	else
		draw_seqs(p);


}

void hscrollnames(BaR bar, int newval, int oldval)
{
	PaneL		p;
        panel_data        data;
 
	p = active_panel.names;
        GetPanelExtra(p, &data);
        data.firstvcol = newval;
        SetPanelExtra(p, &data);
        Select(p);
 
	if (data.vseqs<=0) return;
	draw_names(p);
}

void vscrollnames(BaR bar, int newval, int oldval)
{
	PaneL		p;
        panel_data        data;
 
	p = active_panel.names;
        GetPanelExtra(p, &data);
        data.firstvline = newval;
        SetPanelExtra(p, &data);
        Select(p);
 
	if (data.vseqs<=0) return;
	draw_names(p);
}

void vscrollseqs(BaR bar, int newval, int oldval)
{
	PaneL		p;
        panel_data        data;
        RecT            block,rect;
	int 		l;
 
	p = active_panel.seqs;
        GetPanelExtra(p, &data);
	l=data.firstvline;
        data.firstvline = newval;
        SetPanelExtra(p, &data);
        Select(p);
 
	if (data.vseqs<=0) return;

	if (data.vseqs<3 || data.nseqs-l < data.vseqs)
	{
		fromvscroll=0;
		draw_seqs(p);
		return;
	}

        if (newval == oldval + 1) {
		fromvscroll=1;
                ObjectRect(p, &rect);
		InsetRect(&rect,1,1);
                block.top = rect.top+(data.nhead)*data.lineheight+data.descent+1;
                block.bottom = block.top+(data.vseqs)*data.lineheight;
		block.left=rect.left;
		block.right=rect.right;
                ScrollRect(&block, 0, -data.lineheight);
        } else if (newval == oldval - 1) {
		fromvscroll=-1;
                ObjectRect(p, &rect);
		InsetRect(&rect,1,1);
                block.top = rect.top+(data.nhead)*data.lineheight+data.descent+1;
                block.bottom = block.top+(data.vseqs)*data.lineheight;
		block.left=rect.left;
		block.right=rect.right;
                ScrollRect(&block, 0, data.lineheight);
        } else {
		fromvscroll=0;
        }
	draw_seqs(p);
}

void hscrollseqs(BaR bar, int newval, int oldval)
{
	PaneL p;
        panel_data        data;
        RecT            rect;
 
 
	p = active_panel.seqs;
        GetPanelExtra(p, &data);
        data.firstvcol = newval+data.lockoffset;
        SetPanelExtra(p, &data);
        Select(p);
 
	if (data.vcols<=0) return;

	if (data.vcols<3)
	{
		fromhscroll=0;
		draw_header(p);
		draw_seqs(p);
		draw_footer(p);
		draw_colscores(p);
		return;
	}
        if (newval == oldval + 1) {
		fromhscroll=1;
                ObjectRect(p, &rect);
                InsetRect(&rect,1,1);
                rect.left+=data.charwidth;
                ScrollRect(&rect, -data.charwidth, 0);
        } else if (newval == oldval - 1) {
		fromhscroll=-1;
                ObjectRect(p, &rect);
                InsetRect(&rect,1,1);
                rect.right=rect.left+(data.vcols+1)*data.charwidth;
                ScrollRect(&rect, data.charwidth, 0);
        } else {
		fromhscroll=0;
        }
	draw_header(p);
	draw_seqs(p);
	draw_footer(p);
	draw_colscores(p);
}

void draw_names(PaneL p)
{
	int i,f,l;
	panel_data data;

	UseWindow(mainw);
	Select(p);
	GetPanelExtra(p,&data);
	if(data.lines==NULL) return;
	SelectFont(datafont);
	
	if (fromvscroll==0)
	{
		f=data.firstvline;
		l=data.firstvline+data.vseqs-1;
	}
	else if (fromvscroll==-1)
		f=l=data.firstvline;
	else
		f=l=data.firstvline+data.vseqs-1;
	
	if(l>=data.nseqs) l=data.nseqs-1;
	for(i=f;i<=l;i++)
		if (data.selected[i]==TRUE)
			draw_nameline(p,i,i,HIGHLIGHT);
		else
			draw_nameline(p,i,i,NORMAL);
}

void draw_seqs(PaneL p)
{
	int i,f,l,s,x,y,format;
	int fs,ls;
	panel_data data;
	PoinT pt;
	RecT r,block;

	UseWindow(mainw);
	Select(p);
	GetPanelExtra(p,&data);
	if(data.lines==NULL) return;
	SelectFont(datafont);
	black_on_white();
	if (fromhscroll==-1)
	{
		f=data.firstvcol;
		if ((f>=data.firstsel) && (f<=data.lastsel))
			format=HIGHLIGHT;
		else format=NORMAL; 
		draw_seqcol(p,f,format);
	}
	else if (fromhscroll==1)
	{
		f=data.firstvcol+data.vcols-1;
		if ((f>=data.firstsel) && (f<=data.lastsel))
			format=HIGHLIGHT;
		else format=NORMAL; 
		draw_seqcol(p,f,format);
	}
	else
	{
 		if (fromvscroll==-1)
		{
			f=l=data.firstvline;
		}
		else if (fromvscroll==1)
		{
			f=l=data.firstvline+data.vseqs-1;
		}
		else
		{
			f=data.firstvline;
			l=data.firstvline+data.vseqs-1;
		}
	
		if(l>=data.nseqs) l=data.nseqs-1;
        	s=f-data.firstvline;
        	ObjectRect (p, &r);
        	InsetRect(&r,1,1);
		data_colors();
        	block.top=r.top+((s+data.nhead)*data.lineheight)+data.descent+1;
        	block.bottom=block.top+(l-f+1)*data.lineheight;
        	block.left=r.left;
        	block.right=r.right;
        	EraseRect(&block);
        	if(data.nseqs == 0) return;

		if(data.firstsel != -1)
		{
			if ((data.firstsel>=data.firstvcol && data.firstsel<data.firstvcol+data.vcols)||
	   		(data.lastsel>=data.firstvcol && data.lastsel<data.firstvcol+data.vcols))
			{
				fs=data.firstsel-data.firstvcol;
				if (fs<0) fs=0;
				if (fs>=data.vcols) fs=data.vcols-1;
				ls=data.lastsel-data.firstvcol;
				if (ls<0) ls=0;
				if (ls>=data.vcols) ls=data.vcols-1;
        			block.left=r.left+(fs+1)*data.charwidth;
        			block.right=r.left+(ls+2)*data.charwidth;
				text_colors();
        			EraseRect(&block);
			}
		}
        	x=r.left+data.charwidth;
	 
        	for(i=f;i<=l;i++)
        	{
               		y=block.top+(i-f+1)*data.lineheight-data.descent-1;
			LoadPt(&pt,x,y);
               		draw_seqline(data,i,pt,data.firstvcol,data.firstvcol+data.vcols-1,NORMAL);
        	}
	}

	black_on_white();
	fromvscroll=fromhscroll=0;
}

static void NameClick(PaneL panel, PoinT pt)
{
	int i;
	panel_data data;
	RecT r;

	GetPanelExtra(panel,&data);
	if(data.prf_no==1)
	{
/* revert selected area in profile 2 to normal */
		GetPanelExtra(prf_panel[1].names,&data);
		if(data.nseqs==0)
			draw_seq_pointer(prf_panel[1].names,0,NORMAL);
		for(i=0;i<data.nseqs;i++)
			if (data.selected[i]==TRUE)
				draw_nameline(prf_panel[1].names,i,i,NORMAL);
		SetPanelExtra(prf_panel[1].names,&data);
	}
	else if(data.prf_no==2)
	{
/* revert selected area in profile 1 to normal */
		GetPanelExtra(prf_panel[0].names,&data);
		if(data.nseqs==0)
			draw_seq_pointer(prf_panel[0].names,0,NORMAL);
		for(i=0;i<data.nseqs;i++)
			if (data.selected[i]==TRUE)
				draw_nameline(prf_panel[0].names,i,i,NORMAL);
		SetPanelExtra(prf_panel[0].names,&data);
	}
	GetPanelExtra(panel,&data);
	Select(panel);
	ObjectRect(panel,&r);
	if (!shftKey)
	{
/* revert existing selected area to normal */
		for(i=0;i<data.nseqs;i++)
			if (data.selected[i]==TRUE)
				draw_nameline(panel,i,i,NORMAL);
	}

	selected_seqs.first = (pt.y - r.top-data.lineheight/2)/data.lineheight + data.firstvline-data.nhead;
	if (selected_seqs.first <0) selected_seqs.first=0;
	if (selected_seqs.first >=data.nseqs) selected_seqs.first=data.nseqs-1;
	if (selected_seqs.first==-1 && ncutseqs > 0)
	{
		selected_seqs.last=selected_seqs.first=0;
		draw_seq_pointer(panel,0,HIGHLIGHT);
	}
	else
	{
		selected_seqs.last=selected_seqs.first;
		draw_nameline(panel,selected_seqs.first,selected_seqs.last,HIGHLIGHT);
	}
	black_on_white();

}

static void NameDrag(PaneL panel, PoinT pt)
{
	panel_data data;
	RecT r;
	int s;

	GetPanelExtra(panel,&data);
	Select(panel);
	ObjectRect(panel,&r);
	s = (pt.y - r.top-data.lineheight/2)/data.lineheight + data.firstvline-data.nhead;
	if (s<0) s=0;
	if (s>=data.nseqs) s=data.nseqs-1;
	if (s==selected_seqs.first)
	{
		if (s!=selected_seqs.last)
		{
			draw_nameline(panel,selected_seqs.first,selected_seqs.last,NORMAL);
			draw_nameline(panel,selected_seqs.first,s,HIGHLIGHT);
		}
	}
	else if (s>selected_seqs.first)
	{
		if (s>selected_seqs.last)
			draw_nameline(panel,selected_seqs.last+1,s,HIGHLIGHT);
		else if (s<selected_seqs.last)
			draw_nameline(panel,s+1,selected_seqs.last,NORMAL);
	}
	else
	{
		if (s<selected_seqs.last)
			draw_nameline(panel,s,selected_seqs.last-1,HIGHLIGHT);
		else if (s>selected_seqs.last)
			draw_nameline(panel,selected_seqs.last,s-1,NORMAL);
	}
	selected_seqs.last=s;

	black_on_white();
}

static void NameRelease(PaneL panel, PoinT pt)
{
	int t;
	panel_data data;

	if (selected_seqs.first > selected_seqs.last)
	{
		t=selected_seqs.first;
		selected_seqs.first=selected_seqs.last;
		selected_seqs.last=t;
	}	
	active_panel.names = panel;
	GetPanelExtra(panel,&data);
	active_panel.seqs = data.index;

}

void draw_seq_pointer(PaneL panel,int seq,int format)
{
	RecT r,block;
	panel_data data;

	Select(panel);
	GetPanelExtra(panel,&data);

	ObjectRect(panel,&r);
	InsetRect(&r,1,1);
	block.top=r.top+((seq+data.nhead)*data.lineheight)+data.descent+1;
	block.bottom=block.top+data.lineheight;
	block.left=r.left;
	block.right=r.right;
	if (format==HIGHLIGHT)
		Black();
	else
		White();
	PaintRect(&block);

}

static void SeqClick(PaneL panel, PoinT pt)
{
	int s;
	int f,l;
	panel_data data;
	RecT r;

	GetPanelExtra(panel,&data);
	if(data.prf_no==1)
	{
/* revert selected area in profile 2 to normal */
		GetPanelExtra(prf_panel[1].seqs,&data);
		f=data.firstsel;
		l=data.lastsel;
		data.firstsel=-1;
		data.lastsel=-1;
		SetPanelExtra(prf_panel[1].seqs,&data);
		if (f != -1) highlight_seqrange(prf_panel[1].seqs,f,l,NORMAL);
	}
	else if(data.prf_no==2)
	{
/* revert selected area in profile 1 to normal */
		GetPanelExtra(prf_panel[0].seqs,&data);
		f=data.firstsel;
		l=data.lastsel;
		data.firstsel=-1;
		data.lastsel=-1;
		SetPanelExtra(prf_panel[0].seqs,&data);
		if (f != -1) highlight_seqrange(prf_panel[0].seqs,f,l,NORMAL);
	}
	GetPanelExtra(panel,&data);
	Select(panel);
	ObjectRect(panel,&r);

	s = (pt.x - r.left-data.charwidth)/data.charwidth + data.firstvcol;
	if (s <0) s=0;
	if (s<data.firstvcol) s=data.firstvcol;
	if (s >=data.ncols) s=data.ncols-1;
	if (s >=data.firstvcol+data.vcols) s=data.firstvcol+data.vcols-1;

	if (shftKey && data.firstsel != -1)
	{
		if (s>data.lastsel)
		{
			highlight_seqrange(panel,data.firstsel,s,HIGHLIGHT);
			data.lastsel=s;
		}
		else if (s<data.firstsel)
		{
			highlight_seqrange(panel,s,data.lastsel,HIGHLIGHT);
			data.firstsel=s;
		}
		else
		{
			highlight_seqrange(panel,s+1,data.lastsel,NORMAL);
			highlight_seqrange(panel,data.firstsel,s,HIGHLIGHT);
			data.lastsel=s;
		}
		selected_res.first=data.firstsel;
		selected_res.last=data.lastsel;
	}
	else
	{
/* revert existing selected area to normal */
		f=data.firstsel;
		l=data.lastsel;
		data.firstsel=-1;
		data.lastsel=-1;
		SetPanelExtra(panel,&data);
		if (f != -1) highlight_seqrange(panel,f,l,NORMAL);
		selected_res.first=selected_res.last=s;
		highlight_seqrange(panel,selected_res.first,selected_res.last,HIGHLIGHT);
		data.firstsel=selected_res.first;
		data.lastsel=selected_res.last;
	}

	SetPanelExtra(panel,&data);
	black_on_white();

}

static void SeqDrag(PaneL panel, PoinT pt)
{
	panel_data data;
	RecT r;
	int s;

	GetPanelExtra(panel,&data);
	Select(panel);
	ObjectRect(panel,&r);
	s = (pt.x - r.left-data.charwidth)/data.charwidth + data.firstvcol;
	if (s<0) s=0;
	if (s<data.firstvcol) s=data.firstvcol;
	if (s>=data.ncols) s=data.ncols-1;
	if (s >=data.firstvcol+data.vcols) s=data.firstvcol+data.vcols-1;
	if (s==selected_res.first)
	{
		if (s!=selected_res.last)
		{
			highlight_seqrange(panel,selected_res.first,selected_res.last,NORMAL);
			highlight_seqrange(panel,selected_res.first,s,HIGHLIGHT);
		}
	}
	else if (s>selected_res.first)
	{
		if (s>selected_res.last)
			highlight_seqrange(panel,selected_res.last+1,s,HIGHLIGHT);
		else if (s<selected_res.last)
			highlight_seqrange(panel,s+1,selected_res.last,NORMAL);
	}
	else
	{
		if (s<selected_res.last)
			highlight_seqrange(panel,s,selected_res.last-1,HIGHLIGHT);
		else if (s>selected_res.last)
			highlight_seqrange(panel,selected_res.last,s-1,NORMAL);
	}
	selected_res.last=s;

	black_on_white();
}

static void SeqRelease(PaneL panel, PoinT pt)
{
	int t;
	panel_data data;

        if (selected_res.first > selected_res.last)
        {
                t=selected_res.first;
                selected_res.first=selected_res.last;
                selected_res.last=t;
        }

	active_panel.seqs = panel;
	GetPanelExtra(panel,&data);
	active_panel.names = data.index;
	data.firstsel=selected_res.first;
	data.lastsel=selected_res.last;
	SetPanelExtra(panel,&data);

}

void draw_header(PaneL p)
{ 
	RecT  block,r;
	PoinT pt;
	int i, j, x, y;
	panel_data data;
	char *line;

	UseWindow(mainw);
	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	if(data.header == NULL) return;
	if(data.vlines<data.nhead) return;
	if(data.vcols<=0) return;

	line=(char *)ckalloc((data.vcols+1) * sizeof(char));
	ObjectRect (p, &r);
	InsetRect(&r,1,1);
	block.top=r.top+data.descent/2;
	block.bottom=block.top+(data.nhead*data.lineheight);
	block.left=r.left;
	block.right=r.right;
	text_colors();
	EraseRect(&block);
	if (data.type==NAMES)
        	x=r.left+DNUMBER*data.charwidth;
	else
        	x=r.left+data.charwidth;
        y=r.top+data.lineheight-data.descent/2;
	for(i=0;i<data.nhead;i++)
	{
		for(j=data.firstvcol;j<data.firstvcol+data.vcols && j<data.ncols;j++)
			if(j>=0)
				line[j-data.firstvcol]=data.header[i][j];
			else
				line[j-data.firstvcol]=' ';
		line[j-data.firstvcol]='\0';
		LoadPt(&pt, x, y);
		SetPen(pt);
		PaintString(line);
		y+=data.lineheight;
	}
	black_on_white();
	ckfree(line);
} 

void draw_footer(PaneL p)
{ 
	RecT  block,r;
	PoinT pt;
	int i, j,x, y;
	panel_data data;
	char *line;

	UseWindow(mainw);
	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	if(data.footer == NULL) return;
	if(data.vlines<data.nfoot) return;
	if(data.vcols<=0) return;

	line=(char *)ckalloc((data.vcols+1) * sizeof(char));
	ObjectRect (p, &r);
	InsetRect(&r,1,1);
	block.top=r.top+((data.vlines-data.nfoot)*data.lineheight)+data.descent+data.ascent/2;
	block.bottom=block.top+data.nfoot*data.lineheight;
	block.left=r.left;
	block.right=r.right;
	text_colors();
	EraseRect(&block);
	if(data.type==NAMES)
        	x=block.left+DNUMBER*data.charwidth;
	else
        	x=block.left+data.charwidth;
        y=block.top+data.lineheight-1;
	for(i=0;i<data.nfoot;i++)
	{
		for(j=data.firstvcol;j<data.firstvcol+data.vcols && j<data.ncols;j++)
			if(j>=0)
				line[j-data.firstvcol]=data.footer[i][j];
			else
				line[j-data.firstvcol]=' ';
		line[j-data.firstvcol]='\0';
		LoadPt(&pt, x, y);
		SetPen(pt);
		PaintString(line);
		y+=data.lineheight;
	}
	black_on_white();
	ckfree(line);
} 


void draw_nameline(PaneL p,int fseq,int lseq,int format)
{ 
	RecT  block,r;
	PoinT pt;
	int n,i, j, t, f,l,x, y,ix;
	panel_data data;
	char *line;

	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	
	n=1;
	i=data.nseqs;
	for(;;)
	{
		i/=10;
		if(i==0) break;
		n++;
	}

	line=(char *)ckalloc((data.vcols+1) * sizeof(char));
	if (fseq > lseq)
	{
		t=fseq;
		fseq=lseq;
		lseq=t;
	}	
	if (format==HIGHLIGHT)
		for(i=fseq;i<=lseq;i++) data.selected[i]=TRUE;
	else
		for(i=fseq;i<=lseq;i++) data.selected[i]=FALSE;
	SetPanelExtra(p,&data);
	if (fseq<data.firstvline)
		fseq=data.firstvline;
	if (fseq>=data.firstvline+data.vseqs)
		fseq=data.firstvline+data.vseqs;
	if (lseq<data.firstvline)
		lseq=data.firstvline;
	if (lseq>=data.firstvline+data.vseqs)
		lseq=data.firstvline+data.vseqs-1;
	f=fseq-data.firstvline;
	l=lseq-data.firstvline;
	ObjectRect (p, &r);
	InsetRect(&r,1,1);
	block.top=r.top+((f+data.nhead)*data.lineheight)+data.descent+1;
	block.bottom=block.top+((l-f+1)*data.lineheight);
	block.left=r.left;
	block.right=r.right;
	if (format==HIGHLIGHT)
		white_on_black();
	else
		data_colors();
	EraseRect(&block);
        y=block.top+data.lineheight-data.descent-1;
	for(i=fseq;i<=lseq;i++)
	{
        	x=r.left+data.charwidth;
		sprintf(line,"%*d",n,i+1);
		LoadPt(&pt, x, y);
		SetPen(pt);
		Gray();
		PaintString(line);
		y+=data.lineheight;
	}
        y=block.top+data.lineheight-data.descent-1;
	for(i=fseq;i<=lseq;i++)
	{
		ix=output_index[i+1]-1;
        	x=r.left+DNUMBER*data.charwidth;
		for(j=0;j<data.vcols && j<data.ncols-data.firstvcol;j++)
			line[j]=data.lines[ix][j+data.firstvcol];
		line[j]='\0';
		LoadPt(&pt, x, y);
		SetPen(pt);
		if(format==HIGHLIGHT) White();
		else Black();
		PaintString(line);
		y+=data.lineheight;
	}
	black_on_white();
	ckfree(line);
} 

void draw_seqline(panel_data data,int seq,PoinT pt,int fcol,int lcol,int format)
{ 
	RecT r;
	int i, j, ix;
	char *line[MAXCOLORS+1];

	if(data.nseqs == 0) return;

/* draw colored character on white background */
	for(i=0;i<ncolors;i++)
	{
		line[i]=(char *)ckalloc((data.vcols+1) * sizeof(char));
		for(j=0;j<data.vcols;j++)
		line[i][j]=' ';
		line[i][j]='\0';
	}
	
	ix=output_index[seq+1]-1;
	
	r.top=pt.y-data.lineheight+data.descent+1;
	r.bottom=r.top+data.lineheight;
	for(j=fcol;j<=lcol && j<data.ncols;j++)
	{
		if(j>=0)
		{
		if(segment_exceptions && data.segment_exception[ix][j] > 0)
		{
			r.left=pt.x;
			r.right=r.left+data.charwidth;
			DkGray();
			PaintRect(&r);
			White();
		}
		else if(residue_exceptions && data.residue_exception[ix][j] == TRUE)
		{
			r.left=pt.x;
			r.right=r.left+data.charwidth;
		/*	LtGray(); */
			SelectColor(150,150,150);
			PaintRect(&r);
			White();
		}
		else
		{
                        if(inverted)
                        {
				if(format==HIGHLIGHT || (j>=data.firstsel && j<=data.lastsel))
					Black();
				else
				{
                                	r.left=pt.x;
#ifdef UNIX
                	                r.right=r.left+data.charwidth-1;
#else
                	                r.right=r.left+data.charwidth;
#endif
                                 	SetColor(color_lut[(int)data.colormask[ix][j]].val);
                       	         	PaintRect(&r);
                       	         	Black();
				}
                        }
                        else
                                SetColor(color_lut[(int)data.colormask[ix][j]].val);

		}
		SetPen(pt);
		PaintChar(data.lines[ix][j]);
		}
		pt.x+=data.charwidth;
	}
	for(i=0;i<ncolors;i++)
		ckfree(line[i]);
	Black();
} 

void draw_seqcol(PaneL p,int col,int format)
{ 
	RecT  block,r, r2;
	PoinT pt;
	int totseqs,i, c,x,y,ix;
	panel_data data;

	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	if(data.ncols == 0) return;

	SetPanelExtra(p, &data);
	
	if (col<data.firstvcol)
		col=data.firstvcol;
	if (col>=data.firstvcol+data.vcols)
		col=data.firstvcol+data.vcols-1;
	c=col-data.firstvcol;
	totseqs=data.vseqs;
	if (totseqs>data.nseqs) totseqs=data.nseqs;
	ObjectRect (p, &r);
	InsetRect(&r,1,1);
	block.top=r.top+(data.nhead*data.lineheight)+data.descent+1;
	block.bottom=block.top+(totseqs)*data.lineheight;
	block.left=r.left+(c+1)*data.charwidth;
	block.right=block.left+data.charwidth;
        if (format==HIGHLIGHT)
                text_colors();
        else
                data_colors();
	EraseRect(&block);

	x=r.left+(c+1)*data.charwidth;
       	y=block.top+data.lineheight-data.descent-1;
	r2.left=x;
	r2.right=r2.left+data.charwidth;
	for(i=data.firstvline;i<data.firstvline+data.vseqs && i<data.nseqs;i++)
	{	
		ix=output_index[i+1]-1;
		if(segment_exceptions && data.segment_exception[ix][col] > 0)
		{
			r2.top=y-data.lineheight+data.descent+1;
			r2.bottom=r.top+data.lineheight;
			DkGray();
			PaintRect(&r2);
			White();
		}
		else if(residue_exceptions && data.residue_exception[ix][col] == TRUE)
		{
			r2.top=y-data.lineheight+data.descent+1;
			r2.bottom=r.top+data.lineheight;
		/*	LtGray(); */
			SelectColor(150,150,150);
			PaintRect(&r2);
			White();
		}
		else
		{
                        if(inverted)
                        {
                                r2.top=y-data.lineheight+data.descent+1;
                                r2.bottom=r2.top+data.lineheight;
                                if(format==HIGHLIGHT)
				{
                                        LtGray();
				}
                                else
                                SetColor(color_lut[(int)data.colormask[ix][col]].val);
                                PaintRect(&r2);
                                Black();
                        }
                        else
                                SetColor(color_lut[(int)data.colormask[ix][col]].val);

		}
		LoadPt(&pt,x,y);
		SetPen(pt);
		PaintChar(data.lines[ix][col]);
		y+=data.lineheight;
	}
	Black();
} 

void highlight_seqrange(PaneL p,int fcol,int lcol, int format)
{ 
	RecT  block,r;
	int i,t,x,y;
	int fseq,lseq,s;
	panel_data data;
	PoinT pt;

	Select(p);
	SelectFont(datafont);
	GetPanelExtra(p, &data);
	if(data.nseqs == 0) return;
	if(data.ncols == 0) return;

        if (fcol > lcol)
        {
                t=fcol;
                fcol=lcol;
                lcol=t;
        }

	if ((fcol>=data.firstvcol && fcol<data.firstvcol+data.vcols)||
	   (lcol>=data.firstvcol && lcol<data.firstvcol+data.vcols))
	{
		if (fcol<data.firstvcol) fcol=data.firstvcol;
		if (fcol>=data.firstvcol+data.vcols) fcol=data.firstvcol+data.vcols-1;
		if (lcol<data.firstvcol) lcol=data.firstvcol;
		if (lcol>=data.firstvcol+data.vcols) lcol=data.firstvcol+data.vcols-1;
	}
 
	fseq=data.firstvline;
	lseq=data.firstvline+data.vseqs-1;
	if(lseq>=data.nseqs) lseq=data.nseqs-1;
        s=fseq-data.firstvline;
        ObjectRect (p, &r);
        InsetRect(&r,1,1);
	if(format==HIGHLIGHT)
		text_colors();
	else
		data_colors();
       	block.top=r.top+((s+data.nhead)*data.lineheight)+data.descent+1;
       	block.bottom=block.top+(lseq-fseq+1)*data.lineheight;
       	block.left=r.left+(fcol-data.firstvcol+1)*data.charwidth;
       	block.right=r.left+(lcol-data.firstvcol+2)*data.charwidth;
        EraseRect(&block);

       	x=r.left+(fcol-data.firstvcol+1)*data.charwidth;
	 
       	for(i=fseq;i<=lseq;i++)
       	{
               	y=block.top+(i-fseq+1)*data.lineheight-data.descent-1;
		LoadPt(&pt,x,y);
               	draw_seqline(data,i,pt,fcol,lcol,format);
       	}
	black_on_white();
} 

GrouP make_scroll_area(GrouP w,int prf_no,int nwidth,int swidth,int height,int firstseq,int nseqs,spanel *p)
{
	panel_data ndata,sdata;
        GrouP display;
	RecT rect;
        PoinT pt;
	PaneL names,seqs;
	BaR vscrollbar,hnscrollbar,hsscrollbar;
	BarScrlProc hscrollnameproc, hscrollseqproc, vscrollproc;

	if(prf_no==0)
	{
		hscrollnameproc=HscrollMultiN;
		hscrollseqproc=HscrollMultiS;
		vscrollproc=VscrollMulti;
	}
	else if (prf_no==1)
	{
		hscrollnameproc=HscrollPrf1N;
		hscrollseqproc=HscrollPrf1S;
		vscrollproc=VscrollPrf1;
	}
	else
	{
		hscrollnameproc=HscrollPrf2N;
		hscrollseqproc=HscrollPrf2S;
		vscrollproc=VscrollPrf2;
	}

        display=HiddenGroup(w, 0, 0, NULL);
        SetGroupSpacing(display, 0, 0);
	Hide(display);

        vscrollbar=ScrollBar(display, -1, 1, vscrollproc);

        ObjectRect(vscrollbar, &rect);
        pt.x=rect.right;
        pt.y=rect.top;
        SetNextPosition(display, pt);
        names=make_panel(NAMES,display, nwidth, height, firstseq,nseqs);

        ObjectRect(names, &rect);
        pt.x=rect.right;
        pt.y=rect.top;
        SetNextPosition(display, pt);
        seqs=make_panel(SEQS,display, swidth, height, firstseq,nseqs);

/* horizontal scroll bars */
        ObjectRect(names, &rect);
        pt.x=rect.left;
        pt.y=rect.bottom;
        SetNextPosition(display, pt);
        hnscrollbar=ScrollBar(display, 1, -1, hscrollnameproc);
        ObjectRect(seqs, &rect);
        pt.x=rect.left;
        pt.y=rect.bottom;
        SetNextPosition(display, pt);
        hsscrollbar=ScrollBar(display, 1, -1, hscrollseqproc);
        
	SetRange(hsscrollbar,1,1,0);
	SetRange(hnscrollbar,1,1,0);
	SetRange(vscrollbar,1,1,0);

	GetPanelExtra(names,&ndata);
	ndata.hscrollbar=hnscrollbar;
	ndata.index=seqs;
        ndata.prf_no=prf_no;

	GetPanelExtra(seqs,&sdata);
	sdata.vscrollbar=vscrollbar;
	sdata.hscrollbar=hsscrollbar;
	sdata.index=names;
        sdata.prf_no=prf_no;

	SetPanelClick(names,NameClick, NameDrag, NULL, NameRelease);
	SetPanelClick(seqs,SeqClick, SeqDrag, NULL, SeqRelease);

	p->names = names;
	p->seqs = seqs;

	ndata=alloc_name_data(ndata);
	sdata=alloc_seq_data(sdata);
	SetPanelExtra(names,&ndata);
	SetPanelExtra(seqs,&sdata);

	Show(display);
	return(display);
}


void white_on_black(void)
{
	Black(); InvertColors(); White();
}
void black_on_white(void)
{
	White(); InvertColors(); Black();
}
void text_colors(void)
{
	SelectColor(220,220,220);
	InvertColors();
	Black();
}
void data_colors(void)
{
	White();
	InvertColors();
	Black();
}




void make_ruler(int length, char *name,char *seq)
{

	int i,j;
	char marker[5];
	int marker_len;

	strcpy(name,"ruler");
	seq[0] = '1';
	for (i=1;i<length;i++)
	{
		if ((i+1)%10 > 0)
			seq[i] = '.';
		else
		{
			sprintf(marker,"%d",((i+1)/10)*10);
			marker_len = strlen(marker);
			for (j=0;j<marker_len && i+1+j-marker_len < length;j++)
				seq[i+1+j-marker_len] = marker[j];
		}
	}
	seq[length]='\0';
}

panel_data free_panel_data(panel_data data)
{
	int i;

	if (data.header!=NULL)
	{
		for (i=0;i<mheader;i++)
		{
			if(data.header[i] != NULL) ckfree(data.header[i]);
			data.header[i]=NULL;
		}
		ckfree(data.header);
		data.header=NULL;
	}
	if (data.footer!=NULL)
	{
		for (i=0;i<mfooter;i++)
		{
			if(data.footer[i] != NULL) ckfree(data.footer[i]);
			data.footer[i]=NULL;
		}
		ckfree(data.footer);
		data.footer=NULL;
	}
	if (data.consensus!=NULL)
	{
		ckfree(data.consensus);
		data.consensus=NULL;
	}
	if (data.lines!=NULL)
	{
		for (i=0;i<data.nseqs;i++)
		{
			if(data.lines[i] != NULL) ckfree(data.lines[i]);
			data.lines[i]=NULL;
		}
		ckfree(data.lines);
		data.lines=NULL;
	}
	if (data.colormask!=NULL)
	{
		for (i=0;i<data.nseqs;i++)
		{
			if(data.colormask[i] != NULL) ckfree(data.colormask[i]);
			data.colormask[i]=NULL;
		}
		ckfree(data.colormask);
		data.colormask=NULL;
	}
	if (data.selected!=NULL) ckfree(data.selected);
	data.selected=NULL;

	if (data.seqweight!=NULL) ckfree(data.seqweight);
	data.seqweight=NULL;
	if (data.subgroup!=NULL) ckfree(data.subgroup);
	data.subgroup=NULL;
	if (data.colscore!=NULL) ckfree(data.colscore);
	data.colscore=NULL;
	if (data.residue_exception!=NULL)
	{
		for (i=0;i<data.nseqs;i++)
		{
			if(data.residue_exception[i] != NULL) ckfree(data.residue_exception[i]);
			data.residue_exception[i]=NULL;
		}
		ckfree(data.residue_exception);
		data.residue_exception=NULL;
	}
	if (data.segment_exception!=NULL)
	{
		for (i=0;i<data.nseqs;i++)
		{
			if(data.segment_exception[i] != NULL) ckfree(data.segment_exception[i]);
			data.segment_exception[i]=NULL;
		}
		ckfree(data.segment_exception);
		data.segment_exception=NULL;
	}

	return(data);
}


void make_consensus(panel_data data,char *name,char *seq1)
{
 	char c;
	sint catident1[NUMRES],catident2[NUMRES],ident;
	sint i,j,k,l;


	strcpy(name,"");    
    	for(i=0; i<data.ncols; i++) {
			seq1[i]=' ';
			ident=0;
			for(j=0;res_cat1[j]!=NULL;j++) catident1[j] = 0;
			for(j=0;res_cat2[j]!=NULL;j++) catident2[j] = 0;
			for(j=0;j<data.nseqs;++j) {
				if(isalpha(data.lines[0][i])) {
					if(data.lines[0][i] == data.lines[j][i])
					++ident;
					for(k=0;res_cat1[k]!=NULL;k++) {
					        for(l=0;(c=res_cat1[k][l]);l++) {
					        if (c=='\0') break;
							if (data.lines[j][i]==c)
							{
								catident1[k]++;
								break;
							}
						}
					}
					for(k=0;res_cat2[k]!=NULL;k++) {
					        for(l=0;(c=res_cat2[k][l]);l++) {
					        if (c=='\0') break;
							if (data.lines[j][i]==c)
							{
								catident2[k]++;
								break;
							}
						}
					}
				}
			}
			if(ident==data.nseqs)
				seq1[i]='*';
			else if (!dnaflag) {
				for(k=0;res_cat1[k]!=NULL;k++) {
					if (catident1[k]==data.nseqs) {
						seq1[i]=':';
						break;
					}
				}
				if(seq1[i]==' ')
				for(k=0;res_cat2[k]!=NULL;k++) {
					if (catident2[k]==data.nseqs) {
						seq1[i]='.';
						break;
					}
				}
			}
		}
}

int make_struct_data(int prf_no,int len, char *name,char *seq)
{
	int i,n=0;
	char val;
        char *ss_mask;
 
	seq[0]='\0';
	name[0]='\0';
if (prf_no == 1)
{
        if (struct_penalties1 == SECST && use_ss1 == TRUE) {
		n=1;
		strcpy(name,"Structures");
                ss_mask = (char *)ckalloc((seqlen_array[1]+10) * sizeof(char));
                for (i=0;i<seqlen_array[1];i++)
                        ss_mask[i] = sec_struct_mask1[i];
                print_sec_struct_mask(seqlen_array[1],sec_struct_mask1,ss_mask)
;
                for(i=0; i<len; i++) {
                        val=ss_mask[i];
                        if (val == gap_pos1)
                                seq[i]='-';
                        else
                                seq[i]=val;
                }
                seq[i]=EOS;
        	ckfree(ss_mask);
        }
 
}
else if (prf_no == 2)
{
        if (struct_penalties2 == SECST && use_ss2 == TRUE) {
		n=1;
		strcpy(name,"Structures");
                ss_mask = (char *)ckalloc((seqlen_array[profile1_nseqs+1]+10) *
sizeof(char));
                for (i=0;i<seqlen_array[profile1_nseqs+1];i++)
                        ss_mask[i] = sec_struct_mask2[i];
                print_sec_struct_mask(seqlen_array[profile1_nseqs+1],sec_struct_mask2,ss_mask);
       
                for(i=0; i<len; i++) {
                        val=ss_mask[i];
                        if (val == gap_pos1)
                                seq[i]='-';
                        else
                                seq[i]=val;
                }
                seq[i]=EOS;
       		ckfree(ss_mask);
       }
}
	return(n);
}

int make_gp_data(int prf_no,int len, char *name,char *seq)
{
	int i,n=0;
	char val;

	seq[0]='\0';
	name[0]='\0';
if (prf_no == 1)
{
        if (struct_penalties1 == GMASK && use_ss1 == TRUE) {
		n=1;
		strcpy(name,"Gap Penalties");
                for(i=0; i<len; i++) {
                        val=gap_penalty_mask1[i];
                        if (val == gap_pos1)
                                seq[i]='-';
                        else
                                seq[i]=val;
                }
                seq[i]=EOS;
        }
}
else if (prf_no == 2)
{
        if (struct_penalties2 == GMASK && use_ss2 == TRUE) {
		n=1;
		strcpy(name,"Gap Penalties");
                for(i=0; i<len; i++) {
                        val=gap_penalty_mask2[i];
                        if (val == gap_pos1)
                                seq[i]='-';
                        else
                                seq[i]=val;
                }
                seq[i]=EOS;
        }
}
	return(n);
}

static void VscrollMulti(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=seq_panel;
	vscrollnames(bar, newval, oldval);
	vscrollseqs(bar, newval, oldval);
}

static void HscrollMultiN(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=seq_panel;
	hscrollnames(bar, newval, oldval);
}

static void HscrollMultiS(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=seq_panel;
	hscrollseqs(bar, newval, oldval);
}

static void VscrollPrf1(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=prf_panel[0];
	vscrollnames(bar, newval, oldval);
	vscrollseqs(bar, newval, oldval);
}

static void HscrollPrf1N(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=prf_panel[0];
	hscrollnames(bar, newval, oldval);
}

static void HscrollPrf1S(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=prf_panel[0];
	hscrollseqs(bar, newval, oldval);
	if(fixed_prf_scroll==TRUE)
	{
		active_panel=prf_panel[1];
		hscrollseqs(bar, newval, oldval);
	}
}

static void VscrollPrf2(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=prf_panel[1];
	vscrollnames(bar, newval, oldval);
	vscrollseqs(bar, newval, oldval);
}

static void HscrollPrf2N(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	active_panel=prf_panel[1];
	hscrollnames(bar, newval, oldval);
}

static void HscrollPrf2S(BaR bar, GraphiC p, Nlm_Int2 newval, Nlm_Int2 oldval)
{
	if(fixed_prf_scroll==TRUE)
	{
		active_panel=prf_panel[0];
		hscrollseqs(bar, newval, oldval);
	}
	active_panel=prf_panel[1];
	hscrollseqs(bar, newval, oldval);
}


