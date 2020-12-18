
#include <stdarg.h>
#include <string.h>

#include <vibrant.h>
#include <document.h>

#include "clustalw.h"
#include "xmenu.h"


static void RemoveWin(WindoW w);
static void QuitWinW(WindoW w);
static void QuitWinI(IteM i);
static void QuitHelpW(WindoW w);
static void QuitHelpB(ButtoN b);
static void SearchStrWin (IteM item);
static void SavePSSeqWin (IteM item);
static void SavePSPrf1Win (IteM item);
static void SavePSPrf2Win (IteM item);
static void SaveSeqFileWin (IteM item);
static void SavePrf1FileWin (IteM item);
static void SavePrf2FileWin (IteM item);
static void OpenColorParWin (IteM item);
static void SearchStr(ButtoN but);
static void SavePSSeqFile(ButtoN but);
static void SavePSPrf1File(ButtoN but);
static void SavePSPrf2File(ButtoN but);
static void SaveSeqFile(ButtoN but);
static void SavePrf1File(ButtoN but);
static void SavePrf2File(ButtoN but);
static void SaveScoresWin (IteM item);
static void SaveScores(ButtoN but);
static void OpenColorPar(ButtoN but);
static void CancelWin(ButtoN but);
static void SaveTreeWin (IteM item);
static void CAlignWin (IteM item);
static void RealignSeqsWin (IteM item);
static void RealignSeqRangeWin (IteM item);
static void DrawTreeWin (IteM item);
static void AlignFromTreeWin(IteM item);
static void PrfPrfAlignWin(IteM item);
static void PrfPrfTreeAlignWin(IteM item);
static void SeqPrfAlignWin(IteM item);
static void SeqPrfTreeAlignWin(IteM item);
static void BootstrapTreeWin (IteM item);
static void CreateAlignTree(ButtoN but);
static void CompleteAlign(ButtoN but);
static void RealignSeqs(ButtoN but);
static void RealignSeqRange(ButtoN but);
static void DrawTree(ButtoN but);
static void AlignFromTree(ButtoN but);
static void PrfPrfAlign(ButtoN but);
static void PrfPrfTreeAlign(ButtoN but);
static void SeqPrfAlign(ButtoN but);
static void SeqPrfTreeAlign(ButtoN but);
static void BootstrapTree(ButtoN but);
static void OpenSeqFile (IteM item);
static void AppendSeqFile (IteM item);
static void OpenPrf1File (IteM item);
static void OpenPrf2File (IteM item);
static void ScoreWin(IteM item);
static void SegmentWin(IteM item);
static void ScoreSegments(ButtoN but);
static void PWParameters(IteM item);
static void MultiParameters(IteM item);
static void GapParameters(IteM item);
static void SSParameters(IteM item);
static void OutputParameters(IteM item);
static void OutputTreeParameters(IteM item);
static void HelpProc(IteM item);
static void DefColorPar(IteM item);
static void BlackandWhite(IteM item);
static void set_reset_new_gaps(IteM i);
static void set_reset_all_gaps(IteM i);
static void SearchStringAgain(ButtoN but);

static PopuP make_toggle(GrouP g,CharPtr title,CharPtr true_text, CharPtr false_text,
				 Boolean *value,PupActnProc SetProc);
static PrompT make_scale(GrouP g,CharPtr title,int length,int value,int max,BarScrlProc SetProc);
static PrompT make_prompt(GrouP g,CharPtr title);

static void CutSequences(IteM item);
static void PasteSequences(IteM item);
static void RemoveGaps(IteM item);
static void RemoveGapPos(IteM item);

static void SelectSeqs(IteM item);
static void SelectPrf1(IteM item);
static void SelectPrf2(IteM item);
static void MergeProfiles(IteM item);
static void ClearSeqs(IteM item);
 
static void cut_multiplem(void);
static void cut_profile1(void);
static void cut_profile2(void);
static void ssave(int j);
static void sscpy(int i,int j);
static void sload(int i);
static void clear_seqrange(spanel p);
static void select_seqs(spanel p,Boolean flag);
static void clear_seg_exceptions(spanel p);

static void make_menu_headers(WindoW w);
static void make_help_menu(void);
static void make_score_menu(void);
static void make_file_menu(void);
static void make_edit_menu(void);
static void make_align_menu(void);
static void make_tree_menu(void);
static void make_color_menu(void);

static void save_aln_window(int prf_no,char *title,char *prompt,void save_proc(ButtoN but));
static void save_ps_window(int prf_no,char *prompt,void save_proc(ButtoN but));
static void read_file_window(char *title,char *prompt,char *filename,void read_proc(ButtoN but));
static void do_align_window(WindoW *alignw,TexT *treetext,Boolean treestatus,char *title,void align_proc(ButtoN but));
static void do_palign_window(WindoW *alignw,TexT *tree1text,TexT *tree2test,Boolean treestatus,char *title,void align_proc(ButtoN but));
static Boolean open_aln_files(void);
static void write_file(int fseq,int lseq,int fres,int lres);


Boolean x_menus=FALSE;

int    mheader = 2; /* maximum header lines */
int    mfooter = 1; /* maximum footer lines */
int max_mlines = 20;      /*   multiple align display length */
int min_mlines = 10;      /*   multiple align display length */
int max_plines = 8;     /*   profile align display length */
int min_plines1 = 5;     /*   profile align display length */
int min_plines2 = 3;     /*   profile align display length */

Boolean aln_mode = MULTIPLEM;
Boolean window_displayed = FALSE;

int    save_format = CLUSTAL;
Boolean fixed_prf_scroll = FALSE;
int loffset,boffset,toffset;
int roffset;
int poffset;

int score_cutoff=5;    /* cutoff for residue exceptions */
int score_hwin=5;    /* half window for summing alignment column scores */
int score_scale=5;
int segment_dnascale=5;
int length_cutoff=1;    /* length cutoff for segment exceptions */
Boolean residue_exceptions=FALSE;
Boolean segment_exceptions=FALSE;
int score_matnum=4;
char score_mtrxname[FILENAMELEN];
int segment_matnum=3;
char segment_mtrxname[FILENAMELEN];
int score_dnamatnum=1;
char score_dnamtrxname[FILENAMELEN];
int segment_dnamatnum=1;
char segment_dnamtrxname[FILENAMELEN];

Boolean output_ss;
Boolean output_gp;

extern char     revision_level[];
extern Boolean interactive;

extern char seqname[];
extern char     outfile_name[];
extern char    profile1_name[];
extern char    profile2_name[];
extern char     usermtrxname[], pw_usermtrxname[];
extern char     dnausermtrxname[], pw_dnausermtrxname[];

extern Boolean usemenu;
extern Boolean use_tree_file;
extern Boolean use_tree1_file,use_tree2_file;
extern Boolean  dnaflag;
extern sint     nseqs;
extern sint    profile1_nseqs;
extern sint profile_no;
extern sint     max_aa;
extern sint     *seqlen_array;
extern char     **seq_array;
extern char     **names, **titles;
extern Boolean empty;
extern Boolean profile1_empty, profile2_empty;
extern sint     gap_pos1, gap_pos2;
extern Boolean use_ambiguities;


extern float    gap_open,      gap_extend;
extern float    dna_gap_open,  dna_gap_extend;
extern float    prot_gap_open, prot_gap_extend;
extern float    pw_go_penalty,      pw_ge_penalty;
extern float    dna_pw_go_penalty,  dna_pw_ge_penalty;
extern float    prot_pw_go_penalty, prot_pw_ge_penalty;
extern sint    wind_gap,ktup,window,signif;
extern sint    dna_wind_gap, dna_ktup, dna_window, dna_signif;
extern sint    prot_wind_gap,prot_ktup,prot_window,prot_signif;
extern sint        helix_penalty;
extern sint        strand_penalty;
extern sint        loop_penalty;
extern sint        helix_end_minus;
extern sint        helix_end_plus;
extern sint        strand_end_minus;
extern sint        strand_end_plus;
extern sint        helix_end_penalty;
extern sint        strand_end_penalty;
extern sint     divergence_cutoff;
extern sint     gap_dist;
extern sint boot_ntrials;               /* number of bootstrap trials */
extern unsigned sint boot_ran_seed;     /* random number generator seed */

extern sint        matnum,pw_matnum;
extern char     mtrxname[], pw_mtrxname[];
extern sint        dnamatnum,pw_dnamatnum;
extern char     dnamtrxname[], pw_dnamtrxname[];

extern MatMenu matrix_menu;
extern MatMenu pw_matrix_menu;
extern MatMenu dnamatrix_menu;

extern Boolean  quick_pairalign;
extern sint        matnum,pw_matnum;
extern Boolean  neg_matrix;
extern float    transition_weight;
extern char     hyd_residues[];
extern Boolean  no_var_penalties, no_hyd_penalties, no_pref_penalties;
extern Boolean         use_endgaps;
extern Boolean         endgappenalties;
extern Boolean  output_clustal, output_nbrf, output_phylip, output_gcg, output_gde, output_nexus;
extern Boolean  save_parameters;
extern Boolean  output_tree_clustal, output_tree_phylip, output_tree_distances, output_tree_nexus;
extern Boolean  lowercase; /* Flag for GDE output - set on comm. line*/
extern Boolean  cl_seq_numbers;
extern sint     output_order;
extern sint     *output_index;
extern Boolean  reset_alignments_new;               /* DES */
extern Boolean  reset_alignments_all;               /* DES */

extern FILE     *clustal_outfile, *gcg_outfile, *nbrf_outfile, *phylip_outfile;
extern FILE     *gde_outfile, *nexus_outfile;


extern sint     max_aln_length;

extern Boolean tossgaps;  /* Ignore places in align. where ANY seq. has a gap*/
extern Boolean kimura;    /* Use correction for multiple substitutions */
extern sint    bootstrap_format;      /* bootstrap file format */

extern sint output_struct_penalties;
extern Boolean use_ss1, use_ss2;
extern char *res_cat1[];
extern char *res_cat2[];

extern char     *amino_acid_codes;

PrompT   message;           /* used in temporary message window */

static Char filename[FILENAMELEN]; /* used in temporary file selection window */

Boolean mess_output=TRUE;
Boolean save_log=FALSE;
FILE *save_log_fd=NULL;
static char save_log_filename[FILENAMELEN];
static IteM save_item1,save_item2,exc_item;

spanel  seq_panel;        /* data for multiple alignment area */
spanel  prf_panel[2];       /* data for profile alignment areas */
spanel  active_panel;       /* 'in-use' panel -scrolling,clicking etc. */
static range selected_seqs;           /* sequences selected by clicking on names */
static range selected_res;           /* residues selected by clicking on seqs */
int firstres, lastres;	/* range of alignment for saving as ... */
 
/* data for Search function */

char find_string[MAXFINDSTR]="";
aln_pos find_pos;

/* arrays for storing clustalw data for cut-and-paste sequences */
static sint     *saveseqlen_array=NULL;
static char     **saveseq_array=NULL;
static char     **savenames=NULL, **savetitles=NULL;
sint     ncutseqs=0;

FonT datafont,helpfont;
WindoW mainw=NULL;
WindoW messagew=NULL;
WindoW readfilew=NULL;
WindoW savealnw=NULL;
WindoW savescoresw=NULL;
WindoW savepsw=NULL;
WindoW findw=NULL;
WindoW calignw=NULL;
WindoW ralignw=NULL;
WindoW rralignw=NULL;
WindoW talignw=NULL;
WindoW palignw=NULL;
WindoW salignw=NULL;
WindoW scorew=NULL;
WindoW exceptionw=NULL;
TexT savealntext;
TexT savescorestext;
TexT savepstext;
TexT findtext;
TexT pspartext;
TexT ctreetext;
TexT rtreetext;
TexT rrtreetext;
TexT ttreetext;
TexT ptree1text,ptree2text;
TexT streetext;
TexT readfiletext;
WindoW savetreew=NULL;
TexT savetreetext;
WindoW drawtreew=NULL;
TexT drawnjtreetext;
TexT drawphtreetext;
TexT drawdsttreetext;
TexT drawnxstreetext;
WindoW boottreew=NULL;
TexT bootnjtreetext;
TexT bootphtreetext;
TexT bootnxstreetext;
TexT blocklentext;
PrompT mattext,pwmattext,dnamattext,pwdnamattext,scoremattext,segmentmattext;
PrompT scorednamattext,segmentdnamattext;
GrouP seg_matrix_list,score_matrix_list;
GrouP seg_dnamatrix_list,score_dnamatrix_list;
GrouP matrix_list,pw_matrix_list,dnamatrix_list,pw_dnamatrix_list;

TexT cl_outtext,pir_outtext,msf_outtext,phylip_outtext,gde_outtext,nexus_outtext;
GrouP slow_para,fast_para;
GrouP  seq_display,prf1_display,prf2_display;

MenU   filem,alignm,editm,treem,colorm;
menu_item file_item,align_item,edit_item,tree_item,color_item;
MenU   scorem,helpmenu;
menu_item score_item,help_item;
IteM segment_item;
IteM bw_item,defcol_item,usercol_item;
IteM new_gaps_item,all_gaps_item;
WindoW helpw[MAXHELPW];
int numhelp=0;

PopuP modetext,flisttext;
ButtoN pscrolltext;
PopuP show_seg_toggle;
PrompT residue_cutofftext;
PrompT length_cutofftext;
PrompT scorescaletext;
PrompT segmentdnascaletext;

#define MAXFONTS 6
int nfonts=5;
int av_font[MAXFONTS]={8,10,12,14,18};
int font_size=1;

int ncolors=0;
int ncolor_pars=0;
color color_lut[MAXCOLORS+1];
char def_protpar_file[]="colprot.par";
char def_dnapar_file[]="coldna.par";
char *explicit_par_file = NULL;
char *par_file = NULL;
int    inverted = TRUE;
int usebw=FALSE,usedefcolors=TRUE,useusercolors=FALSE;

char ps_par_file[FILENAMELEN]="colprint.par";
int pagesize=A4;
int orientation=LANDSCAPE;
Boolean ps_header=TRUE;
Boolean ps_ruler=TRUE;
Boolean resize=TRUE;
int first_printres=0,last_printres=0,blocklen;
Boolean ps_curve=TRUE;
Boolean ps_resno=TRUE;
PoinT display_pos;
int namewidth,seqwidth; /* fixed widths of sequence display areas */

Boolean         realign_endgappenalties=TRUE;
Boolean         align_endgappenalties=FALSE;

char helptext[MAXHELPLENGTH];
/* main subroutine called from clustalx.c, initialises windows and enters a
   forever loop monitoring user input */

void x_menu(void)
{
	int i,n;
	char font[30];
	char tstr[30];
	int height;
	PrompT   fsize;
	RecT wr,r,r1;
 

/*  make the pulldown menu bar  */

#ifdef WIN_MAC
	MenU   m;

        m=AppleMenu (NULL);
        DeskAccGroup (m);
	make_menu_headers(NULL);
#endif
#ifndef UNIX
	ProcessUpdatesFirst(FALSE);
#endif

	sprintf(tstr,"Clustal%s",revision_level);
/*#ifdef WIN_MSWIN
	mainw = FixedWindow (-50,-33,-10,-10,tstr,QuitWinW);
#else*/
	mainw = DocumentWindow (-50,-33,-10,-10,tstr,QuitWinW,ResizeWindowProc);
/*#endif*/	SetGroupSpacing(mainw,0,10);
	SetGroupSpacing(mainw,0,10);
 
	x_menus=TRUE;

#ifndef WIN_MAC
	make_menu_headers(mainw);
#endif
/* decide if we're starting in profile or sequence mode */
	if (!profile1_empty) aln_mode=PROFILEM;
	else aln_mode=MULTIPLEM;

	make_file_menu();
	make_edit_menu();
	make_align_menu();
	make_tree_menu();
	make_color_menu();
	make_score_menu();
	make_help_menu();

/*  add a button to switch between multiple and profile alignment modes */

	modetext=PopupList(mainw,TRUE,set_aln_mode);
	PopupItem(modetext,"Multiple Alignment Mode");
	PopupItem(modetext,"Profile Alignment Mode");
	if(aln_mode==MULTIPLEM)
		SetValue(modetext,1);
	else
		SetValue(modetext,2);

	sprintf(font, "%s,%d,%c", "courier", av_font[font_size], 'm');
        datafont=ParseFont(font);

	sprintf(font, "%s,%d,%c", "courier", 10, 'm');
        helpfont=ParseFont(font);

	Advance(mainw);
	shift(mainw,20,0);

/*  add a button to select font size */
	fsize=StaticPrompt(mainw,"Font Size:",0,dialogTextHeight,systemFont,'r');
	Advance(mainw);
	flisttext=PopupList(mainw,TRUE,set_font_size);
	for(i=0;i<nfonts;i++)
	{
		sprintf(tstr,"%d",av_font[i]);
		PopupItem(flisttext,tstr);
	}
	SetValue(flisttext,font_size+1);

	Advance(mainw);
	shift(mainw,20,0);

/*  add a button to switch profile scrolling modes */
        pscrolltext=CheckBox(mainw,"Lock Scroll",set_pscroll_mode);
	if(fixed_prf_scroll) SetStatus(pscrolltext,TRUE);
	Break(mainw);


	selected_seqs.first=selected_seqs.last=-1;
        selected_res.first=selected_res.last=-1;


/*  initialise the multiple alignment display area */

	SelectFont(datafont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
 
	GetNextPosition(mainw,&display_pos);

/* calculate initial pixel width and height of displays */
	namewidth=(DNAMES+DNUMBER+1)*stdCharWidth;
	seqwidth=(DCOLS+2*MARGIN)*stdCharWidth+2;
	n=screenRect.right-screenRect.left;
	if(seqwidth+namewidth>n) seqwidth=n-namewidth;

	height=(max_mlines+mfooter+MARGIN)*stdLineHeight+2+SCOREHEIGHT;
	n=screenRect.bottom-screenRect.top;
	if(height>n) height=n;

	seq_display=make_scroll_area(mainw,0,namewidth,seqwidth,height,1,nseqs,&seq_panel);
	position_scrollbars(seq_panel);

/*  initialise the profile alignment display area */
 
	SetNextPosition(mainw,display_pos);
	height=(max_plines+MARGIN)*stdLineHeight+2+SCOREHEIGHT;
	if(height>n) height=n;
	prf1_display=make_scroll_area(mainw,1,namewidth,seqwidth,height,1,profile1_nseqs,&prf_panel[0]);
	position_scrollbars(prf_panel[0]);

	prf2_display=make_scroll_area(mainw,2,namewidth,seqwidth,height,profile1_nseqs+1,nseqs-profile1_nseqs,&prf_panel[1]);
	position_scrollbars(prf_panel[1]);

/*  add the message line */
	Break(mainw);
	Advance(mainw);
	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	message = StaticPrompt(mainw, "",500, 0,systemFont,'l');

/* save some pixel sizes for future resizing events */
	if(aln_mode==PROFILEM)
	{
		Hide(seq_display);
		profile_no=1;
		Show(prf1_display);
		Show(prf2_display);
		Show(pscrolltext);
		active_panel=prf_panel[0];
		Select(prf1_display);
		load_aln(prf_panel[0],0,profile1_nseqs-1,TRUE);
		load_aln(prf_panel[1],profile1_nseqs,nseqs-1,TRUE);

		Show(mainw);
		ObjectRect(mainw,&wr);
		ObjectRect(prf_panel[0].names,&r);
		ObjectRect(prf_panel[1].names,&r1);
		boffset=wr.bottom-wr.top-r1.bottom;
		loffset=r.left;
		toffset=r.top;
		ObjectRect(prf_panel[0].seqs,&r);
		roffset=wr.right-wr.left-r.right;
	}
	else
	{
		Hide(prf1_display);
		Hide(prf2_display);
		Hide(pscrolltext);
		profile_no=0;
		Show(seq_display);
		active_panel=seq_panel;

		Select(seq_display);
		load_aln(seq_panel,0,nseqs-1,TRUE);

		Show(mainw);
		ObjectRect(mainw,&wr);
		ObjectRect(seq_panel.names,&r);
		boffset=wr.bottom-wr.top-r.bottom;
		loffset=r.left;
		toffset=r.top;
		ObjectRect(seq_panel.seqs,&r);
		roffset=wr.right-wr.left-r.right;
	}
	ObjectRect(prf_panel[0].names,&r);
	ObjectRect(prf_panel[1].names,&r1);
	poffset=r1.top-r.bottom;

/* initialise some variables before we display the window */
        if(orientation==LANDSCAPE)
        {
                if(pagesize==A4) blocklen=150;
                else if (pagesize==A3) blocklen=250;
		else blocklen=150;
        }
        else
        {
                if(pagesize==A4) blocklen=80;
                else if (pagesize==A3) blocklen=150;
		else blocklen=150;
        }

/* ok - Go! */
	window_displayed=TRUE;
	ProcessEvents();

}


static void RemoveWin(WindoW w)
{
	Remove(w);
}


static void QuitWinW(WindoW w)
{
	if(aln_mode == MULTIPLEM)
	{
		if(seq_panel.modified)
			if (Message(MSG_YN,"Alignment has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
	}
	else if(aln_mode == PROFILEM)
	{
		if(prf_panel[0].modified)
			if (Message(MSG_YN,"Profile 1 has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
		if(prf_panel[1].modified)
			if (Message(MSG_YN,"Profile 2 has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
	}
	QuitProgram ();
}

static void SearchStrWin (IteM item)
{
	int i;
	Boolean sel=FALSE;
	GrouP findgr;
	ButtoN find_can,find_ok;
	PopuP ps,or;
	char path[FILENAMELEN];
	char str[FILENAMELEN];
	panel_data data;

	GetPanelExtra(active_panel.names,&data);
	if (data.nseqs==0)
	{
		Message(MSG_OK,"No file loaded.");
		return;
	}
	for (i=0;i<data.nseqs;i++)
		if(data.selected[i]==TRUE)
		{
			sel=TRUE;
			break;
		}
	if(sel==FALSE)
	{
		Message(MSG_OK,"Select sequences by clicking on the names.");
		return;
	}

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	findw=FixedWindow(-50, -33, -10, -10, "SEARCH IN SELECTED SEQUENCES",RemoveWin);
        stdLineHeight=18;
        SelectFont(programFont);
	findtext=DialogText(findw, "", 35, NULL);
	Break(findw);
	find_ok=PushButton(findw, "SEARCH FROM START", SearchStr);
	Break(findw);
	find_ok=PushButton(findw, "SEARCH AGAIN", SearchStringAgain);
	Break(findw);
	find_can=PushButton(findw, "CLOSE", CancelWin);

	Show(findw);
}

static void SavePSSeqWin (IteM item)
{
	if (empty)
	{
		error("No file loaded");
		return;
	}
	save_ps_window(0,"WRITE SEQUENCES TO:",SavePSSeqFile);
}

static void SavePSPrf1Win (IteM item)
{
	if (profile1_empty)
	{
		error("No file loaded");
		return;
	}
	save_ps_window(1,"WRITE PROFILE 1 TO:",SavePSPrf1File);
}

static void SavePSPrf2Win (IteM item)
{
	if (profile2_empty)
	{
		error("No file loaded");
		return;
	}
	save_ps_window(2,"WRITE PROFILE 2 TO:",SavePSPrf2File);
}

static void save_ps_window(int prf_no,char *prompt,void save_proc(ButtoN but))
{
	GrouP savegr;
	ButtoN save_can,save_ok;
	PopuP ps,or;
	char path[FILENAMELEN];
	char str[FILENAMELEN];
	panel_data data;

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	savepsw=FixedWindow(-50, -33, -10, -10, "WRITE POSTSCRIPT FILE",RemoveWin);
	make_prompt(savepsw, prompt);
        stdLineHeight=18;
        SelectFont(programFont);
	savepstext=DialogText(savepsw, "", 35, NULL);
	Break(savepsw);
	make_prompt(savepsw, "PS Colors File :");
	pspartext=DialogText(savepsw, ps_par_file, 35, NULL);
	Break(savepsw);
       	make_prompt(savepsw, "Page Size");
	Advance(savepsw);
	ps=PopupList(savepsw,TRUE,set_pagesize);
	PopupItem(ps,"A4");
	PopupItem(ps,"A3");
	PopupItem(ps,"US Letter");
	if (pagesize == A4)
		SetValue(ps,1);
	else if (pagesize == A3)
		SetValue(ps,2);
	else if (pagesize == USLETTER)
		SetValue(ps,3);
	Break(savepsw);
       	make_prompt(savepsw, "Orientation");
	Advance(savepsw);
	or=PopupList(savepsw,TRUE,set_orientation);
	PopupItem(or,"LANDSCAPE");
	PopupItem(or,"PORTRAIT");
	if (orientation == LANDSCAPE)
		SetValue(or,1);
	else if (orientation == PORTRAIT)
		SetValue(or,2);
	Break(savepsw);
	make_toggle(savepsw,"Print Header :","YES","NO",&ps_header,set_header);
	Advance(savepsw);
	make_toggle(savepsw,"Print Quality Curve :","YES","NO",&ps_curve,set_curve);
	Break(savepsw);
	make_toggle(savepsw,"Print Ruler :","YES","NO",&ps_ruler,set_ruler);
	Advance(savepsw);
	make_toggle(savepsw,"Print Residue Numbers :","YES","NO",&ps_resno,set_resno);
	Break(savepsw);
	make_toggle(savepsw,"Resize to fit page:","YES","NO",&resize,set_resize);
	Break(savepsw);
	first_printres=1;
	if (prf_no==0)
		GetPanelExtra(seq_panel.seqs,&data);
	else if (prf_no==1)
		GetPanelExtra(prf_panel[0].seqs,&data);
	else
		GetPanelExtra(prf_panel[1].seqs,&data);
	last_printres=data.ncols;
        make_prompt(savepsw, "Print from position :");
	Advance(savepsw);
	sprintf(str,"%5d",first_printres);
        DialogText(savepsw, str, 5,set_fpres);
	Advance(savepsw);
        make_prompt(savepsw, "to :");
	Advance(savepsw);
	sprintf(str,"%5d",last_printres);
        DialogText(savepsw, str, 5,set_lpres);
	Break(savepsw);
        make_prompt(savepsw, "Use block length :");
	Advance(savepsw);
	sprintf(str,"%5d",blocklen);
        blocklentext=DialogText(savepsw, str, 5,set_blocklen);
	Break(savepsw);
	savegr=HiddenGroup(savepsw, 2, 0, NULL);
	shift(savegr, 60, 20);
	save_ok=PushButton(savegr, "  OK  ", save_proc);
	shift(savegr, 20,0);
	save_can=PushButton(savegr, "CLOSE", CancelWin);

	if(prf_no==0)
		get_path(seqname,path);
	else if(prf_no==1)
		get_path(profile1_name,path);
	else if(prf_no==2)
		get_path(profile2_name,path);
	strcat(path,"ps");
	SetTitle(savepstext, path);
	Show(savepsw);
}

static void SaveScoresWin (IteM item)
{
	int i;
	Boolean sel=FALSE;
	GrouP scoregr;
	ButtoN score_can,score_ok;
	PopuP ps,or;
	char path[FILENAMELEN];
	char str[FILENAMELEN];
	panel_data data;


	if (empty)
	{
		error("No file loaded");
		return;
	}

	GetPanelExtra(active_panel.names,&data);
	for (i=0;i<data.nseqs;i++)
		if(data.selected[i]==TRUE)
		{
			sel=TRUE;
			break;
		}
	if(sel==FALSE)
	{
		Message(MSG_OK,"Select sequences to be written by clicking on the names.");
		return;
	}

        get_path(seqname,path);
	strcat(path,"qscores");

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	savescoresw=FixedWindow(-50, -33, -10, -10, "SAVE QUALITY SCORES",RemoveWin);
        stdLineHeight=18;
        SelectFont(programFont);
	make_prompt(savescoresw, "SAVE QUALITY SCORES TO:");
        stdLineHeight=18;
        SelectFont(programFont);
	Break(savescoresw);
	savescorestext=DialogText(savescoresw, "", 35, NULL);
	Break(savescoresw);
	scoregr=HiddenGroup(savescoresw, 2, 0, NULL);
	shift(scoregr, 60, 20);
	score_ok=PushButton(scoregr, "  OK  ", SaveScores);
	shift(scoregr, 20,0);
	score_can=PushButton(scoregr, "CANCEL", CancelWin);

	SetTitle(savescorestext, path);
	Show(savescoresw);

	Advance(savescoresw);
	Show(savescoresw);
}

static void SaveScores(ButtoN but)
{
	char c;
	int i,j,val;
	int length=0;
	FILE *outfile;
	panel_data name_data,seq_data;
	Boolean gap;

	GetPanelExtra(active_panel.names,&name_data);
	GetPanelExtra(active_panel.seqs,&seq_data);

	GetTitle(savescorestext, filename, FILENAMELEN);
	stripspace(filename);

	outfile=open_explicit_file(filename); 

/* get the maximum length of the selected sequences */
        for (i=1;i<=nseqs;i++)
           if (name_data.selected[i-1]==TRUE && length < seqlen_array[i]) length = seqlen_array[i];

	for(j=1;j<=length;j++)
	{
/* first check for a column of gaps */
		gap=TRUE;
        	for (i=1;i<=nseqs;i++)
           		if (name_data.selected[i-1]==TRUE)
			{
                                val = seq_array[i][j];
                                if(j<=seqlen_array[i] && (val != gap_pos1) && (val != gap_pos2))
				{
					gap=FALSE;
					break;
				}
			}
		if(gap==FALSE)
		{
        		for (i=1;i<=nseqs;i++)
			{
           			if (name_data.selected[i-1]==TRUE)
				{
                                	val = seq_array[i][j];
                                	if(j>seqlen_array[i] || (val == gap_pos1) || (val == gap_pos2))
                                        	c = '-';
                                	else {
                                        	c = amino_acid_codes[val];
                                	}
	 
					fprintf(outfile,"%c ",c);
				}
			}
			fprintf(outfile,"\t%3d\n",seq_data.colscore[j-1]);
		}

	}
	fclose(outfile);

        if (Visible(savescoresw))
        {
                Remove(savescoresw);
                savescoresw=NULL;
        }



	info("File %s saved",filename);
}

static void SaveSeqFileWin (IteM item)
{
	if (empty)
	{
		error("No file loaded");
		return;
	}
	save_aln_window(0,"SAVE SEQUENCES","SAVE SEQUENCES AS:",SaveSeqFile);
}
static void SavePrf1FileWin (IteM item)
{
	if (profile1_empty)
	{
		error("No file loaded");
		return;
	}
	save_aln_window(1,"SAVE PROFILE","SAVE PROFILE 1 AS:",SavePrf1File);
}
static void SavePrf2FileWin (IteM item)
{
	if (profile2_empty)
	{
		error("No file loaded");
		return;
	}
	save_aln_window(2,"SAVE PROFILE","SAVE PROFILE 2 AS:",SavePrf2File);
}

static void save_aln_window(int prf_no,char *title,char *prompt,void save_proc(ButtoN but))
{
	GrouP savegr;
	ButtoN save_ok, save_can;
	GrouP maing;
        GrouP format_list;
	ButtoN formatb[6];
	PopuP case_toggle,snos_toggle;
	char path[FILENAMELEN];
	char str[FILENAMELEN];

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();

	savealnw=FixedWindow(-50, -33, -10, -10, title,RemoveWin);

        format_list=NormalGroup(savealnw,3,0,"Format",systemFont,set_format);
        formatb[0]=RadioButton(format_list,"CLUSTAL");
        formatb[1]=RadioButton(format_list,"NBRF/PIR");
        formatb[2]=RadioButton(format_list,"GCG/MSF");
        formatb[3]=RadioButton(format_list,"PHYLIP");
        formatb[4]=RadioButton(format_list,"GDE");
        formatb[5]=RadioButton(format_list,"NEXUS");
	if(prf_no==0)
        	get_path(seqname,path);
	else if(prf_no==1)
        	get_path(profile1_name,path);
	else if(prf_no==2)
        	get_path(profile2_name,path);
	if (save_format==CLUSTAL)
	{
        	SetValue(format_list,1);
		strcat(path,"aln");
	}
	else if (save_format==PIR)
	{
        	SetValue(format_list,2);
		strcat(path,"pir");
	}
	else if (save_format==MSF)
	{
        	SetValue(format_list,3);
		strcat(path,"msf");
	}
	else if (save_format==PHYLIP)
	{
        	SetValue(format_list,4);
		strcat(path,"phy");
	}
	else if (save_format==GDE)
	{
        	SetValue(format_list,5);
		strcat(path,"gde");
	}
 	else if (save_format==NEXUS)
	{
        	SetValue(format_list,6);
		strcat(path,"nxs");
	}

 
	maing=HiddenGroup(savealnw,0,0,NULL);
	SetGroupSpacing(maing,0,10);

	case_toggle=make_toggle(maing,"GDE output case :","Lower","Upper",&lowercase,set_case);
	Break(maing);
	snos_toggle=make_toggle(maing,"CLUSTALW sequence numbers :","ON","OFF",&cl_seq_numbers,set_snos);

	Break(maing);
        make_prompt(maing, "Save from residue :");
	Advance(maing);
	sprintf(str,"%5d",firstres);
        DialogText(maing, str, 5,set_fres);
	Advance(maing);
        make_prompt(maing, "to :");
	Advance(maing);
	sprintf(str,"%5d",lastres);
        DialogText(maing, str, 5,set_lres);

	Break(maing);
	shift(savealnw, 0, 20);
	make_prompt(savealnw, prompt);
        stdLineHeight=18;
        SelectFont(programFont);
	Break(savealnw);
	savealntext=DialogText(savealnw, "", 35, NULL);
	Break(savealnw);
	savegr=HiddenGroup(savealnw, 2, 0, NULL);
	shift(savegr, 60, 20);
	save_ok=PushButton(savegr, "  OK  ", save_proc);
	shift(savegr, 20,0);
	save_can=PushButton(savegr, "CANCEL", CancelWin);

	SetTitle(savealntext, path);
	Show(savealnw);
}

static void read_file_window(char *title,char *prompt,char *filename,void read_proc(ButtoN but))
{
	GrouP readgr;
	ButtoN read_ok, read_can;
	GrouP maing;

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	readfilew=FixedWindow(-50, -33, -10, -10, title,RemoveWin);

	maing=HiddenGroup(readfilew,2,0,NULL);
	SetGroupSpacing(maing,0,10);

	shift(readfilew, 0, 20);
	make_prompt(readfilew, prompt);
        stdLineHeight=18;
        SelectFont(programFont);
	Break(readfilew);
	readfiletext=DialogText(readfilew, "", 35, NULL);
	if (filename != NULL) SetTitle(readfiletext, filename);
	Break(readfilew);
	readgr=HiddenGroup(readfilew, 2, 0, NULL);
	shift(readgr, 60, 20);
	read_ok=PushButton(readgr, "  OK  ", read_proc);
	shift(readgr, 20,0);
	read_can=PushButton(readgr, "CANCEL", CancelWin);

	Show(readfilew);
}

static void CancelWin (ButtoN but)
{
	Remove(ParentWindow(but));
}

static void SearchStr(ButtoN but)
{

/* reset the current position */

	find_pos.seq=0;
	find_pos.res=-1;

/* find the next occurrence of the string */
	SearchStringAgain(but);


}

static void SearchStringAgain(ButtoN but)
{
	int i,j,ix,length;
	int seq,res,start_res;
	Boolean in_string,found;
	panel_data ndata,sdata;

	GetTitle(findtext, filename, FILENAMELEN);
	stripspace(filename);

	strncpy(find_string,filename,MAXFINDSTR);
	length=strlen(find_string);
	if(length==0) return;
	for(i=0;i<length;i++)
		find_string[i]=toupper(find_string[i]);

        GetPanelExtra(active_panel.names,&ndata);
        GetPanelExtra(active_panel.seqs,&sdata);

	in_string=FALSE;
	found=FALSE;
	start_res=0;
	ix=0;
	seq=find_pos.seq;
	res=find_pos.res+1;
        while (seq<ndata.nseqs)
	{
                if(ndata.selected[seq]==TRUE)
                {
        		while (res<sdata.ncols)
			{
				if(sdata.lines[seq][res]==find_string[ix])
				{
					if(in_string==FALSE) 
						start_res=res;
					ix++;
					in_string=TRUE;
				}
				else if(in_string==TRUE)
				{
					res=start_res;
					ix=0;
					in_string=FALSE;
				}
				if(ix==length)
				{
					find_pos.seq=seq;
					find_pos.res=start_res;
					found=TRUE;
					break;
				}
				res++;
				while(res<sdata.ncols && sdata.lines[seq][res]=='-')
					res++;
			}
                }
		if(found) break;
		seq++;
		res=0;
	}


	if(found==FALSE)
		info("String %s not found",find_string);
	else
	{
		info("String %s in sequence %s, column %d",find_string,names[find_pos.seq+1],find_pos.res+1);
	}
}

static void SavePSSeqFile(ButtoN but)
{
	char *ps_file;

	GetTitle(savepstext, filename, FILENAMELEN);
	stripspace(filename);

        ps_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	strcpy(ps_file,filename); 

	GetTitle(pspartext, filename, FILENAMELEN);
	stripspace(filename);

	strcpy(ps_par_file,filename); 

	write_ps_file(seq_panel,ps_file,ps_par_file,pagesize,orientation,
		ps_header,ps_ruler,ps_resno,
		resize,first_printres,last_printres,blocklen,ps_curve);

	info("Postscript file %s written",ps_file);
	ckfree(ps_file);

}

static void SavePSPrf1File(ButtoN but)
{
	char *ps_file;
	char *ps_par_file;

	GetTitle(savepstext, filename, FILENAMELEN);
	stripspace(filename);

        ps_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	strcpy(ps_file,filename); 

	GetTitle(pspartext, filename, FILENAMELEN);
	stripspace(filename);

        ps_par_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	strcpy(ps_par_file,filename); 

	write_ps_file(prf_panel[0],ps_file,ps_par_file,pagesize,orientation,
		ps_header,ps_ruler,ps_resno,
		resize,first_printres,last_printres,blocklen,ps_curve);

	info("Postscript file %s written",ps_file);
	ckfree(ps_file);

}

static void SavePSPrf2File(ButtoN but)
{
	char *ps_file;
	char *ps_par_file;

	GetTitle(savepstext, filename, FILENAMELEN);
	stripspace(filename);

        ps_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	strcpy(ps_file,filename); 

	GetTitle(pspartext, filename, FILENAMELEN);
	stripspace(filename);

        ps_par_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	strcpy(ps_par_file,filename); 

	write_ps_file(prf_panel[1],ps_file,ps_par_file,pagesize,orientation,
		ps_header,ps_ruler,ps_resno,
		resize,first_printres,last_printres,blocklen,ps_curve);

	info("Postscript file %s written",ps_file);
	ckfree(ps_file);

}

static void SaveSeqFile(ButtoN but)
{
	write_file(1,nseqs,firstres,lastres);
	seq_panel.modified=FALSE;
	info("File %s saved",filename);
}

static void SavePrf1File(ButtoN but)
{
	write_file(1,profile1_nseqs,firstres,lastres);
	prf_panel[0].modified=FALSE;
	info("File %s saved",filename);
}

static void SavePrf2File(ButtoN but)
{
	write_file(profile1_nseqs+1,nseqs,firstres,lastres);
	prf_panel[1].modified=FALSE;
	info("File %s saved",filename);
}

/* this is equivalent to open_alignment_output(), but uses the window
interface to input file names */

static Boolean open_aln_files(void)
{
	char path[FILENAMELEN];

	if(!output_clustal && !output_nbrf && !output_gcg &&
		 !output_phylip && !output_gde && !output_nexus) {
                error("You must select an alignment output format");
                return FALSE;
        }

	if(output_clustal) {
		GetTitle(cl_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((clustal_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}
	if(output_nbrf) {
		GetTitle(pir_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((nbrf_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}
	if(output_gcg) {
		GetTitle(msf_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((gcg_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}
	if(output_phylip) {
		GetTitle(phylip_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((phylip_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}
	if(output_gde) {
		GetTitle(gde_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((gde_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}
	if(output_nexus) {
		GetTitle(nexus_outtext,filename,FILENAMELEN);
		stripspace(filename);
		if((nexus_outfile = open_explicit_file(
			filename))==NULL) return FALSE;
	}

	if(save_log)
	{
        	get_path(seqname,path);
        	strcpy(save_log_filename,path);
        	strcat(save_log_filename,"log");
		if ((save_log_fd=fopen(save_log_filename,"a"))==NULL)
			error("Cannot open log file %s",save_log_filename);
	}

	return TRUE;
}

static void write_file(int fseq, int lseq, int fres, int lres)
{
	int i,length=0;
	FILE *outfile;

	GetTitle(savealntext, filename, FILENAMELEN);
	stripspace(filename);

	outfile=open_explicit_file(filename); 

        for (i=fseq;i<=lseq;i++)
           if (length < seqlen_array[i]) length = seqlen_array[i];

	if(fres<1) fres=1;
	if(lres<1) lres=length;
	length=lres-fres+1;
 
        if(save_format==CLUSTAL) {
                clustal_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("CLUSTAL format file created  [%s]",filename);
        }
        else if(save_format==PIR)  {
                nbrf_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("NBRF/PIR format file created  [%s]",filename);
        }
        else if(save_format==MSF)  {
                gcg_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("GCG/MSF format file created  [%s]",filename);
        }
        else if(save_format==PHYLIP)  {
                phylip_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("PHYLIP format file created  [%s]",filename);
        }
        else if(save_format==GDE)  {
                gde_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("GDE format file created  [%s]",filename);
        }
        else if(save_format==NEXUS)  {
                nexus_out(outfile, fres, length, fseq, lseq);
                fclose(outfile);
                info("NEXUS format file created  [%s]",filename);
        }


	if (Visible(savealnw))
	{
		Remove(savealnw);
		savealnw=NULL;
	}


}

static void SaveTreeWin (IteM item)
{
	GrouP savegr;
	ButtoN save_ok, save_can;
	char path[FILENAMELEN];

	if (empty)
	{
		error("No file loaded");
		return;
	}
        if (nseqs < 2)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	savetreew=FixedWindow(-50, -33, -10, -10, "CREATE TREE",RemoveWin);
	shift(savetreew, 0, 20);
	make_prompt(savetreew, "SAVE TREE AS :");
	Advance(savetreew);
	shift(savetreew, 0, -10);
	stdLineHeight=18;
	SelectFont(programFont);
	savetreetext=DialogText(savetreew, "", 35, NULL);
	SelectFont(systemFont);
	stdLineHeight=15;
	Break(savetreew);
	savegr=HiddenGroup(savetreew, 2, 0, NULL);
	shift(savegr, 140, 20);
	save_ok=PushButton(savegr, "  OK  ", CreateAlignTree);
	shift(savegr, 20, 0);
	save_can=PushButton(savegr, "CANCEL", CancelWin);

	get_path(seqname,path);
	strcat(path,"dnd");
  
	SetTitle(savetreetext, path);
	Show(savetreew);
}

static void DrawTreeWin (IteM item)
{
	GrouP drawgr;
	GrouP output_list;
	ButtoN draw_ok, draw_can;
	char path[FILENAMELEN];
	char name[FILENAMELEN];

	if (empty)
	{
		error("No file loaded");
		return;
	}
        if (nseqs < 2)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }

	get_path(seqname,path);

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	drawtreew=FixedWindow(-50, -33, -10, -10, "DRAW TREE",RemoveWin);
	output_list=HiddenGroup(drawtreew, 2, 0, NULL);
	if (output_tree_clustal)
	{
		make_prompt(output_list, "SAVE CLUSTAL TREE AS :");
		drawnjtreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"nj");
		SetTitle(drawnjtreetext, name);
		Break(output_list);
	}
	if (output_tree_phylip)
	{
		make_prompt(output_list, "SAVE PHYLIP TREE AS :");
		drawphtreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"ph");
		SetTitle(drawphtreetext, name);
		Break(output_list);
	}
	if (output_tree_distances)
	{
		make_prompt(output_list, "SAVE DISTANCE MATRIX AS :");
		drawdsttreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"dst");
		SetTitle(drawdsttreetext, name);
		Break(output_list);
	}
	if (output_tree_nexus)
	{
		make_prompt(output_list, "SAVE NEXUS TREE AS :");
		drawnxstreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"tre");
		SetTitle(drawnxstreetext, name);
		Break(output_list);
	}
	SelectFont(systemFont);
	stdLineHeight=15;
	Break(drawtreew);
	drawgr=HiddenGroup(drawtreew, 2, 0, NULL);
	shift(drawgr, 140, 20);
	draw_ok=PushButton(drawgr, "  OK  ", DrawTree);
	shift(drawgr, 20, 0);
	draw_can=PushButton(drawgr, "CANCEL", CancelWin);


	Show(drawtreew);
}

static void BootstrapTreeWin (IteM item)
{
	GrouP bootgr;
	ButtoN boot_ok, boot_can;
	TexT seed,ntrials;
	char name[FILENAMELEN];
	char path[FILENAMELEN];
	char str[FILENAMELEN];
	GrouP output_list;

	if (empty)
	{
		error("No file loaded");
		return;
	}
        if (nseqs < 2)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }

	get_path(seqname,path);

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	boottreew=FixedWindow(-50, -33, -10, -10, "BOOTSTRAP TREE",RemoveWin);
        make_prompt(boottreew, "Random number generator seed [1-1000] :");
	Advance(boottreew);
	sprintf(str,"%4d",boot_ran_seed);
        seed=DialogText(boottreew, str, 4,set_ran_seed);
	Break(boottreew);
        make_prompt(boottreew, "Number of bootstrap trials [1-10000] :");
	Advance(boottreew);
	sprintf(str,"%5d",boot_ntrials);
        ntrials=DialogText(boottreew, str, 5,set_ntrials);
	Break(boottreew);

	output_list=HiddenGroup(boottreew, 2, 0, NULL);
	if (output_tree_clustal)
	{
		make_prompt(output_list, "SAVE CLUSTAL TREE AS :");
		bootnjtreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"njb");
		SetTitle(bootnjtreetext, name);
		Break(output_list);
	}
	if (output_tree_phylip)
	{
		make_prompt(output_list, "SAVE PHYLIP TREE AS :");
		bootphtreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"phb");
		SetTitle(bootphtreetext, name);
		Break(output_list);
	}
	if (output_tree_nexus)
	{
		make_prompt(output_list, "SAVE NEXUS TREE AS :");
		bootnxstreetext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
        	strcat(name,"treb");
		SetTitle(bootnxstreetext, name);
		Break(output_list);
	}
	SelectFont(systemFont);
	stdLineHeight=15;
	Break(boottreew);
	bootgr=HiddenGroup(boottreew, 2, 0, NULL);
	shift(bootgr, 140, 20);
	boot_ok=PushButton(bootgr, "  OK  ", BootstrapTree);
	shift(bootgr, 20, 0);
	boot_can=PushButton(bootgr, "CANCEL", CancelWin);


	Show(boottreew);
}

static void CreateAlignTree(ButtoN but)
{
	char path[FILENAMELEN];
	char phylip_name[FILENAMELEN];

	GetTitle(savetreetext, filename, FILENAMELEN);
	strcpy(phylip_name,filename);
	stripspace(filename);

	info("Doing pairwise alignments...");
	if(save_log)
	{
        	get_path(seqname,path);
        	strcpy(save_log_filename,path);
        	strcat(save_log_filename,"log");
		if ((save_log_fd=fopen(save_log_filename,"a"))==NULL)
			error("Cannot open log file %s",save_log_filename);
	}

	WatchCursor();
	if (Visible(savetreew))
	{
		Remove(savetreew);
		savetreew=NULL;
	}
	make_tree(phylip_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
	ArrowCursor();
	info("Tree %s created",filename);
}

static void DrawTree(ButtoN but)
{
	char path[FILENAMELEN];
	char phylip_name[FILENAMELEN];
	char clustal_name[FILENAMELEN];
	char dist_name[FILENAMELEN];
	char nexus_name[FILENAMELEN];

	if(output_tree_clustal)
	{
		GetTitle(drawnjtreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(clustal_name,filename);
	}
	if(output_tree_phylip)
	{
		GetTitle(drawphtreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(phylip_name,filename);
	}
	if(output_tree_distances)
	{
		GetTitle(drawdsttreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(dist_name,filename);
	}
	if(output_tree_nexus)
	{
		GetTitle(drawnxstreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(nexus_name,filename);
	}

	info("Calculating tree...");
	WatchCursor();
	if(save_log)
	{
        	get_path(seqname,path);
        	strcpy(save_log_filename,path);
        	strcat(save_log_filename,"log");
		if ((save_log_fd=fopen(save_log_filename,"a"))==NULL)
			error("Cannot open log file %s",save_log_filename);
	}
	if (Visible(drawtreew))
	{
		Remove(drawtreew);
		drawtreew=NULL;
	}
	phylogenetic_tree(phylip_name,clustal_name,dist_name,nexus_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
	ArrowCursor();
	info("Tree %s created",filename);
}

static void BootstrapTree(ButtoN but)
{
	char phylip_name[FILENAMELEN];
	char clustal_name[FILENAMELEN];
	char nexus_name[FILENAMELEN];
	char path[FILENAMELEN];

	if(output_tree_clustal)
	{
		GetTitle(bootnjtreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(clustal_name,filename);
	}
	if(output_tree_phylip)
	{
		GetTitle(bootphtreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(phylip_name,filename);
	}
	if(output_tree_nexus)
	{
		GetTitle(bootnxstreetext, filename, FILENAMELEN);
		stripspace(filename);
		strcpy(nexus_name,filename);
	}

	info("Bootstrapping tree...");

	WatchCursor();
	if(save_log)
	{
        	get_path(seqname,path);
        	strcpy(save_log_filename,path);
        	strcat(save_log_filename,"log");
		if ((save_log_fd=fopen(save_log_filename,"a"))==NULL)
			warning("Cannot open log file %s",save_log_filename);
	}
	if (Visible(boottreew))
	{
		Remove(boottreew);
		boottreew=NULL;
	}
	bootstrap_tree(phylip_name,clustal_name,nexus_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
	info("Bootstrap tree %s created",filename);
	ArrowCursor();
}

static void OpenSeqFile (IteM item)
{
	int n;
	panel_data data;

	if (nseqs>0)
	{
		if (Message(MSG_YN,"Replace existing sequences ?")==ANS_NO)
			return;
	}

	if (!GetInputFileName (filename,FILENAMELEN,"","")) return;
 
	strcpy(seqname,filename);
	GetPanelExtra(seq_panel.names,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(seq_panel.names,&data);
	GetPanelExtra(seq_panel.seqs,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(seq_panel.seqs,&data);
	
	n=seq_input(FALSE);
	if (n<=0)
	{
		info("File %s not loaded.",seqname);
		return;
	}

	load_aln(seq_panel,0,nseqs-1,TRUE);

	ncutseqs=0;

	info("File %s loaded.",seqname);
}

static void AppendSeqFile (IteM item)
{
	int n;
	panel_data data;

	if (!GetInputFileName (filename,FILENAMELEN,"","")) return;
 
	strcpy(seqname,filename);
	GetPanelExtra(seq_panel.names,&data);
	data.nseqs=0;
	SetPanelExtra(seq_panel.names,&data);
	n=seq_input(TRUE);
	if (n<=0)
	{
		info("File %s not loaded.",seqname);
		return;
	}

	load_aln(seq_panel,0,nseqs-1,FALSE);

	info("File %s appended.",seqname);
}

static void OpenPrf1File (IteM item)
{
	int i,j,n,tmpn=0,tmpfs;
	sint *tmplen_array;
	sint *tmpindex;
	char **tmp_array;
	char **tmpnames;
	char **tmptitles;
	panel_data data;

	if (profile1_nseqs>0)
	{
		if (Message(MSG_YN,"Replace existing sequences ?")==ANS_NO)
			return;
	}

	if (!GetInputFileName (filename,FILENAMELEN,"","")) return;
 
	if(!profile2_empty)
	{
		tmpn=nseqs-profile1_nseqs;
		tmpfs=profile1_nseqs;
		tmpnames=(char **)ckalloc((tmpn+1)*sizeof(char *));
		tmptitles=(char **)ckalloc((tmpn+1)*sizeof(char *));
		tmplen_array=(sint *)ckalloc((tmpn+1)*sizeof(sint));
		tmpindex=(sint *)ckalloc((tmpn+1)*sizeof(sint));
		tmp_array=(char **)ckalloc((tmpn+1)*sizeof(char *));
		for(i=profile1_nseqs+1;i<=nseqs;i++)
		{
			tmpnames[i-profile1_nseqs-1]=(char *)ckalloc((MAXNAMES+2)*sizeof(char));
			tmptitles[i-profile1_nseqs-1]=(char *)ckalloc((MAXTITLES+2)*sizeof(char));
			strcpy(tmpnames[i-profile1_nseqs-1],names[i]);

			strcpy(tmptitles[i-profile1_nseqs-1],titles[i]);
			tmplen_array[i-profile1_nseqs-1]=seqlen_array[i];
			tmpindex[i-profile1_nseqs-1]=output_index[i]-tmpfs+profile1_nseqs;
			tmp_array[i-profile1_nseqs-1]=(char *)ckalloc((seqlen_array[i]+2)*sizeof(char));
			for(j=1;j<=seqlen_array[i];j++)
				tmp_array[i-profile1_nseqs-1][j]=seq_array[i][j];
		}
	}

	strcpy(seqname,filename);
	GetPanelExtra(prf_panel[0].names,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(prf_panel[0].names,&data);
	GetPanelExtra(prf_panel[0].seqs,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(prf_panel[0].seqs,&data);
        profile_no = 1;
        n=profile_input();
	if (n<=0)
	{
		info("File %s not loaded.",seqname);
		return;
	}
	strcpy(profile1_name,seqname);
	load_aln(prf_panel[0],0,profile1_nseqs-1,TRUE);

	if(tmpn!=0)
	{
		nseqs=tmpn+profile1_nseqs;
		realloc_aln(profile1_nseqs+1,nseqs);
		for(i=profile1_nseqs+1;i<=nseqs;i++)
		{
			names[i]=(char *)ckalloc((MAXNAMES+2)*sizeof(char));
			titles[i]=(char *)ckalloc((MAXTITLES+2)*sizeof(char));

			strcpy(names[i],tmpnames[i-profile1_nseqs-1]);
			ckfree(tmpnames[i-profile1_nseqs-1]);
			strcpy(titles[i],tmptitles[i-profile1_nseqs-1]);
			ckfree(tmptitles[i-profile1_nseqs-1]);
			seqlen_array[i]=tmplen_array[i-profile1_nseqs-1];
			output_index[i]=tmpindex[i-profile1_nseqs-1]-tmpfs+profile1_nseqs;
			seq_array[i]=(char *)ckalloc((seqlen_array[i]+2)*sizeof(char));
			for(j=1;j<=seqlen_array[i];j++)
				seq_array[i][j]=tmp_array[i-profile1_nseqs-1][j];
			ckfree(tmp_array[i-profile1_nseqs-1]);
		}
		ckfree(tmpnames);
		ckfree(tmptitles);
		ckfree(tmplen_array);
		ckfree(tmpindex);
		ckfree(tmp_array);
		profile2_empty=FALSE;
	}
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,TRUE);

	ncutseqs=0;

	info("File %s loaded.",profile1_name);
}

static void OpenPrf2File (IteM item)
{
	int n;
	panel_data data;

	if(profile1_empty)
	{
		error("You must load profile 1 first.");
		return;
	}

	if (nseqs>profile1_nseqs)
	{
		if (Message(MSG_YN,"Replace existing sequences ?")==ANS_NO)
			return;
	}

	if (!GetInputFileName (filename,FILENAMELEN,"","")) return;
 
	strcpy(seqname,filename);
	GetPanelExtra(prf_panel[1].names,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(prf_panel[1].names,&data);
	GetPanelExtra(prf_panel[1].seqs,&data);
	data.nseqs=0;
	data.vseqs=0;
	SetPanelExtra(prf_panel[1].seqs,&data);
        profile_no = 2;
        n=profile_input();
	if (n<=0)
	{
		info("File %s not loaded.",seqname);
		return;
	}
	strcpy(profile2_name,seqname);
	ncutseqs=0;
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,TRUE);

	info("File %s loaded.",profile2_name);
}


static void BlackandWhite(IteM item)
{

	ncolors=1;
	
	if (aln_mode == MULTIPLEM)
		color_seqs();
	else
	{
		color_prf1();
		color_prf2();
	}
	usebw=TRUE;
	usedefcolors=FALSE;
	useusercolors=FALSE;
        SetStatus(bw_item,usebw);
        SetStatus(defcol_item,usedefcolors);
        SetStatus(usercol_item,useusercolors);
	info("Done.");
}


static void DefColorPar(IteM item)
{

	if (explicit_par_file != NULL)
		ckfree(explicit_par_file);
	explicit_par_file=NULL;
	if(dnaflag)
		par_file=find_file(def_dnapar_file);
	else
		par_file=find_file(def_protpar_file);
	init_color_parameters(par_file);
	if (aln_mode == MULTIPLEM)
		color_seqs();
	else
	{
		color_prf1();
		color_prf2();
	}
	usebw=FALSE;
	usedefcolors=TRUE;
	useusercolors=FALSE;
        SetStatus(bw_item,usebw);
        SetStatus(defcol_item,usedefcolors);
        SetStatus(usercol_item,useusercolors);
	info("Done.");
}

void set_reset_new_gaps(IteM i)
{
        reset_alignments_new=GetStatus(i);
	if(reset_alignments_new==TRUE)
	{
		reset_alignments_all=FALSE;
        	SetStatus(all_gaps_item,reset_alignments_all);
	}
}
void set_reset_all_gaps(IteM i)
{
        reset_alignments_all=GetStatus(i);
	if(reset_alignments_all==TRUE)
	{
		reset_alignments_new=FALSE;
        	SetStatus(new_gaps_item,reset_alignments_new);
	}
}
 

static void OpenColorParWin(IteM item)
{

	read_file_window("Input Color File","COLOR PARAMETER FILE NAME:",explicit_par_file,OpenColorPar);
}

static void OpenColorPar(ButtoN but)
{
	GetTitle(readfiletext, filename, FILENAMELEN);
	stripspace(filename);
 
	if (explicit_par_file != NULL)
		ckfree(explicit_par_file);
	explicit_par_file=(char *)ckalloc(FILENAMELEN*sizeof(char));
	if (par_file != NULL)
		ckfree(par_file);
	par_file=(char *)ckalloc(FILENAMELEN*sizeof(char));

	strcpy(explicit_par_file,filename);
	strcpy(par_file,filename);
        info("Loading color file: %s\n",par_file);
	init_color_parameters(par_file);
	if (Visible(readfilew))
	{
		Remove(readfilew);
		readfilew=NULL;
	}
	if (aln_mode == MULTIPLEM)
		color_seqs();
	else
	{
		color_prf1();
		color_prf2();
	}
	usebw=FALSE;
	usedefcolors=FALSE;
	useusercolors=TRUE;
        SetStatus(bw_item,usebw);
        SetStatus(defcol_item,usedefcolors);
        SetStatus(usercol_item,useusercolors);
	info("Done.");
}

static void RemoveGapPos(IteM item)
{
	int i,j,sl;
	Boolean sel=FALSE;

	if (nseqs==0)
	{
		Message(MSG_OK,"No file loaded.");
		return;
	}

	if (Message(MSG_YN,"Remove positions that contain gaps in all sequences ?")==ANS_NO)
		return;

	if(aln_mode==MULTIPLEM)
	{
		remove_gap_pos(1,nseqs,0);
		load_aln(seq_panel,0,nseqs-1,FALSE);
	}
	else
	{
		remove_gap_pos(1,profile1_nseqs,1);
		load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
		remove_gap_pos(profile1_nseqs+1,nseqs,2);
		load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	}
	info("Gap positions removed.");
}


static void RemoveGaps(IteM item)
{
	int i,j,sl;
	panel_data data;
	Boolean sel=FALSE;

	if (nseqs==0)
	{
		Message(MSG_OK,"No file loaded.");
		return;
	}
	GetPanelExtra(active_panel.names,&data);
	for (i=0;i<data.nseqs;i++)
		if(data.selected[i]==TRUE)
		{
			sel=TRUE;
			break;
		}
	if(sel==FALSE)
	{
		Message(MSG_OK,"Select sequences by clicking on the names.");
		return;
	}

	if (Message(MSG_YN,"Remove gaps from selected sequences ?")==ANS_NO)
		return;

	for (i=data.firstseq+1;i<=data.firstseq+data.nseqs;i++)
		if(data.selected[i-data.firstseq-1]==TRUE)
		{
                	sl=0;
                	for(j=1;j<=seqlen_array[i];++j) {
                        	if((seq_array[i][j] == gap_pos1) ||
                        	   (seq_array[i][j] == gap_pos2)) continue;
                        	++sl;
                        	seq_array[i][sl]=seq_array[i][j];
                	}
                        seq_array[i][sl+1]=-3;
                	seqlen_array[i]=sl;
		}
	load_aln(active_panel,data.firstseq,data.firstseq+data.nseqs-1,FALSE);
	active_panel.modified=TRUE;
	info("Gaps in selected sequences removed.");
}

static void CutSequences(IteM item)
{
	int i,pos;
	Boolean sel=FALSE;
	panel_data data;

	if (nseqs==0)
	{
		Message(MSG_OK,"No file loaded.");
		return;
	}

	GetPanelExtra(active_panel.names,&data);
	for (i=0;i<data.nseqs;i++)
		if(data.selected[i]==TRUE)
		{
			sel=TRUE;
			pos=i;
			break;
		}
	if(sel==FALSE)
	{
		Message(MSG_OK,"Select sequences to be cut by clicking on the names.");
		return;
	}

	if(ncutseqs>0)
	{
		if (Message(MSG_YN,"The previously cut sequences will be lost.\nDo you want to continue?")==ANS_NO) return;
	}

	if (saveseqlen_array!=NULL) ckfree(saveseqlen_array);
	if (saveseq_array!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (saveseq_array[i]!=NULL) ckfree(saveseq_array[i]);
		}
		ckfree(saveseq_array);
	}
	if (savetitles!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (savetitles[i]!=NULL) ckfree(savetitles[i]);
		}
		ckfree(savetitles);
	}
	if (savenames!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (savenames[i]!=NULL) ckfree(savenames[i]);
		}
		ckfree(savenames);
	}
	ncutseqs=0;

	savenames=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	savetitles=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	saveseq_array=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	saveseqlen_array=(sint *)ckalloc((data.nseqs+1) * sizeof(sint));
	for(i=0;i<data.nseqs;i++)
	{
		savenames[i]=NULL;
		savetitles[i]=NULL;
		saveseq_array[i]=NULL;
	}
	if (data.prf_no == 0)
		cut_multiplem();
	else if (data.prf_no == 1)
		cut_profile1();
	else if (data.prf_no == 2)
		cut_profile2();

	GetPanelExtra(active_panel.names,&data);
	if(pos>=data.nseqs) pos=data.nseqs-1;
	if(data.nseqs>0)
		data.selected[pos]=TRUE;
	SetPanelExtra(active_panel.names,&data);
	DrawPanel(active_panel.names);

	active_panel.modified=TRUE;
	info("Cut %d sequences.",ncutseqs);
}

static void cut_multiplem(void)
{
	int i,j;
	panel_data data;

        GetPanelExtra(active_panel.names,&data);
	for (i=data.nseqs;i>0;i--)
	{
		if(data.selected[i-1]==TRUE)
		{
			ssave(i);
			for(j=i;j<data.nseqs;j++)
				sscpy(j,j+1);
		}
	}
        nseqs-=ncutseqs;
        if (nseqs<=0) empty=TRUE;
	if (ncutseqs>0)
		if(nseqs<=data.vseqs)
			load_aln(active_panel,0,nseqs-1,TRUE);
		else
			load_aln(active_panel,0,nseqs-1,FALSE);
}

static void cut_profile1(void)
{
	int i,j;
	panel_data data;

        GetPanelExtra(active_panel.names,&data);
	for (i=data.nseqs;i>0;i--)
	{
		if(data.selected[i-1]==TRUE)
		{
			ssave(i);
			for(j=i;j<nseqs;j++)
				sscpy(j,j+1);
		}
	}
        profile1_nseqs-=ncutseqs;
	nseqs-=ncutseqs;
        if (profile1_nseqs<=0) profile1_empty=TRUE;
        if (nseqs<=0) empty=TRUE;
	if (ncutseqs>0)
	{
		if(profile1_nseqs<=data.vseqs)
			load_aln(active_panel,0,profile1_nseqs-1,TRUE);
		else
			load_aln(active_panel,0,profile1_nseqs-1,FALSE);
		if (!profile2_empty)
			load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	}
}

static void cut_profile2(void)
{
	int i,j;
	panel_data data;

        GetPanelExtra(active_panel.names,&data);
	for (i=data.nseqs;i>0;i--)
	{
		if(data.selected[i-1]==TRUE)
		{
			ssave(i+profile1_nseqs);
			for(j=i+profile1_nseqs;j<nseqs;j++)
				sscpy(j,j+1);
		}
	}
        nseqs-=ncutseqs;
        if (nseqs-profile1_nseqs<=0) profile2_empty=TRUE;
        if (nseqs<=0) empty=TRUE;
	if (ncutseqs>0)
		if(nseqs-profile1_nseqs<=data.vseqs)
			load_aln(active_panel,profile1_nseqs,nseqs-1,FALSE);
		else
			load_aln(active_panel,profile1_nseqs,nseqs-1,TRUE);
}

static void PasteSequences(IteM item)
{
	int insert;
	int i,n;
	panel_data data;

	if (ncutseqs<=0)
	{
		Message(MSG_OK,"No sequences available for pasting.\n"
                  " Cut selected sequences first.");
		 return;
	}

	GetPanelExtra(active_panel.names,&data);
	n=ncutseqs;
	insert=-1;
        if (data.nseqs>0)
        {
                for(i=data.nseqs-1;i>=0;i--)
                        if(data.selected[i]==TRUE)
                        {
                                insert=i;
                                break;
                        }
                if (insert==-1)
                {
                        Message(MSG_OK,"Select a sequence by clicking on the name.\n"
                        " Cut sequences will be pasted after this one.");
			return;
                }
        }

	if (data.prf_no == 2)
	{
		insert += profile1_nseqs;
		for(i=profile1_nseqs+data.nseqs;i>insert+1;i--)
			sscpy(i+ncutseqs,i);
		for(i=1;ncutseqs>0;i++)
			sload(insert+i+1);
	}


       	else
	{
		for(i=nseqs;i>insert+1;i--)
			sscpy(i+ncutseqs,i);
		for(i=1;ncutseqs>0;i++)
			sload(insert+i+1);
	}

        if(data.prf_no==0)
        {
                nseqs=data.nseqs+n;
                if (nseqs>0) empty=FALSE;
		load_aln(seq_panel,0,nseqs-1,FALSE);
        }
        else if(data.prf_no==1)
        {
                profile1_nseqs=data.nseqs+n;
                nseqs+=n;
                if (profile1_nseqs>0) profile1_empty=FALSE;
		load_aln(active_panel,0,profile1_nseqs-1,FALSE);
		if (!profile2_empty)
			load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
        }
        else if(data.prf_no==2)
        {
                nseqs=profile1_nseqs+data.nseqs+n;
                if (profile1_nseqs<nseqs)
		{
			profile2_empty=FALSE;
			empty=FALSE;
		}
/*
		load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
*/
		load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
        }

	active_panel.modified=TRUE;
	info("Pasted %d sequences.",n);
}

/* copies a sequence from clustal arrays position j to temp arrays */
static void ssave(int j)
{
	int k;

	if (saveseq_array[ncutseqs] != NULL) ckfree(saveseq_array[ncutseqs]);
	if (savenames[ncutseqs] != NULL) ckfree(savenames[ncutseqs]);
	if (savetitles[ncutseqs] != NULL) ckfree(savetitles[ncutseqs]);
	savenames[ncutseqs]=(char *)ckalloc((MAXNAMES+2)*sizeof(char));
	savetitles[ncutseqs]=(char *)ckalloc((MAXTITLES+2)*sizeof(char));

	strcpy(savenames[ncutseqs],names[j]);
	strcpy(savetitles[ncutseqs],titles[j]);
	saveseqlen_array[ncutseqs]=seqlen_array[j];
	saveseq_array[ncutseqs]=(char *)ckalloc((seqlen_array[j]+2)*sizeof(char));
	for(k=1;k<=seqlen_array[j];k++)
		saveseq_array[ncutseqs][k]=seq_array[j][k];
	saveseq_array[ncutseqs][k]= -3;
	ncutseqs++;
}

/* copies a sequence from clustal arrays position i to position j */
static void sscpy(int i,int j)
{
	int k;


	strcpy(names[i],names[j]);
	strcpy(titles[i],titles[j]);
	seqlen_array[i]=seqlen_array[j];
	realloc_seq(i,seqlen_array[i]);

	for(k=1;k<=seqlen_array[j];k++)
		seq_array[i][k]=seq_array[j][k];
	seq_array[i][k]= -3;
}

/* copies last sequence in temp arrays to clustal arrays after entry i */
static void sload(int i)
{
	int k;

	if (ncutseqs<1) return;

	ncutseqs--;
	strcpy(names[i],savenames[ncutseqs]);
	strcpy(titles[i],savetitles[ncutseqs]);
	seqlen_array[i]=saveseqlen_array[ncutseqs];
	realloc_seq(i,seqlen_array[i]);
	for(k=1;k<=seqlen_array[i];k++)
		seq_array[i][k]=saveseq_array[ncutseqs][k];
	seq_array[i][k]= -3;
}

static void SelectSeqs(IteM item)
{
	select_seqs(seq_panel,TRUE);
}

static void SelectPrf1(IteM item)
{
	select_seqs(prf_panel[0],TRUE);
}

static void SelectPrf2(IteM item)
{
	select_seqs(prf_panel[1],TRUE);
}

static void MergeProfiles(IteM item)
{
        if (profile2_empty)
        {
                error("Profile 2 not loaded");
                return;
        }
	profile_no=1;
	profile1_nseqs=nseqs;
	profile2_empty=TRUE;
	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	active_panel=prf_panel[0];

	info("Added Profile 2 to Profile 1.");
}


static void ClearSeqRange(IteM item)
{
	if(aln_mode==MULTIPLEM)
		clear_seqrange(seq_panel);
	else
	{
		clear_seqrange(prf_panel[1]);
		clear_seqrange(prf_panel[0]);
	}
}

static void ClearSeqs(IteM item)
{
	if(aln_mode==MULTIPLEM)
		select_seqs(seq_panel,FALSE);
	else
	{
		select_seqs(prf_panel[1],FALSE);
		select_seqs(prf_panel[0],FALSE);
	}
}

static void clear_seqrange(spanel p)
{
	int f,l;
	panel_data data;

	GetPanelExtra(p.seqs,&data);
	f=data.firstsel;
	l=data.lastsel;
	data.firstsel=data.lastsel=-1;
	SetPanelExtra(p.seqs,&data);
	highlight_seqrange(p.seqs,f,l,NORMAL);
}

static void select_seqs(spanel p,Boolean flag)
{
	int i;
	panel_data data;

	GetPanelExtra(p.names,&data);
	if (data.nseqs == 0) return;

	for (i=0;i<data.nseqs;i++)
		data.selected[i]=flag;

	SetPanelExtra(p.names,&data);
	draw_names(p.names);
	if(flag==TRUE) active_panel=p;
}

static void CAlignWin (IteM item)
{
        if (empty)
        {
                error("No sequences loaded");
                return;
        }
        if (nseqs <= 1)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }
	do_align_window(&calignw,&ctreetext,NEW,"Complete Alignment",CompleteAlign);
}

void CompleteAlign(ButtoN but)
{
	char phylip_name[FILENAMELEN];

	GetTitle(ctreetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);

	if (!open_aln_files()) return;
 
	WatchCursor();
	if (Visible(calignw))
	{
		Remove(calignw);
		calignw=NULL;
	}
        align(phylip_name);

	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
/* reload the sequences from the output file (so that the sequence order is
correct - either INPUT or ALIGNED , don't output messages */
	reload_alignment();

        load_aln(seq_panel,0,nseqs-1,FALSE);
	ArrowCursor();
}

static void RealignSeqsWin (IteM item)
{
	int i;
	Boolean sel=FALSE;
	panel_data data;

        if (empty)
        {
                error("No sequences loaded");
                return;
        }
        if (nseqs <= 1)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }

/* check some sequences have been selected */

	GetPanelExtra(seq_panel.names,&data);
	for (i=0;i<data.nseqs;i++)
		if(data.selected[i]==TRUE)
		{
			sel=TRUE;
			break;
		}
	if(sel==FALSE)
	{
		Message(MSG_OK,"Select sequences to be realigned\n"
                               "by clicking on the names.");
		return;
	}

	do_align_window(&ralignw,&rtreetext,NEW,"Realign Sequences",RealignSeqs);
}

static void RealignSeqs(ButtoN but)
{
	int insert;
	int i,j,n;
	panel_data data;
	char phylip_name[FILENAMELEN];

	GetTitle(rtreetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);

	if (!open_aln_files()) return;

/* cut selected sequences */

	GetPanelExtra(seq_panel.names,&data);

	if (saveseqlen_array!=NULL) ckfree(saveseqlen_array);
	if (saveseq_array!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (saveseq_array[i]!=NULL) ckfree(saveseq_array[i]);
		}
		ckfree(saveseq_array);
	}
	if (savetitles!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (savetitles[i]!=NULL) ckfree(savetitles[i]);
		}
		ckfree(savetitles);
	}
	if (savenames!=NULL)
	{
		for(i=0;i<ncutseqs;i++)
		{
			if (savenames[i]!=NULL) ckfree(savenames[i]);
		}
		ckfree(savenames);
	}
	ncutseqs=0;

	savenames=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	savetitles=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	saveseq_array=(char **)ckalloc((data.nseqs+1) * sizeof(char *));
	saveseqlen_array=(sint *)ckalloc((data.nseqs+1) * sizeof(sint));
	for(i=0;i<data.nseqs;i++)
	{
		savenames[i]=NULL;
		savetitles[i]=NULL;
		saveseq_array[i]=NULL;
	}
	for (i=data.nseqs;i>0;i--)
	{
		if(data.selected[i-1]==TRUE)
		{
			ssave(i);
			for(j=i;j<data.nseqs;j++)
				sscpy(j,j+1);
		}
	}
        nseqs=data.nseqs-ncutseqs;
        if (nseqs<=0) empty=TRUE;

/* paste selected sequences at the end */
	n=ncutseqs;
	profile1_nseqs=nseqs;
	insert=profile1_nseqs-1;

	for(i=nseqs;i>insert+1;i--)
		sscpy(i+ncutseqs,i);
	for(i=1;ncutseqs>0;i++)
		sload(insert+i+1);

        nseqs=profile1_nseqs+n;

/* align profile 2 sequences to profile 1 */
	WatchCursor();
	if (Visible(ralignw))
	{
		Remove(ralignw);
		ralignw=NULL;
	}
        new_sequence_align(phylip_name);

	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
/* reload the sequences from the output file (so that the sequence order is
correct - either INPUT or ALIGNED */
	reload_alignment();
	load_aln(seq_panel,0,nseqs-1,FALSE);

	GetPanelExtra(seq_panel.names,&data);
	for (i=0;i<profile1_nseqs;i++)
		data.selected[i]=FALSE;
	for (i=profile1_nseqs;i<nseqs;i++)
		data.selected[i]=TRUE;
	SetPanelExtra(seq_panel.names,&data);
	draw_names(seq_panel.names);
	ArrowCursor();
	info("Selected sequences realigned.");
}

static void RealignSeqRangeWin (IteM item)
{
	panel_data data;
	GrouP aligngr;
	GrouP output_list;
	ButtoN align_ok, align_can;
	GrouP maing;
	PopuP end_gap_toggle;
	char name[FILENAMELEN+1];
	char path[FILENAMELEN+1];

        if (empty)
        {
                error("No sequences loaded");
                return;
        }
        if (nseqs <= 1)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }

/* check a range has been selected */

	GetPanelExtra(seq_panel.seqs,&data);
	if(data.firstsel==-1)
	{
		Message(MSG_OK,"Select residue range to be realigned\n"
                               "by clicking in the sequence display area.");
		return;
	}


	get_path(seqname,path);
	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	rralignw=FixedWindow(-50, -33, -10, -10,"Realign Residue Range",RemoveWin);

	maing=HiddenGroup(rralignw,2,0,NULL);
	SetGroupSpacing(maing,0,10);

	make_prompt(rralignw, "Output Guide Tree File:");
        stdLineHeight=18;
        SelectFont(programFont);
	Break(rralignw);
	rrtreetext=DialogText(rralignw, "", 35, NULL);
	strcpy(name,path);
	strcat(name,"dnd");
	SetTitle(rrtreetext, name);
	Break(rralignw);

	make_prompt(rralignw, "Output Alignment Files:");
	output_list=HiddenGroup(rralignw, 2, 0, NULL);
	if(output_clustal) {
		make_prompt(output_list,"Clustal: ");
		cl_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"aln");
		SetTitle(cl_outtext, name);
		Break(output_list);
	}
	if(output_nbrf) {
		make_prompt(output_list,"NBRF/PIR: ");
		pir_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"pir");
		SetTitle(pir_outtext, name);
		Break(output_list);
	}
	if(output_gcg) {
		make_prompt(output_list,"GCG/MSF: ");
		msf_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"msf");
		SetTitle(msf_outtext, name);
		Break(output_list);
	}
	if(output_phylip) {
		make_prompt(output_list,"Phylip: ");
		phylip_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"phy");
		SetTitle(phylip_outtext, name);
		Break(output_list);
	}
	if(output_gde) {
		make_prompt(output_list,"GDE: ");
		gde_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"gde");
		SetTitle(gde_outtext, name);
		Break(output_list);
	}
	if(output_nexus) {
		make_prompt(output_list,"Nexus: ");
		nexus_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"nxs");
		SetTitle(nexus_outtext, name);
		Break(output_list);
	}

	Break(rralignw);
	end_gap_toggle=make_toggle(rralignw,"Realign Segment End Gap Penalties","ON","OFF",&realign_endgappenalties,set_realign_endgappenalties);
	Break(rralignw);

	aligngr=HiddenGroup(rralignw, 2, 0, NULL);
	shift(aligngr, 60, 20);
	align_ok=PushButton(aligngr, " ALIGN ", RealignSeqRange);
	shift(aligngr, 20,0);
	align_can=PushButton(aligngr, "CANCEL", CancelWin);

	Show(rralignw);
}

static void RealignSeqRange(ButtoN but)
{
	int i,j;
	int fs,save_order,length,length1,length2;
	panel_data data;
	sint *tmplen_array;
	char **tmp_array;
	sint *newlen_array;
	char **new_array;
	char phylip_name[FILENAMELEN];

	GetTitle(rrtreetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);

	if (!open_aln_files()) return;

	WatchCursor();
	if (Visible(rralignw))
	{
		Remove(rralignw);
		rralignw=NULL;
	}
/* save the alignment into a temporary area */

	GetPanelExtra(seq_panel.seqs,&data);

	tmplen_array=(sint *)ckalloc((data.nseqs+2) * sizeof(sint));
	tmp_array=(char **)ckalloc((data.nseqs+2) * sizeof(char *));
	for (i=1;i<=data.nseqs;i++)
	{
		tmplen_array[i]=seqlen_array[i];
		tmp_array[i]=(char *)ckalloc((data.ncols+2) * sizeof(char));
		for(j=1;j<=seqlen_array[i];j++)
			tmp_array[i][j]=seq_array[i][j];
		for(j=seqlen_array[i]+1;j<=data.ncols;j++)
			tmp_array[i][j]=gap_pos2;
	}

/* copy the selected residue range to the clustal alignment arrays */

	fs=data.firstsel;
	length=data.lastsel-data.firstsel+1;
	max_aln_length=2*length;
	for (i=1;i<=data.nseqs;i++)
	{
		seqlen_array[i]=length;
		realloc_seq(i,length);
		for(j=data.firstsel;j<=data.lastsel;j++)
			seq_array[i][j-data.firstsel+1]=tmp_array[i][j+1];
		seq_array[i][j-data.firstsel+1]=-3;
	}
/* temporarily set the output order to be the same as the input */
	save_order=output_order;
	output_order=INPUT;
/* set the end gaps penalties */
	endgappenalties=realign_endgappenalties;

/* align the residue range */
        align(phylip_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}

	output_order=save_order;
/* reset the end gaps penalties */
	endgappenalties=align_endgappenalties;

/* remove positions that contain just gaps */
	remove_gap_pos(1,nseqs,0);


/* save the new alignment into another temporary area */
	newlen_array=(sint *)ckalloc((data.nseqs+2) * sizeof(sint));
	new_array=(char **)ckalloc((data.nseqs+2) * sizeof(char *));
	for (i=1;i<=data.nseqs;i++)
	{
		newlen_array[i]=seqlen_array[i];
		new_array[i]=(char *)ckalloc((seqlen_array[i]+2) * sizeof(char));
		for(j=1;j<=seqlen_array[i];j++)
			new_array[i][j]=seq_array[i][j];
	}

/* paste the realigned range back into the alignment */
	max_aln_length=0;
	length1=length2=0;
	for (i=1;i<=data.nseqs;i++)
	{
		length1=tmplen_array[i]-length+newlen_array[i];
		if(length1>max_aln_length) max_aln_length=length1;
		length2=newlen_array[i];
		seqlen_array[i]=length1;
		realloc_seq(i,length1);
		for(j=1;j<=data.firstsel;j++)
			seq_array[i][j]=tmp_array[i][j];
		for(j=data.firstsel+1;j<=data.firstsel+length2;j++)
			seq_array[i][j]=new_array[i][j-data.firstsel];
		for(j=data.firstsel+length2+1;j<=length1;j++)
			seq_array[i][j]=tmp_array[i][data.lastsel+j-data.firstsel-length2+1];
	}
	max_aln_length*=2;
	ckfree(tmplen_array);
	for(i=1;i<=data.nseqs;i++)
		ckfree(tmp_array[i]);
	ckfree(tmp_array);
	ckfree(newlen_array);
	for(i=1;i<=data.nseqs;i++)
		ckfree(new_array[i]);
	ckfree(new_array);

	if (open_aln_files())
        	create_alignment_output(1,data.nseqs);

	load_aln(seq_panel,0,nseqs-1,FALSE);
	GetPanelExtra(seq_panel.seqs,&data);
	data.firstsel=fs;
	data.lastsel=data.firstsel+length2-1;
	SetPanelExtra(seq_panel.seqs,&data);
	highlight_seqrange(seq_panel.seqs,data.firstsel,data.lastsel,HIGHLIGHT);
	ArrowCursor();
	info("Selected sequence range realigned.");
}


void AlignFromTreeWin(IteM item)
{
        if (empty)
        {
                error("No sequences loaded");
                return;
        }
        if (nseqs < 2)
	{
                error("Alignment has only %d sequences",nseqs);
                return;
        }
	do_align_window(&talignw,&ttreetext,OLD,"Alignment from Guide Tree",AlignFromTree);
}

static void do_align_window(WindoW *ralignw,TexT *rtreetext,Boolean treestatus,char *title,void align_proc(ButtoN but))
{
	WindoW alignw;
	TexT treetext;
	GrouP aligngr;
	GrouP output_list;
	ButtoN align_ok, align_can;
	GrouP maing;
	char name[FILENAMELEN+1];
	char path[FILENAMELEN+1];

	get_path(seqname,path);
	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	alignw=FixedWindow(-50, -33, -10, -10,title,RemoveWin);

	maing=HiddenGroup(alignw,2,0,NULL);
	SetGroupSpacing(maing,0,10);

	if(treestatus==NEW)
		make_prompt(alignw, "Output Guide Tree File:");
	else 
		make_prompt(alignw, "Input Guide Tree File:");
       	stdLineHeight=18;
       	SelectFont(programFont);
	Break(alignw);
	treetext=DialogText(alignw, "", 35, NULL);
	strcpy(name,path);
	strcat(name,"dnd");
	SetTitle(treetext, name);
	Break(alignw);

	make_prompt(alignw, "Output Alignment Files:");
	output_list=HiddenGroup(alignw, 2, 0, NULL);
	if(output_clustal) {
		make_prompt(output_list,"Clustal: ");
		cl_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"aln");
		SetTitle(cl_outtext, name);
		Break(output_list);
	}
	if(output_nbrf) {
		make_prompt(output_list,"NBRF/PIR: ");
		pir_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"pir");
		SetTitle(pir_outtext, name);
		Break(output_list);
	}
	if(output_gcg) {
		make_prompt(output_list,"GCG/MSF: ");
		msf_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"msf");
		SetTitle(msf_outtext, name);
		Break(output_list);
	}
	if(output_phylip) {
		make_prompt(output_list,"Phylip: ");
		phylip_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"phy");
		SetTitle(phylip_outtext, name);
		Break(output_list);
	}
	if(output_gde) {
		make_prompt(output_list,"GDE: ");
		gde_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"gde");
		SetTitle(gde_outtext, name);
		Break(output_list);
	}
	if(output_nexus) {
		make_prompt(output_list,"Nexus: ");
		nexus_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"nxs");
		SetTitle(nexus_outtext, name);
		Break(output_list);
	}

	Break(alignw);
	aligngr=HiddenGroup(alignw, 2, 0, NULL);
	shift(aligngr, 60, 20);
	align_ok=PushButton(aligngr, " ALIGN ", align_proc);
	shift(aligngr, 20,0);
	align_can=PushButton(aligngr, "CANCEL", CancelWin);

	*ralignw=alignw;
	*rtreetext=treetext;
	Show(alignw);
}


static void do_palign_window(WindoW *ralignw,TexT *rtree1text,TexT *rtree2text,Boolean treestatus,char *title,void align_proc(ButtoN but))
{
	Boolean istree=FALSE;
	WindoW alignw;
	TexT tree1text,tree2text;
	GrouP aligngr;
	GrouP output_list;
	ButtoN align_ok, align_can;
	GrouP maing;
	char name[FILENAMELEN+1];
	char path[FILENAMELEN+1];

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	alignw=FixedWindow(-50, -33, -10, -10,title,RemoveWin);

	maing=HiddenGroup(alignw,2,0,NULL);
	SetGroupSpacing(maing,0,10);

	if(treestatus==NEW) 
		make_prompt(alignw, "Output Guide Tree Files:");
	else 
		make_prompt(alignw, "Input Guide Tree Files:");
        stdLineHeight=18;
        SelectFont(programFont);
	Break(alignw);
	tree1text=DialogText(alignw, "", 35, NULL);
	get_path(profile1_name,path);
	strcpy(name,path);
	strcat(name,"dnd");
	SetTitle(tree1text, name);
	Break(alignw);
	tree2text=DialogText(alignw, "", 35, NULL);
	get_path(profile2_name,path);
	strcpy(name,path);
	strcat(name,"dnd");
	SetTitle(tree2text, name);
	Break(alignw);

	make_prompt(alignw, "Output Alignment Files:");
	output_list=HiddenGroup(alignw, 2, 0, NULL);
	if(output_clustal) {
		make_prompt(output_list,"Clustal: ");
		cl_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"aln");
		SetTitle(cl_outtext, name);
		Break(output_list);
	}
	if(output_nbrf) {
		make_prompt(output_list,"NBRF/PIR: ");
		pir_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"pir");
		SetTitle(pir_outtext, name);
		Break(output_list);
	}
	if(output_gcg) {
		make_prompt(output_list,"GCG/MSF: ");
		msf_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"msf");
		SetTitle(msf_outtext, name);
		Break(output_list);
	}
	if(output_phylip) {
		make_prompt(output_list,"Phylip: ");
		phylip_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"phy");
		SetTitle(phylip_outtext, name);
		Break(output_list);
	}
	if(output_gde) {
		make_prompt(output_list,"GDE: ");
		gde_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"gde");
		SetTitle(gde_outtext, name);
		Break(output_list);
	}
	if(output_nexus) {
		make_prompt(output_list,"Nexus: ");
		nexus_outtext=DialogText(output_list, "", 35, NULL);
		strcpy(name,path);
		strcat(name,"nxs");
		SetTitle(nexus_outtext, name);
		Break(output_list);
	}

	Break(alignw);
	aligngr=HiddenGroup(alignw, 2, 0, NULL);
	shift(aligngr, 60, 20);
	align_ok=PushButton(aligngr, " ALIGN ", align_proc);
	shift(aligngr, 20,0);
	align_can=PushButton(aligngr, "CANCEL", CancelWin);

	*ralignw=alignw;
	*rtree1text=tree1text;
	*rtree2text=tree2text;
	Show(alignw);
}


void AlignFromTree(ButtoN but)
{
	FILE *tree;
	char phylip_name[FILENAMELEN];

	GetTitle(ttreetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);
#ifdef VMS
        if((tree=fopen(phylip_name,"r","rat=cr","rfm=var"))==NULL) {
#else
        if((tree=fopen(phylip_name,"r"))==NULL) {
#endif
                error("Cannot open tree file [%s]",phylip_name);
                return;
        }

	if (!open_aln_files()) return;

	WatchCursor();
	info("Doing alignments from guide tree...");
	if (Visible(talignw))
	{
		Remove(talignw);
		talignw=NULL;
	}
	get_tree(phylip_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
/* reload the sequences from the output file (so that the sequence order is
correct - either INPUT or ALIGNED */
	reload_alignment();
        load_aln(seq_panel,0,nseqs-1,FALSE);
	ArrowCursor();
	info("Done.");
}

static void PrfPrfAlignWin (IteM item)
{
        if (profile1_empty)
        {
                error("Profile 1 not loaded");
                return;
        }
        if (profile2_empty)
        {
                error("Profile 2 not loaded");
                return;
        }
	do_palign_window(&palignw,&ptree1text,&ptree2text,NEW,"Profile to Profile Alignment",PrfPrfAlign);
}

static void PrfPrfTreeAlignWin (IteM item)
{
        if (profile1_empty)
        {
                error("Profile 1 not loaded");
                return;
        }
        if (profile2_empty)
        {
                error("Profile 2 not loaded");
                return;
        }
	do_palign_window(&palignw,&ptree1text,&ptree2text,OLD,"Profile Alignment from Tree",PrfPrfTreeAlign);
}

static void SeqPrfAlignWin (IteM item)
{
        if (profile1_empty)
        {
                error("Profile 1 not loaded");
                return;
        }
        if (profile2_empty)
        {
                error("Profile 2 not loaded");
                return;
        }
	do_align_window(&salignw,&streetext,NEW,"Sequence to Profile Alignment",SeqPrfAlign);
}

static void SeqPrfTreeAlignWin (IteM item)
{
        if (profile1_empty)
        {
                error("Profile 1 not loaded");
                return;
        }
        if (profile2_empty)
        {
                error("Profile 2 not loaded");
                return;
        }
	do_align_window(&salignw,&streetext,OLD,"Sequence to Profile Alignment from Tree",SeqPrfTreeAlign);
}

static void PrfPrfAlign(ButtoN but)
{
	char p1_tree_name[FILENAMELEN];
	char p2_tree_name[FILENAMELEN];

	GetTitle(ptree1text, filename, FILENAMELEN);
	stripspace(filename);
	use_tree1_file=FALSE;
       	strcpy(p1_tree_name,filename);

	GetTitle(ptree2text, filename, FILENAMELEN);
	stripspace(filename);
	use_tree2_file=FALSE;
       	strcpy(p2_tree_name,filename);

	if (!open_aln_files()) return;

	WatchCursor();
	if (Visible(palignw))
	{
		Remove(palignw);
		palignw=NULL;
	}
        profile_align(p1_tree_name,p2_tree_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	ArrowCursor();
}

static void PrfPrfTreeAlign(ButtoN but)
{
	char p1_tree_name[FILENAMELEN];
	char p2_tree_name[FILENAMELEN];

	GetTitle(ptree1text, filename, FILENAMELEN);
	stripspace(filename);
	if(filename[0]!=EOS) use_tree1_file=TRUE;
       	strcpy(p1_tree_name,filename);

	GetTitle(ptree2text, filename, FILENAMELEN);
	stripspace(filename);
	if(filename[0]!=EOS) use_tree2_file=TRUE;
       	strcpy(p2_tree_name,filename);

	if (!open_aln_files()) return;

	WatchCursor();
	if (Visible(palignw))
	{
		Remove(palignw);
		palignw=NULL;
	}
        profile_align(p1_tree_name,p2_tree_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}
	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	ArrowCursor();
}

static void SeqPrfAlign(ButtoN but)
{
	char phylip_name[FILENAMELEN];

	GetTitle(streetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);
	use_tree_file=FALSE;

	if (!open_aln_files()) return;

	WatchCursor();
	if (Visible(salignw))
	{
		Remove(salignw);
		salignw=NULL;
	}
        new_sequence_align(phylip_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}

/* reload the sequences from the output file (so that the sequence order is
correct - either INPUT or ALIGNED */
	reload_alignment();

	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	ArrowCursor();
}

static void SeqPrfTreeAlign(ButtoN but)
{
	char phylip_name[FILENAMELEN];

	GetTitle(streetext, filename, FILENAMELEN);
	stripspace(filename);

        strcpy(phylip_name,filename);
	use_tree_file=TRUE;

	if (!open_aln_files()) return;

	WatchCursor();
	if (Visible(salignw))
	{
		Remove(salignw);
		salignw=NULL;
	}
        new_sequence_align(phylip_name);
	if(save_log && save_log_fd!=NULL)
	{
		fclose(save_log_fd);
		save_log_fd=NULL;
	}

/* reload the sequences from the output file (so that the sequence order is
correct - either INPUT or ALIGNED */
	reload_alignment();

	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
	ArrowCursor();
}
void reload_alignment(void)
{
	int i,k;
	sint     *sseqlen_array;
	char     **sseq_array;
	char     **snames, **stitles;

	if (nseqs==0) return;
	if (output_order == INPUT) return;


	snames=(char **)ckalloc((nseqs+2) * sizeof(char *));
	stitles=(char **)ckalloc((nseqs+2) * sizeof(char *));
	sseq_array=(char **)ckalloc((nseqs+2) * sizeof(char *));
	sseqlen_array=(sint *)ckalloc((nseqs+2) * sizeof(sint));
	for (i=1;i<=nseqs;i++)
	{
		snames[i]=(char *)ckalloc((MAXNAMES+2)*sizeof(char));
		stitles[i]=(char *)ckalloc((MAXTITLES+2)*sizeof(char));
		sseq_array[i]=(char *)ckalloc((seqlen_array[output_index[i]]+2)*sizeof(char));
		strcpy(snames[i],names[output_index[i]]);
		strcpy(stitles[i],titles[output_index[i]]);
		sseqlen_array[i]=seqlen_array[output_index[i]];
		for(k=1;k<=seqlen_array[output_index[i]];k++)
			sseq_array[i][k]=seq_array[output_index[i]][k];
	}
	for (i=1;i<=nseqs;i++)
	{
		strcpy(names[i],snames[i]);
		strcpy(titles[i],stitles[i]);
		seqlen_array[i]=sseqlen_array[i];
		realloc_seq(i,seqlen_array[i]);
		for(k=1;k<=seqlen_array[i];k++)
			seq_array[i][k]=sseq_array[i][k];
		output_index[i]=i;
	}

	ckfree(sseqlen_array);
	for(i=1;i<=nseqs;i++)
		ckfree(sseq_array[i]);
	ckfree(sseq_array);
	for(i=1;i<=nseqs;i++)
		ckfree(stitles[i]);
	ckfree(stitles);
	for(i=1;i<=nseqs;i++)
		ckfree(snames[i]);
	ckfree(snames);
	ncutseqs=0;
}

static void SegmentWin(IteM item)
{
	WindoW w;
	GrouP maing;
	ButtoN closeb;
	GrouP mat_list;
	ButtoN matrixb[5];

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	w=FixedWindow(-50, -33, -10, -10, "Low-Scoring Segment Parameters",RemoveWin);
	maing=HiddenGroup(w,0,0,NULL);
	SetGroupSpacing(maing,20,10);
	closeb=PushButton(maing, "CLOSE", CancelWin);
	Break(maing);

	/*PushButton(maing, "Calculate Low-Scoring Segments", calc_segment_exceptions);
	Break(maing);*/

	length_cutofftext=make_scale(maing,"Minimum Length of Segments:",9,length_cutoff,19,set_lengthcutoff);
	Break(maing);

	segmentdnascaletext=make_scale(maing,"DNA Marking Scale:",9,segment_dnascale,9,set_segment_dnascale);
	if(!dnaflag) Disable(segmentdnascaletext);
	Break(maing);


	mat_list=NormalGroup(maing,4,0,"Protein Weight Matrix",systemFont,set_segment_matrix);
	matrixb[0]=RadioButton(mat_list,"Gonnet PAM 80");
	matrixb[1]=RadioButton(mat_list,"Gonnet PAM 120");
	matrixb[2]=RadioButton(mat_list,"Gonnet PAM 250");
	matrixb[3]=RadioButton(mat_list,"Gonnet PAM 350");
	matrixb[4]=RadioButton(mat_list,"User defined");
	SetValue(mat_list,segment_matnum);
	seg_matrix_list=mat_list;
	Break(maing);
	PushButton(maing, "Load protein matrix: ", set_segment_user_matrix);
	Advance(maing);
	segmentmattext=StaticPrompt(maing,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(segmentmattext,segment_mtrxname);
	Break(maing);

	mat_list=NormalGroup(maing,4,0,"DNA Weight Matrix",systemFont,set_segment_dnamatrix);
	matrixb[0]=RadioButton(mat_list,"IUB");
	matrixb[1]=RadioButton(mat_list,"CLUSTALW(1.6)");
	matrixb[2]=RadioButton(mat_list,"User defined");
	SetValue(mat_list,segment_dnamatnum);
	seg_dnamatrix_list=matrix_list;
	Break(maing);
	PushButton(maing, "Load DNA matrix: ", set_segment_user_dnamatrix);
	Advance(maing);
	segmentdnamattext=StaticPrompt(maing,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(segmentdnamattext,segment_dnamtrxname);
	Break(maing);


	Show(w);
}

static void ScoreWin(IteM item)
{
	WindoW w;
	GrouP maing;
	ButtoN closeb;
	GrouP mat_list;
	PopuP show_exceptions;
	ButtoN matrixb[5];

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	w=FixedWindow(-50, -33, -10, -10, "Score Parameters",RemoveWin);
	maing=HiddenGroup(w,0,0,NULL);
	SetGroupSpacing(maing,20,10);
	closeb=PushButton(maing, "CLOSE", CancelWin);
	Break(maing);


/* add a scale to set the scaling value for the alignment scoring function */
	scorescaletext=make_scale(maing,"Score Plot Scale:",9,score_scale,9,set_scorescale);
	Break(maing);

	residue_cutofftext=make_scale(maing,"Residue Exception Cutoff:",9,score_cutoff,9,set_scorecutoff);
	Break(maing);

	mat_list=NormalGroup(maing,4,0,"Protein Weight Matrix",systemFont,set_score_matrix);
	matrixb[0]=RadioButton(mat_list,"Identity");
	matrixb[1]=RadioButton(mat_list,"Gonnet PAM 80");
	matrixb[2]=RadioButton(mat_list,"Gonnet PAM 120");
	matrixb[3]=RadioButton(mat_list,"Gonnet PAM 250");
	matrixb[4]=RadioButton(mat_list,"Gonnet PAM 350");
	matrixb[5]=RadioButton(mat_list,"User defined");
	SetValue(mat_list,score_matnum);
	score_matrix_list=mat_list;
	Break(maing);
	PushButton(maing, "Load protein matrix: ", set_score_user_matrix);
	Advance(maing);
	scoremattext=StaticPrompt(maing,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(scoremattext,score_mtrxname);
	Break(maing);

	mat_list=NormalGroup(maing,4,0,"DNA Weight Matrix",systemFont,set_score_dnamatrix);
	matrixb[0]=RadioButton(mat_list,"IUB");
	matrixb[1]=RadioButton(mat_list,"CLUSTALW(1.6)");
	matrixb[2]=RadioButton(mat_list,"User defined");
	SetValue(mat_list,score_dnamatnum);
	score_dnamatrix_list=mat_list;
	Break(maing);
	PushButton(maing, "Load DNA matrix: ", set_score_user_dnamatrix);
	Advance(maing);
	scorednamattext=StaticPrompt(maing,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(scorednamattext,score_dnamtrxname);
	Break(maing);

	Show (w);
}


static void PWParameters(IteM item)
{
	int i;
	WindoW w;
	PoinT pt;
	GrouP maing;
	ButtoN closeb;
	TexT go_scale,ge_scale;
	TexT gp_scale,ktuple_scale,topdiags_scale,window_scale;
	PopuP fs_toggle;
	GrouP mat_list;
	ButtoN matrixb[5];
	char str[FILENAMELEN];

        if(dnaflag) {
                gap_open   = dna_gap_open;
                gap_extend = dna_gap_extend;
                pw_go_penalty     = dna_pw_go_penalty;
                pw_ge_penalty     = dna_pw_ge_penalty;
                ktup       = dna_ktup;
                window     = dna_window;
                signif     = dna_signif;
                wind_gap   = dna_wind_gap;
 
        }
        else {
                gap_open   = prot_gap_open;
                gap_extend = prot_gap_extend;
                pw_go_penalty     = prot_pw_go_penalty;
                pw_ge_penalty     = prot_pw_ge_penalty;
                ktup       = prot_ktup;
                window     = prot_window;
                signif     = prot_signif;
                wind_gap   = prot_wind_gap;
 
        }

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	w=FixedWindow(-50, -33, -10, -10, "Pairwise Parameters",RemoveWin);
	maing=HiddenGroup(w,0,0,NULL);
	SetGroupSpacing(maing,0,10);
	closeb=PushButton(maing, "CLOSE", CancelWin);
	Break(maing);
	fs_toggle=make_toggle(maing,"Pairwise Alignments :","Fast-Approximate","Slow-Accurate",&quick_pairalign,set_fs_toggle);
	Break(maing);

	GetNextPosition(maing,&pt);
	slow_para=NormalGroup(maing,0,0,"Pairwise Parameters",systemFont,NULL);
	SetGroupSpacing(slow_para,0,10);

	make_prompt(slow_para, "Gap Opening [0-100] :");
	Advance(slow_para);
	sprintf(str,"%.2f",pw_go_penalty);
	go_scale=DialogText(slow_para, str, 5, set_pw_go_penalty);
	Break(slow_para);

	make_prompt(slow_para, "Gap Extension [0-100] :");
	Advance(slow_para);
	sprintf(str,"%.2f",pw_ge_penalty);
	ge_scale=DialogText(slow_para, str, 5, set_pw_ge_penalty);
	Break(slow_para);
	mat_list=NormalGroup(slow_para,4,0,"Protein Weight Matrix",systemFont,set_pw_matrix);
	for(i=0;i<pw_matrix_menu.noptions;i++)
		matrixb[i]=RadioButton(mat_list,pw_matrix_menu.opt[i].title);
	SetValue(mat_list,pw_matnum);
	pw_matrix_list=mat_list;
	Break(slow_para);
	PushButton(slow_para, "Load protein matrix: ", set_pw_user_matrix);
	Advance(slow_para);
	pwmattext=StaticPrompt(slow_para,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(pwmattext,pw_usermtrxname);

	Break(slow_para);
	mat_list=NormalGroup(slow_para,4,0,"DNA Weight Matrix",systemFont,set_pw_dnamatrix);
	for(i=0;i<dnamatrix_menu.noptions;i++)
		matrixb[i]=RadioButton(mat_list,dnamatrix_menu.opt[i].title);
	SetValue(mat_list,pw_dnamatnum);
	pw_dnamatrix_list=mat_list;
	Break(slow_para);
	PushButton(slow_para, "Load DNA matrix: ", set_pw_user_dnamatrix);
	Advance(slow_para);
	pwdnamattext=StaticPrompt(slow_para,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(pwdnamattext,pw_dnausermtrxname);

	Break(slow_para);


/* fast parameters */
	SetNextPosition(maing,pt);
	fast_para=NormalGroup(maing,2,0,"Pairwise Parameters",systemFont,NULL);
	SetGroupSpacing(fast_para,0,10);
	make_prompt(fast_para, "Gap Penalty [1-500]:");
	sprintf(str,"%d",wind_gap);
	gp_scale=DialogText(fast_para, str, 3, set_gp);
	make_prompt(fast_para, "K-Tuple Size [1-2]:");
	sprintf(str,"%d",ktup);
	ktuple_scale=DialogText(fast_para, str, 1, set_ktuple);
	make_prompt(fast_para, "Top Diagonals [1-50]:");
	sprintf(str,"%d",signif);
	topdiags_scale=DialogText(fast_para, str, 2, set_topdiags);
	make_prompt(fast_para, "Window Size [1-50]:");
	sprintf(str,"%d",window);
	window_scale=DialogText(fast_para, str, 2, set_window);

	if (quick_pairalign)
	{
		Hide(slow_para);
		Show(fast_para);
	}
	else
	{
		Hide(fast_para);
		Show(slow_para);
	}

	Break(maing);

	Show (w);
}


static void MultiParameters(IteM item)
{
	int i;
	WindoW w;
	GrouP maing;
	ButtoN closeb;
	TexT go_scale,ge_scale;
	GrouP mat_list;
	ButtoN matrixb[5];
	GrouP multi_para;
	TexT div_seq;
	TexT transitions;
	PopuP neg_mat_toggle;
	PopuP end_gap_toggle;
	char str[FILENAMELEN];

        if(dnaflag) {
                gap_open   = dna_gap_open;
                gap_extend = dna_gap_extend;
                pw_go_penalty     = dna_pw_go_penalty;
                pw_ge_penalty     = dna_pw_ge_penalty;
                ktup       = dna_ktup;
                window     = dna_window;
                signif     = dna_signif;
                wind_gap   = dna_wind_gap;
 
        }
        else {
                gap_open   = prot_gap_open;
                gap_extend = prot_gap_extend;
                pw_go_penalty     = prot_pw_go_penalty;
                pw_ge_penalty     = prot_pw_ge_penalty;
                ktup       = prot_ktup;
                window     = prot_window;
                signif     = prot_signif;
                wind_gap   = prot_wind_gap;
 
        }

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	w=FixedWindow(-50, -33, -10, -10, "Alignment Parameters",RemoveWin);
	maing=HiddenGroup(w,0,0,NULL);
	SetGroupSpacing(maing,0,10);
	closeb=PushButton(maing, "CLOSE", CancelWin);
	Break(maing);
 
/* multiple alignment parameters */

	multi_para=NormalGroup(maing,0,0,"Multiple Parameters",systemFont,NULL);
	SetGroupSpacing(multi_para,0,10);
        make_prompt(multi_para, "Gap Opening [0-100] :");
	Advance(multi_para);
	sprintf(str,"%.2f",gap_open);
        go_scale=DialogText(multi_para, str, 5, set_go_penalty);
	Advance(multi_para);
        make_prompt(multi_para, "Gap Extention [0-100] :");
	Advance(multi_para);
	sprintf(str,"%.2f",gap_extend);
        ge_scale=DialogText(multi_para, str, 5, set_ge_penalty);
	Break(multi_para);
        make_prompt(multi_para, "Delay Divergent Sequences (%) :");
	Advance(multi_para);
	sprintf(str,"%d",divergence_cutoff);
        div_seq=DialogText(multi_para, str, 3, set_div_seq);
	Break(multi_para);
 
        make_prompt(multi_para, "DNA Transition Weight [0-1] :");
        Advance(multi_para);
        sprintf(str,"%.2f",transition_weight);
        transitions=DialogText(multi_para, str, 5, set_transitions);
        Break(multi_para);


	neg_mat_toggle=make_toggle(multi_para,"Use Negative Matrix","ON","OFF",&neg_matrix,set_neg_matrix);

	Break(multi_para);
        mat_list=NormalGroup(multi_para,2,0,"Protein Weight Matrix",systemFont,set_matrix);
	for(i=0;i<matrix_menu.noptions;i++)
		matrixb[i]=RadioButton(mat_list,matrix_menu.opt[i].title);
        SetValue(mat_list,matnum);
	matrix_list=mat_list;
	Break(multi_para);
	PushButton(multi_para, "Load protein matrix: ", set_user_matrix);
	Advance(multi_para);
	mattext=StaticPrompt(multi_para,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(mattext,usermtrxname);

	Break(multi_para);
        mat_list=NormalGroup(multi_para,2,0,"DNA Weight Matrix",systemFont,set_dnamatrix);
	for(i=0;i<dnamatrix_menu.noptions;i++)
		matrixb[i]=RadioButton(mat_list,dnamatrix_menu.opt[i].title);
        SetValue(mat_list,dnamatnum);
	dnamatrix_list=mat_list;
	Break(multi_para);
	PushButton(multi_para, "Load DNA: ", set_user_dnamatrix);
	Advance(multi_para);
	dnamattext=StaticPrompt(multi_para,"", MAXPROMPTLEN, dialogTextHeight, systemFont, 'l');
	SetTitle(dnamattext,dnausermtrxname);
	Show (w);
}

static void GapParameters(IteM item)
{
	WindoW gapparaw;
	GrouP maing;
	ButtoN closeb;
	PopuP rp_toggle,vp_toggle,hp_toggle,end_gap_toggle;
	TexT gdist,hyd_text;
	char str[80];
 
	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	gapparaw=FixedWindow(-50, -33, -10, -10, "Protein Gap Parameters",RemoveWin);
	maing=HiddenGroup(gapparaw,0,0,NULL);
	SetGroupSpacing(maing,0,10);
	closeb=PushButton(maing, "CLOSE", CancelWin);
	Break(maing);
	rp_toggle=make_toggle(maing,"Residue-specific Penalties","OFF","ON",&no_pref_penalties,set_pref_penalties);
	Break(maing);
	hp_toggle=make_toggle(maing,"Hydrophilic Penalties","OFF","ON",&no_hyd_penalties,set_hyd_penalties);
	Break(maing);
        make_prompt(maing, "Hydrophilic Residues :");
	Advance(maing);
        hyd_text=DialogText(maing, hyd_residues, 20, set_hyd_res);
	Break(maing);
        make_prompt(maing, "Gap Separation Distance [0-100] :");
	Advance(maing);
	sprintf(str,"%d",gap_dist);
        gdist=DialogText(maing, str, 3, set_gap_dist);
	Break(maing);
	end_gap_toggle=make_toggle(maing,"End Gap Separation","ON","OFF",&use_endgaps,set_endgaps);

	Show (gapparaw);
}

static void SSParameters(IteM item)
{
	WindoW ssparaw;
	GrouP maing;
	ButtoN closeb;
	PopuP use_p1,use_p2;
	TexT helix_gp,strand_gp,loop_gp,terminal_gp,helix_minus,helix_plus;
	TexT strand_minus,strand_plus;
	GrouP output_list;
	ButtoN outputb[4]; 
	char str[80];
 
	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	ssparaw=FixedWindow(-50, -33, -10, -10, "Secondary Structure Options",RemoveWin);
	closeb=PushButton(ssparaw, "CLOSE", CancelWin);
	Break(ssparaw);

	use_p1=make_toggle(ssparaw,"Use profile 1 secondary structure / penalty mask","YES","NO",&use_ss1,set_use_ss1);
	Break(ssparaw);
	use_p2=make_toggle(ssparaw,"Use profile 2 secondary structure / penalty mask","YES","NO",&use_ss2,set_use_ss2);

	Break(ssparaw);
        output_list=NormalGroup(ssparaw,2,0,"Output ",systemFont,NULL)
;
        outputb[0]=CheckBox(output_list,"Secondary Structure",set_ss_output);
	if(output_struct_penalties==0 || output_struct_penalties==2) 
		output_ss=TRUE;
	else
		output_ss=FALSE;
	SetStatus(outputb[0],output_ss);
        outputb[1]=CheckBox(output_list,"Gap Penalty Mask",set_gp_output);
	if(output_struct_penalties==1 || output_struct_penalties==2) 
		output_gp=TRUE;
	else
		output_gp=FALSE;
	SetStatus(outputb[1],output_gp);

	Break(ssparaw);
	maing=HiddenGroup(ssparaw,2,0,NULL);
	SetGroupSpacing(maing,0,10);
        make_prompt(maing, "Helix Gap Penalty [0-9] :");
	sprintf(str,"%d",helix_penalty);
        helix_gp=DialogText(maing, str, 1, set_helix_gp);
        make_prompt(maing, "Strand Gap Penalty [0-9] :");
	sprintf(str,"%d",strand_penalty);
        strand_gp=DialogText(maing, str, 1, set_strand_gp);
        make_prompt(maing, "Loop Gap Penalty [0-9] :");
	sprintf(str,"%d",loop_penalty);
        loop_gp=DialogText(maing, str, 1, set_loop_gp);
        make_prompt(maing, "Secondary Structure Terminal Penalty [0-9] :");
	sprintf(str,"%d",helix_end_penalty);
        terminal_gp=DialogText(maing, str, 1, set_terminal_gp);
        make_prompt(maing, "Helix Terminal Positions [0-3]            within:");
	sprintf(str,"%d",helix_end_minus);
        helix_minus=DialogText(maing, str, 1, set_helix_minus);
        make_prompt(maing, "outside:");
	sprintf(str,"%d",helix_end_plus);
        helix_plus=DialogText(maing, str, 1, set_helix_plus);
        make_prompt(maing, "Strand Terminal Penalty [0-3]             within:");
	sprintf(str,"%d",strand_end_minus);
        strand_minus=DialogText(maing, str, 1, set_strand_minus);
        make_prompt(maing, "outside:");
	sprintf(str,"%d",strand_end_plus);
        strand_plus=DialogText(maing, str, 1, set_strand_plus);


	Show (ssparaw);
}

static void OutputParameters(IteM item)
{
	WindoW outputparaw;
	GrouP maing;
	ButtoN closeb;
	GrouP output_list;
	ButtoN outputb[6];
	PopuP order_toggle,para_toggle;
	PopuP case_toggle,snos_toggle;

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	outputparaw=FixedWindow(-50, -33, -10, -10, "Output Format Options",RemoveWin);
	closeb=PushButton(outputparaw, "CLOSE", CancelWin);
	Break(outputparaw);
        output_list=NormalGroup(outputparaw,2,0,"Output Files",systemFont,NULL);

	SelectFont(systemFont);
        outputb[0]=CheckBox(output_list,"CLUSTAL format",set_output_clustal);
	if(output_clustal) SetStatus(outputb[0],TRUE);

        outputb[1]=CheckBox(output_list,"NBRF/PIR format",set_output_nbrf);
	if(output_nbrf) SetStatus(outputb[1],TRUE);

        outputb[2]=CheckBox(output_list,"GCG/MSF format",set_output_gcg);
	if(output_gcg) SetStatus(outputb[2],TRUE);

        outputb[3]=CheckBox(output_list,"PHYLIP format",set_output_phylip);
	if(output_phylip) SetStatus(outputb[3],TRUE);

        outputb[4]=CheckBox(output_list,"GDE format",set_output_gde);
	if(output_gde) SetStatus(outputb[4],TRUE);

        outputb[5]=CheckBox(output_list,"NEXUS format",set_output_nexus);
	if(output_nexus) SetStatus(outputb[5],TRUE);

	maing=HiddenGroup(outputparaw,2,0,NULL);
	SetGroupSpacing(maing,0,10);


	case_toggle=make_toggle(maing,"GDE output case :","Lower","Upper",&lowercase,set_case);
	snos_toggle=make_toggle(maing,"CLUSTALW sequence numbers :","ON","OFF",&cl_seq_numbers,set_snos);
       	make_prompt(maing, "Output order");
	order_toggle=PopupList(maing,TRUE,set_output_order);
	PopupItem(order_toggle,"INPUT");
	PopupItem(order_toggle,"ALIGNED");
	if (output_order == INPUT)
		SetValue(order_toggle,1);
	else if (output_order == ALIGNED)
		SetValue(order_toggle,2);
	para_toggle=make_toggle(maing,"Parameter output","ON","OFF",&save_parameters,set_save_paras);

	Show (outputparaw);
}


static void OutputTreeParameters(IteM item)
{
	WindoW outputtreeparaw;
	ButtoN closeb;
	GrouP output_list;
	ButtoN outputb[3];
	PopuP boot_format;

	SelectFont(systemFont);
	stdCharWidth=CharWidth('A');
	stdLineHeight=LineHeight();
	outputtreeparaw=FixedWindow(-50, -33, -10, -10, "Output Tree Format Options",RemoveWin);
	closeb=PushButton(outputtreeparaw, "CLOSE", CancelWin);
	Break(outputtreeparaw);
        output_list=NormalGroup(outputtreeparaw,2,0,"Output Files",systemFont,NULL);
        outputb[0]=CheckBox(output_list,"CLUSTAL format tree",set_output_tree_clustal);
	if(output_tree_clustal) SetStatus(outputb[0],TRUE);
        outputb[1]=CheckBox(output_list,"Phylip format tree",set_output_tree_phylip);
	if(output_tree_phylip) SetStatus(outputb[1],TRUE);
        outputb[2]=CheckBox(output_list,"Phylip distance matrix",set_output_tree_distances);
	if(output_tree_distances) SetStatus(outputb[2],TRUE);
        outputb[3]=CheckBox(output_list,"Nexus format tree",set_output_tree_nexus);
	if(output_tree_nexus) SetStatus(outputb[3],TRUE);

	Break(outputtreeparaw);
       	make_prompt(outputtreeparaw, "Bootstrap labels on:");
	Advance(outputtreeparaw);
	boot_format=PopupList(outputtreeparaw,TRUE,set_boot_format);
	PopupItem(boot_format ,"NODE");
	PopupItem(boot_format ,"BRANCH");
	if (bootstrap_format == BS_NODE_LABELS)
		SetValue(boot_format ,1);
	if (bootstrap_format == BS_BRANCH_LABELS)
		SetValue(boot_format ,2);
	Show (outputtreeparaw);
}

static PrompT make_prompt(GrouP g,CharPtr title)
{
	PrompT p=NULL;

	if (title != NULL)
        	p=StaticPrompt(g, title, 0, dialogTextHeight, systemFont, 'l');

	return p;
}

static PrompT make_scale(GrouP g,CharPtr title,int length, int value,int max,BarScrlProc SetProc)
{
	char str[FILENAMELEN];
	BaR scale;
	PrompT t;

	sprintf(str,"%s   %3d",title,value);
	t=make_prompt(g,str);
	Advance(g);
        scale=ScrollBar(g, length, -1, SetProc);
       	CorrectBarPage(scale,(Int4)1,(Int4)1);
        CorrectBarMax(scale,(Int4)max);
        CorrectBarValue(scale,(Int4)value);
	return t;
}

static PopuP make_toggle(GrouP g,CharPtr title,CharPtr true_text, CharPtr false_text, Boolean *value,PupActnProc SetProc)
{
	PopuP p;

	if (title != NULL)
        	make_prompt(g, title);
	Advance(g);
	p=PopupList(g,TRUE,SetProc);
	PopupItem(p,true_text);
	PopupItem(p,false_text);
	if (*value)
		SetValue(p,1);
	else
		SetValue(p,2);
	return p;
}


void switch_mode(void)
{
	char path[FILENAMELEN];

	if(aln_mode==MULTIPLEM)
	{
		Hide(prf1_display);
		Hide(prf2_display);
		Hide(pscrolltext);
		SetValue(modetext,1);
		resize_multi_window();
		profile_no=0;
		check_menus(file_item,PROFILEM);
		check_menus(align_item,PROFILEM);
		check_menus(edit_item,PROFILEM);
		check_menus(tree_item,PROFILEM);
		check_menus(color_item,PROFILEM);
		active_panel=seq_panel;
		fix_gaps();
		load_aln_data(seq_panel,0,nseqs-1,TRUE);
		Show(seq_display);
	}
	else if(aln_mode==PROFILEM)
	{
		Hide(seq_display);
		resize_prf_window(nseqs,0);
		SetValue(modetext,2);
		profile_no=1;
		profile1_nseqs=nseqs;
		if (profile1_nseqs > 0) profile1_empty = FALSE;
		profile2_empty = TRUE;

		check_menus(file_item,MULTIPLEM);
		check_menus(align_item,MULTIPLEM);
		check_menus(edit_item,MULTIPLEM);
		check_menus(tree_item,MULTIPLEM);
		check_menus(color_item,MULTIPLEM);
		active_panel=prf_panel[0];
		get_path(seqname,path);
		strcpy(profile1_name,path);
		strcat(profile1_name,"1.");
		strcpy(profile2_name,path);
		strcat(profile2_name,"2.");
		fix_gaps();
		load_aln_data(prf_panel[0],0,profile1_nseqs-1,TRUE);
		load_aln_data(prf_panel[1],profile1_nseqs,nseqs-1,TRUE);
		Show(prf1_display);
		Show(prf2_display);
		Show(pscrolltext);
	}
		
}


static void make_menu_headers(WindoW w)
{
        filem = PulldownMenu (w,"File");
        editm = PulldownMenu (w,"Edit");
        alignm = PulldownMenu (w,"Alignment");
	treem = PulldownMenu(w,"Trees");
	colorm = PulldownMenu(w,"Colors");
	scorem = PulldownMenu(w,"Quality");
	helpmenu = PulldownMenu(w,"Help");

}

static void make_file_menu(void)
{
	int n=0;

	file_item.mode[n] = MULTIPLEM;
       	file_item.i[n] = CommandItem (filem,"Load Sequences", OpenSeqFile); n++;

	file_item.mode[n] = MULTIPLEM;
       	file_item.i[n] = CommandItem (filem,"Append Sequences", AppendSeqFile); n++;
	file_item.mode[n] = MULTIPLEM;
       	file_item.i[n] = CommandItem (filem,"Save Sequences as...", SaveSeqFileWin); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Load Profile 1", OpenPrf1File); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Load Profile 2", OpenPrf2File); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Save Profile 1 as...", SavePrf1FileWin); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Save Profile 2 as...", SavePrf2FileWin); n++;
	file_item.mode[n] = MULTIPLEM;
       	file_item.i[n] = CommandItem (filem,"Write Alignment as PostScript", SavePSSeqWin); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Write Profile 1 as PostScript", SavePSPrf1Win); n++;
	file_item.mode[n] = PROFILEM;
       	file_item.i[n] = CommandItem (filem,"Write Profile 2 as PostScript", SavePSPrf2Win); n++;
       	file_item.i[n] = CommandItem (filem,"Quit", QuitWinI); n++;
	file_item.num = n;
	if(aln_mode==MULTIPLEM)
		check_menus(file_item,PROFILEM);
	else
		check_menus(file_item,MULTIPLEM);
}

static void QuitWinI (IteM i)
{
	if(aln_mode == MULTIPLEM)
	{
		if(seq_panel.modified)
			if (Message(MSG_YN,"Alignment has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
	}
	else if(aln_mode == PROFILEM)
	{
		if(prf_panel[0].modified)
			if (Message(MSG_YN,"Profile 1 has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
		if(prf_panel[1].modified)
			if (Message(MSG_YN,"Profile 2 has not been saved.\n"
			"Quit program anyway?")==ANS_NO) return;
	}
	QuitProgram ();
}

static void make_score_menu(void)
{
	int n=0;

	score_item.i[n] = CommandItem (scorem,"Calculate Low-Scoring Segments", calc_segment_exceptions); n++;
	segment_item=score_item.i[n]=StatusItem(scorem, "Show Low-Scoring Segments", set_show_segments);
        SetStatus(score_item.i[n],segment_exceptions); n++;
	score_item.i[n]=StatusItem(scorem, "Show Exceptional Residues", set_residue_exceptions);
        SetStatus(score_item.i[n],residue_exceptions); n++;
       	score_item.i[n] = CommandItem (scorem,"Low-Scoring Segment Parameters",SegmentWin); n++;
	score_item.i[n] = CommandItem (scorem,"Column Score Parameters",ScoreWin); n++;
	score_item.mode[n] = MULTIPLEM;
	score_item.i[n]=CommandItem(scorem, "Save Column Scores to File", SaveScoresWin); n++;


	score_item.num = n;
}

static void make_help_menu(void)
{
	int n=0;

	help_item.ptr[n] = 'G';
       	help_item.i[n] = CommandItem (helpmenu,"General",HelpProc); n++;
	help_item.ptr[n] = 'F';
       	help_item.i[n] = CommandItem (helpmenu,"Input & Output Files",HelpProc); n++;
	help_item.ptr[n] = 'E';
       	help_item.i[n] = CommandItem (helpmenu,"Editing Alignments",HelpProc); n++;
	help_item.ptr[n] = 'M';
       	help_item.i[n] = CommandItem (helpmenu,"Multiple Alignments",HelpProc); n++;
	help_item.ptr[n] = 'P';
       	help_item.i[n] = CommandItem (helpmenu,"Profile Alignments",HelpProc); n++;
	help_item.ptr[n] = 'B';
       	help_item.i[n] = CommandItem (helpmenu,"Secondary Structures",HelpProc); n++;
	help_item.ptr[n] = 'T';
       	help_item.i[n] = CommandItem (helpmenu,"Trees",HelpProc); n++;
	help_item.ptr[n] = 'C';
       	help_item.i[n] = CommandItem (helpmenu,"Colors",HelpProc); n++;
	help_item.ptr[n] = 'Q';
       	help_item.i[n] = CommandItem (helpmenu,"Alignment Quality",HelpProc); n++;
	help_item.ptr[n] = '9';
       	help_item.i[n] = CommandItem (helpmenu,"Command Line Parameters",HelpProc); n++;
	help_item.ptr[n] = 'R';
       	help_item.i[n] = CommandItem (helpmenu,"References",HelpProc); n++;

	help_item.num = n;
}

static void HelpProc(IteM item)
{
	int n,index=-1;
	FILE *fd;
	int  i, number, nlines;
	Boolean found_help;
	char temp[MAXLINE+1];
	char token;
	char *digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	char *help_marker    = ">>HELP";

	TexT htext;
	char *help_file = NULL;

	extern char *help_file_name;

	helptext[0]='\0';

	for(n=0;n<help_item.num;n++)
		if (item==help_item.i[n])
		{
			index=n;
			break;
		}

	if (index==-1)
	{
		error("Problem with HELP routines\n");
		return;
	}


	help_file=find_file(help_file_name);
	if(help_file==NULL) {
            error("Cannot find help file");
            return;
        }
	
        if((fd=fopen(help_file,"r"))==NULL) {
            error("Cannot open help file [%s]",help_file);
            return;
        }
	nlines = 0;
	number = -1;
	found_help = FALSE;

	while(TRUE) {
		if(fgets(temp,MAXLINE+1,fd) == NULL) {
			if(!found_help)
				error("No help found in help file");
			fclose(fd);
			return;
		}
		if(strstr(temp,help_marker)) {
                        token = ' ';
			for(i=strlen(help_marker); i<8; i++)
				if(strchr(digits, temp[i])) {
					token = temp[i];
					break;
				}
			if(token == help_item.ptr[index]) {
				found_help = TRUE;
				while(fgets(temp,MAXLINE+1,fd)) {
					if(strstr(temp, help_marker)) break;
					if(strlen(helptext)+strlen(temp) >MAXHELPLENGTH)
						break;
						for(i=strlen(temp)-1;i>=0;i--)
							if(iscntrl(temp[i])||isspace(temp[i]))
								temp[i]='\0';
							else break;
/* ignore lines starting with < - these are for html output processing */
					if(temp[0]!='<') {
			       			strcat(helptext,temp);
#ifdef WIN_MAC
			       			strcat(helptext,"\r");
#else
#ifdef WIN32
			       			strcat(helptext,"\r\n");
#else
			       			strcat(helptext,"\n");
#endif
#endif
			       			++nlines;
					}
				}
				fclose(fd);
#ifdef WIN_MAC
				if(numhelp>=1)
#else
				if(numhelp>=MAXHELPW)
#endif
				{
					error("Too many help windows");
					return;
				}
				numhelp++;
 
				helpw[numhelp]=FixedWindow(-50, -33, -10, -10, "", QuitHelpW);
				SelectFont(helpfont);
#ifdef WIN_MAC
				htext=ScrollText(helpw[numhelp], 60, 20, helpfont, TRUE, NULL);
#else
				htext=ScrollText(helpw[numhelp], 80, 30, helpfont, TRUE, NULL);
#endif
				Break(helpw[numhelp]);
				PushButton(helpw[numhelp], "OK", QuitHelpB);
				SetTitle(htext, helptext);
				Show(helpw[numhelp]);
				return;
			}
		}
	}
}

void QuitHelpB(ButtoN b)
{
 
	Remove(ParentWindow(b));
	numhelp--;
}
 
void QuitHelpW(WindoW w)
{
 
	Remove(w);
	numhelp--;
}
static void make_edit_menu(void)
{
	int n=0;

       	edit_item.i[n] = CommandItem (editm,"Cut Sequences", CutSequences); n++;
       	edit_item.i[n] = CommandItem (editm,"Paste Sequences", PasteSequences); n++;
	edit_item.mode[n] = MULTIPLEM;
       	edit_item.i[n] = CommandItem (editm,"Select All Sequences", SelectSeqs); n++;
	edit_item.mode[n] = PROFILEM;
       	edit_item.i[n] = CommandItem (editm,"Select Profile 1", SelectPrf1); n++;
	edit_item.mode[n] = PROFILEM;
       	edit_item.i[n] = CommandItem (editm,"Select Profile 2", SelectPrf2); n++;
	edit_item.mode[n] = PROFILEM;
       	edit_item.i[n] = CommandItem (editm,"Add Profile 2 to Profile 1", MergeProfiles); n++;
       	edit_item.i[n] = CommandItem (editm,"Clear Sequence Selection",ClearSeqs); n++;
       	edit_item.i[n] = CommandItem (editm,"Clear Range Selection",ClearSeqRange); n++;
	SeparatorItem(editm);
       	edit_item.i[n] = CommandItem (editm,"Search for String", SearchStrWin); n++;
	SeparatorItem(editm);
       	edit_item.i[n] = CommandItem (editm,"Remove All Gaps", RemoveGaps); n++;
       	edit_item.i[n] = CommandItem (editm,"Remove Gap-Only Columns", RemoveGapPos); n++;
	edit_item.num = n;
	if(aln_mode==MULTIPLEM)
		check_menus(edit_item,PROFILEM);
	else
		check_menus(edit_item,MULTIPLEM);
}

static void make_align_menu(void)
{
	MenU parasm;
	int n=0;

	align_item.mode[n] = MULTIPLEM;
       	align_item.i[n] = CommandItem (alignm,"Do Complete Alignment",CAlignWin); n++;
	align_item.mode[n] = MULTIPLEM;
       	align_item.i[n] = CommandItem (alignm,"Produce Guide Tree Only",SaveTreeWin); n++;
	align_item.mode[n] = MULTIPLEM;
       	align_item.i[n] = CommandItem (alignm,"Do Alignment from Guide Tree",AlignFromTreeWin); n++;
	SeparatorItem(alignm);
	align_item.mode[n] = MULTIPLEM;
       	align_item.i[n] = CommandItem (alignm,"Realign Selected Sequences",RealignSeqsWin); n++;
	align_item.mode[n] = MULTIPLEM;
       	align_item.i[n] = CommandItem (alignm,"Realign Selected Residue Range",RealignSeqRangeWin); n++;
	align_item.mode[n] = PROFILEM;
       	align_item.i[n] = CommandItem (alignm,"Align Profile 2 to Profile 1",PrfPrfAlignWin); n++;
	align_item.mode[n] = PROFILEM;
       	align_item.i[n] = CommandItem (alignm,"Align Profiles from Guide Trees",PrfPrfTreeAlignWin); n++;
	align_item.mode[n] = PROFILEM;
       	align_item.i[n] = CommandItem (alignm,"Align Sequences to Profile 1",SeqPrfAlignWin); n++;
	align_item.mode[n] = PROFILEM;
       	align_item.i[n] = CommandItem (alignm,"Align Sequences to Profile 1 from Tree",SeqPrfTreeAlignWin); n++;

       	SeparatorItem(alignm);
	parasm=SubMenu(alignm,"Alignment Parameters");
	new_gaps_item=align_item.i[n]=StatusItem(parasm, "Reset New Gaps before Alignment", set_reset_new_gaps);
	SetStatus(align_item.i[n],reset_alignments_new); n++;
	all_gaps_item=align_item.i[n]=StatusItem(parasm, "Reset All Gaps before Alignment", set_reset_all_gaps);
	SetStatus(align_item.i[n],reset_alignments_all); n++;
       	align_item.i[n] = CommandItem (parasm,"Pairwise Alignment Parameters",PWParameters); n++;
       	align_item.i[n] = CommandItem (parasm,"Multiple Alignment Parameters",MultiParameters); n++;
       	align_item.i[n] = CommandItem (parasm,"Protein Gap Parameters",GapParameters); n++;
	align_item.mode[n] = PROFILEM;
       	align_item.i[n] = CommandItem (parasm,"Secondary Structure Parameters",SSParameters); n++;
	align_item.i[n]=StatusItem(alignm, "Save Log File", set_save_log);
	save_item1=align_item.i[n];
	SetStatus(save_item1,save_log); n++;
       	align_item.i[n] = CommandItem (alignm,"Output Format Options",OutputParameters); n++;
	align_item.num = n;
	if(aln_mode==MULTIPLEM)
		check_menus(align_item,PROFILEM);
	else
		check_menus(align_item,MULTIPLEM);

}

void set_save_log(IteM i)
{
        save_log=GetStatus(i);
	SetStatus(save_item1,save_log);
	SetStatus(save_item2,save_log);
}

static void make_tree_menu(void)
{
	int n=0;

	tree_item.mode[n] = MULTIPLEM;
        tree_item.i[n] = CommandItem (treem,"Draw N-J Tree",DrawTreeWin); n++;
	tree_item.mode[n] = MULTIPLEM;
        tree_item.i[n] = CommandItem (treem,"Bootstrap N-J Tree",BootstrapTreeWin); n++;
        SeparatorItem(treem);
	tree_item.mode[n] = MULTIPLEM;
	tree_item.i[n]=StatusItem(treem, "Exclude Positions with Gaps", set_tossgaps);
	SetStatus(tree_item.i[n],tossgaps); n++;
	tree_item.mode[n] = MULTIPLEM;
	tree_item.i[n]=StatusItem(treem, "Correct for Multiple Substitutions", set_kimura);
	SetStatus(tree_item.i[n],kimura); n++;
        SeparatorItem(treem);
	tree_item.mode[n] = MULTIPLEM;
	tree_item.i[n]=StatusItem(treem, "Save Log File", set_save_log);
	save_item2=tree_item.i[n];
	SetStatus(save_item2,save_log); n++;
	tree_item.mode[n] = MULTIPLEM;
        tree_item.i[n] = CommandItem (treem,"Output Format Options",OutputTreeParameters); n++;
	tree_item.mode[n] = MULTIPLEM;
	tree_item.num = n;
	if(aln_mode==MULTIPLEM)
		check_menus(tree_item,PROFILEM);
	else
		check_menus(tree_item,MULTIPLEM);

}

static void make_color_menu(void)
{
	int n=0;

        color_item.i[n]=StatusItem(colorm, "Background Coloring", set_inverted);
        SetStatus(color_item.i[n],inverted); n++;
       	SeparatorItem(colorm);
        bw_item=color_item.i[n] = StatusItem (colorm,"Black and White",BlackandWhite);
        SetStatus(color_item.i[n],usebw); n++;
        defcol_item=color_item.i[n] = StatusItem (colorm,"Default Colors",DefColorPar);
        SetStatus(color_item.i[n],usedefcolors); n++;
        usercol_item=color_item.i[n] = StatusItem (colorm,"Load Color Parameter File",OpenColorParWin);
        SetStatus(color_item.i[n],useusercolors); n++;
	color_item.num = n;

}

void check_menus(menu_item m,int mode)
{
	int i;

	for (i=0;i<m.num;i++)
		if (m.mode[i] == mode)
			Disable(m.i[i]);
		else
			Enable(m.i[i]);
}


