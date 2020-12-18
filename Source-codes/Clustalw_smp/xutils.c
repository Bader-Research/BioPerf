#include <stdio.h>
#include <stdarg.h>
#include <string.h>
 
#include <vibrant.h>
 
#include "clustalw.h"
#include "xmenu.h"

static int get_series_matrixname(char *usermtrxname, short *usermat,short *aa_xref,int usermatnum,int *matnum,PrompT mattext);

extern Boolean x_menus;
extern WindoW mainw;
extern GrouP matrix_list,pw_matrix_list;
extern GrouP dnamatrix_list,pw_dnamatrix_list;
extern GrouP seg_matrix_list,seg_dnamatrix_list;
extern GrouP score_matrix_list,score_dnamatrix_list;
extern Boolean  interactive;
extern Boolean  dnaflag;
extern char     hyd_residues[];
extern sint     gap_dist;
extern Boolean  no_var_penalties, no_hyd_penalties, no_pref_penalties;
extern Boolean  use_endgaps;
extern Boolean  realign_endgappenalties;
extern Boolean  align_endgappenalties;
extern sint     divergence_cutoff;
extern Boolean  lowercase; /* Flag for GDE output - set on comm. line*/
extern Boolean  cl_seq_numbers;
extern sint     output_order;
extern Boolean  save_log;
extern Boolean  quick_pairalign;
extern Boolean  neg_matrix;
extern Boolean  output_clustal, output_nbrf, output_phylip, output_gcg, output_gde, output_nexus;
extern Boolean  save_parameters;
extern Boolean  output_tree_clustal, output_tree_phylip, output_tree_distances, output_tree_nexus;
extern char	seqname[];
extern float	transition_weight;
extern float    gap_open,      gap_extend;
extern float    dna_gap_open,  dna_gap_extend;
extern float    prot_gap_open, prot_gap_extend;
extern float    pw_go_penalty,      pw_ge_penalty;
extern float    dna_pw_go_penalty,  dna_pw_ge_penalty;
extern float    prot_pw_go_penalty, prot_pw_ge_penalty;
extern sint    wind_gap,ktup,window,signif;
extern sint    dna_wind_gap, dna_ktup, dna_window, dna_signif;
extern sint    prot_wind_gap,prot_ktup,prot_window,prot_signif;
extern Boolean tossgaps;  /* Ignore places in align. where ANY seq. has a gap*/
extern Boolean kimura;    /* Use correction for multiple substitutions */
extern sint boot_ntrials;               /* number of bootstrap trials */
extern unsigned sint boot_ran_seed;     /* random number generator seed */
extern sint    bootstrap_format;
extern sint struct_penalties,struct_penalties1,struct_penalties2;
extern sint output_struct_penalties;
extern sint    profile1_nseqs;
extern sint     nseqs;
extern Boolean use_ss1, use_ss2;
extern int inverted;
extern char     mtrxname[], pw_mtrxname[];
extern char     usermtrxname[], pw_usermtrxname[];
extern sint        matnum,pw_matnum;
extern short    usermat[], pw_usermat[];
extern short    aa_xref[], pw_aa_xref[];
extern char     dnamtrxname[], pw_dnamtrxname[];
extern char     dnausermtrxname[], pw_dnausermtrxname[];
extern sint        dnamatnum,pw_dnamatnum;
extern short    userdnamat[], pw_userdnamat[];
extern short    dna_xref[], pw_dna_xref[];
extern Boolean  use_ambiguities;

extern MatMenu matrix_menu;
extern MatMenu dnamatrix_menu;
extern MatMenu pw_matrix_menu;

extern sint        helix_penalty;
extern sint        strand_penalty;
extern sint        loop_penalty;
extern sint        helix_end_minus;
extern sint        helix_end_plus;
extern sint        strand_end_minus;
extern sint        strand_end_plus;
extern sint        helix_end_penalty;
extern sint        strand_end_penalty;

extern TexT savealntext;
extern GrouP slow_para,fast_para;

extern PrompT   message;           /* used in temporary message window */
extern Boolean mess_output;
extern FILE *save_log_fd;
extern color color_lut[];
extern spanel  seq_panel;        /* data for multiple alignment area */
extern spanel  prf_panel[];       /* data for profile alignment areas */
extern Boolean aln_mode;
extern Boolean fixed_prf_scroll;
extern Boolean output_ss;
extern Boolean output_gp;
extern PrompT mattext,pwmattext;
extern PrompT dnamattext,pwdnamattext;
extern int    save_format;
extern Boolean residue_exceptions;
extern Boolean segment_exceptions;
extern int font_size;
extern FonT datafont;
extern int av_font[];
extern TexT blocklentext;
extern IteM segment_item;

extern int      pagesize;
extern int      orientation;
extern Boolean  ps_ruler,ps_header,resize,ps_curve,ps_resno;
extern int      first_printres,last_printres,blocklen;
extern int      firstres,lastres;



void set_go_penalty(TexT t)
{
	char str[10];
	float temp;

	GetTitle(t,str,10);
	temp = atof(str);
	if (temp < 0 || temp > 100)
		return;
	gap_open=temp;

        if(dnaflag)
         	dna_gap_open     = gap_open;
	else
         	prot_gap_open     = gap_open;
}

void set_ge_penalty(TexT t)
{
	char str[10];
	float temp;

	GetTitle(t,str,10);
	temp = atof(str);
	if (temp < 0 || temp > 100)
		return;
	gap_extend=temp;

        if(dnaflag)
         	dna_gap_extend     = gap_extend;
	else
         	prot_gap_extend     = gap_extend;
}

void set_gap_dist(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 100)
		return;
	gap_dist = temp;

}

void set_ntrials(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 10000)
                return;
        boot_ntrials = temp;
}

void set_ran_seed(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
        temp = atoi(str);
        if (temp < 0 || temp > 1000)
                return;
        boot_ran_seed = temp;
}

void set_div_seq(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 100)
		return;
	divergence_cutoff = temp;
}

void set_pw_go_penalty(TexT t)
{
        char str[10];
        float temp;
 
        GetTitle(t,str,10);
        temp = atof(str);
        if (temp < 0 || temp > 100)
                return;

	pw_go_penalty = temp;
        if(dnaflag)
         	dna_pw_go_penalty     = pw_go_penalty;
	else
         	prot_pw_go_penalty     = pw_go_penalty;
}

void set_pw_ge_penalty(TexT t)
{
        char str[10];
        float temp;
 
        GetTitle(t,str,10);
        temp = atof(str);
        if (temp < 0 || temp > 100)
                return;

	pw_ge_penalty = temp;
        if(dnaflag)
         	dna_pw_ge_penalty     = pw_ge_penalty;
	else
         	prot_pw_ge_penalty     = pw_ge_penalty;
}

void set_gp(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
        temp = atoi(str);
        if (temp < 0 || temp > 100)
                return;

	wind_gap = temp;
        if(dnaflag)
         	dna_wind_gap       = wind_gap;
	else
         	prot_wind_gap       = wind_gap;
}

void set_ktuple(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
        temp = atoi(str);
        if (temp < 0 || temp > 100)
                return;

	ktup = temp;
        if(dnaflag)
         	dna_ktup       = ktup;
	else
         	prot_ktup       = ktup;
}

void set_topdiags(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
        temp = atoi(str);
        if (temp < 0 || temp > 100)
                return;

	signif = temp;
        if(dnaflag)
         	dna_signif     = signif;
	else
         	prot_signif     = signif;
}

void set_window(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
        temp = atoi(str);
        if (temp < 0 || temp > 100)
                return;

	window = temp;
        if(dnaflag)
         	dna_window     = window;
	else
         	prot_window     = window;
}

void set_hyd_res(TexT t)
{
	int i,j;
	char tstr[27];

	GetTitle(t,tstr,27);
	for (i=0,j=0;i<strlen(hyd_residues) && i<27;i++)
	{
		if (isalpha(tstr[i]))
			hyd_residues[j++] = tstr[i];
	}
	hyd_residues[j]='\0';
}

void set_button(ButtoN l,Boolean *value)
{
	int tmp;

	tmp = GetStatus(l);
	if (tmp == TRUE)
		*value = TRUE;
	else
		*value = FALSE;
}

void set_toggle(PopuP l,Boolean *value)
{
	int tmp;

	tmp = GetValue(l);
	if (tmp == 1)
		*value = TRUE;
	else
		*value = FALSE;
}

void set_pref_penalties(PopuP l)
{
	set_toggle(l,&no_pref_penalties);
}

void set_hyd_penalties(PopuP l)
{
	set_toggle(l,&no_hyd_penalties);
}
void set_var_penalties(PopuP l)
{
	set_toggle(l,&no_var_penalties);
}
void set_endgaps(PopuP l)
{
	set_toggle(l,&use_endgaps);
}
void set_align_endgappenalties(PopuP l)
{
	set_toggle(l,&align_endgappenalties);
}
void set_realign_endgappenalties(PopuP l)
{
	set_toggle(l,&realign_endgappenalties);
}
void set_case(PopuP l)
{
	set_toggle(l,&lowercase);
}
void set_snos(PopuP l)
{
	set_toggle(l,&cl_seq_numbers);
}
void set_save_paras(PopuP l)
{
	set_toggle(l,&save_parameters);
}
void set_transitions(TexT t)
{
        char str[10];
        float temp;
 
        GetTitle(t,str,10);
        temp = atof(str);
        if (temp < 0 || temp > 100)
                return;
 
        transition_weight = temp;
}

void set_ambiguities(PopuP l)
{
	set_toggle(l,&use_ambiguities);
}

void set_neg_matrix(PopuP l)
{
	set_toggle(l,&neg_matrix);
}

void set_output_nbrf(ButtoN l)
{
	set_button(l,&output_nbrf);
}
void set_output_phylip(ButtoN l)
{
	set_button(l,&output_phylip);
}
void set_output_gcg(ButtoN l)
{
	set_button(l,&output_gcg);
}

void set_output_order(PopuP g)
{
	int tmp;
	tmp = GetValue(g);
	if (tmp == 1)
		output_order=INPUT;
	else
		output_order=ALIGNED;
}

void set_pagesize(PopuP g)
{
	int tmp;
	char tstr[10];

	tmp = GetValue(g);
	if (tmp == 1)
		pagesize=A4;
	else if (tmp == 2)
		pagesize=A3;
	else
		pagesize=USLETTER;
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
	sprintf(tstr,"%d",blocklen);
	SetTitle(blocklentext,tstr);
}
void set_orientation(PopuP g)
{
	int tmp;
	char tstr[10];

	tmp = GetValue(g);
	if (tmp == 1)
		orientation=LANDSCAPE;
	else
		orientation=PORTRAIT;

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
	sprintf(tstr,"%d",blocklen);
	SetTitle(blocklentext,tstr);
}
void set_resno(PopuP l)
{
	set_toggle(l,&ps_resno);
}
void set_curve(PopuP l)
{
	set_toggle(l,&ps_curve);
}
void set_ruler(PopuP l)
{
	set_toggle(l,&ps_ruler);
}
void set_header(PopuP l)
{
	set_toggle(l,&ps_header);
}
void set_resize(PopuP l)
{
	set_toggle(l,&resize);
}
void set_fres(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 100000)
                return;
        firstres = temp;
}
void set_lres(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 100000)
                return;
        lastres = temp;
}
void set_fpres(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 10000)
                return;
        first_printres = temp;
}
void set_lpres(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 10000)
                return;
        last_printres = temp;
}
void set_blocklen(TexT t)
{
        char str[10];
        int temp;
 
        GetTitle(t,str,10);
	if (str == NULL) return;
	temp = atoi(str);
        if (temp < 0 || temp > 10000)
                return;
        blocklen = temp;
}
void set_output_tree_nexus(ButtoN l)
{
	set_button(l,&output_tree_nexus);
}
void set_output_tree_clustal(ButtoN l)
{
	set_button(l,&output_tree_clustal);
}
void set_output_tree_phylip(ButtoN l)
{
	set_button(l,&output_tree_phylip);
}
void set_output_tree_distances(ButtoN l)
{
	set_button(l,&output_tree_distances);
}
void set_tossgaps(IteM i)
{
	tossgaps=GetStatus(i);
}
void set_kimura(IteM i)
{
	kimura=GetStatus(i);
}
void set_boot_format(PopuP g)
{
	int tmp;
	tmp = GetValue(g);
	if (tmp == 1)
		bootstrap_format=BS_NODE_LABELS;
	else
		bootstrap_format=BS_BRANCH_LABELS;
}

char prompt_for_yes_no(char *title,char *prompt)
{
        char lin2[MESSLENGTH*MESSLINES];
 
	if(!x_menus) return;

	strcpy(lin2,title);
	strcat(lin2,".\n");
	strcat(lin2,prompt);
	strcat(lin2,"?");
	if (Message(MSG_YN,lin2)==ANS_NO)
                return('n');
        else
                return('y');
 
}
 

/*
*	fatal()
*
*	Prints error msg and exits.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void fatal( char *msg,...)
{
	va_list ap;
	char istr[MESSLENGTH*MESSLINES] = "FATAL ERROR: ";
	char vstr[1000];

	
	va_start(ap,msg);
	vsprintf(vstr,msg,ap);
	va_end(ap);
	strncat(istr,vstr,MESSLENGTH*MESSLINES-20);
	Message(MSG_FATAL,istr);
}

/*
*	error()
*
*	Prints error msg.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void error( char *msg,...)
{
	va_list ap;
	char istr[MESSLENGTH*MESSLINES] = "ERROR: ";
	char vstr[1000];

	
	va_start(ap,msg);
	vsprintf(vstr,msg,ap);
	va_end(ap);
	strncat(istr,vstr,MESSLENGTH*MESSLINES-10);
	if (!interactive)
		fprintf(stdout,"%s",istr);
	else
		Message(MSG_ERROR,istr);
}

/*
*	warning()
*
*	Prints warning msg.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void warning( char *msg,...)
{
	va_list ap;
	char istr[MESSLENGTH*MESSLINES] = "WARNING: ";
	char vstr[1000];

	
	va_start(ap,msg);
	vsprintf(vstr,msg,ap);
	va_end(ap);
	strncat(istr,vstr,MESSLENGTH*MESSLINES-10);
	if (!interactive)
		fprintf(stdout,"%s",istr);
	else
		Message(MSG_ERROR,istr);
}

/*
*	info()
*
*	Prints info msg.
*	Variadic parameter list can be passed.
*
*	Return values:
*		none
*/

void info( char *msg,...)
{
	va_list ap;
	char istr[MESSLENGTH+10] = "";
	char vstr[1000];
    
	if (!mess_output) return;

	va_start(ap,msg);
	vsprintf(vstr,msg,ap);
	va_end(ap);
	strncat(istr,vstr,MESSLENGTH);
	if (!interactive)
		fprintf(stdout,"%s\n",istr);
	else
	{
		UseWindow(mainw);
		SelectFont(systemFont);	
		SetTitle(message,istr);
		if(save_log && save_log_fd!=NULL)
			fprintf(save_log_fd,"%s\n",istr);
		Update();
	}

}


void set_helix_gp(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	helix_penalty = temp;

}

void set_strand_gp(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	strand_penalty = temp;

}

void set_loop_gp(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	loop_penalty = temp;

}

void set_terminal_gp(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	helix_end_penalty = temp;

}

void set_helix_minus(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	helix_end_minus = temp;

}

void set_helix_plus(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	helix_end_plus = temp;
}

void set_strand_plus(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	strand_end_plus = temp;
}

void set_strand_minus(TexT t)
{
	char str[10];
	int temp;

	GetTitle(t,str,10);
	temp = atoi(str);
	if (temp < 0 || temp > 9)
		return;
	strand_end_minus = temp;
}

 
void set_inverted(IteM i)
{
        inverted=GetStatus(i);
        if (inverted==FALSE)
        {
                strcpy(color_lut[0].name,"BLACK");
                color_lut[0].r=0.4;
                color_lut[0].g=0.4;
                color_lut[0].b=0.4;
                SelectColor(color_lut[0].r*255, color_lut[0].g*255, color_lut[0].b*255);
                color_lut[0].val=GetColor();
        }
        else
        {
                strcpy(color_lut[0].name,"WHITE");
                color_lut[0].r=1.0;
                color_lut[0].g=1.0;
                color_lut[0].b=1.0;
                SelectColor(color_lut[0].r*255, color_lut[0].g*255, color_lut[0].b*255);
                color_lut[0].val=GetColor();
        }
 
        if(aln_mode==MULTIPLEM)
                DrawPanel(seq_panel.seqs);
        else
        {
                DrawPanel(prf_panel[0].seqs);
                DrawPanel(prf_panel[1].seqs);
        }
 
}

void set_ss_output(ButtoN b)
{
        int tmp;
 
        tmp = GetStatus(b);
        if (tmp) output_ss = TRUE;
	else output_ss = FALSE;

	if (output_ss && output_gp)
		output_struct_penalties=2;
	else if (output_ss)
		output_struct_penalties=0;
	else if (output_gp)
		output_struct_penalties=1;
	else
		output_struct_penalties=3;
}
void set_gp_output(ButtoN b)
{
        int tmp;
 
        tmp = GetStatus(b);
        if (tmp) output_gp = TRUE;
	else output_gp = FALSE;

	if (output_ss && output_gp)
		output_struct_penalties=2;
	else if (output_ss)
		output_struct_penalties=0;
	else if (output_gp)
		output_struct_penalties=1;
	else
		output_struct_penalties=3;
}

void set_user_matrix(ButtoN but)
{
	if(get_series_matrixname(usermtrxname,usermat,aa_xref,5,&matnum,mattext))
		strcpy(mtrxname,usermtrxname);
	SetValue(matrix_list,matnum);
}


void set_pw_user_matrix(ButtoN but)
{
	if(get_user_matrixname(pw_usermtrxname,pw_usermat,pw_aa_xref,5,&pw_matnum,pwmattext))
		strcpy(pw_mtrxname,pw_usermtrxname);
	SetValue(pw_matrix_list,pw_matnum);
}

void set_pw_matrix(GrouP g)
{
        int tmp;
 
        tmp = GetValue(g);
        if (tmp>0 && tmp<pw_matrix_menu.noptions) 
        {
                pw_matnum = tmp;
		strcpy(pw_mtrxname,pw_matrix_menu.opt[tmp-1].string);
        }
        else if(pw_usermtrxname[0]=='\0')
	{
		if(get_user_matrixname(pw_usermtrxname,pw_usermat,pw_aa_xref,pw_matrix_menu.noptions,&pw_matnum,pwmattext))
			strcpy(pw_mtrxname,pw_usermtrxname);
	}
	else
		pw_matnum=pw_matrix_menu.noptions;
	SetValue(pw_matrix_list,pw_matnum);
}
void set_matrix(GrouP g)
{
	int tmp;
	int status;
 
        tmp = GetValue(g);
        if (tmp>0 && tmp<matrix_menu.noptions) 
        {
        	matnum = tmp;
		strcpy(mtrxname,matrix_menu.opt[tmp-1].string);
        }
        else if(usermtrxname[0]=='\0')
	{
		if(get_series_matrixname(usermtrxname,usermat,aa_xref,matrix_menu.noptions,&matnum,mattext))
			strcpy(mtrxname,usermtrxname);
	}
	else matnum=matrix_menu.noptions;

	SetValue(matrix_list,matnum);
}

static int get_series_matrixname(char *usermtrxname, short *usermat,short *aa_xref,int usermatnum,int *matnum,PrompT mattext)
{
	int ret=0;
	static Char filename[FILENAMELEN];

        if (GetInputFileName(filename,FILENAMELEN,"","")) 
	{
        	if(user_mat_series(filename, usermat, aa_xref))
        	{
                	strcpy(usermtrxname,filename);
                	*matnum=usermatnum;
			SetTitle(mattext,usermtrxname);
			ret=1;
        	}
	} 

	return ret;
}
int get_user_matrixname(char *usermtrxname, short *usermat,short *aa_xref,int usermatnum,int *matnum,PrompT mattext)
{
	int ret=0;
	static Char filename[FILENAMELEN];

        if (GetInputFileName(filename,FILENAMELEN,"","")) 
	{
        	if(user_mat(filename, usermat, aa_xref))
        	{
                	strcpy(usermtrxname,filename);
                	*matnum=usermatnum;
			SetTitle(mattext,usermtrxname);
			ret=1;
        	}
	} 

	return ret;
}

void set_user_dnamatrix(ButtoN but)
{
	if(get_user_matrixname(dnausermtrxname,userdnamat,dna_xref,3,&dnamatnum,dnamattext))
		strcpy(dnamtrxname,dnausermtrxname);
	SetValue(dnamatrix_list,dnamatnum);
}


void set_pw_user_dnamatrix(ButtoN but)
{
	if(get_user_matrixname(pw_dnausermtrxname,pw_userdnamat,pw_dna_xref,3,&pw_dnamatnum,pwdnamattext))
		strcpy(pw_dnamtrxname,pw_dnausermtrxname);
	SetValue(pw_dnamatrix_list,pw_dnamatnum);
}

void set_pw_dnamatrix(GrouP g)
{
        int tmp;
 
        tmp = GetValue(g);
        if (tmp>0 && tmp<dnamatrix_menu.noptions) 
        {
                pw_dnamatnum = tmp;
		strcpy(pw_dnamtrxname,dnamatrix_menu.opt[tmp-1].string);
        }
        else if(pw_dnausermtrxname[0]=='\0')
	{
		if(get_user_matrixname(pw_dnausermtrxname,pw_userdnamat,pw_dna_xref,dnamatrix_menu.noptions,&pw_dnamatnum,pwdnamattext))
			strcpy(pw_dnamtrxname,pw_dnausermtrxname);
	}
	else pw_dnamatnum=dnamatrix_menu.noptions;
	SetValue(pw_dnamatrix_list,pw_dnamatnum);
}
void set_dnamatrix(GrouP g)
{
	int tmp;
 
        tmp = GetValue(g);
        if (tmp>0 && tmp<dnamatrix_menu.noptions) 
        {
        	dnamatnum = tmp;
		strcpy(dnamtrxname,dnamatrix_menu.opt[tmp-1].string);
        }
        else if(dnausermtrxname[0]=='\0')
	{
		if(get_user_matrixname(dnausermtrxname,userdnamat,dna_xref,dnamatrix_menu.noptions,&dnamatnum,dnamattext))
			strcpy(dnamtrxname,dnausermtrxname);
	}
	else dnamatnum=dnamatrix_menu.noptions;
	SetValue(dnamatrix_list,dnamatnum);
}

FILE *  open_input_file(char *file_name)
{
        FILE * file_handle;
 
        if (*file_name == EOS) {
                error("Bad input file [%s]",file_name);
                return NULL;
        }
#ifdef VMS
        if((file_handle=fopen(file_name,"r","rat=cr","rfm=var"))==NULL) {
#else
        if((file_handle=fopen(file_name,"r"))==NULL) {
#endif
                error("Cannot open input file [%s]",file_name);
                return NULL;
        }
        return file_handle;
}

 
void set_use_ss1(PopuP l)
{
        set_toggle(l,&use_ss1);
	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
}
void set_use_ss2(PopuP l)
{
        set_toggle(l,&use_ss2);
	load_aln(prf_panel[0],0,profile1_nseqs-1,FALSE);
	load_aln(prf_panel[1],profile1_nseqs,nseqs-1,FALSE);
}


void set_output_clustal(ButtoN l)
{
        set_button(l,&output_clustal);
}
void set_output_gde(ButtoN l)
{
        set_button(l,&output_gde);
}
void set_output_nexus(ButtoN l)
{
        set_button(l,&output_nexus);
}

void set_format(GrouP g)
{
        int i;
	char path[FILENAMELEN];

	get_path(seqname,path);
	GetTitle(savealntext, path,FILENAMELEN);
/* remove the current extension */
	for(i=strlen(path)-1;i>=0;i--)
		if(path[i]=='.')
		{
			path[i]='\0';
			break;
		}
			
        i = GetValue(g);
        if (i==1)
	{
                save_format=CLUSTAL;
		strcat(path,".aln");
	}
        else if (i==2)
	{
                save_format=PIR;
		strcat(path,".pir");
	}
        else if (i==3)
	{
                save_format=MSF;
		strcat(path,".msf");
	}
        else if (i==4)
	{
                save_format=PHYLIP;
		strcat(path,".phy");
	}
        else if (i==5)
	{
                save_format=GDE;
		strcat(path,".gde");
	}
        else if (i==6)
	{
                save_format=NEXUS;
		strcat(path,".nxs");
	}

	SetTitle(savealntext, path);
}

void set_residue_exceptions(IteM i)
{
	if (residue_exceptions==FALSE)
		residue_exceptions=TRUE;
	else
		residue_exceptions=FALSE;
	if (aln_mode==MULTIPLEM)
		DrawPanel(seq_panel.seqs);
	else
	{
		DrawPanel(prf_panel[0].seqs);
		DrawPanel(prf_panel[1].seqs);
	}
}
 
 
void set_fs_toggle(PopuP l)
{
        set_toggle(l,&quick_pairalign);
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
}

void set_font_size(PopuP g)
{
        int tmp;
	char font[30];

        tmp = GetValue(g);
	font_size=tmp-1;
	sprintf(font, "%s,%d,%c", "courier", av_font[font_size], 'm');
        datafont=ParseFont(font);

	if (aln_mode==MULTIPLEM)
	{
		DrawPanel(seq_panel.names);
		DrawPanel(seq_panel.seqs);
	}
	else
	{
		DrawPanel(prf_panel[0].names);
		DrawPanel(prf_panel[0].seqs);

		DrawPanel(prf_panel[1].names);
		DrawPanel(prf_panel[1].seqs);
	}
	correct_name_bars(FALSE);
	correct_seq_bars(FALSE);
}
 
void set_pscroll_mode(ButtoN l)
{
	panel_data data;

        set_button(l,&fixed_prf_scroll);
	GetPanelExtra(prf_panel[0].seqs,&data);
	if(fixed_prf_scroll)
		data.lockoffset=data.firstvcol;
	else
		data.lockoffset=0;
	SetPanelExtra(prf_panel[0].seqs,&data);
	GetPanelExtra(prf_panel[1].seqs,&data);
	if(fixed_prf_scroll)
		data.lockoffset=data.firstvcol;
	else
		data.lockoffset=0;
	SetPanelExtra(prf_panel[1].seqs,&data);
	correct_seq_bars(FALSE);
}
 
void set_aln_mode(PopuP g)
{
        int tmp;
        tmp = GetValue(g);
        if (tmp == 1)
                aln_mode = MULTIPLEM;
        else
                aln_mode = PROFILEM;
        switch_mode();
}
 
void set_show_segments(IteM l)
{
        if (segment_exceptions==FALSE)
                segment_exceptions=TRUE;
        else
                segment_exceptions=FALSE;
	calc_seg_exceptions();
	SetStatus(segment_item,segment_exceptions);
        show_segment_exceptions();
}
 

void shift(Handle a, int dx, int dy)
{
	PoinT pt;
 
	GetNextPosition (a, &pt);
	pt.x+=dx;
	pt.y+=dy;
	SetNextPosition(a, pt);
}

void stripspace(char *str)
{
	register int i,j,p;
	char *tstr;
 
#ifndef UNIX
	return;
#endif
        p = strlen(str) - 1;
 
        while ( isspace(str[p]) )
                p--;
 
        str[p + 1] = EOS;

	tstr=(char *)ckalloc((p+2)*sizeof(char));
 
	for(i=0,j=0;i<=p;i++)
		if(!isspace(str[i]))
			tstr[j++]=str[i];
	tstr[j] = EOS;
	strcpy(str,tstr);
	ckfree(tstr);

}

