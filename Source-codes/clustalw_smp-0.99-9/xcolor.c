#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include <vibrant.h>

#include "clustalw.h"
#include "xmenu.h"

#define SIMPLE 1
#define COMPOUND 2

#define LEFTMARGIN 20
#define SEPARATION 2
#define CHARHEIGHT 10
#define CHARWIDTH 6
#define A4X 564
#define A4Y 800
#define A3X 832
#define A3Y 1159
#define USLETTERX 564
#define USLETTERY 750
#define SCOREY 3
#define HEADER 7
#define NOHEADER 0
#define MAXRESNO 6

#define MAXPARLEN 10
#define MAXPAR 100
 
static void print_ps_info(FILE *fd,int pagesize);
static void print_page_header(FILE *fd,int ps_rotation,int maxx,int maxy,
int page,int numpages,Boolean header,char *str_time,
char *ps_file,int ps_xtrans,int ps_ytrans,float ps_scale);
static void print_header_line(FILE *fd,panel_data name_data, panel_data seq_data,
int ix,int fr,int lr);
static void print_footer_line(FILE *fd,panel_data name_data, panel_data seq_data,
int ix,int fr,int lr);
static void print_quality_curve(FILE *fd,panel_data seq_data
,int fr,int lr,int score_height);
static void print_seq_line(FILE *fd,panel_data name_data, panel_data seq_data,
int row,int seq,int fr,int lr,int res_number);


typedef struct consensus_parameters
{
char consensus;
int cutoff;
int length;
char cutoff_list[20];
} consensus_para;
 
typedef struct color_parameters
{
int type;
char residue;
int color;
int length;
char cons_list[20];
} color_para;

static void init_color_lut(FILE *fd);
static int init_printer_lut(char *filename);
static char *init_consensus(panel_data data);
static int SaveColPara(char word[MAXPAR][MAXPARLEN],int num_words,int count);
static int SaveConPara(char word[MAXPAR][MAXPARLEN],int num_words,int count);
static int get_line(char *sinline,char word[MAXPAR][MAXPARLEN]);
static int residue_color(char res,char consensus);
static Boolean commentline(char *line);

#define DEF_NCOLORS 4
#define MAX_NCOLORS 8
#define DEFAULT_COLOR 0

typedef struct rgb_color {
	char name[20];
	float r,g,b;
} rgb_color;

rgb_color def_color_lut[MAX_NCOLORS]={
	"RED"          ,0.9, 0.1, 0.1,
	"BLUE"         ,0.1, 0.1, 0.7,
	"GREEN"        ,0.1, 0.9, 0.1,
	"ORANGE"       ,0.9, 0.6, 0.3,
	"CYAN"         ,0.1, 0.9, 0.9,
	"PINK"         ,0.9, 0.5, 0.5,
	"MAGENTA"      ,0.9, 0.1, 0.9,
	"YELLOW"       ,0.9, 0.9, 0.0,
};

char def_aacolor[MAX_NCOLORS][26]={"krh",
				"fwy",
				"ilmv",
				"gpst"};

char def_dnacolor[MAX_NCOLORS][26]={"a",
				"c",
				"tu",
				"g"};

extern char revision_level[];

extern int max_names;
 
extern int ncolors;
extern int ncolor_pars;
extern color color_lut[];
extern int inverted;
extern Boolean residue_exceptions;
extern Boolean segment_exceptions;
extern Boolean dnaflag;

int NumColParas;
int NumConParas;

color_para Col_Par[100];
consensus_para Con_Par[100];


void make_colormask(panel_data data)
{
	int i,j;

	for(i=0;i<data.nseqs;i++)
		for(j=0;j<data.ncols;j++)
			data.colormask[i][j] = DEFAULT_COLOR;

	if (ncolors > 1)
	{
        	data.consensus=init_consensus(data);

		for(i=0;i<data.nseqs;i++)
			for(j=0;j<data.ncols;j++)
				data.colormask[i][j] = residue_color(data.lines[i][j],data.consensus[j]);

	}
}

static void init_color_lut(FILE *fd)
{ 
	char sinline[1025];
	char *args[10];
	int i,numargs;
	Boolean found=FALSE;

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

	ncolors=1;
	if (fd != NULL)
	{
		for (;fgets(sinline,1024,fd)!=NULL;)
		{
			sinline[strlen(sinline)-1] = '\0';
			if (strcmp(sinline,"@rgbindex")==0) 
			{
				found = TRUE;
				break;
			}
		}
	}
	if (found == TRUE)
	{
		for (;fgets(sinline,1024,fd)!=NULL;)
		{
			if (commentline(sinline)) continue;
			if (sinline[0]=='@') break;
			numargs = getargs(sinline, args, 4);
			if (numargs != 4)
			{
				error("Problem in color rgb index - line %d\n",ncolors+1);
				break;
			}
			else
			{
				strcpy(color_lut[ncolors].name, args[0]);
				color_lut[ncolors].r=atof(args[1]);
				color_lut[ncolors].g=atof(args[2]);
				color_lut[ncolors].b=atof(args[3]);
				SelectColor(color_lut[ncolors].r*255, color_lut[ncolors].g*255, color_lut[ncolors].b*255);
				color_lut[ncolors].val=GetColor();
				ncolors++;
				if (ncolors>=MAXCOLORS) 
				{
					warning("Only using first %d colors in rgb index.",MAXCOLORS);
					break;
				}
			}
		}

	}

/* if we can't find a table, use the hard-coded colors */
        if (ncolors==1)
        {
		ncolors=MAX_NCOLORS+1;
		for(i=1;i<ncolors;i++)
		{
			strcpy(color_lut[i].name,def_color_lut[i-1].name);
			color_lut[i].r=def_color_lut[i-1].r;
			color_lut[i].g=def_color_lut[i-1].g;
			color_lut[i].b=def_color_lut[i-1].b;
			SelectColor(color_lut[i].r*255, color_lut[i].g*255, color_lut[i].b*255);
			color_lut[i].val=GetColor();
		}
	}

}

void init_color_parameters(char *par_file)
{

	int i,j,err;
	char sinline[1025];
	int maxparas = 50;
	char inword[MAXPAR][MAXPARLEN];
	int num_words;
	int in_consensus=FALSE,in_color=FALSE;
	int consensus_found=FALSE,color_found=FALSE;
	FILE *par_fd=NULL;

	if(par_file!=NULL)
		par_fd=fopen(par_file,"r");
	if(par_fd==NULL)
	{
                info("No color file found - using defaults");
		ncolor_pars=0;
	}

	init_color_lut(par_fd);
	if (par_fd != NULL) rewind(par_fd);
	if (ncolors==0) return;

	NumColParas=0;
	NumConParas=0;
	if (par_fd != NULL)
	{
		for(;fgets(sinline,1024,par_fd) != NULL;)
		{
			sinline[strlen(sinline)-1] = '\0';
			if (commentline(sinline)) continue;
			switch(sinline[0])
			{	
				case '\0':
					break;
				case '@':
					if (strcmp((char*)(sinline+1),"consensus")==0) 
					{
						in_consensus = TRUE;
						in_color = FALSE;
						consensus_found = TRUE;
					}
					else if (strcmp((char*)(sinline+1),"color")==0)
					{
						in_consensus = FALSE;
						in_color = TRUE;
						color_found = TRUE;
					}
					break;
				default:
					num_words = get_line(sinline,inword);
					if (in_consensus == TRUE) 
					{
						err = SaveConPara(inword,num_words,NumConParas);
						if (err == 0) NumConParas++;
					}
					else if (in_color == TRUE)
					{
						err = SaveColPara(inword,num_words,NumColParas);
						if (err == 0) NumColParas++;
					}
	
					if((NumColParas>maxparas) || (NumConParas>maxparas))
				     	error("Too many parameters in color file");
	
			}
		}
		if (color_found == FALSE)
		{
			error("@color not found in parameter file - using defaults\n");
			ncolor_pars=0;
		}
		fclose(par_fd);
	}
	ncolor_pars=NumColParas;

/* if no color parameters found, use the default aa groupings */
	if(ncolor_pars==0)
	{
		if (dnaflag)
		{
			for(i=0;i<DEF_NCOLORS;i++)
			{
				for(j=0;j<strlen(def_dnacolor[i]);j++)
				{
					Col_Par[ncolor_pars].type=SIMPLE;
					Col_Par[ncolor_pars].residue=def_dnacolor[i][j];
					Col_Par[ncolor_pars].color=i+1;
					ncolor_pars++;
				}
			}
		}
		else
		{
			for(i=0;i<DEF_NCOLORS;i++)
			{
				for(j=0;j<strlen(def_aacolor[i]);j++)
				{
					Col_Par[ncolor_pars].type=SIMPLE;
					Col_Par[ncolor_pars].residue=def_aacolor[i][j];
					Col_Par[ncolor_pars].color=i+1;
					ncolor_pars++;
				}
			}
		}
	}
	NumColParas=ncolor_pars;
}

char *find_file(char *def_file)
{
	char filename[FILENAMELEN];
	char *retname;
	FILE *fd;
	Boolean found=FALSE;
#ifdef UNIX
        char *path, *path1, *deb, *fin;
        sint lf, ltot;
        char *home;
#endif


        strcpy(filename,def_file);
        fd = fopen(filename,"r");
	if (fd != NULL)
		found=TRUE;
#ifdef UNIX
        if (found == FALSE)
        {
                home = getenv("HOME");
		if (home != NULL)
		{
                	sprintf(filename,"%s/%s",home,def_file);
                	fd = fopen(filename,"r");
			if (fd != NULL)
				found=TRUE;
		}
                if (found == FALSE)
                {
			path=getenv("PATH");/* get the list of path directories,
                        			separated by : */
			/* added for File System Standards  - Francois */
			path1=(char *)ckalloc((strlen(path)+64)*sizeof(char));
			strcpy(path1,path);
			strcat(path1,"/usr/share/clustalx:/usr/local/share/clustalx"); 

        		lf=(sint)strlen(def_file);
        		deb=path1;
        		do
                	{
                		fin=strchr(deb,':');
                		if(fin!=NULL)
                        	{
					strncpy(filename,deb,fin-deb);
					ltot=fin-deb;
				}
                		else
                        	{
					strcpy(filename,deb);
					ltot=(sint)strlen(filename);
				}
                		/* now one directory is in filename */
                		if( ltot + lf + 1 <= FILENAMELEN)
                        	{
                        		filename[ltot]='/';
                        		strcpy(filename+ltot+1,def_file); /* now dir is appended with filename */
                        		if( (fd = fopen(filename,"r") ) != NULL)
					{
						found=TRUE;
						break;
                        		}
                        	}
                		else found = FALSE;
                		deb=fin+1;
                	}
        		while (fin != NULL);
                }
        }
#endif
	if (found == TRUE)
	{
		fclose(fd);
		retname=(char *)ckalloc((strlen(filename)+1)*sizeof(char));
		strcpy(retname,filename);
	}
	else
		retname=NULL;
	return(retname);
}
 
static char *init_consensus(panel_data data)
{
	char *cons_data;
        int num_res,seq,res,par,cons_total,i;
        char residue;
 
	cons_data=(char *)ckalloc((data.ncols+1)*sizeof(char));

        for (res=0;res<data.ncols;res++)
        {
                cons_data[res] = '.';
                for (par=0;par<NumConParas;par++)
                {
                        cons_total = num_res = 0;
                        for (seq=0;seq<data.nseqs;seq++)
                        {
				residue=tolower(data.lines[seq][res]);
				if (isalpha(residue))
					num_res++;
                                for (i=0;i<Con_Par[par].length;i++)
                                        if (residue==tolower(Con_Par[par].cutoff_list[i]))
                                                cons_total++;
                        }
                        if (num_res != 0)
                                if (((cons_total*100)/num_res) >= Con_Par[par].cutoff)
                                        cons_data[res] = Con_Par[par].consensus;
                }
        }

	return(cons_data);
}

static int SaveColPara(char word[MAXPAR][MAXPARLEN],int num_words,int count)
{

	int i;

	if (num_words < 3)
	{
		error("Wrong format in color list");
		return(1);
	}

	if (word[1][0] != '=')
	{
		error("Wrong format in color list");
		return(2);
	}

	if (num_words == 3)
	{
		Col_Par[count].type = SIMPLE;
		Col_Par[count].residue = word[0][0];
		Col_Par[count].color = -1;
		for (i=0;i<ncolors;i++)
			if (strcmp(word[2],color_lut[i].name)==0) Col_Par[count].color = i;
		if (Col_Par[count].color == -1)
		{
			error("%s not found in rgb index - using %s",word[2],color_lut[0].name);
			Col_Par[count].color = 0;
		}
	}
	else
	{
		if (strcmp(word[3],"if")==0)
		{
			Col_Par[count].type = COMPOUND;
			Col_Par[count].residue = word[0][0];
			Col_Par[count].color = -1;
			for (i=0;i<ncolors;i++)
				if (strcmp(word[2],color_lut[i].name)==0) Col_Par[count].color = i;
			if (Col_Par[count].color == -1)
			{
				error("%s not found in rgb index - using %s",word[2],color_lut[0].name);
				Col_Par[count].color = 0;
			}
			Col_Par[count].length = num_words - 4;
			for (i=4;i<num_words;i++)
				Col_Par[count].cons_list[i-4] = word[i][0];
		}
		else
		{
			error("Wrong format in color list");
			return(3);
		}
	}

	return(0);
		
}


static int SaveConPara(char word[MAXPAR][MAXPARLEN],int num_words,int count)
{

	int i;

	if (num_words < 3)
	{
		error("Wrong format in consensus list");
		return(1);
	}

	if (word[1][0] != '=')
	{
		error("Wrong format in consensus list");
		return(2);
	}

	Con_Par[count].consensus = word[0][0];
	for (i=0;i<MAXPARLEN-1;i++)
	{
		if(word[2][i]=='%') word[2][i] = '\0';
	}
	Con_Par[count].cutoff = atoi(word[2]);
	Con_Par[count].length = num_words - 3;
	for (i=3;i<num_words;i++)
	{
		Con_Par[count].cutoff_list[i-3] = word[i][0];
	}

	return(0);
		
}

static int get_line(char *sinline,char word[MAXPAR][MAXPARLEN])
{
	int i=0, j, word_count=0, char_count=0;
	int in_word=FALSE;

	for(i=0;i<MAXPAR-1;i++)
		for(j=0;j<MAXPARLEN-1;j++)
			word[i][j]='\0';

	for (i=0;i<=strlen(sinline);i++)
	{
		switch (sinline[i])
		{
			case ' ':
			case '\t':
			case '\0':
			case ':':
				if (in_word)
				{
					word[word_count][char_count] = '\0';
					word_count++;
					char_count = 0;
					in_word = FALSE;
				}
				break;
			default:
				in_word = TRUE;
				word[word_count][char_count] = sinline[i];
				char_count++;
				break;
		}		

	}
	return(word_count);
}

static int residue_color(char res,char consensus)
{
	int i,j;

        for (i=0;i<NumColParas;i++)
        {
                if (tolower(res) == tolower(Col_Par[i].residue))
                {
                        switch (Col_Par[i].type)
                        {
                        case SIMPLE:
                                return(Col_Par[i].color);
                        case COMPOUND:
                                for (j=0;j<Col_Par[i].length;j++)
                                {
                                        if (consensus == Col_Par[i].cons_list[j]
)
                                                return(Col_Par[i].color);
                                }
                                break;
                        default:
                                return(DEFAULT_COLOR);
                        }
                }
        }
        return(DEFAULT_COLOR);
}

static Boolean commentline(char *line)
{
        int i;
 
	if (line[0] == '#') return TRUE;
        for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
                if( !isspace(line[i]) )
                        return FALSE;
        }
        return TRUE;
}

int block_height,block_left,block_top;
int header_top,seq_top,footer_top,curve_top;

void write_ps_file(spanel p,char *ps_file,char *par_file,int pagesize,
int orientation,Boolean header, Boolean ruler, Boolean resno, Boolean resize,
int first_printres,int last_printres,
int blength,Boolean show_curve)
{
	int i,j,bn,seq,numseqs;
	int err;
	int blocklen,numpages;
	int fr,lr;
	int page,row;
	int ps_rotation=0,ps_xtrans=0,ps_ytrans=0;
	float ps_scale,hscale,wscale;
	int maxseq;
	int maxx=0,maxy=0;
	int score_height=0;
	int main_header=0;
	int numelines,numecols;
	int nhead,nfoot;
	int ppix_width;    /* width of the page in pixels */
	int pchar_height;    /* height of the page in chars for sequences */
	int ppix_height;    /* height of the page in pixels for sequences */
        int blocksperpage,numblocks;
	int *res_number;
	panel_data name_data,seq_data;
	FILE *fd;

	time_t *tptr=NULL,ttime;
	char *str_time;

/* open the output file */
	if ((fd=fopen(ps_file,"w"))==NULL)
	{
		error("Cannot open file %s",ps_file);
		return;
	}

/* check for printer-specific rgb values */
	err=init_printer_lut(par_file);
	if(err>0) warning("No PS Colors file: using default colors\n");

/* get the page size parameters */

	if (pagesize==A4)
	{
		if (orientation==PORTRAIT)
		{
			maxx=A4X;
			maxy=A4Y;
			ps_rotation=0;
		}
		else
		{
			maxx=A4Y;
			maxy=A4X;
			ps_rotation=-90;
		}
	}
	else if (pagesize==A3)
	{
		if (orientation==PORTRAIT)
		{
			maxx=A3X;
			maxy=A3Y;
			ps_rotation=0;
		}
		else
		{
			maxx=A3Y;
			maxy=A3X;
			ps_rotation=-90;
		}
	}
	else if (pagesize==USLETTER)
	{
		if (orientation==PORTRAIT)
		{
			maxx=USLETTERX;
			maxy=USLETTERY;
			ps_rotation=0;
		}
		else
		{
			maxx=USLETTERY;
			maxy=USLETTERX;
			ps_rotation=-90;
		}
	}
	if(show_curve) score_height=SCOREY;
	if(header) main_header=HEADER;
	else main_header=NOHEADER;
        ppix_width=maxx-LEFTMARGIN*2;
        ppix_height=maxy-main_header*CHARHEIGHT;

/* get the name data */
	GetPanelExtra(p.names,&name_data);

/* get the sequence data */
	GetPanelExtra(p.seqs,&seq_data);
	numseqs=seq_data.nseqs;
	nhead=seq_data.nhead;
	if(ruler)
		nfoot=seq_data.nfoot;
	else
		nfoot=seq_data.nfoot-1;
	numelines=nhead+nfoot+score_height+SEPARATION;

/* check the block length, residue range parameters */
	if(first_printres<=0)
		first_printres=1;
	if((last_printres<=0) || (last_printres>seq_data.ncols))
		last_printres=seq_data.ncols;
	if(first_printres>last_printres)
	{
		error("Bad residue range - cannot write postscript");
		return;
	}
	if (blength==0 || last_printres-first_printres+1<blength) 
		blocklen=last_printres-first_printres+1;
	else
		blocklen=blength;

	res_number=(int *)ckalloc((name_data.nseqs+1)*sizeof(int));
	for(i=0;i<name_data.nseqs;i++)
	{
		res_number[i]=0;
		for(j=0;j<first_printres-1;j++)
			if(isalpha(seq_data.lines[i][j])) res_number[i]++;
	}
	if(resno)
		numecols=MAXRESNO+1+max_names;
	else
		numecols=1+max_names;

/* print out the PS revision level etc. */
	ttime = time(tptr);
	str_time = ctime(&ttime);
	print_ps_info(fd,pagesize);

/* calculate scaling factors, block sizes to fit the page etc. */

        if (resize==FALSE || blocklen==last_printres-first_printres+1)
        {
/* split the alignment into blocks of sequences. If the blocks are too long
for the page - tough! */
		if(resize==FALSE)
                	ps_scale=1.0;
		else
			ps_scale=(float)ppix_width/(float)((blocklen+numecols)*CHARWIDTH);
		ps_xtrans= LEFTMARGIN * (1-ps_scale);
		ps_ytrans= ppix_height * (1-ps_scale);
		if (pagesize!=A3 && orientation==LANDSCAPE)
			ps_xtrans-=LEFTMARGIN;

        	pchar_height=((maxy/CHARHEIGHT)-main_header)/ps_scale;
        	maxseq=pchar_height-numelines;
		block_height = (maxseq+numelines) * CHARHEIGHT;
		numpages = (numseqs/maxseq) + 1;
		seq=0;
		for (page=0;page<numpages;page++)
		{
/* print the top of page header */
			print_page_header(fd,ps_rotation,maxx,maxy,
			   page,numpages,header,str_time,
			   ps_file,ps_xtrans,ps_ytrans,ps_scale);

			block_top = maxy - main_header*CHARHEIGHT;
			block_left = LEFTMARGIN + (1+max_names)*CHARWIDTH; 
			header_top = block_top;

			fr=first_printres-1;
			lr=last_printres-1;
/*  show the header lines */
			for (i=0;i<nhead;i++)
				print_header_line(fd,name_data,seq_data,i,fr,lr);

			seq_top = block_top-nhead*CHARHEIGHT;
/*  show the sequence lines */
			for (row=0;row<maxseq ;row++)
			{
				if(resno)
				{
					for(i=fr;i<=lr;i++)
						if(isalpha(seq_data.lines[seq][i]))
							res_number[seq]++;
				}
				print_seq_line(fd,name_data,seq_data,row,seq,fr,lr,res_number[seq]);
				seq++;
				if(seq>=numseqs)
				{
					row++;
					break;
				}
			}

			footer_top = seq_top-row*CHARHEIGHT;
/*  show the footer lines */
			for (i=0;i<nfoot;i++)
				print_footer_line(fd,name_data,seq_data,i,fr,lr);

			curve_top = footer_top-nfoot*CHARHEIGHT;
/* show the quality curve */
			if(show_curve)
				print_quality_curve(fd,seq_data,fr,lr,score_height);

			fprintf(fd,"\nshowpage\n");
			fprintf(fd,"restore\n");
		}
        }
        else
        {
/* split the alignment into blocks of residues, and scale the blocks to fit the page */
        	maxseq=ppix_height/CHARHEIGHT-numelines-main_header;
		hscale=(float)maxseq/(float)numseqs;
		wscale=(float)ppix_width/(float)((blocklen+numecols)*CHARWIDTH);
                ps_scale=MIN(hscale,wscale);
		ps_xtrans= LEFTMARGIN * (1-ps_scale);
		ps_ytrans= ppix_height * (1-ps_scale);
		if (pagesize!=A3 && orientation==LANDSCAPE)
			ps_xtrans-=LEFTMARGIN;

        	pchar_height=((maxy/CHARHEIGHT)-main_header)/ps_scale;
        	maxseq=pchar_height-numelines;
		block_height = (numseqs+numelines) * CHARHEIGHT;
		blocksperpage = pchar_height/(numseqs+numelines);
		if (blocksperpage==0)
		{
			error("illegal combination of print parameters");
			return;
		}
		numblocks = (last_printres-first_printres) / blocklen + 1;
        	if (numblocks % blocksperpage == 0)
			numpages = numblocks / blocksperpage;
        	else
			numpages = numblocks / blocksperpage + 1;

		for (bn=0;bn<numblocks;bn++)
		{
			page = bn / blocksperpage;
/* print the top of page header */
			if (bn % blocksperpage == 0)
				print_page_header(fd,ps_rotation,maxx,maxy,
			   	page,numpages,header,str_time,
			   	ps_file,ps_xtrans,ps_ytrans,ps_scale);

			block_top = maxy - main_header*CHARHEIGHT-block_height*(bn%blocksperpage);
			block_left = LEFTMARGIN + (1+max_names)*CHARWIDTH; 
			header_top = block_top;
			seq_top = block_top-nhead*CHARHEIGHT;
			footer_top = block_top-(nhead+numseqs)*CHARHEIGHT;
			curve_top = block_top-(nhead+numseqs+nfoot)*CHARHEIGHT;

			fr=first_printres-1 + blocklen*bn;
			lr=fr+blocklen-1;
			if(lr>=last_printres) lr=last_printres-1;
/*  show the header lines */
			for (i=0;i<nhead;i++)
				print_header_line(fd,name_data,seq_data,i,fr,lr);

/*  show the sequence lines */
			for (i=0;i<numseqs;i++)
			{
				row = i % maxseq;
				if(resno)
				{
					for(j=fr;j<=lr;j++)
						if(isalpha(seq_data.lines[i][j]))
							res_number[i]++;
				}
				print_seq_line(fd,name_data,seq_data,row,i,fr,lr,res_number[i]);
			}
/*  show the footer lines */
			for (i=0;i<nfoot;i++)
				print_footer_line(fd,name_data,seq_data,i,fr,lr);

/* show the quality curve */
			if(show_curve)
				print_quality_curve(fd,seq_data,fr,lr,score_height);

			if ((bn == (numblocks-1)) || ((bn % blocksperpage == blocksperpage-1)))
			{
				fprintf(fd,"\nshowpage\n");
				fprintf(fd,"restore\n");
			}
		}
	}
	fclose(fd);
	return;
}

static int init_printer_lut(char *filename)
{ 
	FILE *fd;
	char sinline[1025];
	char *args[10];
	char name[20];
	int i,numargs;
	Boolean found=FALSE;
	char *par_file=NULL;

/* reset the printer rgb colors to the color file rgb values */
	for(i=0;i<ncolors;i++)
	{
		color_lut[i].pr=color_lut[i].r;
		color_lut[i].pg=color_lut[i].g;
		color_lut[i].pb=color_lut[i].b;
	}

/* search for the printer color file */
	if(filename[0]==EOS) return 1;
	par_file=find_file(filename);
	if(par_file==NULL)
	{
		error("Cannot find printer file %s",filename);
		return 1;
	}
	if ((fd=fopen(par_file,"r"))==NULL)
	{
		error("Cannot open printer file %s",par_file);
		return 1;
	}

	for (;fgets(sinline,1024,fd)!=NULL;)
	{
		if (commentline(sinline)) continue;
		numargs = getargs(sinline, args, 4);
		if (numargs != 4)
		{
			error("Problem in parameter file - line %d\n",ncolors+1);
			break;
		}
		else
		{
/* we've found a color - find the index the color lut */
			strcpy(name, args[0]);
			for(i=0;i<ncolors;i++)
			{
				if(strcmp(name,color_lut[i].name)==0)
				{
					color_lut[i].pr=atof(args[1]);
					color_lut[i].pg=atof(args[2]);
					color_lut[i].pb=atof(args[3]);
				}
			}
		}
	}
	ckfree(par_file);
	return 0;
}

static void print_ps_info(FILE *fd,int pagesize)
{
	fprintf(fd,"%%!PS-Adobe-1.0\n");
	fprintf(fd,"%%%%Creator: Julie Thompson\n");
	fprintf(fd,"%%%%Title:ClustalX Alignment\n");
	fprintf(fd,"%%%%EndComments\n");
	fprintf(fd,"/box { newpath\n");
	fprintf(fd,"\t-0 -3 moveto\n");
	fprintf(fd,"\t-0 %d lineto\n",CHARHEIGHT-3);
	fprintf(fd,"\t%d %d lineto\n",CHARWIDTH,CHARHEIGHT-3);
	fprintf(fd,"\t%d -3 lineto\n",CHARWIDTH);
	fprintf(fd,"\tclosepath\n");
	fprintf(fd,"      } def\n\n");
	
	fprintf(fd,"/color_char { gsave\n");
	fprintf(fd,"\tsetrgbcolor\n");
	fprintf(fd,"\tmoveto\n");
	fprintf(fd,"\tshow\n");
	fprintf(fd,"\tgrestore\n");
	fprintf(fd,"      } def\n\n");
	
	fprintf(fd,"/cbox { gsave\n");
	fprintf(fd,"\ttranslate\n");
	fprintf(fd,"\tnewpath\n");
	fprintf(fd,"\t0 0 moveto\n");
	fprintf(fd,"\tlineto\n");
	fprintf(fd,"\tlineto\n");
	fprintf(fd,"\tlineto\n");
	fprintf(fd,"\tclosepath\n");
	fprintf(fd,"\tfill\n");
	fprintf(fd,"\tgrestore\n");
	fprintf(fd,"      } def\n\n");

	fprintf(fd,"/color_inv { gsave\n");
	fprintf(fd,"\tsetrgbcolor\n");
	fprintf(fd,"\ttranslate\n");
	fprintf(fd,"\tbox fill\n");
	fprintf(fd,"\tgrestore\n");
	fprintf(fd,"\tmoveto\n");
	fprintf(fd,"\tshow\n");
	fprintf(fd,"      } def\n\n");

        fprintf(fd,"/white_inv { gsave\n");
        fprintf(fd,"\tsetrgbcolor\n");
        fprintf(fd,"\ttranslate\n");
        fprintf(fd,"\tbox fill\n");
        fprintf(fd,"\tgrestore\n");
        fprintf(fd,"\tgsave\n");
        fprintf(fd,"\tsetrgbcolor\n");
        fprintf(fd,"\tmoveto\n");
        fprintf(fd,"\tshow\n");
        fprintf(fd,"\tgrestore\n");
        fprintf(fd,"      } def\n\n");

	if (pagesize==A3)
		fprintf(fd,"statusdict begin a3 end\n\n");
/* For canon color printer, use a3tray instead of a3!! */
}

static void print_page_header(FILE *fd,int ps_rotation,int maxx,int maxy,
int page,int numpages,Boolean header,char *str_time,
char *ps_file,int ps_xtrans,int ps_ytrans,float ps_scale)
{
	int ps_x,ps_y;
	char tstr[50];

	fprintf(fd,"%%%%Page: P%d\n",page);
	fprintf(fd,"save\n\n");

	if (ps_rotation==-90)
	{
		fprintf(fd,"0 %d translate\n",maxx);
		fprintf(fd,"%d rotate\n",ps_rotation);
	}

	if (header)
	{
		sprintf(tstr,"CLUSTAL %s MULTIPLE SEQUENCE ALIGNMENT",revision_level);
		ps_x = (maxx-strlen(tstr)*10)/2;
		ps_y = maxy - 2*CHARHEIGHT;
		fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
		fprintf(fd,"/Times-Bold findfont 14 scalefont setfont\n");
		fprintf(fd,"(%s) show\n\n",tstr);

		ps_x = 20;
		ps_y = maxy - 4*CHARHEIGHT;
		fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
		fprintf(fd,"(File: %s) show\n\n",ps_file);

		sprintf(tstr,"Date: %s",str_time);
		ps_x = maxx-strlen(tstr)*8-20;
		ps_y = maxy - 4*CHARHEIGHT;
		fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
		fprintf(fd,"(%s) show\n\n",tstr);

		sprintf(tstr,"Page %d of %d",page+1,numpages);
		ps_x = 20;
		ps_y = maxy - 5*CHARHEIGHT-4;
		fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
		fprintf(fd,"(%s) show\n\n",tstr);
	}	
	fprintf(fd,"%d %d translate\n",ps_xtrans,ps_ytrans);
	fprintf(fd,"%#3.2f %#3.2f scale\n",ps_scale,ps_scale);
	fprintf(fd,"/Courier-Bold findfont 10 scalefont setfont\n");
}

static void print_header_line(FILE *fd,panel_data name_data, panel_data seq_data,
int ix,int fr,int lr)
{
	int i;
	int ps_x,ps_y;

	ps_x = LEFTMARGIN;
	ps_y = header_top - (ix * CHARHEIGHT);
	fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
	fprintf(fd,"(%*s ) show\n",max_names,name_data.header[ix]);
	for(i=fr;i<=lr;i++)
	{
		ps_x = block_left + (i-fr) * CHARWIDTH; 
		fprintf(fd,"(");
		fprintf(fd,"%c",seq_data.header[ix][i]);
		fprintf(fd,") ");
		fprintf(fd,"%d %d %d %d 1.0 1.0 1.0 color_inv\n",ps_x,ps_y,ps_x,ps_y);
	}
	fprintf(fd,"\n");
}

static void print_footer_line(FILE *fd,panel_data name_data, panel_data seq_data,
int ix,int fr,int lr)
{
	int i;
	int ps_x,ps_y;

	ps_x = LEFTMARGIN;
	ps_y = footer_top - (ix * CHARHEIGHT);
	fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
	fprintf(fd,"(%*s ) show\n",max_names,name_data.footer[ix]);
	for(i=fr;i<=lr;i++)
	{
		ps_x = block_left + (i-fr) * CHARWIDTH; 
		fprintf(fd,"(");
		fprintf(fd,"%c",seq_data.footer[ix][i]);
		fprintf(fd,") ");
		fprintf(fd,"%d %d %d %d 1.0 1.0 1.0 color_inv\n",ps_x,ps_y,ps_x,ps_y);
	}
	fprintf(fd,"\n");
}

static void print_quality_curve(FILE *fd,panel_data seq_data,
int fr,int lr,int score_height)
{
	int i,w,h;
	int ps_x,ps_y,curve_bottom;

	w=CHARWIDTH;
	ps_x = block_left;
       	curve_bottom=curve_top-score_height*CHARHEIGHT;
	fprintf(fd,"0.3 0.3 0.3 setrgbcolor\n");
	for(i=fr;i<=lr;i++)
	{
		fprintf(fd,"%d %d moveto\n",ps_x,curve_bottom);
       		h=score_height*CHARHEIGHT*((float)seq_data.colscore[i]/100.0);
		if(h<1) h=1;
		fprintf(fd,"%d 0 %d %d 0 %d %d %d cbox\n",w,w,h,h,ps_x,curve_bottom);
		ps_x+=CHARWIDTH;
	}
	fprintf(fd,"0.0 0.0 0.0 setrgbcolor\n");
}

static void print_seq_line(FILE *fd,panel_data name_data, panel_data seq_data,
int row,int seq,int fr,int lr,int res_number)
{
	int i,color;
	int ps_x,ps_y;
	float red, green, blue;

	ps_x = LEFTMARGIN;
	ps_y = seq_top - (row * CHARHEIGHT);
	fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
	fprintf(fd,"(%*s ) show\n",max_names,name_data.lines[seq]);
	for(i=fr;i<=lr;i++)
	{
		color = seq_data.colormask[seq][i];
		red = color_lut[color].pr;
		green = color_lut[color].pg;
		blue = color_lut[color].pb;
		ps_x = block_left + (i-fr) * CHARWIDTH; 
		fprintf(fd,"(");
		fprintf(fd,"%c",seq_data.lines[seq][i]);
		fprintf(fd,") ");
                if(segment_exceptions && seq_data.segment_exception[seq][i] > 0)
                {
                       fprintf(fd,"%d %d %1.1f %1.1f %1.1f %d %d %1.1f %1.1f %1.1f white_inv\n",
				ps_x,ps_y,1.0,1.0,1.0,ps_x,ps_y,0.1,0.1,0.1);
                }
                else if(residue_exceptions && seq_data.residue_exception[seq][i] == TRUE)
                {
                       fprintf(fd,"%d %d %1.1f %1.1f %1.1f %d %d %1.1f %1.1f %1.1f white_inv\n",
				ps_x,ps_y,1.0,1.0,1.0,ps_x,ps_y,0.4,0.4,0.4);
                }
                else
                {
                       if(inverted)
                              fprintf(fd,"%d %d %d %d %1.1f %1.1f %1.1f color_inv\n",
				ps_x,ps_y,ps_x,ps_y,red,green,blue);
                       else
                              fprintf(fd,"%d %d %1.1f %1.1f %1.1f color_char\n",
				ps_x,ps_y,red,green,blue);
                }
	}

	if(res_number>0)
	{
		ps_x = block_left + (lr-fr+1) * CHARWIDTH; 
		ps_y = seq_top - (row * CHARHEIGHT);
		fprintf(fd,"%d %d moveto\n",ps_x,ps_y);
		fprintf(fd,"(%*d) show\n",MAXRESNO,res_number);
	}
	fprintf(fd,"\n");
}
