#include "prd.h"
#include "prdprot.h"

extern char   *AMINOACIDS;
extern char   *ALNSYMBOLS;
extern char   *STRUCTURES;

char *OneToThree(char One)
{

  if( One == 'A' ) return("ALA");
  else if( One == 'R' ) return("ARG");
  else if( One == 'N' ) return("ASN");
  else if( One == 'D' ) return("ASP");
  else if( One == 'B' ) return("ASX");
  else if( One == 'C' ) return("CYS");
  else if( One == 'Q' ) return("GLN");
  else if( One == 'E' ) return("GLU");
  else if( One == 'Z' ) return("GLX");
  else if( One == 'G' ) return("GLY");
  else if( One == 'H' ) return("HIS");
  else if( One == 'I' ) return("ILE");
  else if( One == 'L' ) return("LEU");
  else if( One == 'K' ) return("LYS");
  else if( One == 'M' ) return("MET");
  else if( One == 'P' ) return("PRO");
  else if( One == 'F' ) return("PHE");
  else if( One == 'S' ) return("SER");
  else if( One == 'T' ) return("THR");
  else if( One == 'W' ) return("TRP");
  else if( One == 'Y' ) return("TYR");
  else if( One == 'V' ) return("VAL");
  else if( One == 'X' ) return("UNK");
  else return("***");
}

void Char2Int(char *String, int *IString, int Len)
{

  register int i;

  for( i=0; i<Len; i++ ) {
    IString[i] = (int)(strchr(AMINOACIDS,String[i])-AMINOACIDS);
    if( IString[i] == 22 )
      IString[i] = 2;
    else
    if( IString[i] == 21 )
      IString[i] = 5;
    else
    if( IString[i] == 20 )
      IString[i] = 7;
  }
}


void Int2Char(int *IString, char *String, int Len)
{
  register int i;

  for( i=0; i<Len; i++ ) {
    if( IString[i] == 22 )
      IString[i] = 2;
    else
    if( IString[i] == 21 )
      IString[i] = 5;
    else
    if( IString[i] == 20 )
      IString[i] = 7;
    String[i] = AMINOACIDS[IString[i]];
  }
}


void Str2Int(char *String, int *IString, int Len)
{

  register int i;

  for( i=0; i<Len; i++ )
    IString[i] = (int)(strchr(STRUCTURES,_2C(String[i]))-STRUCTURES);
}

void CorrectAsn(char *Asn, int Length, char SecStrType, char EditChar, int MaxLength)
{

  int NStr = 0, Res, Flag = 0, Bound[MAX_ASSIGN][2], i;

  for( Res=0; Res<Length; Res++ ) {
    if( Asn[Res] == SecStrType && Flag == 0 ) {
      Flag = 1; 
      Bound[NStr][0] = Res;
    }
    else
    if( Asn[Res] != SecStrType && Flag == 1 ) {
      Flag = 0; 
      Bound[NStr++][1] = Res-1;
    }
  }

  for( i=0; i<NStr; i++ )
    if( Bound[i][1]-Bound[i][0]+1 <= MaxLength )
      for( Res=Bound[i][0]; Res<=Bound[i][1]; Res++ ) 
	Asn[Res] = EditChar;
}

