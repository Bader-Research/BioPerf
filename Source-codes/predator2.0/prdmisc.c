#include "prd.h"
#include "prdprot.h"

/* Split a char string into text fields */

int SplitString(char *Buffer, char **Fields, int MaxField)
{
  int FieldCnt, SymbCnt, FieldFlag, BuffLen;
  static char LocalBuffer[BUFSZ];


  FieldCnt =0; FieldFlag = 0;
  BuffLen = (int)strlen(Buffer) - 1;

  strcpy(LocalBuffer,Buffer);

  for(SymbCnt=0; SymbCnt<BuffLen; SymbCnt++) {
    if( (isspace(LocalBuffer[SymbCnt])) && FieldFlag == 0 && SymbCnt != BuffLen-1 ) continue;
    if( (!isspace(LocalBuffer[SymbCnt])) && FieldFlag == 1 && SymbCnt == BuffLen-1 ) {
      LocalBuffer[SymbCnt+1] = '\0';
      return(FieldCnt);
    }
    else
      if( (isspace(LocalBuffer[SymbCnt])) && FieldFlag == 1 ) {
        LocalBuffer[SymbCnt] = '\0';
        FieldFlag = 0;
        if( FieldCnt == MaxField ) return(FieldCnt);
      }
      else
        if( (!isspace(LocalBuffer[SymbCnt])) && FieldFlag == 0 ) {
          FieldFlag = 1;
          Fields[FieldCnt] = LocalBuffer+SymbCnt;
          FieldCnt++;
        }
  }

  return(FieldCnt);
}

void die(char *format, ... ) {
  void exit(int return_code);
  va_list ptr;
  va_start(ptr,format);
  vfprintf(stderr,format,ptr);
  exit(1);
  va_end(ptr);
}

#ifndef __ICM__
bool wrer(char *format, ... ) {
va_list ptr;

va_start(ptr,format);
vfprintf(stderr,format,ptr);
va_end(ptr);
return(0);
}
#endif
void NullMatrix(float **Matrix, int N, int M)
{

  register int i, j;
  
  for( i=0; i<N; i++ )
    for( j=0; j<M; j++ )
      Matrix[i][j] = 0.0;
}

