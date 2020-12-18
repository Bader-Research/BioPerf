/*
 
Copyright (c)  1995-2000   The Regents of the University of California
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this PDB
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this PDB Software into commercial products
or use for commercial purposes should contact the Technology Transfer
Office, University of California, San Diego, 9500 Gilman Drive, La Jolla,
CA 92093-0910, Ph: (619) 534-5815, FAX: (619) 534-7345.
 
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS PDB SOFTWARE, EVEN IF THE
UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
 
THE PDB SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE PDB SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
 
*/
////////////////////////////////////////////////////////////////////
//     POM V2.0 (C) 1996-2000  I.Shindyalov, H.Weissig, P.Bourne  //
////////////////////////////////////////////////////////////////////
#ifndef H_POM
#define H_POM
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <dirent.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <unistd.h>
////////////////////////////////////////////////////////////////////
#ifndef DEFAULT_CACHE
#define DEFAULT_CACHE 1000
#endif

#define PATH_SEPARATOR               "/"

//////////////////////////////////////////////////////////////////// 
typedef char int1;
#ifndef T3E
typedef int int4;
#else
typedef short int4;
#endif
#ifndef INT2_OFF
typedef short int2;
#else
typedef int4 int2;
#endif
typedef float flt4;
typedef double flt8;
///////////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
enum DataType {D_INT1, D_INT2, D_INT4, D_FLT4, D_FLT8};
#else
enum DataType {D_INT1, D_INT2=2, D_INT4=2, D_FLT4, D_FLT8};
#endif
enum CopyArray {COPY_ARRAY_OFF, COPY_ARRAY_ON};
enum LoadInMemory {LOAD_IN_MEMORY_OFF, LOAD_IN_MEMORY_ON};
///////////////////////////////////////////////////////////////////////////
class IAtom;
class IEntity;
class ISubentity;
class IMonomer;
class IBond;
class Monomers;
//////////////////////////////////////////////////////////////////
struct XYZ {
public:
  double X, Y, Z;
void operator +=(const XYZ& xyz)
  { this->X+=xyz.X; this->Y+=xyz.Y; this->Z+=xyz.Z; };
void operator -=(const XYZ& xyz)
  { this->X-=xyz.X; this->Y-=xyz.Y; this->Z-=xyz.Z; };
void operator /=(int i)
  { this->X/=(double)i; this->Y/=(double)i; this->Z/=(double)i; };
void operator /=(double i)
  { this->X/=i; this->Y/=i; this->Z/=i; };
void operator *=(double i) { this->X*=i; this->Y*=i; this->Z*=i; };
void operator =(double i) { this->X=i; this->Y=i; this->Z=i; };
void XYZ2() { X=X*X; Y=Y*Y; Z=Z*Z; };
double dist(XYZ xyz) {XYZ xyz_=xyz; xyz_-=*this; xyz_.XYZ2();
		    double sum=xyz_.X+xyz_.Y+xyz_.Z;
		    return(sum>0.0?sqrt(sum):0.0); };
double dist2(XYZ xyz) {XYZ xyz_=xyz; xyz_-=*this; xyz_.XYZ2();
		    double sum=xyz_.X+xyz_.Y+xyz_.Z; return(sum); };
};
//////////////////////////////////////////////////////////////////
class AminoacidProperty {
public:
  char *code1;
  char *code3;
  int number;
  double molWeight;
  double volume;
  double polarity;
  double isoelectricPoint;
  double hydrophobicity;
  double meanExposure;
  double ChouFasmanAlpha;
  double ChouFasmanBeta;
};
///////////////////////////////////////////////////////////////////////////
class DB {
 public:
  char **path;
  int nPath;
  Monomers *monomers;

  DB(char *path_=NULL);
  ~DB();
  void setDefault();
  void setPath(char *path_);
  void addPath(char *path_);
  char* findPath(char *file_);
  int openFile(char *file_);
  void readFile(void *array, int4 nArray, int4 adr, int4 size, int fd);
  void writeFile(void *array, int4 nArray, int4& adr, int4 size, int fd);
  void closeFile(int );
  int testFile(char *file);
  static int testPathFile(char *path_, char *file_);
  int createFile(char *file);
  void testAdr(int4 adr);
  void addAdr(int4 adr1, int4 nBytes, int4& adr2);
};
////////////////////////////////////////////////////////////////////
class Property {
 public:
  int fd; // file descriptor
  int load_in_memory, load_in_memory_orig;
  static DB *db_default;
  DB *db;
  int writeFlag;
  char *name;

  int4 nItems;
  int4 nOptItems;
  DataType type; 
  int4 size;    // 0-var, otherwise size;
  int4 adrProp;
  int4 adrRef;

  int1 *prop1;
#ifndef INT2_OFF
  int2 *prop2;
#endif
  int4 *prop4;
  flt4 *propf;
  flt8 *propd;
  int4 nProp;

  int4 *ref;

  int4 iObj, iPos;

  Property(DB *db_=db_default);
  Property(char *name_, DB *db_=db_default);  
  Property(char *name_, int load_in_memory_, DB *db_=db_default);  
  ~Property();
  void create(char *name_, DataType type_, int4 size_, 
	      int load_in_memory_=LOAD_IN_MEMORY_ON, DB *db_=db_default, int cacheSize = DEFAULT_CACHE);
  void clear();
  void open(char *name_, DB *db_=db_default);
  void open(char *name_, int load_in_memory_, DB *db_=db_default,
	    int cacheSize = DEFAULT_CACHE);
  int isOpen();
  void saveHeader(int with_load_mode=0);
  void close(int export_array=0);
  void addItem(int1 item);
#ifndef INT2_OFF
  void addItem(int2 item);
#endif
  void addItem(int4 item);
  void addItem(flt4 item);
  void addItem(flt8 item);
  void addItem(int1* item);
#ifndef INT2_OFF
  void addItem(int2* item);
#endif
  void addItem(int4* item);
  void addItem(flt4* item);
  void addItem(flt8* item);
  void addItem(int1* item, int4 nValues);
#ifndef INT2_OFF
  void addItem(int2* item, int4 nValues);
#endif
  void addItem(int4* item, int4 nValues);
  void addItem(flt4* item, int4 nValues);
  void addItem(flt8* item, int4 nValues);
  void addItem();

  void setItem(int4 index, int1 item);
#ifndef INT2_OFF
  void setItem(int4 index, int2 item);
#endif
  void setItem(int4 index, int4 item);
  void setItem(int4 index, flt4 item);
  void setItem(int4 index, flt8 item);
  void setItem(int4 index, int1* item);
#ifndef INT2_OFF
  void setItem(int4 index, int2* item);
#endif
  void setItem(int4 index, int4* item);
  void setItem(int4 index, flt4* item);
  void setItem(int4 index, flt8* item);
  void setItem(int4 index, int1* item, int4 nValues);
#ifndef INT2_OFF
  void setItem(int4 index, int2* item, int4 nValues);
#endif
  void setItem(int4 index, int4* item, int4 nValues);
  void setItem(int4 index, flt4* item, int4 nValues);
  void setItem(int4 index, flt8* item, int4 nValues);

  void extendItem(int4 index, int1* item);
  void extendItem(int4 index, int1* item, int4 nValues);
  void extendItem(int4 index, int4 item);
  void extendItem(int4 index, int4* item, int4 nValues);

  int4 getItemSize(int4);
  int4 getSize() {return(size);};
  DataType getType() {return(type);};
  int4 getObjectSize() {return(nItems);};
  int4 getPropertySize() {return(nProp);};

  int1* item1(int4 index, int copy_array=COPY_ARRAY_OFF); // >1
#ifndef INT2_OFF
  int2* item2(int4 index, int copy_array=COPY_ARRAY_OFF); // >1
#else
  int2* item2(int4 index, int copy_array=COPY_ARRAY_OFF) 
    {return(item4(index, copy_array));};
#endif
  int4* item4(int4 index, int copy_array=COPY_ARRAY_OFF); // >1 
  flt4* itemf(int4 index, int copy_array=COPY_ARRAY_OFF); // >1 
  flt8* itemd(int4 index, int copy_array=COPY_ARRAY_OFF); // >1 
  int1* item1n(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF);// 0
#ifndef INT2_OFF
  int2* item2n(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF); // 0
#else
  int2* item2n(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF)
    {return(item4n(index, nValues, copy_array));}; 
#endif
  int4* item4n(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF); // 0
  flt4* itemfn(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF); // 0
  flt8* itemdn(int4 index, int4& nValues, int copy_array=COPY_ARRAY_OFF); // 0

  void array(int1 *);
#ifndef INT2_OFF
  void array(int2 *);
#endif
  void array(int4 *);
  void array(flt4 *);
  void array(flt8 *);
  int4 find(int1 *text, int mismatches=0);
  int4 find(int1 *text, int mismatches, int4& pos) {
            int next = 0; return(find(text, mismatches, pos, next));};
  int4 find(int1 *text, int mismatches, int4& pos, int&);
  int4 find(int1 t1, int1 t2, int& next);
  int4 find(int1 t1, int1 t2=-1){int next=0; return(find(t1, t2, next));};
#ifndef INT2_OFF
  int4 find(int2 t1, int2 t2, int& next);
  int4 find(int2 t1, int2 t2=-1){int next=0; return(find(t1, t2, next));};
#endif
  int4 find(int4 t1, int4 t2, int& next);
  int4 find(int4 t1, int4 t2=-1){int next=0; return(find(t1, t2, next));};
  int4 find(flt4 t1, flt4 t2, int& next);
  int4 find(flt4 t1, flt4 t2=-1){int next=0; return(find(t1, t2, next));};
  int4 find(flt8 t1, flt8 t2, int& next);
  int4 find(flt8 t1, flt8 t2=-1){int next=0; return(find(t1, t2, next));};

private:
  void updateRef(int4 index, int4 sizeDiff);
};
/////////////////////////////////////////////////////////////////////
class Collection {
 public:
  int fd; // file descriptor
  static DB *db_default;
  DB *db;
  int writeFlag;
  char *name;
  
  int4 nClasses;
  DataType type; 
  int4 adrColl;
  int4 adrNColl;
  int load_in_memory; 
  
#ifndef INT2_OFF
  int2 *coll2t;
#endif
  int4 *coll4t;
  int4 *nCollt;
  int4 index_;
  
#ifndef INT2_OFF
  int2 **coll2;
#endif
  int4 **coll4;
  int4 *nColl;
  
  Collection(DB *db_=db_default);
  Collection(char *name_, DB *db_=db_default);  
  ~Collection();
  void create(char *name_, DataType type_, int4 nClasses_, DB *db_=db_default);
  void clear();
  void open(char *name_, DB *db_=db_default);
  void open(char *name_, int load_in_memory, DB *db_=db_default, 
	    int cacheSize = DEFAULT_CACHE);
  int isOpen();
  void save();
  void close();
#ifndef INT2_OFF
  void add(int4 index, int2 item);
#endif
  void add(int4 index, int4 item);
#ifndef INT2_OFF
  void add(int4 index, int2* item, int4 nValues);
#endif
  void add(int4 index, int4* item, int4 nValues);
  
  static void expandCol(Collection *, char *, int);
  
  DataType getType() {return(type);};
  int4 getCollectionSize() {return(nClasses);};
  int4 getClassSize(int4 index) {return(nColl[index]);};
  int4* getClassSize() {return(nColl);};
  
#ifndef INT2_OFF
  int2 itemValue2(int4 index, int4 cindex);
  int2* item2(int4 index);
#else
  int2 itemValue2(int4 index, int4 cindex){return(itemValue2(index, cindex));};
  int2* item2(int4 index){return(item4(index));};
#endif
  int4 itemValue4(int4 index, int4 cindex);
  int4 itemValue(int4 index, int4 cindex);
  int4* item4(int4 index); 
#ifndef INT2_OFF
  int2* item2n(int4 index, int4& nValues); 
#else
  int2* item2n(int4 index, int4& nValues){return(item4n(index,nValues));}; 
#endif
  int4* item4n(int4 index, int4& nValues); 

#ifndef INT2_OFF
  int2** array2(){return(coll2);};
#else
  int2** array2(){return(array4());};
#endif
  int4** array4(){return(coll4);};
};
/////////////////////////////////////////////////////////////////////
class Monomers {
 public:
  static DB *db_default;
  DB *db;
  int nMonomers;
 private:
  Property *monCode3;
  Property *monCode1;
  Property *monNbond;
  Property *monBond;
  Property *monNatom;
  Property *monAtom;
  Property *monPrev;
  Property *monNext;
  Property *monType;

  Property *comID;
  
 public:
  Monomers(DB *db_=db_default);
  ~Monomers();
  int findCom(char *text);
  int1* code3(int iM);
  int1 code1(int iM);
  int nBond(int iM);
  int2* bonds(int iM);
  int2 bond1(int iM, int iB);
  int2 bond2(int iM, int iB);
  int2 nAtom(int iM);
  int1* atoms(int iM);
  int1* atom(int iM, int iA);
  int type(int iM);
  int findAtom(int iM, char *);
  int prev(int iM);
  int next(int iM);
};
////////////////////////////////////////////////////////////////////
class Entities {
  friend class IEntity;
  friend class ISubentity;
  friend class IAtom;
  friend class IMonomer;
  friend class IBond;
 public:
  static DB *db_default;
  DB *db;
  Monomers *monomers;

 private:

// Mandatory Loading

  int4 nEnc;
  int1 **encName;
  int4 *comNum;
  int4 *encNum;
  int4 *enpNum;
  int2 *encNSE;

// Optional Loading: ENC properties
  int2 **encSE;
  int1 ***encSEnum;

  int4 *encNXYZ;
  flt4 **encXYZ;
  flt4 **encBfac;
  int4 **enc_se_xyz;
  int4 **enc_xyz_se;

// Optional Loading: ENP properties
  int1 **encpSeq;
  int1 **encpKS;

  flt4 **encpExp;
  flt4 **encpPol;
  flt4 **encpBfac;

  flt4 **encpCa;
  int1 **encpMonType;

 public:
  Entities(DB *db_=db_default);
  ~Entities();
  void addCom(char *com_id, int isAdd=0);
  void addCom(int4 com_ind, int isAdd=0);
  void addEnc(char *enc_id, int isAdd=0);
  void addEnc(int4 enc_ind, int isAdd=0);
  void addEnp(char *enp_id, int isAdd=0);
  void addEnp(int4 enp_ind, int isAdd=0);

  void load_encSE(int);
  void load_encSEnum(int);
  void load_encNXYZ(int);
  void load_encXYZ(int);
  void load_encBfac(int);
  void load_enc_se_xyz(int);
  void load_enc_xyz_se(int);
  void load_encpSeq(int);
  void load_encpKS(int);
  void load_encpExp(int);
  void load_encpPol(int);
  void load_encpBfac(int);
  void load_encpCa(int);
  void load_monType(int);

  int test_encp(int iE){return(enpNum[iE]>0);};

  void reread();
};
////////////////////////////////////////////////////////////////////
class IEntity {
  int iE;
 public:
  Entities *entities;

  IEntity() {entities=NULL; iE=-1;};
  IEntity(Entities *entities_) {entities=entities_; iE=0;};
  operator void*() {return((void*)(iE>=0 && iE<entities->nEnc));};
  void operator++() {if(iE>=0) iE++;};
  void operator++(int) {if(iE>=0) iE++;};
  int operator() () {return(iE);};
  void reset() {if(iE>0) iE=0;};
  int1* name() {return(entities->encName[iE]);};
  int1* atoms(int is);
  int1* atom(int is, int ia);
  int NSE() {return(entities->encNSE[iE]);};
  void set(int iE_) {iE=iE_;};
};
////////////////////////////////////////////////////////////////////
class ISubentity {
  friend class IMonomer;
  int iE;
  int iS;
  int iM;
 public:
  Entities *entities;
  Monomers *monomers;

  ISubentity() {entities=NULL; monomers=NULL; iE=-1; iS=-1; iM=-1;};
  ISubentity(IEntity);
  void set(IEntity ie) {entities=ie.entities; 
                        monomers=entities->monomers; iE=ie(); iS=0; iM=-1;};
  void set(IEntity ie, int is) {entities=ie.entities;   
                        monomers=entities->monomers; iE=ie(); iS=is; iM=-1;};
  void set(int is) {iS=is;iM=-1;};
  operator void*() {return((void*)(iS>=0 && iS<entities->encNSE[iE]));};
  void operator++() {if(iS>=0) {iS++;iM=-1;}};
  void operator++(int) {if(iS>=0) {iS++;iM=-1;}};
  void operator+=(int s) {if(iS>=0) {iS+=s;iM=-1;}};
  void reset() {if(iS>0) {iS=0;iM=-1;}};
  int operator() () {return(iS);};
  int myEntityIndex() {return(iE);};
  ISubentity neighbor(int nb);

  // entity info
  int1* pdbNum() {if(!entities->encSEnum[iE]) entities->load_encSEnum(iE);
    return(*(entities->encSEnum[iE]+iS));};
  int1 ks() {if(entities->enpNum[iE]<0) return(' ');
	     if(!entities->encpKS[iE]) entities->load_encpKS(iE);
	     return(*(entities->encpKS[iE]+iS));};
  flt4 exp() {if(entities->enpNum[iE]<0) return(-1.0);
	      if(!entities->encpExp[iE]) entities->load_encpExp(iE);
	      return(*(entities->encpExp[iE]+iS));};
  flt4 pol() {if(entities->enpNum[iE]<0) return(-1.0);
	      if(!entities->encpPol[iE]) entities->load_encpPol(iE);
	      return(*(entities->encpPol[iE]+iS));};
  int monType() {if(entities->enpNum[iE]<0) return(0);
		 if(!entities->encpMonType[iE]) entities->load_monType(iE);
		 return(*(entities->encpMonType[iE]+iS));};
  flt4 bfac() {if(entities->enpNum[iE]<0) return(-1.0);
	      if(!entities->encpBfac[iE]) entities->load_encpBfac(iE);
	      return(*(entities->encpBfac[iE]+iS));};

  flt4 ca_x() {if(!entities->encpCa[iE]) entities->load_encpCa(iE);
		 return(*(entities->encpCa[iE]+iS*3));};
  flt4 ca_y() {if(!entities->encpCa[iE]) entities->load_encpCa(iE);
		 return(*(entities->encpCa[iE]+iS*3+1));};
  flt4 ca_z() {if(!entities->encpCa[iE]) entities->load_encpCa(iE);
		 return(*(entities->encpCa[iE]+iS*3+2));};
  XYZ ca_xyz() {if(!entities->encpCa[iE]) entities->load_encpCa(iE);
		XYZ s; flt4 *p=entities->encpCa[iE]+iS*3; s.X=*p;
		s.Y=*(p+1); s.Z=*(p+2); return(s);};

  // monomer info
  int getMon();
  int iFindAtom(char *text) {if(iM==-1)iM=getMon();
                              return(monomers->findAtom(iM,text));};
  IAtom findAtom(char *text);
  int isAtom(int1 *text);
  int1 code1(){if(iM==-1)iM=getMon();return(monomers->code1(iM));};
  int1* code3(){if(iM==-1)iM=getMon();return(monomers->code3(iM));};
  int code();
  int nBond(){if(iM==-1)iM=getMon();return(monomers->nBond(iM));};
  int2* bonds(){if(iM==-1)iM=getMon();return(monomers->bonds(iM));};
  int2 bond1(int iB){if(iM==-1)iM=getMon();return(monomers->bond1(iM,iB));};
  int2 bond2(int iB){if(iM==-1)iM=getMon();return(monomers->bond2(iM,iB));};
  int nAtom(){if(iM==-1)iM=getMon();return(monomers->nAtom(iM));};
  int1* atoms(){if(iM==-1)iM=getMon();return(monomers->atoms(iM));};
  int1* atom(int iA){if(iM==-1)iM=getMon();return(monomers->atom(iM,iA));};
  int4 prev(int iM){if(iM==-1)iM=getMon();return(monomers->prev(iM));};
  int4 next(int iM){if(iM==-1)iM=getMon();return(monomers->next(iM));};
};
////////////////////////////////////////////////////////////////////
class IAtom {
  friend class ISubentity;
  friend class IMonomer;
  friend class IBond;
  int iE;
  int iS;
  int iA;
  char _name[5];
 public:
  Entities *entities;

  IAtom() {entities=NULL; iE=-1; iS=-1; iA=-1;};
  IAtom(IEntity);
  IAtom(ISubentity);
  operator void*() {if(!entities->encNXYZ[iE]) entities->load_encNXYZ(iE);
		    if(!entities->enc_xyz_se[iE]) 
		      entities->load_enc_xyz_se(iE);
		    if(iA<0 || iA>=entities->encNXYZ[iE]) return((void*)0);
		    if(iS>=0 && iS!=*(entities->enc_xyz_se[iE]+iA)) 
		      return((void*)0); return ((void*)1);};
  void operator++() {if(iA>=0) iA++;};
  void operator++(int) {if(iA>=0) iA++;};
  void reset() {if(!entities->enc_se_xyz[iE]) entities->load_enc_se_xyz(iE);
    if(iA>0) iA=(iS>=0?*(entities->enc_se_xyz[iE]+iS):0);};
  void set(int ia) {iA=ia+(*(entities->enc_se_xyz[iE]+iS));};
  int operator() () {return(iA-(*(entities->enc_se_xyz[iE]+iS)));};
  int myEntityIndex() {return(iE);};
  int mySubentityIndex() {return(iS);};
  flt4 x() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	      return(*(entities->encXYZ[iE]+iA*3));};
  flt4 y() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	      return(*(entities->encXYZ[iE]+iA*3+1));};
  flt4 z() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	      return(*(entities->encXYZ[iE]+iA*3+2));};
  flt4 bfac() {if(!entities->encBfac[iE]) entities->load_encBfac(iE);
		 return(*(entities->encBfac[iE]+iA));};
  int ix() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	    return((int)*(entities->encXYZ[iE]+iA*3)*1000);};
  int iy() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	    return((int)*(entities->encXYZ[iE]+iA*3+1)*1000);};
  int iz() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	    return((int)*(entities->encXYZ[iE]+iA*3+2)*1000);};
  XYZ xyz() {if(!entities->encXYZ[iE]) entities->load_encXYZ(iE);
	     XYZ s; flt4 *p=entities->encXYZ[iE]+iA*3; s.X=*p;
	     s.Y=*(p+1); s.Z=*(p+2); return(s);};
  int1* name();
};
/////////////////////////////////////////////////////////////////////
class IMonomer {
  int iM;
  char _name[5];
 public:
  Monomers *monomers;

  IMonomer() {monomers=NULL; iM=-1;};
  IMonomer(ISubentity);
  operator void*(){return((void*)(iM<monomers->nMonomers));};
  void operator++(){if(iM>=0) iM++;};
  void operator++(int){if(iM>=0) iM++;};
  int operator() () {return(iM);};
  void reset() {if(iM>0) iM=0;};
  int1 code1() {return(monomers->code1(iM));};
  int1* code3() {return(monomers->code3(iM));};
  int nBond() {return(monomers->nBond(iM));};
  int nAtom() {return(monomers->nAtom(iM));};
  int1* atomName(IAtom);
};
//////////////////////////////////////////////////////////////////
class IBond {
  int iB;
  int iM;
 public:
  Monomers *monomers;

  IBond() {monomers=NULL; iB=-1; iM=-1;};
  IBond(IMonomer mon) {monomers=mon.monomers; iB=0; iM=mon();};
  IBond(IAtom);
  operator void*() {return((void*)(iB<monomers->nBond(iM)));};
  void operator++() {if(iB>0) iB++;};
  void operator++(int) {if(iB>0) iB++;};
  int operator() () {return(iB);};
  void reset() {if(iB>0) iB=0;};
  int atom1() {return(monomers->bond1(iM, iB));};
  int atom2() {return(monomers->bond2(iM, iB));};
};
//////////////////////////////////////////////////////////////////
#endif
