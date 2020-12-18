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
#include "pom.h"
#include "miscutil.h"

int _getSpaceNeeded(int nArray, int cacheSize = DEFAULT_CACHE);
int _getSpaceNeeded(int nArray, int cacheSize) {
    return (nArray + (cacheSize - leftOver(nArray, cacheSize)));
} 
////////////////////////////////////////////////////////////////////
template <class A> void InsertItems(A **array,
				    int4 n,
				    A *vals,
				    int4 nVals,
				    int4 *nObj,
				    int4 diff,
				    int cacheSize = DEFAULT_CACHE) {

   (*nObj) += diff;
   int space = _getSpaceNeeded((*nObj), cacheSize);
   
#ifdef DEBUG
   fprintf(stderr, "InsertItems(array, %d, vals, %d, %d, %d) allocated %d\n",  n, nVals, (*nObj)-diff, diff, space);
#endif

   A *newArray = new A[space];
   int i, old, nw = 0;

   for (old = 0; old < n; old++, nw++) newArray[nw] = *(*array + old);

   for (i = 0; i < nVals; i++, nw++) newArray[nw] = vals[i];

   old = n + nVals - diff;
   for (; nw < (*nObj); nw++, old++) newArray[nw] = *(*array + old);

   if (((*nObj) - diff)) delete [] (*array);
   (*array) = newArray;
};
////////////////////////////////////////////////////////////////////
template <class A> void AddToArray(A ** array,
				   A value,
				   int Narray,
				   int cacheSize = DEFAULT_CACHE) {

	int space;

	if (leftOver(Narray, cacheSize)) {
		*(*array + Narray) = value;
	} else {
		A * array_;

		space = _getSpaceNeeded(Narray+1, cacheSize);
#ifdef DEBUG
    fprintf(stderr, "AddToArray allocated %d\n", space);
#endif
		array_=new A[space];
		    
		for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
		array_[Narray]=value;
		if(Narray) delete [] (*array);
		(*array)=array_;
	}
};
///////////////////////////////////////////////////////////////////////////
template <class A> void AddToArrayN(A ** array,
				    A * values,
				    int Narray,
				    int Nvalues,
				    int cacheSize = DEFAULT_CACHE) {
	int space;
	if (Nvalues <= 0) return;

	if (leftOver(Narray, cacheSize) &&
	    Nvalues < (cacheSize - leftOver(Narray, cacheSize))) {
	    for (int i = 0; i < Nvalues; i++) *(*array+Narray+i) = values[i];
	} else {
	    A * array_; int i;
	    space = _getSpaceNeeded(Narray+Nvalues+1, cacheSize);
#ifdef DEBUG
    fprintf(stderr, "AddToArrayN allocated %d\n", space);
#endif
	    array_=new A[space];
	    for(i=0; i<Narray; i++) array_[i]=*(*array+i);
	    for(i=0; i<Nvalues; i++) array_[Narray+i]=values[i];
	    if(Narray) delete [] (*array);
	    (*array)=array_;
	}
};
//////////////////////////////////////////////////////////////////
template <class A> void ExtendArray(A ** array,
				    int NarrayOld,
				    int NarrayNew) {
	A * array_;
	array_=new A[NarrayNew];
	for(int i=0; i<NarrayOld; i++) array_[i]=*(*array+i);
	if(NarrayOld) delete [] (*array);
	(*array)=array_;
};
//////////////////////////////////////////////////////////////////
template <class A> void ExtendArrayValue(A ** array,
					 int NarrayOld,
					 int NarrayNew,
					 A value) {
	int i;
	A * array_;
	array_=new A[NarrayNew];
	for(i=0; i<NarrayOld; i++) array_[i]=*(*array+i);
	for(i=NarrayOld; i<NarrayNew; i++) array_[i]=value;
	if(NarrayOld) delete [] (*array);
	(*array)=array_;
};
//////////////////////////////////////////////////////////////////
template <class A> void DeleteArray(A ** array, int Narray) {
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};
//////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
template <class A> void DeleteArray(A ** array, int2 Narray) {
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};
#endif
//////////////////////////////////////////////////////////////////
DB * Property::db_default=NULL;
DB * Collection::db_default=NULL;
DB * Monomers::db_default=NULL;
DB * Entities::db_default=NULL;

#ifdef READ_WRITE
	const int FMODE=O_RDWR;
#else
	const int FMODE=O_RDONLY;
#endif

////////////////////////////////////////////////////////////////////
AminoacidProperty properties[]={
  "A", "ALA",  1,   71.09, 15.9,  25.9,  39.2,  23.1,  37.4, 1.33, 0.79,
  "V", "VAL",  2,   99.15, 47.7,   8.6,  38.5,  49.6,  19.6, 0.93, 1.50,
  "L", "LEU",  3,  113.17, 63.6,   0.0,  38.6,  57.6,  10.1, 1.20, 1.09,
  "I", "ILE",  4,  113.17, 63.6,   0.0,  39.2,  83.6,   7.5, 1.02, 1.44,
  "C", "CYS",  5,  103.15, 28.0,   7.4,  26.3,  40.3,   7.4, 0.98, 1.11,
  "M", "MET",  6,  131.21, 62.8,   4.9,  35.7,  44.3,   3.9, 1.25, 0.99,
  "P", "PRO",  7,   97.13, 41.0,  21.0,  40.2,  73.5,  66.2, 0.55, 0.73,
  "F", "PHE",  8,  147.19, 77.2,   1.2,  38.6,  76.1,   5.5, 1.00, 1.20,
  "Y", "TYR",  9,  163.19, 78.5,   9.9,  34.4,  70.8,  30.1, 0.75, 1.44,
  "W", "TRP", 10,  186.23, 100.0,  4.9,  37.7, 100.0,  13.8, 1.03, 1.17,
  "D", "ASP", 11,  115.10, 31.3, 100.0,   0.0,  17.5,  45.0, 1.04, 0.68,
  "N", "ASN", 12,  114.12, 35.4,  63.0,  31.3,   2.4,  46.1, 0.96, 0.80,
  "E", "GLU", 13,  129.13, 47.2,  93.8,   3.2,  17.8,  48.6, 1.41, 0.59,
  "Q", "GLN", 14,  128.15, 51.3,  45.7,  34.4,   0.0,  43.6, 1.07, 0.91,
  "H", "HIS", 15,  137.16, 49.2,  43.2,  59.2,  23.1,  28.1, 1.18, 0.94,
  "S", "SER", 16,   87.09, 18.1,  32.1,  34.8,   1.9,  40.5, 0.81, 1.00,
  "T", "THR", 17,  101.12, 34.0,  21.0,  45.7,   1.9,  35.3, 0.81, 1.15,
  "R", "ARG", 18,  156.20, 70.8,  51.9, 100.0,  22.6,  50.1, 0.97, 1.08,
  "K", "LYS", 19,  128.19, 68.0,  64.2,  86.9,  43.5,  54.3, 1.15, 0.81,
  "G", "GLY", 20,   57.07,  0.0,  37.0,  38.5,   2.7,  54.0, 0.67, 0.96};
///////////////////////////////////////////////////////////////////////////
int aa_decode[25]={1, 0, 5, 11, 13, 8, 20, 15, 4, 0, 19, 3, 6, 12, 0, 7, 14, 
		   18, 16, 17, 0, 2, 10, 0, 9};
char aa_code[21]={'A', 'V', 'L', 'I', 'C', 'M', 'P', 'F', 'Y', 'W',
		  'D', 'N', 'E', 'Q', 'H', 'S', 'T', 'R', 'K', 'G', 'X'};
///////////////////////////////////////////////////////////////////////////
DB::DB(char *path_) {
  // setting number of files - begin
  // this is important for sun platform
#ifndef T3E
  rlimit rlp;
  getrlimit(RLIMIT_NOFILE, &rlp);
  rlp.rlim_cur=600;
  setrlimit(RLIMIT_NOFILE, &rlp);
#endif
  // setting number of files - end

  nPath = 0;
  setPath(path_);
  nPath= path_ ? 1 : 0;
  setDefault();
  monomers=NULL;
}
////////////////////////////////////////////////////////////////////
DB::~DB() {
  if(nPath) DeleteArray(&path, nPath); 
}
////////////////////////////////////////////////////////////////////
void DB::setDefault() {
  Property::db_default=this;
  Collection::db_default=this;  
  Monomers::db_default=this;
  Entities::db_default=this;
}
////////////////////////////////////////////////////////////////////
void DB::addPath(char *path_) {
  int lPath=strlen(path_);
  if(path_) {
    AddBuffer(&path, path_, nPath); nPath++;
    if(*(path[nPath-1]+lPath-1)!=*PATH_SEPARATOR) 
      strcat(path[nPath-1], PATH_SEPARATOR);
  }
}
////////////////////////////////////////////////////////////////////
void DB::setPath(char *path_) {
  if(!path_) return;
  int lPath=strlen(path_);
  if(nPath) {
    DeleteArray(&path, nPath); nPath=0;
  }
  if(path_) {
    AddBuffer(&path, path_, nPath); nPath++;
    if(*(path[nPath-1]+lPath-1)!=*PATH_SEPARATOR) 
      addText(&path[nPath-1], PATH_SEPARATOR);
  }
}
////////////////////////////////////////////////////////////////////
char* DB::findPath(char *file_) {
  if(nPath==0) {
    printf("DB::findPath -- no db path is set\n");
    exit(0);
  }
  if(nPath==1) {
    return(path[0]);
  }
  else {
    for(int i=0; i<nPath; i++)
      if(testPathFile(path[i], file_)) return(path[i]);
  }
  printf("DB::findPath -- no path is found for %s in:\n", file_);
  for(int i=0; i<nPath; i++) printf("  (%d) '%s'\n", i+1, path[i]);
  
  exit(0);
}
////////////////////////////////////////////////////////////////////
int DB::openFile(char *file_) {
  char *buffer=NULL;
  int fd;

  addText(&buffer, findPath(file_));
  addText(&buffer, file_);
  if((fd=open(buffer, FMODE))==-1) {
    char *tmp = new char[500];
    tmp[0] = '\0';
    strcpy(tmp, "DB::openFile error - can't open ");
    strcat(tmp, buffer);
    perror(tmp);
    exit(0);
  }
  delete [] buffer;
  return(fd);
}
////////////////////////////////////////////////////////////////////
void DB::readFile(void *array, int4 nArray, int4 adr, int4 size, int fd) {
  testAdr(adr);
  
  int4 lbyte=size*nArray;
  
  lseek (fd, adr, SEEK_SET);
  read (fd, array, lbyte);
}
////////////////////////////////////////////////////////////////////
void DB::writeFile(void *array, int4 nArray, int4& adr, int4 size, int fd) {
  int4 lbyte=nArray*size; 
  
  if(lseek (fd, adr, SEEK_SET)==-1) {
    perror("DB::writeFile - lseek error");
    exit(0);
  }
  
  if(nArray>0)
    if(write (fd, array, lbyte)==-1) {
      perror("DB::writeFile - write error");
      exit(0);
    }
  addAdr(adr, lbyte, adr);
}
////////////////////////////////////////////////////////////////////
void DB::closeFile(int fd) {
  close(fd);
}
////////////////////////////////////////////////////////////////////
int DB::testFile(char *file_) {
  DIR *handleDir;
  dirent *d;
  
  if(path) {
    for(int i=0; i<nPath; i++) {
      char *searchDir=NULL;
      addText(&searchDir, path[i]);
      searchDir[strlen(searchDir)-1]='\0';
      handleDir=opendir(searchDir);
      if(searchDir) delete [] searchDir;
      
      while((d=readdir(handleDir)))
	if(!strcmp(d->d_name, file_)) {
	  closedir(handleDir); return(1);
	}
      closedir(handleDir);
    }
  }
  else {
    handleDir=opendir("."); 
    while((d=readdir(handleDir)))
      if(!strcmp(d->d_name, file_)) {
	closedir(handleDir); return(1);
      }
    closedir(handleDir);
  }
  return(0);
}
////////////////////////////////////////////////////////////////////
int DB::testPathFile(char *path_, char *file_) {
  DIR *handleDir;
  dirent *d;
  
  handleDir=opendir(path_);
  if(handleDir==NULL) return(0);
  while((d=readdir(handleDir)))
    if(!strcmp(d->d_name, file_)) {
      closedir(handleDir); return(1);
    }
  closedir(handleDir);
  return(0);
}
////////////////////////////////////////////////////////////////////
int DB::createFile(char *file_) {
  int fd=-1;
  char *buffer;

  setText(&buffer, path[nPath-1]);
  addText(&buffer, file_);
  
  if(FMODE!=O_RDWR) {
    printf("DB::createFile error - creation of files is not allowed\n");
    exit(0);
  }

  if((fd=open(buffer, O_RDWR | O_CREAT, 0666))==-1) {
    char *tmp = new char[500];
    tmp[0] = '\0';
    strcpy(tmp, "DB::createFile error - can't open ");
    strcat(tmp, buffer);
    perror(tmp);
    exit(0);
  }
  delete [] buffer;
  return(fd);
}
////////////////////////////////////////////////////////////////////
void DB::testAdr(int4 adr) {
  if (adr<0 || adr>2140000000l) {
    printf("DB::testAdr error - adr=%ld\n", adr);
    exit (0);
  }
}
////////////////////////////////////////////////////////////////////
void DB::addAdr(int4 adr1, int4 nBytes, int4& adr2) {
  testAdr(adr1);
  adr2=adr1+nBytes;
}
////////////////////////////////////////////////////////////////////
Property::Property(DB *db_) {
  db=db_;
  prop1=NULL;
#ifndef INT2_OFF
  prop2=NULL;
#endif
  prop4=NULL;
  propf=NULL;
  propd=NULL;
  name=NULL;
  writeFlag=0;
}
////////////////////////////////////////////////////////////////////
Property::Property(char *name_, DB *db_) {
  db=db_;
  prop1=NULL;
#ifndef INT2_OFF
  prop2=NULL;
#endif
  prop4=NULL;
  propf=NULL;
  propd=NULL;
  name=NULL;
  writeFlag=0;
  open(name_, db_);
}
////////////////////////////////////////////////////////////////////
Property::Property(char *name_, int load_in_memory_, DB *db_) {
  db=db_;
  prop1=NULL;
#ifndef INT2_OFF
  prop2=NULL;
#endif
  prop4=NULL;
  propf=NULL;
  propd=NULL;
  name=NULL;
  writeFlag=0;
  open(name_, load_in_memory_, db_);
}
////////////////////////////////////////////////////////////////////
Property::~Property() {
  close();
}
////////////////////////////////////////////////////////////////////
void Property::create(char *name_, DataType type_, int4 size_, 
		      int load_in_memory_, DB *db_, int cacheSize) {
  db=db_;
  if(!db) {
    printf("Property::create %s - no db\n", name_);
    exit(0);
  }

  if(isOpen()) close();
  
  prop1=NULL;
#ifndef INT2_OFF
  prop2=NULL;
#endif
  prop4=NULL;
  propf=NULL;
  propd=NULL;
  writeFlag=0;

  if(name_) {
    if(name) delete [] name; 
    name=new char [strlen(name_)+1];
    strcpy(name, name_);
  }
  else {
      printf("Property::create - name_ == NULL: %s\n", this->name);
      exit(0);
  }
  
  if(!db->testFile(name)) {
    fd=db->createFile(name);
  }
  else {
    fd=db->openFile(name);
  }

  nItems=0;                     // 0
  nOptItems=0;                    // 1
  type=type_;                     // 2
  size=size_;                     // 3
  adrProp=(8+nOptItems)*4;        // 4
  nProp=0;                        // 5
  adrRef=(8+nOptItems)*4;         // 6
  load_in_memory=load_in_memory_; // 7
  
  if(!size)
    {
      ref=new int4 [_getSpaceNeeded(nItems+1, cacheSize)];
      ref[nItems]=nProp;
    }

  saveHeader(1);
}
////////////////////////////////////////////////////////////////////
void Property::clear() {
  char *name_tmp=new char [strlen(name)+1];
  DataType type_tmp=type;
  int4 size_tmp=size;
  int load_in_memory_tmp=load_in_memory_orig;
  strcpy(name_tmp, name);
  
  close();
  create(name_tmp, type_tmp, size_tmp, load_in_memory_tmp);
  delete [] name_tmp;
}
////////////////////////////////////////////////////////////////////
void Property::open(char *name_, DB *db_) {
  open(name_, -1, db_, DEFAULT_CACHE);
}
////////////////////////////////////////////////////////////////////
void Property::open(char *name_, int load_in_memory_, DB *db_, int cacheSize) {
  int space;

  db=db_;
  if(!db) {
    printf("Property::open %s - no db\n", name_);
    exit(0);
  }

  close();
  
  if(name_) {
    if(name) delete [] name; 
    name=new char [strlen(name_)+1];
    strcpy(name, name_);
  }
  int4 hdr[8], adr=0;
  fd=db->openFile(name);
  db->readFile(hdr, 8, adr, 4, fd);

  nItems=hdr[0]; 
  nOptItems=hdr[1]; 
  type=(DataType)hdr[2]; 
  size=hdr[3]; 
  adrProp=hdr[4];
  nProp=hdr[5]; 
  adrRef=hdr[6];
  load_in_memory=(load_in_memory_==-1?hdr[7]:load_in_memory_);
  load_in_memory_orig = hdr[7];
  writeFlag=0;

  /* fprintf(stderr, "%-20s %8d %4d %4d %8d %1d\n", name, nItems, type, size, nProp, hdr[7]);
  fprintf(stderr, "%-20s %8d %4d %4d %8d %1d\n", name, nItems, type, size, nProp, load_in_memory);
  fflush(stderr); */

  if(load_in_memory && nProp) {

    space = _getSpaceNeeded(nProp, cacheSize);

#ifdef DEBUG
	fprintf(stderr, "opened %s with %d allocated space\n", name, space);
#endif

    if(type==D_INT1) {
      prop1=new int1 [space];
      db->readFile(prop1, nProp, adrProp, 1, fd);
    }
      
#ifndef INT2_OFF
    if(type==D_INT2) {
      prop2=new int2 [space];
      db->readFile(prop2, nProp, adrProp, 2, fd);
    }
#endif   
 
    if(type==D_INT4) {
      prop4=new int4 [space];
      db->readFile(prop4, nProp, adrProp, 4, fd);
    }
    
    if(type==D_FLT4) {
      propf=new flt4 [space];
      db->readFile(propf, nProp, adrProp, 4, fd);
    }

    if(type==D_FLT8) {
      propd=new flt8 [space];
      db->readFile(propf, nProp, adrProp, 8, fd);
    }
  }

  if(!size) {
    space = _getSpaceNeeded(nItems+1, cacheSize);
    ref=new int4 [space];
    db->readFile(ref, nItems+1, adrRef, 4, fd);
  }
  
}
////////////////////////////////////////////////////////////////////
int Property::isOpen() {
  return(name==NULL?0:1);
}
////////////////////////////////////////////////////////////////////
void Property::saveHeader(int with_load_mode) {
  int4 hdr[8];
  hdr[0]=nItems; 
  hdr[1]=nOptItems; 
  hdr[2]=type; 
  hdr[3]=size; 
  hdr[4]=adrProp; 
  hdr[5]=nProp; 
  hdr[6]=adrRef;
  hdr[7]=load_in_memory;
  writeFlag=0;
  
  int4 adrHdr=0, adrRef_=adrRef;
  
  db->writeFile(hdr, with_load_mode?8:7, adrHdr, 4, fd);
  if(!size) db->writeFile(ref, nItems+1, adrRef_, 4, fd);
}
////////////////////////////////////////////////////////////////////
void Property::close(int export_array) {
  if(name) {

#ifdef DEBUG
    fprintf(stderr, "%s closed\n", name);
    fflush(stderr);
#endif

    if(writeFlag) saveHeader();
    db->closeFile(fd);
    delete [] name;
    name=NULL;
    if(!size) delete [] ref;
    if(!nItems) return;
    if(nProp && !export_array) {
      if(type==D_INT1 && prop1) {
	delete [] prop1; prop1=NULL;
      }
      
#ifndef INT2_OFF
      if(type==D_INT2 && prop2) {
	delete [] prop2; prop2=NULL;
      }
#endif
	
      if(type==D_INT4 && prop4) {
	delete [] prop4; prop4=NULL;
      }
      
      if(type==D_FLT4 && propf) {
	delete [] propf; propf=NULL;
      }

      if(type==D_FLT8 && propd) {
	delete [] propd; propd=NULL;
      }
    }
  }
}
////////////////////////////////////////////////////////////////////
int4 Property::getItemSize(int4 index) {
	return (size ? size : ref[index+1] - ref[index]);
}
////////////////////////////////////////////////////////////////////
void Property::updateRef(int4 index, int4 sizeDiff) {
	for (int i = index; i < nItems + 1; i++) {
		ref[i] = ref[i] + sizeDiff;
	}
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int1 item) {
	setItem(index, &item, 1);
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::setItem(int4 index, int2 item) {
	setItem(index, &item, 1);
}
#endif
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int4 item) {
	setItem(index, &item, 1);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt4 item) {
	setItem(index, &item, 1);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt8 item) {
	setItem(index, &item, 1);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int1 *item) {
  int nValues;
  if (!size) {
    nValues = (item ? strlen(item) : 0) + 1;
  } 
  else {
    nValues=size;
  }
  setItem(index, item, nValues);
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::setItem(int4 index, int2 *item) {
	if (!size) {
		printf("Property::setItem - %s size is 0\n", name);
		exit(0);
	}
	setItem(index, item, size);
}
#endif
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int4 *item) {
	if (!size) {
		printf("Property::setItem - %s size is 0\n", name);
		exit(0);
	}
	setItem(index, item, size);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt4 *item) {
	if (!size) {
		printf("Property::setItem - %s size is 0\n", name);
		exit(0);
	}
	setItem(index, item, size);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt8 *item) {
	if (!size) {
		printf("Property::setItem - %s size is 0\n", name);
		exit(0);
	}
	setItem(index, item, size);
}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int1 *item, int4 nValues) {
	if (type != D_INT1) {
	   printf("Property::setItem - %s type mismatch %d vs %d\n", 
		  name, type, D_INT1);
	   exit(0);
	}

	if (!size && !load_in_memory) {
	   printf("Property::setItem - %s is not editable\n", name);
	   exit(0);
	}

	if (size && size != nValues) {
		printf("Property::setItem - %s size mismatch: %d %d\n",
		       name, size, nValues);
		exit(0);
	}

	if (index > nItems) {
	   printf("Property::setItem - %s Array index out of bounds: %d\n", name, index);
	   exit(0);
	}

	int4 nPrePro  = (!size ? ref[index] : index*size);
	int4 sizeDiff = (!size ? (nValues-(ref[index+1]-ref[index])) : 0);

	if (size && load_in_memory) {
		for (int i = 0; i < nValues; i++) prop1[nPrePro+i] = item[i];
	} else {
		InsertItems(&prop1, nPrePro, item, nValues, &nProp, sizeDiff);
	}

	int4 addr = adrProp + nPrePro;

	if (!size && sizeDiff) {
		updateRef(index+1, sizeDiff);
		adrRef += sizeDiff;
		db->writeFile(&prop1[nPrePro], nProp-nPrePro, addr, 1, fd);
		writeFlag = 1;
	} else {
		db->writeFile(item, nValues, addr, 1, fd);
	}

}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::setItem(int4 index, int2 *item, int4 nValues) {
	if (type != D_INT2) {
	   printf("Property::setItem - %s type mismatch %d vs %d\n", 
		  name, type, D_INT2);
	   exit(0);
	}

	if (!size && !load_in_memory) {
	   printf("Property::setItem - %s is not editable\n", name);
	   exit(0);
	}

	if (size && size != nValues) {
		printf("Property::setItem - %s size mismatch: %d %d\n",
		       name, size, nValues);
		exit(0);
	}

	if (index > nItems) {
	   printf("Property::setItem - %s Array index out of bounds: %d\n", name, index);
	   exit(0);
	}

	int4 nPrePro  = (!size ? ref[index] : index*size);
	int4 sizeDiff = (!size ? (nValues-(ref[index+1]-ref[index])) : 0);

	if ((size && load_in_memory) || (!size && !sizeDiff)) {
		for (int i = 0; i < nValues; i++) prop2[nPrePro+i] = item[i];
	} else if (!size && sizeDiff) {
		InsertItems(&prop2, nPrePro, item, nValues, &nProp, sizeDiff);
	}

	int4 addr = adrProp + nPrePro*2;

	if (!size && sizeDiff) {
		updateRef(index+1, sizeDiff);
		adrRef += sizeDiff*2;
		db->writeFile(&prop2[nPrePro], nProp-nPrePro, addr, 2, fd);
		writeFlag = 1;
	} else {
		db->writeFile(item, nValues, addr, 2, fd);
	}

}
#endif
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, int4 *item, int4 nValues) {
	if (type != D_INT4) {
	   printf("Property::setItem - %s type mismatch %d vs %d\n", 
		  name, type, D_INT4);
	   exit(0);
	}

	if (!size && !load_in_memory) {
	   printf("Property::setItem - %s is not editable\n", name);
	   exit(0);
	}

	if (size && size != nValues) {
		printf("Property::setItem - %s size mismatch: %d %d\n",
		       name, size, nValues);
		exit(0);
	}

	if (index > nItems) {
	   printf("Property::setItem - %s Array index out of bounds: %d\n", name, index);
	   exit(0);
	}

	int4 nPrePro  = (!size ? ref[index] : index*size);
	int4 sizeDiff = (!size ? (nValues-(ref[index+1]-ref[index])) : 0);

#ifdef DEBUG
	if (sizeDiff < 0) {
		int i;
		fprintf(stderr, "shrinking %d in %s\n", index, this->name); 
		fprintf(stderr, "was: ");
		for (i = 0; i < (ref[index+1]-ref[index]); i++)
			fprintf(stderr, "%d, ", prop4[i]);
		fprintf(stderr, "\b\b\nis: ");
		for (i = 0; i < nValues; i++)
			fprintf(stderr, "%d, ", item[i]);
		fprintf(stderr, "\b\b\n");
	}
#endif

	if (size && load_in_memory) {
		for (int i = 0; i < nValues; i++) prop4[nPrePro+i] = item[i];
	} else {
		InsertItems(&prop4, nPrePro, item, nValues, &nProp, sizeDiff);
	}

	int4 addr = adrProp + nPrePro*4;

	if (!size && sizeDiff) {
		updateRef(index+1, sizeDiff);
		adrRef += sizeDiff*4;
		db->writeFile(&prop4[nPrePro], nProp-nPrePro, addr, 4, fd);
		writeFlag = 1;
	} else {
		db->writeFile(item, nValues, addr, 4, fd);
	}

}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt4 *item, int4 nValues) {
	if (type != D_FLT4) {
	   printf("Property::setItem - %s type mismatch %d vs %d\n", 
		  name, type, D_FLT4);
	   exit(0);
	}

	if (!size && !load_in_memory) {
	   printf("Property::setItem - %s is not editable\n", name);
	   exit(0);
	}

	if (size && size != nValues) {
		printf("Property::setItem - %s size mismatch: %d %d\n",
		       name, size, nValues);
		exit(0);
	}

	if (index > nItems) {
	   printf("Property::setItem - %s Array index out of bounds: %d\n", 
		  name, index);
	   exit(0);
	}

	int4 nPrePro  = (!size ? ref[index] : index*size);
	int4 sizeDiff = (!size ? (nValues-(ref[index+1]-ref[index])) : 0);

	if (size && load_in_memory) {
		for (int i = 0; i < nValues; i++) propf[nPrePro+i] = item[i];
	} else {
		InsertItems(&propf, nPrePro, item, nValues, &nProp, sizeDiff);
	}

	int4 addr = adrProp + nPrePro*4;

	if (!size && sizeDiff) {
		updateRef(index+1, sizeDiff);
		adrRef += sizeDiff*4;
		db->writeFile(&propf[nPrePro], nProp-nPrePro, addr, 4, fd);
		writeFlag = 1;
	} else {
		db->writeFile(item, nValues, addr, 4, fd);
	}

}
////////////////////////////////////////////////////////////////////
void Property::setItem(int4 index, flt8 *item, int4 nValues) {
	if (type != D_FLT8) {
	   printf("Property::setItem - %s type mismatch %d vs %d\n", 
		  name, type, D_FLT8);
	   exit(0);
	}

	if (!size && !load_in_memory) {
	   printf("Property::setItem - %s is not editable\n", name);
	   exit(0);
	}

	if (size && size != nValues) {
		printf("Property::setItem - %s size mismatch: %d %d\n",
		       name, size, nValues);
		exit(0);
	}

	if (index > nItems) {
	   printf("Property::setItem - %s Array index out of bounds: %d\n", name,
									index);
	   exit(0);
	}

	int4 nPrePro  = (!size ? ref[index] : index*size);
	int4 sizeDiff = (!size ? (nValues-(ref[index+1]-ref[index])) : 0);

	if (size && load_in_memory) {
		for (int i = 0; i < nValues; i++) propd[nPrePro+i] = item[i];
	} else {
		InsertItems(&propd, nPrePro, item, nValues, &nProp, sizeDiff);
	}

	int4 addr = adrProp + nPrePro*8;

	if (!size && sizeDiff) {
		updateRef(index+1, sizeDiff);
		adrRef += sizeDiff*8;
		db->writeFile(&propd[nPrePro], nProp-nPrePro, addr, 8, fd);
		writeFlag = 1;
	} else {
		db->writeFile(item, nValues, addr, 8, fd);
	}

}
////////////////////////////////////////////////////////////////////
void Property::extendItem(int4 index, int1* item) {
  if(item) {
    int nItems;
    if(!size) {
      nItems=(item?strlen(item):0)+1;
    }
    else {
      nItems=size;
    }
    extendItem(index, item, nItems);
  }
  else {
    extendItem(index, "", 1);
  }
}
////////////////////////////////////////////////////////////////////
void Property::extendItem(int4 index, int1* item, int4 nValues) {
  int4 nValues_;
  int1 *item_=item1n(index, nValues_);
  int1 *_item=new int1 [nValues_+nValues];
  int i;
  for(i=0; i<nValues_; i++) _item[i]=item_[i];
  for(i=0; i<nValues; i++) _item[nValues_+i]=item[i];
  setItem(index, _item, nValues_+nValues);
  delete [] _item;
}
////////////////////////////////////////////////////////////////////
void Property::extendItem(int4 index, int4 item) {
  extendItem(index, &item, 1);
}
////////////////////////////////////////////////////////////////////
void Property::extendItem(int4 index, int4* item, int4 nValues) {
  int4 nValues_;
  int4 *item_=item4n(index, nValues_);
  int4 *_item=new int4 [nValues_+nValues];
  int i;
  for(i=0; i<nValues_; i++) _item[i]=item_[i];
  for(i=0; i<nValues; i++) _item[nValues_+i]=item[i];
  setItem(index, _item, nValues_+nValues);
  delete [] _item;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(int1 item) {
  if(type!=D_INT1) {
    printf("Property::addItem - %s type mismatch %d vs %d\n", name, type, D_INT1); 
    exit(0);
  }
  
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  
  if(load_in_memory) AddToArray(&prop1, item, nProp);
  db->writeFile(&item, 1, adrRef, 1, fd);
  nProp++;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::addItem(int2 item) {
  if(type!=D_INT2) {
    printf("Property::addItem - %s type mismatch %d vs %d\n", name, type, D_INT2); 
    exit(0);
  }
  
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  
  if(load_in_memory) AddToArray(&prop2, item, nProp);
  db->writeFile(&item, 1, adrRef, 2, fd);
  nProp++;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
#endif
////////////////////////////////////////////////////////////////////
void Property::addItem(int4 item) {
  if(type!=D_INT4) {
    printf("Property::addItem - %s type mismatch %d vs %d\n", name, type, D_INT4); 
    exit(0);
  }
  
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  
  if(load_in_memory) AddToArray(&prop4, item, nProp);
  db->writeFile(&item, 1, adrRef, 4, fd);
  nProp++;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt4 item) {
  if(type!=D_FLT4) {
    printf("Property::addItem - %s type mismatch %d vs %d\n", name, type, D_FLT4); 
    exit(0);
  }
  
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  
  if(load_in_memory) AddToArray(&propf, item, nProp);
  db->writeFile(&item, 1, adrRef, 4, fd);
  nProp++;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt8 item) {
  if(type!=D_FLT8) {
    printf("Property::addItem - %s type mismatch %d vs %d\n", name, type, D_FLT8); 
    exit(0);
  }
  
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  
  if(load_in_memory) AddToArray(&propd, item, nProp);
  db->writeFile(&item, 1, adrRef, 8, fd);
  nProp++;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(int1* item) {
  if(item) {
    int nItems;
    if(!size) {
      nItems=(item?strlen(item):0)+1;
    }
    else {
      nItems=size;
    }
    addItem(item, nItems);
  }
  else {
    addItem("", 1);
  }
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::addItem(int2* item) {
  if(!size) {
    printf("Property::addItem - %s size is 0\n", name); 
    exit(0);
  }
  addItem(item, size);
}
#endif
////////////////////////////////////////////////////////////////////
void Property::addItem(int4* item) {
  if(!size) {
    printf("Property::addItem - %s size is 0\n", name); 
    exit(0);
  }
  addItem(item, size);
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt4* item) {
  if(!size) {
    printf("Property::addItem - %s size is 0\n", name); 
    exit(0);
  }
  addItem(item, size);
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt8* item) {
  if(!size) {
    printf("Property::addItem - %s size is 0\n", name); 
    exit(0);
  }
  addItem(item, size);
}
////////////////////////////////////////////////////////////////////
void Property::addItem(int1* item, int4 nValues) {
  if(!size) AddToArray(&ref, nProp, nItems+1);

  nItems++;
  
  if(load_in_memory) AddToArrayN(&prop1, item, nProp, nValues);
  db->writeFile(item, nValues, adrRef, 1, fd);
  nProp+=nValues;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::addItem(int2* item, int4 nValues) {
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;

  if(load_in_memory) AddToArrayN(&prop2, item, nProp, nValues);
  db->writeFile(item, nValues, adrRef, 2, fd);
  nProp+=nValues;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
#endif
////////////////////////////////////////////////////////////////////
void Property::addItem(int4* item, int4 nValues) {
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;

  if(load_in_memory) AddToArrayN(&prop4, item, nProp, nValues);
  db->writeFile(item, nValues, adrRef, 4, fd);
  nProp+=nValues;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt4* item, int4 nValues) {
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;

  if(load_in_memory) AddToArrayN(&propf, item, nProp, nValues);
  db->writeFile(item, nValues, adrRef, 4, fd);
  nProp+=nValues;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem(flt8* item, int4 nValues) {
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;

  if(load_in_memory) AddToArrayN(&propd, item, nProp, nValues);
  db->writeFile(item, nValues, adrRef, 8, fd);
  nProp+=nValues;
  if(!size) ref[nItems]=nProp;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
void Property::addItem() {
  if(!size) AddToArray(&ref, nProp, nItems+1);
  nItems++;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
int1* Property::item1(int4 index, int copy_array) {
  if(type!=D_INT1) {
    printf("Property::item1 - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  int i_tmp=size?size*index:ref[index];
  if(load_in_memory && !copy_array) return(&prop1[i_tmp]);
  
  int i, l_tmp=size?size:ref[index+1]-ref[index];
  int1 *tmp=new int1 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp, 1, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop1[i_tmp+i];
    
  if(!load_in_memory && !copy_array) {
    if(prop1) delete [] prop1;
    prop1=tmp;
  }

  return(tmp);
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int2* Property::item2(int4 index, int copy_array) {
  if(type!=D_INT2) {
    printf("Property::item2 - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  int i_tmp=size?size*index:ref[index];
  if(load_in_memory && !copy_array) return(&prop2[i_tmp]);
  
  int i, l_tmp=size?size:ref[index+1]-ref[index]; 
  int2 *tmp=new int2 [l_tmp];
  
  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*2, 2, fd);
  
  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop2[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(prop2) delete [] prop2;
    prop2=tmp;
  }
  
  return(tmp);
}
#endif
////////////////////////////////////////////////////////////////////
int4* Property::item4(int4 index, int copy_array) {
  if(type!=D_INT4) {
    printf("Property::item4 - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  int i_tmp=size?size*index:ref[index];
  if(load_in_memory && !copy_array) return(&prop4[i_tmp]);
  
  int i, l_tmp=size?size:ref[index+1]-ref[index]; 
  int4 *tmp=new int4 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*4, 4, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop4[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(prop4) delete [] prop4;
    prop4=tmp;
  }
  
  return(tmp);
}
////////////////////////////////////////////////////////////////////
flt4* Property::itemf(int4 index, int copy_array) {
  if(type!=D_FLT4) {
    printf("Property::itemf - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  int i_tmp=size?size*index:ref[index];
  if(load_in_memory && !copy_array) return(&propf[i_tmp]);
  
  int i, l_tmp=size?size:ref[index+1]-ref[index];
  flt4 *tmp=new flt4 [l_tmp];
  
  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*4, 4, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=propf[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(propf) delete [] propf;
      propf=tmp;
  }

  return(tmp);
}
////////////////////////////////////////////////////////////////////
flt8* Property::itemd(int4 index, int copy_array) {
  if(type!=D_FLT8) {
    printf("Property::itemf - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  int i_tmp=size?size*index:ref[index];
  if(load_in_memory && !copy_array) return(&propd[i_tmp]);
  
  int i, l_tmp=size?size:ref[index+1]-ref[index];
  flt8 *tmp=new flt8 [l_tmp];
  
  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*8, 8, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=propd[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(propd) delete [] propd;
      propd=tmp;
  }

  return(tmp);
}
////////////////////////////////////////////////////////////////////
int1* Property::item1n(int4 index, int4& l_tmp, int copy_array) {
  if(type!=D_INT1) {
    printf("Property::item1 - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  l_tmp=size?size:ref[index+1]-ref[index]; 
  int i_tmp=size?size*index:ref[index];

  if(load_in_memory && !copy_array) return(&prop1[i_tmp]);

  int i; 
  int1 *tmp=new int1 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp, 1, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop1[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(prop1) delete [] prop1;
    prop1=tmp;
  }
  
  return(tmp);
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int2* Property::item2n(int4 index, int4& l_tmp, int copy_array) {
  if(type!=D_INT2) {
    printf("Property::item2 - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  l_tmp=size?size:ref[index+1]-ref[index]; 
  int i_tmp=size?size*index:ref[index];
  
  if(load_in_memory && !copy_array) return(&prop2[i_tmp]);
  
  int i;
  int2 *tmp=new int2 [l_tmp];
  
  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*2, 2, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop2[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(prop2) delete [] prop2;
    prop2=tmp;
  }
  
  return(tmp);
}
#endif
////////////////////////////////////////////////////////////////////
int4* Property::item4n(int4 index, int4& l_tmp, int copy_array) {
  if(type!=D_INT4) {
    printf("Property::item4n - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  l_tmp=size?size:ref[index+1]-ref[index];
  int i_tmp=size?size*index:ref[index];

  if(load_in_memory && !copy_array) return(&prop4[i_tmp]);

  int i; 
  int4 *tmp=new int4 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*4, 4, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=prop4[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(prop4) delete [] prop4;
    prop4=tmp;
  }
  
  return(tmp);
}
////////////////////////////////////////////////////////////////////
flt4* Property::itemfn(int4 index, int4& l_tmp, int copy_array) {
  if(type!=D_FLT4){
    printf("Property::itemf - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  l_tmp=size?size:ref[index+1]-ref[index];
  int i_tmp=size?size*index:ref[index];

  if(load_in_memory && !copy_array) return(&propf[i_tmp]);

  int i; 
  flt4 *tmp=new flt4 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*4, 4, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=propf[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(propf) delete [] propf;
    propf=tmp;
  }

  return(tmp);
}
////////////////////////////////////////////////////////////////////
flt8* Property::itemdn(int4 index, int4& l_tmp, int copy_array) {
  if(type!=D_FLT8){
    printf("Property::itemf - '%s' access error size=%d type=%d\n", 
	   name, size, type);
    exit(0);
  }
  l_tmp=size?size:ref[index+1]-ref[index];
  int i_tmp=size?size*index:ref[index];

  if(load_in_memory && !copy_array) return(&propd[i_tmp]);

  int i; 
  flt8 *tmp=new flt8 [l_tmp];

  if(!load_in_memory) db->readFile(tmp, l_tmp, adrProp+i_tmp*8, 8, fd);

  if(load_in_memory && copy_array) 
    for(i=0; i<l_tmp; i++) tmp[i]=propd[i_tmp+i];

  if(!load_in_memory && !copy_array) {
    if(propd) delete [] propd;
    propd=tmp;
  }

  return(tmp);
}
////////////////////////////////////////////////////////////////////
void Property::array(int1 *prop) {
  prop=prop1;
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Property::array(int2 *prop) {
  prop=prop2;
}
#endif
////////////////////////////////////////////////////////////////////
void Property::array(int4 *prop) {
  prop=prop4;
}
////////////////////////////////////////////////////////////////////
void Property::array(flt4 *prop) {
  prop=propf;
}
////////////////////////////////////////////////////////////////////
void Property::array(flt8 *prop) {
  prop=propd;
}
////////////////////////////////////////////////////////////////////
int4 Property::find(char *text, int mismatches) {
	if(type != D_INT1 || !load_in_memory) {
	   printf("Property::find - '%s' error type=%d load_in_memory=%d\n",
		  name, type, load_in_memory);
	   exit(0);
	}

	int iObj = 0, j, k, l = strlen(text), mm, kk;

	if (!size) {
		for(; iObj<nItems; iObj++) {
		  for(j=ref[iObj]; j<ref[iObj+1]-l; j++) {
		      mm=0; kk=0;
			for(k=j; k<j+l; k++, kk++) {
			  if(prop1[k]!=text[kk]) {
			      mm++; if(mm>mismatches) goto no_match1;
			  }
			}
		      return(iObj);
		    no_match1: ;
		  }
		}
	} else {
		for(; iObj<nItems; iObj++) {
		  for(j=iObj*size; j<(iObj+1)*size-l; j++) {
		      mm=0; kk=0;
			for(k=j; k<j+l; k++, kk++) {
			  if(prop1[k]!=text[kk]) {
			      mm++; if(mm>mismatches) goto no_match2;
			  }
			}
		      return(iObj);
		    no_match2: ;
		  }
		}
	}
	return(-1);
}
////////////////////////////////////////////////////////////////////
int4 Property::find(char *text, int mismatches, int4& pos, int& next) {
  if(type!=D_INT1 || !load_in_memory) {
    printf("Property::find - '%s' access error type=%d load_in_memory=%d\n", 
	   name, type, load_in_memory);
    exit(0);
  }

  if(!next) {
    iObj=0; iPos=0; next=1;
  }
  else {
    if(next==1) iPos++;
    if(next==2) {
      iObj++; iPos=0;
    }
  }
  int j, k, l=strlen(text), mm, kk;

  if(!size) {
    for(; iObj<nItems; iObj++) {
      for(j=ref[iObj]+iPos; j<ref[iObj+1]-l; j++) {
	mm=0; kk=0;
	for(k=j; k<j+l; k++, kk++) {
	  if(prop1[k]!=text[kk]) {
	    mm++; if(mm>mismatches) goto no_match1;
	  }
	}
	iPos=j-ref[iObj];
	pos=iPos;
	return(iObj);
      no_match1: ;
      }
      iPos=0;
    }
  }
  else {
    for(; iObj<nItems; iObj++) {
      for(j=iObj*size+iPos; j<(iObj+1)*size-l; j++) {
	mm=0; kk=0;
	for(k=j; k<j+l; k++, kk++) {
	  if(prop1[k]!=text[kk]) {
	    mm++; if(mm>mismatches) goto no_match2;
	  }
	}
	iPos=j-iObj*size;
	pos=iPos;
	return(iObj);
      no_match2: ;
      }
      iPos=0;
    }
  }
  return(-1);
}
////////////////////////////////////////////////////////////////////
int4 Property::find(int1 t1, int1 t2, int& next) {
  return(find((int4)t1, (int4)t2, next));
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int4 Property::find(int2 t1, int2 t2, int& next) {
  return(find((int4)t1, (int4)t2, next));
}
#endif
////////////////////////////////////////////////////////////////////
int4 Property::find(int4 t1, int4 t2, int& next) {
  if(size!=1  || (type!=D_INT1 && type!=D_INT2 && type!=D_INT4) || 
     !load_in_memory) {
    printf("Property::find - '%s' access error size=%d type=%d load_in_memory=%d\n", name, size, type, load_in_memory);
    exit(0);
  }

  if(!next) {
    iObj=0; next=1;
  }
  else {
    iObj++;
  }
  
  if(type==D_INT1) {
    for(; iObj<nItems; iObj++) {
      if(t2==-1) {
	if(prop1[iObj]==t1) return(iObj);
      }
      else {
	if(prop1[iObj]>=t1 && prop1[iObj]<=t2) return(iObj);
      }
    }
    return(-1);
  }

#ifndef INT2_OFF
  if(type==D_INT2) {
    for(; iObj<nItems; iObj++) {
      if(t2==-1) {
	if(prop2[iObj]==t1) return(iObj);
      }
      else {
	if(prop2[iObj]>=t1 && prop2[iObj]<=t2) return(iObj);
      }
    }
    return(-1);
  }
#endif

  if(type==D_INT4) {
    for(; iObj<nItems; iObj++) {
      if(t2==-1) {
	if(prop4[iObj]==t1) return(iObj);
      }
      else {
	if(prop4[iObj]>=t1 && prop4[iObj]<=t2) return(iObj);
      }
    }
    return(-1);
  }
}
////////////////////////////////////////////////////////////////////
int4 Property::find(flt4 t1, flt4 t2, int& next) {
  if(size!=1 || type!=D_FLT4 || !load_in_memory) {
    printf("Property::find - '%s' access error size=%d type=%d load_in_memory=%d\n", name, size, type, load_in_memory);
    exit(0);
  }

  if(!next) {
    iObj=0; next=1;
  }
  else {
    iObj++;
  }

  for(; iObj<nItems; iObj++) {
    if(t2==-1) {
      if(propf[iObj]==t1) return(iObj);
    }
    else {
      if(propf[iObj]>=t1 && propf[iObj]<=t2) return(iObj);
    }
  }
  return(-1);

}
////////////////////////////////////////////////////////////////////
Collection::Collection(DB *db_) {
  db=db_;
#ifndef INT2_OFF
  coll2=NULL;
#endif
  coll4=NULL;
#ifndef INT2_OFF
  coll2t=NULL;
#endif
  coll4t=NULL;
  name=NULL;
  writeFlag=0;
}
////////////////////////////////////////////////////////////////////
Collection::Collection(char *name_, DB *db_) {
  db=db_;
#ifndef INT2_OFF
  coll2=NULL;
#endif
  coll4=NULL;
#ifndef INT2_OFF
  coll2t=NULL;
#endif
  coll4t=NULL;
  name[0]=NULL;
  writeFlag=0;
  open(name_, db_);
}
////////////////////////////////////////////////////////////////////
Collection::~Collection() {
  close();
}
////////////////////////////////////////////////////////////////////
void Collection::create(char *name_, DataType type_, int4 nClasses_, DB *db_) {
  db=db_;
  if(!db) {
    printf("Collection::create %s - no db\n", name_);
    exit(0);
  }

  close();

  if(name_) {
    if(name) delete [] name;
    name=new char [strlen(name_)+1];
    strcpy(name, name_);
  }
  else {
    printf("Collection::create - name_ == NULL: %s\n", this->name);
    exit(0);
  }

  if(!db->testFile(name)) fd=db->createFile(name);
  else fd=db->openFile(name);

  nClasses=nClasses_;             // 0
  type=type_;                     // 1
  adrColl=16;                     // 2
  adrNColl=16;                    // 3
  writeFlag=0;

  nColl=new int4 [nClasses];
  int i;
  for(i=0; i<nClasses; i++) nColl[i]=0;

#ifndef INT2_OFF
  coll2t=NULL;
#endif
  coll4t=NULL;
  index_=-1;
  load_in_memory=1;
#ifndef INT2_OFF
  if(type==D_INT2) {
    coll2=new int2* [nClasses];  
    for(i=0; i<nClasses; i++) coll2[i]=NULL;
  }
#endif
  if(type==D_INT4) {
    coll4=new int4* [nClasses];
    for(i=0; i<nClasses; i++) coll4[i]=NULL;
  }

  save();
}
////////////////////////////////////////////////////////////////////
void Collection::clear() {
  char *name_tmp=new char [strlen(name)+1];
  DataType type_tmp=type;
  int4 nClasses_tmp=nClasses;
  strcpy(name_tmp, name);

  close();
  create(name_tmp, type_tmp, nClasses_tmp);
  delete [] name_tmp;
}
////////////////////////////////////////////////////////////////////
void Collection::open(char *name_, DB *db_) {
  open(name_, 1);
}
////////////////////////////////////////////////////////////////////
void Collection::open(char *name_, int load_in_memory_, DB *db_, 
		      int cacheSize) {
  db=db_;
  if(!db) {
    printf("Collection::open %s - no db\n", name_);
    exit(0);
  }
  
  close();
  
  if(name_) {
    if(name) delete [] name;
    name=new char [strlen(name_)+1];
    strcpy(name, name_);
  }
  int4 hdr[4], adr=0;
  fd=db->openFile(name);
  db->readFile(hdr, 4, adr, 4, fd);
  load_in_memory=load_in_memory_;
    
  nClasses=hdr[0];                // 0
  type=(DataType)hdr[1];          // 1
  adrColl=hdr[2];                 // 2
  adrNColl=hdr[3];                // 3

  writeFlag=0;
  nColl=new int4 [nClasses];
  db->readFile(nColl, nClasses, adrNColl, 4, fd);

  int i;
  adr=adrColl;
  
  if(load_in_memory) {
#ifndef INT2_OFF
    if(type==D_INT2) coll2=new int2* [nClasses];
#endif
    if(type==D_INT4) coll4=new int4* [nClasses];
    for(i=0; i<nClasses; i++) {
      if(nColl[i]) {
#ifndef INT2_OFF
	if(type==D_INT2) {
	  coll2[i]=new int2 [nColl[i]+cacheSize];
	  db->readFile(coll2[i], nColl[i], adr, 2, fd);
	  db->addAdr(adr, nColl[i]*2, adr);
	}
#endif
	if(type==D_INT4) {
	  coll4[i]=new int4 [nColl[i]+cacheSize];
	  db->readFile(coll4[i], nColl[i], adr, 4, fd);
	  db->addAdr(adr, nColl[i]*4, adr);
	}
      }
      else {
#ifndef INT2_OFF
	if(type==D_INT2) coll2[i]=NULL;
#endif
	if(type==D_INT4) coll4[i]=NULL;
      }
    }
  }
  else {
    int sum=0;
    nCollt=new int4[nClasses];
    for(i=0; i<nClasses; i++) {
      nCollt[i]=sum;
      sum+=nColl[i];
    }
#ifndef INT2_OFF
    coll2t=NULL;
#endif
    coll4t=NULL;
    index_=-1;
  }

}
////////////////////////////////////////////////////////////////////
int Collection::isOpen() {
  return(name!=NULL);
}
////////////////////////////////////////////////////////////////////
void Collection::save() {
  if(!load_in_memory) {
    printf("Collection::save - not loaded in memory: %s\n", this->name);
    exit(0);
  }

  int4 hdr[4], adr=16;

  int i;
  for(i=0; i<nClasses; i++) {
  //  if(nColl[i]<0 || nColl[i]>10000) 
  //    printf("nColl error %d %d\n", i, nColl[i]);
    if(nColl[i]) {
//	  printf("Collection::save - writing %d out of %d size=%d\n", 
//		 i, nClasses, nColl[i]);
#ifndef INT2_OFF
      if(type==D_INT2) {
	db->writeFile(coll2[i], nColl[i], adr, 2, fd);
      }
#endif
      if(type==D_INT4) {
	db->writeFile(coll4[i], nColl[i], adr, 4, fd);
      }
    }
  }
  
  adrNColl=adr;
  db->writeFile(nColl, nClasses, adr, 4, fd);
  writeFlag=0;

  adr=0;
  hdr[0]=nClasses;        // 0
  hdr[1]=type;            // 1
  hdr[2]=adrColl;         // 2
  hdr[3]=adrNColl;        // 3

  db->writeFile(hdr, 4, adr, 4, fd);
}
////////////////////////////////////////////////////////////////////
void Collection::close() {
  if(!isOpen()) return;
  if(writeFlag) save();
  db->closeFile(fd);
  if(load_in_memory) {
#ifndef INT2_OFF
    if(type==D_INT2) DeleteArray(&coll2, nClasses);
#endif
    if(type==D_INT4) DeleteArray(&coll4, nClasses);
  }
  else {
#ifndef INT2_OFF
    if(coll2t) delete [] coll2t;
#endif
    if(coll4t) delete [] coll4t;
#ifndef INT2_OFF
    coll2t=NULL;
#endif
    coll4t=NULL;
    index_=-1;
    delete [] nCollt;
  }
  delete [] nColl;
  if(name) {
    delete [] name;
    name=NULL;
  }
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Collection::add(int4 index, int2 item) {
  if(!load_in_memory) {
    printf("Collection::addItem - not loaded in memory: %s\n", this->name);
    exit(0);
  }

  if(type!=D_INT2) {
    printf("Collection::Adding wrong type int2: %s\n", this->name);
    exit(0);
  }
  
  AddToArray(&coll2[index], item, nColl[index]);
  nColl[index]++;
  writeFlag=1;
}
#endif
////////////////////////////////////////////////////////////////////
void Collection::add(int4 index, int4 item) {
  if(!load_in_memory) {
    printf("Collection::AddItem - not loaded in memory: %s\n", this->name);
    exit(0);
  }
  
  if(type!=D_INT4) {
    printf("Collection::Adding wrong type int4: %s\n", this->name);
    exit(0);
  }
  
  AddToArray(&coll4[index], item, nColl[index]);
  nColl[index]++;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
void Collection::add(int4 index, int2* item, int4 nValues) {
  if(!load_in_memory) {
    printf("Collection::addItem - not loaded in memory: %s\n", this->name);
    exit(0);
  }
  
  if(type!=D_INT2) {
    printf("Collection::Adding wrong type int2\n");
    exit(0);
  }
  
  AddToArrayN(&coll2[index], item, nColl[index], nValues);
  nColl[index]+=nValues;
  writeFlag=1;
}
#endif
////////////////////////////////////////////////////////////////////
void Collection::add(int4 index, int4* item, int4 nValues) {
  if(!load_in_memory) {
    printf("Collection::addItem - not loaded in memory: %s\n", this->name);
    exit(0);
  }

  if(type!=D_INT4) {
    printf("Collection::Adding wrong type int4: %s: %s\n", this->name);
    exit(0);
  }
  
  AddToArrayN(&coll4[index], item, nColl[index], nValues);
  nColl[index]+=nValues;
  writeFlag=1;
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int2 Collection::itemValue2(int4 index, int4 cindex) {
  if(load_in_memory) {
    return(*(coll2[index]+cindex));
  }
  else {
    if(index_!=index) {
      if(coll2t) delete [] coll2t;
      db->readFile(coll2, nColl[index], adrColl+nCollt[index]*2, 2, fd);
      index_=index;
    }
    return(coll2t[cindex]);
  }
} 
#endif
////////////////////////////////////////////////////////////////////
int4 Collection::itemValue4(int4 index, int4 cindex) {
  if(load_in_memory) {
    return(*(coll4[index]+cindex));
  }
  else {
    if(index_!=index) {
      if(coll4t) delete [] coll4t;
      db->readFile(coll4, nColl[index], adrColl+nCollt[index]*4, 4, fd);
      index_=index;
    }
    return(coll4t[cindex]);
  }
} 
////////////////////////////////////////////////////////////////////
int4  Collection::itemValue(int4 index, int4 cindex) {
#ifndef INT2_OFF
  if(type==D_INT2) return((int)itemValue2(index, cindex));
#endif
  if(type==D_INT4) return((int)itemValue4(index, cindex));
}
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int2* Collection::item2(int4 index) {
  int n;
  return(item2n(index, n));
}
#endif
////////////////////////////////////////////////////////////////////
#ifndef INT2_OFF
int2* Collection::item2n(int4 index, int4& nValues) {
  nValues=nColl[index];
  if(load_in_memory) {
    return(coll2[index]);
  }
  else {
    if(index_!=index) {
      if(coll2t) {
	delete [] coll2t; coll2t=NULL;
      }
      if(nColl[index]) {
	coll2t=new int2 [nColl[index]];
	db->readFile(coll2t, nColl[index], adrColl+nCollt[index]*2, 2, fd);
	index_=index;
      }
    }
    return(coll2t);
  }
} 
#endif
////////////////////////////////////////////////////////////////////
int4* Collection::item4(int4 index) {
  int4 n;
  return(item4n(index, n));
}
////////////////////////////////////////////////////////////////////
int4* Collection::item4n(int4 index, int4& nValues) {
  nValues=nColl[index];
  if(load_in_memory) {
      
    return(coll4[index]);
  }
  else {
    if(index_!=index) {
      if(coll4t) {
	delete [] coll4t; coll4t=NULL;
      }
      if(nColl[index]) {
	coll4t=new int4 [nColl[index]];
	db->readFile(coll4t, nColl[index], adrColl+nCollt[index]*4, 4, fd);
	index_=index;
      }
    }
    return(coll4t);
  }
}
////////////////////////////////////////////////////////////////////
void Collection::expandCol(Collection *col, char *name, int n) {
	int4 **array, *nObjs, k = 0;
	DataType type = (*col).getType();
	if (type == D_INT4) {
		array = (*col).array4();
	} else array = (int4 **) (*col).array2();

	int nClasses = (*col).getCollectionSize();

	nObjs = new int4 [nClasses];

	int i;

	for (i = 0; i < nClasses; i++) nObjs[i] = (*col).getClassSize(i);

	(*col).create(name, type, nClasses+n);

	for (i = 0; i < nClasses; i++) {
		for (int j = 0; j < nObjs[i]; j++) {
			if (type == D_INT4) {
				(*col).add(i, *(*array+k));
			} 
#ifndef INT2_OFF
			else (*col).add(i, (int2) *(*array+k));
#endif
			k++;
		}
	}

	(*col).save();

	delete [] nObjs;
}
////////////////////////////////////////////////////////////////////
Monomers::Monomers(DB *db_) {
  db=db_;
  db->monomers=this;

  monCode1=NULL;
  monNbond=NULL;
  monBond=NULL;
  monNatom=NULL;
  monAtom=NULL;
  monPrev=NULL;
  monNext=NULL;
  monType=NULL;

  monCode3=new Property("code3.mon", db);
  nMonomers=monCode3->getObjectSize();

  comID=new Property("id.com", db);
}
//////////////////////////////////////////////////////////////////
Monomers::~Monomers() {

  delete monCode3;
  delete monCode1;
  delete monNbond;
  delete monBond;
  delete monNatom;
  delete monAtom;
  delete monPrev;
  delete monNext;
  delete monType;

  delete comID;

  db->monomers=NULL;
}
//////////////////////////////////////////////////////////////////
int Monomers::findCom(char *text) {
  int4 pos;
  return((int)(comID->find(text, (int4)0, pos)));
}
//////////////////////////////////////////////////////////////////
char* Monomers::code3(int iM) {
  return(monCode3->item1(iM));
}
//////////////////////////////////////////////////////////////////
char Monomers::code1(int iM) {
  if(!monCode1) monCode1=new Property("code1.mon", db);
  return((char)(*monCode1->item1(iM)));
}
//////////////////////////////////////////////////////////////////
int Monomers::nBond(int iM) {
  if(!monNbond) monNbond=new Property("n_bond.mon", db);
  return(*monNbond->item2(iM));
}
//////////////////////////////////////////////////////////////////
int2* Monomers::bonds(int iM) {
  if(!monBond) monBond=new Property("bond.mon", db);
  return(monBond->item2(iM));
}
//////////////////////////////////////////////////////////////////
int2 Monomers::bond1(int iM, int iB) {
  if(!monBond) monBond=new Property("bond.mon", db);
  int2 *bonds=monBond->item2(iM);
  return(*(bonds+iB*2));
}
//////////////////////////////////////////////////////////////////
int2 Monomers::bond2(int iM, int iB) {
  if(!monBond) monBond=new Property("bond.mon", db);
  int2 *bonds=monBond->item2(iM);
  return(*(bonds+iB*2+1));
}
//////////////////////////////////////////////////////////////////
int2 Monomers::nAtom(int iM) {
  if(!monNatom) monNatom=new Property("n_atom.mon", db);
  return(*monNatom->item2(iM));
}
//////////////////////////////////////////////////////////////////
char* Monomers::atoms(int iM) {
  if(!monAtom) monAtom=new Property("atom.mon", db);
  return(monAtom->item1(iM));
}
//////////////////////////////////////////////////////////////////
char* Monomers::atom(int iM, int iA) {
  if(!monAtom) monAtom=new Property("atom.mon", db);
  char *atoms=monAtom->item1(iM);
  return(&atoms[iA*4]);
}
//////////////////////////////////////////////////////////////////
int Monomers::type(int iM) {
  if(!monType) monType=new Property("type.mon", db);
  return(*monType->item2(iM));
}
//////////////////////////////////////////////////////////////////
int Monomers::findAtom(int iM, char *text) {
  if(!monAtom) monAtom=new Property("atom.mon", db);
  if(!monNatom) monNatom=new Property("n_atom.mon", db);
  char *atoms=monAtom->item1(iM);
  for(int i=0; i<*monNatom->item2(iM); i++)
    if(!strncmp(atoms+i*4, text, 4)) return(i);
  return(-1);
}
//////////////////////////////////////////////////////////////////
int Monomers::prev(int iM) {
  if(!monPrev) monPrev=new Property("prev.mon", db);
  return(*monPrev->item2(iM));
}
//////////////////////////////////////////////////////////////////
int Monomers::next(int iM) {
  if(!monNext) monNext=new Property("next.mon", db);
  return(*monNext->item2(iM));
}
//////////////////////////////////////////////////////////////////
Entities::~Entities() {
  if(nEnc) {
    DeleteArray(&encName, nEnc);
    delete [] comNum;
    delete [] encNum;
    delete [] enpNum;
    delete [] encNSE;
    
    delete [] encNXYZ;
    DeleteArray(&encSE, nEnc);
    for(int i=0; i<nEnc; i++) 
      if(encSEnum[i]) DeleteArray(&encSEnum[i], encNSE[i]);
    DeleteArray(&encSEnum, nEnc);
    DeleteArray(&encXYZ, nEnc);
    DeleteArray(&encBfac, nEnc);
    DeleteArray(&enc_se_xyz, nEnc);
    DeleteArray(&enc_xyz_se, nEnc);
    
    DeleteArray(&encpSeq, nEnc);
    DeleteArray(&encpKS, nEnc);
    
    DeleteArray(&encpExp, nEnc);
    DeleteArray(&encpPol, nEnc);
    DeleteArray(&encpBfac, nEnc);
    DeleteArray(&encpCa, nEnc);
    DeleteArray(&encpMonType, nEnc);
  }

}
//////////////////////////////////////////////////////////////////
Entities::Entities(DB *db_) {
  db=db_;
  if(db->monomers==NULL) db->monomers=new Monomers(db_);
  monomers=db->monomers;

  nEnc=0;
}
//////////////////////////////////////////////////////////////////
void Entities::addCom(char *text, int addMode) {
  int iCom=monomers->findCom(text);
  if(iCom==-1) {
    printf("Compound %s is not found\n", text);
    exit(0);
  }
  addCom(iCom, addMode);
}
//////////////////////////////////////////////////////////////////
void Entities::addCom(int4 iCom, int addMode) {
  int err, i, j, k, iView; int2 nEncAdd, nEnpAdd, iXYZ, nAtom;
  char *comEnt_; int4 nEncOld;
  int4 iEnc;

  if(!addMode && nEnc) {
    DeleteArray(&encName, nEnc);
    delete [] comNum;
    delete [] encNum;
    delete [] enpNum;
    delete [] encNSE;
    
    delete [] encNXYZ;
    DeleteArray(&encSE, nEnc);
    for(i=0; i<nEnc; i++) 
      if(encSEnum[i]) DeleteArray(&encSEnum[i], encNSE[i]);
    DeleteArray(&encSEnum, nEnc);
    DeleteArray(&encXYZ, nEnc);
    DeleteArray(&encBfac, nEnc);
    DeleteArray(&enc_se_xyz, nEnc);
    DeleteArray(&enc_xyz_se, nEnc);
    
    DeleteArray(&encpSeq, nEnc);
    DeleteArray(&encpKS, nEnc);
    
    DeleteArray(&encpExp, nEnc);
    DeleteArray(&encpPol, nEnc);
    DeleteArray(&encpBfac, nEnc);
    DeleteArray(&encpCa, nEnc);
    DeleteArray(&encpMonType, nEnc);
    
    nEnc=0;
  }
  nEncOld=nEnc;
  
  static Property n_enc_com("n_enc.com", db);

  nEncAdd=*n_enc_com.item2(iCom);

  nEnc+=nEncAdd;
  ExtendArray(&encName, nEncOld, nEnc);
  ExtendArray(&comNum, nEncOld, nEnc);
  ExtendArray(&encNum, nEncOld, nEnc);
  ExtendArray(&enpNum, nEncOld, nEnc);
  ExtendArray(&encNSE, nEncOld, nEnc);

  ExtendArrayValue(&encNXYZ, nEncOld, nEnc, (int4) (-1));
  ExtendArrayValue(&encSE, nEncOld, nEnc, (int2*) NULL);
  ExtendArrayValue(&encSEnum, nEncOld, nEnc, (char**) NULL);

  ExtendArrayValue(&encXYZ, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encBfac, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&enc_se_xyz, nEncOld, nEnc, (int4*) NULL);
  ExtendArrayValue(&enc_xyz_se, nEncOld, nEnc, (int4*) NULL);

  ExtendArrayValue(&encpSeq, nEncOld, nEnc, (char*) NULL);
  ExtendArrayValue(&encpKS, nEncOld, nEnc, (char*) NULL);

  ExtendArrayValue(&encpExp, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpPol, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpBfac, nEncOld, nEnc, (flt4*) NULL);

  ExtendArrayValue(&encpCa, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpMonType, nEncOld, nEnc, (int1*) NULL);

  static Property i_enc_com("i_enc.com", db);
  static Property i_enp_enc("i_enp.enc", db);
  int b_enc_com=*i_enc_com.item4(iCom);

  static Property name_enc("name.enc", 0, db);
  static Property n_se_enc("n_se.enc", 0, db);

  for(i=0; i<nEncAdd; i++) {
    iView=i+nEncOld;
    iEnc=b_enc_com+i;
    encName[iView]=name_enc.item1(iEnc, COPY_ARRAY_ON);
    
    comNum[iView]=iCom;
    encNum[iView]=iEnc;
    enpNum[iView]=*i_enp_enc.item4(iEnc);
    encNSE[iView]=*n_se_enc.item2(iEnc);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::addEnc(char *enc_id, int isAdd) {
  static Property name_enc("name.enc", db);
  int4 iEnc=name_enc.find(enc_id);
  if(iEnc==-1) {
    printf("Entity %s is not found\n", enc_id);
    exit(0);
  }
  addEnc(iEnc, isAdd);
}
//////////////////////////////////////////////////////////////////
void Entities::addEnc(int4 iEnc, int isAdd) {
  int err, i, j, k, iView; int2 nEncAdd, nEnpAdd, iXYZ, nAtom;
  char *comEnt_; int4 nEncOld;

  if(!isAdd && nEnc) {
    DeleteArray(&encName, nEnc);
    delete [] comNum;
    delete [] encNum;
    delete [] enpNum;
    delete [] encNSE;
    
    delete [] encNXYZ;
    DeleteArray(&encSE, nEnc);
    for(i=0; i<nEnc; i++) 
      if(encSEnum[i]) DeleteArray(&encSEnum[i], encNSE[i]);
    DeleteArray(&encSEnum, nEnc);
    DeleteArray(&encXYZ, nEnc);
    DeleteArray(&encBfac, nEnc);
    DeleteArray(&enc_se_xyz, nEnc);
    DeleteArray(&enc_xyz_se, nEnc);
    
    DeleteArray(&encpSeq, nEnc);
    DeleteArray(&encpKS, nEnc);
    
    DeleteArray(&encpExp, nEnc);
    DeleteArray(&encpPol, nEnc);
    DeleteArray(&encpBfac, nEnc);
    DeleteArray(&encpCa, nEnc);
    DeleteArray(&encpMonType, nEnc);
    
    nEnc=0;
  }
  nEncOld=nEnc;
  

  nEnc+=1;
  ExtendArray(&encName, nEncOld, nEnc);
  ExtendArray(&comNum, nEncOld, nEnc);
  ExtendArray(&encNum, nEncOld, nEnc);
  ExtendArray(&enpNum, nEncOld, nEnc);
  ExtendArray(&encNSE, nEncOld, nEnc);

  ExtendArrayValue(&encNXYZ, nEncOld, nEnc, (int4) (-1));
  ExtendArrayValue(&encSE, nEncOld, nEnc, (int2*) NULL);
  ExtendArrayValue(&encSEnum, nEncOld, nEnc, (char**) NULL);

  ExtendArrayValue(&encXYZ, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encBfac, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&enc_se_xyz, nEncOld, nEnc, (int4*) NULL);
  ExtendArrayValue(&enc_xyz_se, nEncOld, nEnc, (int4*) NULL);

  ExtendArrayValue(&encpSeq, nEncOld, nEnc, (char*) NULL);
  ExtendArrayValue(&encpKS, nEncOld, nEnc, (char*) NULL);

  ExtendArrayValue(&encpExp, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpPol, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpBfac, nEncOld, nEnc, (flt4*) NULL);

  ExtendArrayValue(&encpCa, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpMonType, nEncOld, nEnc, (int1*) NULL);

  static Property i_com_enc("i_com.enc", db);
  static Property i_enp_enc("i_enp.enc", db);

  static Property name_enc("name.enc", 0, db);
  static Property n_se_enc("n_se.enc", 0, db);

  int iCom=*i_com_enc.item4(iEnc);
  iView=nEncOld;
  encName[iView]=name_enc.item1(iEnc, COPY_ARRAY_ON);
    
  comNum[iView]=iCom;
  encNum[iView]=iEnc;
  enpNum[iView]=*i_enp_enc.item4(iEnc);
  encNSE[iView]=*n_se_enc.item2(iEnc);
}
//////////////////////////////////////////////////////////////////
void Entities::addEnp(char *enp_id, int isAdd) {
  static Property name_enp("name.enp", db);
  int4 iEnp=name_enp.find(enp_id);
  if(iEnp==-1) {
    printf("Entity %s is not found\n", enp_id);
    exit(0);
  }
  addEnp(iEnp, isAdd);
}
//////////////////////////////////////////////////////////////////
void Entities::addEnp(int4 iEnp, int isAdd) {
  int err, i, j, k, iView; int2 nEncAdd, nEnpAdd, iXYZ, nAtom;
  char *comEnt_; int4 nEncOld;

  static Property i_enc_enp("i_enc.enp", db);
  int iEnc=*i_enc_enp.item4(iEnp);
  
  if(!isAdd && nEnc) {
    DeleteArray(&encName, nEnc);
    delete [] comNum;
    delete [] encNum;
    delete [] enpNum;
    delete [] encNSE;
    
    delete [] encNXYZ;
    DeleteArray(&encSE, nEnc);
    for(i=0; i<nEnc; i++) 
      if(encSEnum[i]) DeleteArray(&encSEnum[i], encNSE[i]);
    DeleteArray(&encSEnum, nEnc);
    DeleteArray(&encXYZ, nEnc);
    DeleteArray(&encBfac, nEnc);
    DeleteArray(&enc_se_xyz, nEnc);
    DeleteArray(&enc_xyz_se, nEnc);
    
    DeleteArray(&encpSeq, nEnc);
    DeleteArray(&encpKS, nEnc);
    
    DeleteArray(&encpExp, nEnc);
    DeleteArray(&encpPol, nEnc);
    DeleteArray(&encpBfac, nEnc);
    DeleteArray(&encpCa, nEnc);
    DeleteArray(&encpMonType, nEnc);
    
    nEnc=0;
  }
  nEncOld=nEnc;
  
  nEnc+=1;
  ExtendArray(&encName, nEncOld, nEnc);
  ExtendArray(&comNum, nEncOld, nEnc);
  ExtendArray(&encNum, nEncOld, nEnc);
  ExtendArray(&enpNum, nEncOld, nEnc);
  ExtendArray(&encNSE, nEncOld, nEnc);

  ExtendArrayValue(&encNXYZ, nEncOld, nEnc, (int4) (-1));
  ExtendArrayValue(&encSE, nEncOld, nEnc, (int2*) NULL);
  ExtendArrayValue(&encSEnum, nEncOld, nEnc, (char**) NULL);

  ExtendArrayValue(&encXYZ, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encBfac, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&enc_se_xyz, nEncOld, nEnc, (int4*) NULL);
  ExtendArrayValue(&enc_xyz_se, nEncOld, nEnc, (int4*) NULL);

  ExtendArrayValue(&encpSeq, nEncOld, nEnc, (char*) NULL);
  ExtendArrayValue(&encpKS, nEncOld, nEnc, (char*) NULL);

  ExtendArrayValue(&encpExp, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpPol, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpBfac, nEncOld, nEnc, (flt4*) NULL);

  ExtendArrayValue(&encpCa, nEncOld, nEnc, (flt4*) NULL);
  ExtendArrayValue(&encpMonType, nEncOld, nEnc, (int1*) NULL);

  static Property i_com_enc("i_com.enc", db);
  static Property name_enc("name.enc", 0, db);
  static Property n_se_enc("n_se.enc", 0, db);

  int iCom=*i_com_enc.item4(iEnc);
  iView=nEncOld;
  encName[iView]=name_enc.item1(iEnc, COPY_ARRAY_ON);
    
  comNum[iView]=iCom;
  encNum[iView]=iEnc;
  enpNum[iView]=iEnp;
  encNSE[iView]=*n_se_enc.item2(iEnc);
}
//////////////////////////////////////////////////////////////////
void Entities::load_encSE(int iE) {
  if(!encSE[iE]) {
    static Property prop("se.enc", 0, db);
    encSE[iE]=prop.item2(encNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encSEnum(int iE) {
  if(!encSEnum[iE]) {
    static Property prop("sen_pdb.enc", 0, db);
    char *tmp=prop.item1(encNum[iE], 1);
    BufferToStringArray(encNSE[iE], &tmp, &encSEnum[iE]);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encNXYZ(int iE) {
  if(encNXYZ[iE]<0) {
    static Property prop("n_xyz.enc", 0, db);
    encNXYZ[iE]=*prop.item4(encNum[iE]);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encXYZ(int iE) {
  load_encNXYZ(iE);
  if(encNXYZ[iE] && !encXYZ[iE]) {
    static Property prop("xyz.enc", 0, db);
    encXYZ[iE]=prop.itemf(encNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encBfac(int iE) {
  load_encNXYZ(iE);
  if(encNXYZ[iE] && !encBfac[iE]) {
    static Property prop("bfac.enc", 0, db);
    encBfac[iE]=prop.itemf(encNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_enc_se_xyz(int iE) {
  if(!enc_se_xyz[iE]) {
    static Property prop("se_xyz.enc", 0, db);
    enc_se_xyz[iE]=prop.item4(encNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_enc_xyz_se(int iE) {
  load_encNXYZ(iE);
  if(encNXYZ[iE] && !enc_xyz_se[iE]) {
    static Property prop("xyz_se.enc", 0, db);
    enc_xyz_se[iE]=prop.item4(encNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpSeq(int iE) {
  if(!encpSeq[iE] && enpNum[iE]>=0) {
    static Property prop("seq.enp", 0, db);
    encpSeq[iE]=prop.item1(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpKS(int iE) {
  if(!encpKS[iE] && enpNum[iE]>=0) {
    static Property prop("k_s.enp", 0, db);
    encpKS[iE]=prop.item1(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpExp(int iE) {
  if(!encpExp[iE] && enpNum[iE]>=0) {
    static Property prop("exp.enp", 0, db);
    encpExp[iE]=prop.itemf(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpPol(int iE) {
  if(!encpPol[iE] && enpNum[iE]>=0) {
    static Property prop("pol.enp", 0, db);
    encpPol[iE]=prop.itemf(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpBfac(int iE) {
  if(!encpBfac[iE] && enpNum[iE]>=0) {
    static Property prop("bfac_flt.enp", 0, db);
    encpBfac[iE]=prop.itemf(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_encpCa(int iE) {
  if(!encpCa[iE] && enpNum[iE]>=0) {
    static Property prop("c_a.enp", 0, db);
    encpCa[iE]=prop.itemf(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::load_monType(int iE) {
  if(!encpMonType[iE] && enpNum[iE]>=0) {
    static Property prop("se_type.enp", 0, db);
    encpMonType[iE]=prop.item1(enpNum[iE], 1);
  }
}
//////////////////////////////////////////////////////////////////
void Entities::reread() {
  for(int i=0; i<nEnc; i++) addEnc(encNum[i], i?1:0);
}
//////////////////////////////////////////////////////////////////
char* IEntity::atoms(int is) {
  if(!entities->encSE[iE]) entities->load_encSE(iE);
  return(entities->monomers->atoms(*(entities->encSE[iE]+is)));
}
//////////////////////////////////////////////////////////////////
char* IEntity::atom(int is, int ia) {
  if(!entities->enc_se_xyz[iE]) entities->load_enc_se_xyz(iE);
  if(!entities->encSE[iE]) entities->load_encSE(iE);

  return(entities->monomers->atom(*(entities->encSE[iE]+is), ia));
}
//////////////////////////////////////////////////////////////////
ISubentity::ISubentity(IEntity ientity) {
  entities=ientity.entities; 
  monomers=entities->monomers; 
  iE=ientity(); 
  iS=0;
  iM=-1;
}
//////////////////////////////////////////////////////////////////
ISubentity ISubentity::neighbor(int shift) {
  int iS_=iS+shift; 
  ISubentity se=*this;
  if(iS_<0 || iS_>=entities->encNSE[iE]) iS_=-1;
  se.iS=iS_;
  iM=-1;
  return(se);
}
//////////////////////////////////////////////////////////////////
int ISubentity::getMon() {
  IMonomer mon(*this);
  return(mon());
}
//////////////////////////////////////////////////////////////////
int ISubentity::code() {
  int c=(int)code1()-(int)'A';
  if(c<0 || c>24) return(0);
  return(aa_decode[c]);
}
//////////////////////////////////////////////////////////////////
IAtom ISubentity::findAtom(char *text) { 
  if(iM==-1) iM=getMon();  
  IAtom atom(*this);
  
  int i=entities->monomers->findAtom(iM, text);

  if(i>=0) {
    if(!entities->enc_se_xyz[iE]) entities->load_enc_se_xyz(iE);
    atom.iA=*(entities->enc_se_xyz[iE]+iS)+i;
    return(atom);
  }
  atom.iA=-1;
  return(atom);
}
//////////////////////////////////////////////////////////////////
int ISubentity::isAtom(char *text) { 
  if(iM==-1) iM=getMon();  
  return(entities->monomers->findAtom(iM, text)>=0?1:0);
}
//////////////////////////////////////////////////////////////////
IAtom::IAtom(IEntity ientity) {
  entities=ientity.entities; 
  iE=ientity(); 
  iS=-1;
  iA=0;
}
//////////////////////////////////////////////////////////////////
IAtom::IAtom(ISubentity isubentity) {
  entities=isubentity.entities; 
  iE=isubentity.myEntityIndex(); 
  iS=isubentity();
  if(!entities->enc_se_xyz[iE]) entities->load_enc_se_xyz(iE);
  iA=*(entities->enc_se_xyz[iE]+iS);
  _name[4]='\0';
}
//////////////////////////////////////////////////////////////////
char* IAtom::name() {
  if(!entities->enc_xyz_se[iE]) entities->load_enc_xyz_se(iE);
  if(!entities->enc_se_xyz[iE]) entities->load_enc_se_xyz(iE);
  if(!entities->encSE[iE]) entities->load_encSE(iE);
  
  int iS_=iS>=0?iS:*(entities->enc_xyz_se[iE]+iA);
  strncpy(_name, entities->monomers->atom(*(entities->encSE[iE]+iS_),
	 iA-*(entities->enc_se_xyz[iE]+iS_)), 4);
  return(_name);
}
//////////////////////////////////////////////////////////////////
IMonomer::IMonomer(ISubentity isubentity) {
  if(!isubentity.entities->encSE[isubentity.iE]) 
    isubentity.entities->load_encSE(isubentity.iE);
  monomers=isubentity.entities->monomers; 
  iM=*(isubentity.entities->encSE[isubentity.myEntityIndex()]+
       isubentity());
  _name[4]='\0';
}
//////////////////////////////////////////////////////////////////
char* IMonomer::atomName(IAtom atom) {
  if(!atom.entities->enc_xyz_se[atom.iE]) 
    atom.entities->load_enc_xyz_se(atom.iE);
  if(!atom.entities->enc_se_xyz[atom.iE]) 
    atom.entities->load_enc_se_xyz(atom.iE);
  int iS_=atom.mySubentityIndex()>=0?atom.mySubentityIndex():
  *(atom.entities->enc_xyz_se[atom.myEntityIndex()]+atom());
  strncpy(_name, monomers->atom
	  (iM, atom()-*(atom.entities->enc_se_xyz[atom.myEntityIndex()]+
			atom.mySubentityIndex())), 4);
  return(_name);
}  
//////////////////////////////////////////////////////////////////
IBond::IBond(IAtom atom){
  if(!atom.entities->enc_xyz_se[atom.iE]) 
    atom.entities->load_enc_xyz_se(atom.iE);
  if(!atom.entities->encSE[atom.iE]) atom.entities->load_encSE(atom.iE);
  monomers=atom.entities->monomers;
  int iS_=atom.mySubentityIndex()>=0?atom.mySubentityIndex():
  *(atom.entities->enc_xyz_se[atom.myEntityIndex()]+atom());
  iM=*(atom.entities->encSE[atom.myEntityIndex()]+atom.mySubentityIndex())-1;
  iB=0;
}
//////////////////////////////////////////////////////////////////
