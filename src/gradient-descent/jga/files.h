#ifndef jga_files         // ************************************************
#define jga_files 1       // ***  jga/files.h                 17-ix-1999  ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  File handing:                                                     ***
// ***                                                                    ***
// ***     - open and close files with possible error messages            ***
// ***     - counting number of rows, columns, chars, etc.                ***
// ***                                                                    ***
// ***  ----------------------------------------- ROUTINES -------------  ***
// ***                                                                    ***
// ***  --- open files and display error message if it cannot be opened   ***
// ***                                                                    ***
// ***         FILE *open_file(fout,name,mode);                           ***
// ***         FILE *open_file(name,mode);                                ***
// ***                                                                    ***
// ***      the file name is 'name', the opening mode is 'mode' (e.g.     ***
// ***      "w"), and the message is written to file fout or stdout;      ***
// ***      no action if name=="none", in which case it returns NULL;     ***
// ***      use                                                           ***
// ***                                                                    ***
// ***         close_file(ff);                                            ***
// ***                                                                    ***
// ***      to close files keeping in mind that ff can be NULL            ***
// ***                                                                    ***
// ***  --- skip(fin,ls);                                                 ***
// ***                                                                    ***
// ***      skip ls lines in file fin                                     ***
// ***                                                                    ***
// ***  --- number_of_columns(ls, *name);                                 ***
// ***  --- number_of_words(ls, *name);                                   ***
// ***  --- number_of_rows(ls, *name);                                    ***
// ***  --- number_of_chars(ls, *name);                                   ***
// ***                                                                    ***
// ***      return number of columns, words, rows, and chars in file of   ***
// ***      name "name", skipping ls lines at the begining                ***
// ***                                                                    ***
// ***  --- number_of_columns(*name);                                     ***
// ***  --- number_of_words(*name);                                       ***
// ***  --- number_of_rows(*name);                                        ***
// ***  --- number_of_chars(*name);                                       ***
// ***                                                                    ***
// ***      ident., without skipping lines                                ***
// ***                                                                    ***
// ***  --- limit_values(ls, *name, &nc, *vmin, *vmax, *vave);            ***
// ***                                                                    ***
// ***      on output, vmin,vmax,vave are arrays of nc positions (number  ***
// ***      of columns in file "name") containing the minimum, maximum    ***
// ***      and average value of each column; ls lines are skipped        ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** open and close files                                               ***
// **************************************************************************

FILE *open_file(FILE *fout, char *name, const char *mode)
{
  FILE *ff;  ff=NULL;

  if(!strcmpC(name,(char *)"inline"))   return stdin;   // stdin if name=="inline"
  if(!strcmpC(name,(char *)"terminal")) return stdout;  // stdout if name=="terminal"
  if(strcmpC(name,(char *)"none"))                      // no action if name=="none"
  if((ff=fopen(name,mode))==NULL)
    on_error(fout, "open_file", "cannot open file", name);

  return ff;
}

// --------------------------------------------------------------------------

FILE *open_file(char *name, const char *mode)
  {return open_file(foutput,name,mode);}

// --------------------------------------------------------------------------

void close_file(FILE *ff)
  {if(ff!=NULL && ff!=stdout && ff!=stdin) fclose(ff);}


// **************************************************************************
// *** skip lines at the begining of file                                 ***
// **************************************************************************

void skip(FILE *fin, int ls)
{
  int i;  char c;
  for(i=0; i<ls; i++)  do fscanf(fin,"%c",&c); while(c!='\n' && !feof(fin));
}


// **************************************************************************
// *** count columns, rows, words, and chars                              ***
// **************************************************************************

int number_of_columns(int ls, char *name)
{
  FILE *fin;
  int  nc=0;
  char c;

  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin)) {
               do fscanf(fin,"%c",&c); while(c==' '||c=='\t'||c=='\n');
    do {nc++;  do fscanf(fin,"%c",&c); while(c!=' '&&c!='\t'&&c!='\n');
      if(c!='\n') do fscanf(fin,"%c",&c); while(c==' '||c=='\t');
    } while(c!='\n');
  }
  close_file(fin);

  return nc;
}

// --------------------------------------------------------------------------

int number_of_words(int ls, char *name)
{
  FILE *fin;
  char palabra[80];
  int  nw=-1;       // -1 to avoid counting the end of file as an extra word

  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin))  do {fscanf(fin,"%s",palabra);  nw++;} while(!feof(fin));
  close_file(fin);

  return nw;
}

// --------------------------------------------------------------------------

int number_of_rows(int ls, char *name)
{
  int nc=number_of_columns(ls,name);
  int nw=number_of_words(ls,name);

  if(nw%nc!=0)  on_warning("number_of_rows", "unstructured file");

  return  nw/nc;
}

// --------------------------------------------------------------------------

int number_of_chars(int ls, char *name)
{
  FILE *fin;
  int  nc=0;
  char c;

  fin=open_file(name,"r");
  skip(fin,ls);
  if(!feof(fin))  do {fscanf(fin,"%c",&c);  nc++;} while(!feof(fin));
  close_file(fin);

  return nc;
}

// --------------------------------------------------------------------------

int number_of_columns(char *name) {return number_of_columns(0,name);}
int number_of_words  (char *name) {return number_of_words  (0,name);}
int number_of_rows   (char *name) {return number_of_rows   (0,name);}
int number_of_chars  (char *name) {return number_of_chars  (0,name);}


// **************************************************************************
// *** skip lines at the begining of file                                 ***
// **************************************************************************

void limit_values(int ls, char *name, int &nc,
                  numero *vmin, numero *vmax, numero *vave)
{
  FILE *fin;
  int  nr=number_of_rows(ls,name);  nc=number_of_columns(ls,name);
  int  i,j;
  float  v;

  for(j=0; j<nc; j++)  {vmin[j]=infinity; vmax[j]=-infinity; vave[j]=0;}

  fin=open_file(name,"r");
  skip(fin,ls);
  for(i=0; i<nr; i++)
  for(j=0; j<nc; j++) {
     fscanf(fin,"%f",&v);  vave[j]+=v;
     if(v>vmax[j])  vmax[j]=v;
     if(v<vmin[j])  vmin[j]=v; 
  }
  close_file(fin);

  for(j=0; j<nc; j++)  vave[j]=vave[j]/nr;
}

#endif  // ******************************************************************
