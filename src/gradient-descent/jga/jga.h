#ifndef jga_jga           // ************************************************
#define jga_jga 1         // ***  jga/jga.h                               ***
                          // ***                                          ***
#include <stdio.h>        // ***  <-- standard libraries of common use    ***
#include <stdlib.h>       // ***                                          ***
#include <string.h>       // ***                                          ***
#include <time.h>         // ***                                          ***
#include <math.h>         // ***                                          ***
                          // ***                                          ***
#define numero double     // ***  <-- definition of 'numero'              ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***  Routines of general use in the jga library                        ***
// ***                                                                    ***
// ***  ----------------------------------------- CONSTANTS ------------  ***
//                                                                        ***
        numero                                        //                  ***
                                                      //                  ***
           pi      = 3.141592653589793238462643,      //                  ***
                                                      //                  ***
           au_eV   = 27.2113834,         // 1 a.u. of enery in eV         ***
           a0_au   = 0.529177208,        // Bohr radius a0 in Angstroms   ***
           nm      = 10/a0_au,           // 1 nm in a.u.                  ***
           c_au    = 137.03599971;       // 1/alfa=speed of light (a.u.)  ***
                                         //                               ***
           #define infinite      1.0e20  // effective infinity            ***
           #define infinity      1.0e20  // effective infinity            ***
           #define infinite_int  10000   // effective infinite integer    ***
           #define infinitesimal 1.0e-40 // effective epsilon             ***
//                                                                        ***
// ***  ----------------------------------------- VARIABLES ------------  ***
// ***                                                                    ***
// ***  --- fact[n]      =  n!           (init_fact()  required)          ***
// ***  --- fact2[n]     = (n!)^(1/2)    (init_fact()  required)          ***
// ***  --- ffact[n]     =  n!!          (init_ffact() required)          ***
// ***  --- ffact2[n]    =  Gamma(n/2)   (init_ffact() required)          ***
// ***                                                                    ***
// ***      notice that init_ffact() also calls init_fact()               ***
// ***                                                                    ***
// ***  ----------------------------------------- ROUTINES -------------  ***
// ***                                                                    ***
// ***  --- n!, etc.                                                      ***
// ***                                                                    ***
// ***         numero factorial(int n);      n!                           ***
// ***         numero combination(n,k);      n!/[k!(n-k)!]                ***
// ***                                                                    ***
// ***  --- elementary function                                           ***
// ***                                                                    ***
// ***         sqr(x);                       x^2 for int's and numeros    ***
// ***         ABS(x);                       |x| for int's and numeros    ***
// ***         sign_1l(l)                    (-1)^l                       ***
// ***         SIGN(x)                     x<0 -> -1, x>0 -> 1, x=0 -> 0  ***
// ***                                                                    ***
// ***  --- *copy(numero *x, int n) returns a pointer to a numero array   ***
// ***      with n elements assigned to the first n values of x[i]        ***
// ***                                                                    ***
// ***  --- *copy(int *x, int n), like *copy(numero *x, int n) for ints.  ***
// ***                                                                    ***
// ***  --- cartesian_to_spherical(x,y,z,r,th,fi);                        ***
// ***                                                                    ***
// ***      on output, (r,th,fi) are the spherical coord. of (x,y,z)      ***
// ***                                                                    ***
// ***  --- cartesian_to_spherical_mu(x,y,z,r,cth,fi);                    ***
// ***                                                                    ***
// ***      like cartesian_to_spherical(...), but cth=cos(th)             ***
// ***                                                                    ***
// ***  --- normalize_angles(th,fi);                                      ***
// ***                                                                    ***
// ***      on output, 0<=th<=pi and -pi<fi<=pi                           ***
// ***                                                                    ***
// ***  --- strcmpC(str1,str2)                                            ***
// ***                                                                    ***
// ***      like strcmp(...) with additional features                     ***
// ***                                                                    ***
// ***  --- input strings, int's, and numeros disregarding comments       ***
// ***      with the same format as comments in C++; also, 'infinity'     ***
// ***      and 'zero' are understood as a large number and 0, respect.   ***
// ***      (see routine also for special role of 'all')                  ***
// ***                                                                    ***
// ***         read_name(f,name);            read 'name' from file f      ***
// ***         int read_int(f);              read int from file f         ***
// ***         numero read_numero(f);        read numero from file f      ***
// ***         int read_int(name);           read int from string name    ***
// ***         numero read_numero(name);     read numero from str. name   ***
// ***         read_name(name);              read 'name' from stdin       ***
// ***         int read_int();               read int from stdin          ***
// ***         numero read_numero();         read numero from stdin       ***
// ***                                                                    ***
// ***  --- char *ftostr(x,d);    // ftostr(x,d,1) forces d figures       ***
// ***                                                                    ***
// ***      real x -> string with d figures after the decimal point       ***
// ***                                                                    ***
// ***  --- error and warning messages; on_error abort the program after  ***
// ***      displaying an error message, and on_warning continues the     ***
// ***      program execution; their format is                            ***
// ***                                                                    ***
// ***         on_error  (where, comment);                                ***
// ***         on_warning(where, comment);                                ***
// ***         on_error  (where, comment1, comment2);                     ***
// ***         on_warning(where, comment1, comment2);                     ***
// ***         on_error  (where, comment1, x, comment2);                  ***
// ***         on_warning(where, comment1, x, comment2);                  ***
// ***         on_error  (fout, where, comment);                          ***
// ***         on_warning(fout, where, comment);                          ***
// ***         on_error  (fout, where, comment1, comment2);               ***
// ***         on_warning(fout, where, comment1, comment2);               ***
// ***         on_error  (fout, where, comment1, x, comment2);            ***
// ***         on_warning(fout, where, comment1, x, comment2);            ***
// ***                                                                    ***
// ***      and the error/warning/message displayed is                    ***
// ***                                                                    ***
// ***         *** error/warning in where: comment                        ***
// ***         *** error/warning in where: comment1 comment2              ***
// ***         *** error/warning in where: comment1 x comment2            ***
// ***                                                                    ***
// ***      either in stdout or in file fout, respectively                ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** miscellaneous                                                      ***
// **************************************************************************

int  verbose=0;              // general purpose flag for logs

int  relativistic_flag=1;
void relativistic_on (void)  {c_au = 137.03599976;  relativistic_flag=1;}
void relativistic_off(void)  {c_au = 137035.99976;  relativistic_flag=0;}

typedef enum {electrones, fotones, sonido} particle_type_def;
particle_type_def  particle_type=fotones;

FILE *foutput=stdout;        // file(s) normally used as the standard output

// **************************************************************************
// *** factorial and combinations                                         ***
// **************************************************************************

#define fact_max    170               // maxium n stored in fact[n]=n!
#define fact2_max   300               // maxium n stored in fact2[n]
numero  *fact,                        // fact[n] = n!
        *fact2,                       // fact2[n]=(n!)^(1/2)
        *ffact,                       // ffact[n]= n!!
        *ffact2;                      // ffact2[n]= Gamma(n/2)
int     fact_initialized=0,           // flag to show if fact has been init.
        ffact_initialized=0;          // flag to show if ffact has been init.

// --------------------------------------------------------------------------

numero factorial(int n)               // n!
{
  numero p=1;
  if(n<0) {printf("fact: negative number"); return 1;}
  while(n) p*=n--;
  return p;
}

// --------------------------------------------------------------------------

int init_fact(void)                   // init fact (n!) and fact2 (n!^0.5)
{
  if(fact_initialized)  return 0;

  int n;

  fact   = new numero [fact_max+1];       fact_initialized=1;
  fact2  = new numero [fact2_max+1];

  fact[0]=fact2[0]=1;
  for(n=1; n<=fact_max; n++) {
    fact[n]=n*fact[n-1];  fact2[n]=sqrt(n)*fact2[n-1];
  }
  for(n=fact_max+1; n<=fact2_max; n++)  fact2[n]=sqrt(n)*fact2[n-1];

  return 0;
}

// --------------------------------------------------------------------------

int init_ffact(void)                  // init ffact (n!!)
{                                     // and ffact2 (gamma(n/2))
  if(ffact_initialized)  return 0;

  int n;    init_fact();

  ffact   = new numero [fact2_max+1];     ffact_initialized=1;
  ffact2  = new numero [fact2_max+1];

  ffact[0]=ffact2[0]=ffact[1]=1;  ffact2[1]=sqrt(pi);
  for(n=2; n<=fact2_max; n++)  {
     ffact[n]=n*ffact[n-2];
     if(n%2==0)  ffact2[n]=fact[n/2-1];
     else        ffact2[n]=(n-2)*(ffact2[n-2]/2);
  }

  return 0;
}

// --------------------------------------------------------------------------

numero combination(int n, int k)      // n!/[k!(n-k)!]
{
  if(fact_initialized)  return (fact[n]     /(fact[k]     *fact[n-k]     ));
  else                  return (factorial(n)/(factorial(k)*factorial(n-k)));
}


// **************************************************************************
// *** elementary functions for int's and numeros                         ***
// **************************************************************************

numero sqr(numero a) {return a*a;}                        // a^2
int    sqr(int    i) {return i*i;}                        // i^2

// --------------------------------------------------------------------------

numero ABS(numero a) {return (a<0)?-a:a;}                 // |a|
int    ABS(int    i) {return (i<0)?-i:i;}                 // |i|

// --------------------------------------------------------------------------

numero SIGN(numero a) {return (a<0)?-1:((a>0)?1:0);}      // sign(a)
int    SIGN(int    i) {return (i<0)?-1:((i>0)?1:0);}      // sign(i)

// --------------------------------------------------------------------------

int sign_1l(int l)  {if(l%2)  return -1;  return 1;}      // (-1)^l


// **************************************************************************
// *** numerical arrays manipulations                                     ***
// **************************************************************************

numero *copy(numero *x, int n)
{
  numero *val;  val=new numero [n];
  int i;

  for(i=0; i<n; i++)  val[i]=x[i];

  return val;
}

// --------------------------------------------------------------------------

int *copy(int *x, int n)
{
  int *val;  val=new int [n];
  int i;

  for(i=0; i<n; i++)  val[i]=x[i];

  return val;
}


// **************************************************************************
// *** transform cartesion coordinates into spherical                     ***
// **************************************************************************

void cartesian_to_spherical_mu(numero x, numero y, numero z,
                               numero &d, numero &cth, numero &fi)
{
  d=sqrt(x*x+y*y+z*z);                        // On output, (d,th,fi) are
  if(d>0)  cth=z/d;  else  cth=1;             // spherical coord. of (x,y,z)
  if(!x&&!y) fi=0; else fi=atan2(y,x);        // and cth=cos(th)
}

// --------------------------------------------------------------------------

void cartesian_to_spherical(numero x, numero y, numero z,
                            numero &d, numero &th, numero &fi)
{
  cartesian_to_spherical_mu(x,y,z, d,th,fi);  // th is returned here
  if(th>1)  th=0;  else  if(th<-1)  th=pi;    // rather than cos(th)
  th=acos(th);
}


// **************************************************************************
// *** normalize polar angles to be   0<=th<=pi  and -pi<fi<=pi           ***
// **************************************************************************

void normalize_angles(numero &th, numero &fi)        // sets  0<=th<=pi
{                                                    // and  -pi<fi<=pi
  while(th>pi) th-=2*pi;   while(th<-pi) th+=2*pi;
  if(th<0) {th=-th; fi=pi+fi;}
  while(fi>pi) fi-=2*pi;   while(fi<=-pi) fi+=2*pi;
}


// **************************************************************************
// *** comparare two strings with independence of the case                ***
// **************************************************************************

int strcmpC(char *str1, char *str2)
{
  return strcmp(str1,str2);           // implement!!!  (jga)
}


// **************************************************************************
// *** read input string/int/float and disregard comments                 ***
// **************************************************************************

char read_command[100];                          // internal variable

// --------------------------------------------------------------------------

void read_name(FILE *f, char *command)
{
  fscanf(f,"%s",command);  char c;  int i;       // comments start with /*
                                                 // and end with '*/';
  if(strlen(command)>1) {                        // also, '//' like in C++
    if(command[0]=='/' && command[1]=='*') {
      c=1;
      do {
        fscanf(f,"%s",command);  i=strlen(command);
        if(i>1)  if(command[i-2]=='*')  if(command[i-1]=='/')  c=0;
      } while(c);
      read_name(f, command);
    }  else {
      if(command[0]=='/' && command[1]=='/') {
        do fscanf(f,"%c",&c); while(c!='\n');
        read_name(f,command);
} } } }

// --------------------------------------------------------------------------

numero read_numero(char *command)
{
  if(!strcmpC(command,(char *)"infinite"))   return infinity;
  if(!strcmpC(command,(char *)"infinity"))   return infinity;
  if(!strcmpC(command,(char *)"zero"))       return 0;
                                             return atof(command);
}

// --------------------------------------------------------------------------

int read_int(char *command)
{
  if(!strcmpC(command,(char *)"infinite"))   return 3*infinite_int;
  if(!strcmpC(command,(char *)"infinity"))   return 3*infinite_int;
  if(!strcmpC(command,(char *)"all"))        return 2*infinite_int;
  if(!strcmpC(command,(char *)"zero"))       return 0;
  if(!strcmpC(command,(char *)"off"))        return 0;
  if(!strcmpC(command,(char *)"on"))         return 1;
  if(!strcmpC(command,(char *)"false"))      return 0;
  if(!strcmpC(command,(char *)"true"))       return 1;
                                             return atoi(command);
}

// --------------------------------------------------------------------------

numero read_numero(FILE *f)
{
  read_name(f,read_command);
  return read_numero(read_command);
}

// --------------------------------------------------------------------------

int read_int(FILE *f)
{
  read_name(f,read_command);
  return read_int(read_command);
}

// --------------------------------------------------------------------------

numero read_numero(void)  {return read_numero(stdin);}
int    read_int(void)     {return read_int(stdin);}


// **************************************************************************
// *** transform real number x into string with d figures after the point ***
// **************************************************************************

char ftostr_data[50];

char *ftostr(numero x, int d, int opt)
{
  if(d>0)  if(x-floor(x)==0 && opt==0)  return ftostr(x,0,0);

  char *val,*v;  val=ftostr_data;  v=val;
  numero dec=1e20, eps=1e-12;
  int i=0, j;

  if(x<0) {val[i]='-';  i++;  x=ABS(x);}

  while(floor(x/dec+eps)==0 && dec>1)  dec=dec/10;

  while(dec>=1) {
    j=((int) floor(x/dec+eps));
    x=x-j*dec;  val[i]='0'+j;  i++;  dec=dec/10;
  }

  if(ABS(x)>eps || opt) {
    val[i]='.';  i++;
    while(d) {
      j=int(floor(x/dec+eps));   x=x-j*dec;
      if(d==1 && j<9) if(x/dec>=0.5) j++;
      val[i]='0'+j;  i++;  dec=dec/10;
      d--;
  } }

  val[i]='\0';

  return v;
}

char *ftostr(numero x, int d) {return ftostr(x,d,0);}


// **************************************************************************
// *** error and warning messages                                         ***
// **************************************************************************

void on_warning(FILE *fout, const char *where, const char *com1, numero x,
                const char *com2)
  {fprintf(fout,"*** warning in '%s': %s %g %s\n",where,com1,x,com2);}

void on_error(FILE *fout, const char *where, const char *com1, numero x,
                const char *com2)
  {fprintf(fout,"*** error in '%s': %s %g %s\n",where,com1,x,com2); exit(1);}

void on_warning(FILE *fout, const char *where, const char *com1,
                const char *com2)
  {fprintf(fout,"*** warning in '%s': %s %s\n",where,com1,com2);}

void on_error(FILE *fout, const char *where, const char *com1,
                const char *com2)
  {fprintf(fout,"*** error in '%s': %s %s\n",where,com1,com2); exit(1);}

void on_warning(FILE *fout, const char *where, const char *com)
  {fprintf(fout,"*** warning in '%s': %s\n",where,com);}

void on_error(FILE *fout, const char *where, const char *com)
  {fprintf(fout,"*** error in '%s': %s\n",where,com); exit(1);}

void on_warning(const char *where, const char *com1, numero x,
                const char *com2)
  {on_warning(foutput, where,com1,x,com2);}

void on_error(const char *where, const char *com1, numero x, const char *com2)
  {on_error(foutput, where,com1,x,com2);}

void on_warning(const char *where, const char *com1, const char *com2)
  {on_warning(foutput, where,com1,com2);}

void on_error(const char *where, const char *com1, const char *com2)
  {on_error(foutput, where,com1,com2);}

void on_warning(const char *where, const char *com)
  {on_warning(foutput, where,com);}

void on_error(const char *where, const char *com)
  {on_error(foutput, where,com);}

void on_syntax_error(const char *where, const char *com)
  {fprintf(foutput,"*** syntax error in '%s': %s\n",where,com);  exit(1);}

void on_syntax_error(const char *where, const char *c1, const char *c2)
  {fprintf(foutput,"*** syntax error in '%s': %s%s\n",where,c1,c2);  exit(1);}

#endif  // ******************************************************************
