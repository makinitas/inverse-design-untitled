#ifndef jga_complex       // ************************************************
#define jga_complex 1     // ***  jga/complex.h                           ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Complex number algebra                                            ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Operations:                                                       ***
// ***                                                                    ***
// ***    a+b  a+x  x+a  a-b  a-x  x-a  a*b  a*x  x*a  a/b  a/x  x/a      ***
// ***    a=b  a=x  -a  a+=b  a+=x  a-=b  a-=x  a*=b  a*=x  a/=b  a/=x    ***
// ***                                                                    ***
// ***  Functions:                                                        ***
// ***                                                                    ***
// ***    real(a) = Re(a)     imag(a) = Im(a)                             ***
// ***    mod(a) = |a|        mod2(a) = |a|^2                             ***
// ***    conj(a) = complex conjugate of a                                ***
// ***    exp(a) = e^a        sqrt(a) = a^0.5 (Re>=0)      sqr(a) = a^2   ***
// ***    log(a) = ln(a)      sin(a)  cos(a)  sinh(a)  cosh(a)            ***
// ***    llog(a) = ln(a) for very small arguments                        ***
// ***    complex(a,b) = a + b i                                          ***
// ***    euler(x,y) = x exp(i y)                                         ***
// ***    arg(a) = argument of complex a, in between -pi and pi           ***
// ***    pow(a,x) = a^x                                                  ***
// ***                                                                    ***
// ***  Constants:                                                        ***
// ***                                                                    ***
// ***   i_c = (0,1)                                                      ***
// ***                                                                    ***
// ***  Miscellaneous:                                                    ***
// ***                                                                    ***
// ***    i_l(l) = i^l,  read_complex(file)                               ***
// ***                                                                    ***
// ***  Here a and b are 'complex' numbers and x is a number              ***
// ***  of type 'numero'.                                                 ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** class definition                                                   ***
// **************************************************************************

class complex {
  friend complex operator+(numero   x, const complex &v);
  friend complex operator+(const complex &u, numero   x);
  friend complex operator+(const complex &u, const complex &v);
  friend complex operator-(const complex &v);
  friend complex operator-(numero   x, const complex &v);
  friend complex operator-(const complex &u, numero   x);
  friend complex operator-(const complex &u, const complex &v);
  friend complex operator*(numero   x, const complex &v);
  friend complex operator*(const complex &u, numero   x);
  friend complex operator*(const complex &u, const complex &v);
  friend complex operator/(numero   x, const complex &v);
  friend complex operator/(const complex &u, numero   x);
  friend complex operator/(const complex &u, const complex &v);
  friend numero  real(const complex &v);
  friend numero  imag(const complex &v);
  friend numero  mod(const complex  &v);
  friend numero  abs(const complex  &v);
  friend numero  mod2(const complex &v);
  friend complex conj(const complex &v);
  friend complex exp(const complex  &v);
  friend complex pow(const complex  &v, int m);
  friend complex pow(const complex  &a, numero x);
  friend complex pow(const complex  &a, const complex &b);
  friend complex sqrt(const complex &v);
  friend complex sqr(const complex &v);
  friend complex log(const complex  &v);
  friend complex llog(const complex  &v);
  friend complex sin(const complex  &v);
  friend complex cos(const complex  &v);
  friend complex atan(const complex  &v);
  friend complex sinh(const complex &v);
  friend complex cosh(const complex &v);
  friend complex euler(numero a, numero b);
  friend numero  arg(complex a);
public:
  numero &real(void);
  numero &imag(void);
  complex(void) {}
  complex(numero re);
  complex(int re);
  complex(numero re, numero im);
  complex(const complex &v);
  complex &operator=(numero    x);
  complex &operator=(int       x);
  complex &operator=(const complex  &v);
  complex &operator+=(numero   x);
  complex &operator+=(const complex &v);
  complex &operator-=(numero   x);
  complex &operator-=(const complex &v);
  complex &operator*=(numero   x);
  complex &operator*=(const complex &v);
  complex &operator/=(numero   x);
  complex &operator/=(const complex &v);
private:
  numero re, im;
};


// **************************************************************************
// *** implementation                                                     ***
// **************************************************************************

inline numero &complex::real(void)
       {return re;}
inline numero &complex::imag(void)
       {return im;}
inline complex::complex(numero r)
       {re=r; im=0.0;}
inline complex::complex(int    r)
       {re=r; im=0.0;}
inline complex::complex(numero r, numero i)
       {re=r; im=i;}
inline complex::complex(const complex &v)
       {re=v.re; im=v.im;}
inline complex &complex::operator=(numero    x)
       {re=x; im=0.0; return *this;}
inline complex &complex::operator=(int       x)
       {re=x; im=0.0; return *this;}
inline complex &complex::operator=(const complex  &v)
       {re=v.re; im=v.im; return *this;}
inline complex &complex::operator+=(numero   x)
       {re+=x; return *this;}
inline complex &complex::operator+=(const complex &v)
       {re+=v.re; im+=v.im; return *this;}
inline complex &complex::operator-=(numero   x)
       {re-=x; return *this;}
inline complex &complex::operator-=(const complex &v)
       {re-=v.re; im-=v.im; return *this;}
inline complex &complex::operator*=(numero   x)
       {re*=x; im*=x; return *this;}
inline complex &complex::operator*=(const complex &v)
       {numero r=re;  re=re*v.re-im*v.im; im=r*v.im+im*v.re; return *this;}
inline complex &complex::operator/=(numero   x)
       {re/=x; im/=x; return *this;}
inline complex &complex::operator/=(const complex &v)
       {numero m=v.re*v.re+v.im*v.im, r=re;
        re=(re*v.re+im*v.im)/m;  im=(im*v.re-r*v.im)/m; return *this;}

inline complex operator+(numero   x, const complex &v)
       {return complex(x+v.re,v.im);}
inline complex operator+(const complex &u, numero   x)
       {return complex(u.re+x,u.im);}
inline complex operator+(const complex &u, const complex &v)
       {return complex(u.re+v.re,u.im+v.im);}
inline complex operator-(const complex &v)
       {return complex(-v.re,-v.im);}
inline complex operator-(numero   x, const complex &v)
       {return complex(x-v.re,-v.im);}
inline complex operator-(const complex &u, numero   x)
       {return complex(u.re-x,u.im);}
inline complex operator-(const complex &u, const complex &v)
       {return complex(u.re-v.re,u.im-v.im);}
inline complex operator*(numero   x, const complex &v)
       {return complex(x*v.re,x*v.im);}
inline complex operator*(const complex &u, numero   x)
       {return complex(u.re*x,u.im*x);}
inline complex operator*(const complex &u, const complex &v)
       {return complex(u.re*v.re-u.im*v.im,u.re*v.im+u.im*v.re);}
inline complex operator/(numero   x, const complex &v)
       {numero m=v.re*v.re+v.im*v.im;
        return complex(x*v.re/m,-x*v.im/m);}
inline complex operator/(const complex &u, numero   x)
       {return complex(u.re/x,u.im/x);}
inline complex operator/(const complex &u, const complex &v)
       {numero m=v.re*v.re+v.im*v.im;
        return complex((u.re*v.re+u.im*v.im)/m,(u.im*v.re-u.re*v.im)/m);}
inline numero real(const complex &v)
       {return v.re;}
inline numero imag(const complex &v)
       {return v.im;}
inline numero abs(const complex  &v)
       {return sqrt(v.re*v.re+v.im*v.im);}
inline numero mod(const complex  &v)
       {return sqrt(v.re*v.re+v.im*v.im);}
inline numero mod2(const complex &v)
       {return v.re*v.re+v.im*v.im;}
inline complex conj(const complex &v)
       {return complex(v.re,-v.im);}
inline complex exp(const complex  &v)
       {return exp(v.re)*complex(cos(v.im),sin(v.im));}
inline complex sqr(const complex &v)
       {return complex(v.re*v.re-v.im*v.im, 2*(v.re*v.im));}
inline complex euler(numero a, numero b)
       {return complex(a*cos(b), a*sin(b));}
inline complex sin(const complex &v)
       {return complex(sin(v.re)*cosh(v.im), cos(v.re)*sinh(v.im));}
inline complex cos(const complex &v)
       {return complex(cos(v.re)*cosh(v.im),-sin(v.re)*sinh(v.im));}
inline complex sinh(const complex &v)
       {return (exp(v)-exp(-v))/2;}
inline complex cosh(const complex &v)
       {return (exp(v)+exp(-v))/2;}

complex pow(const complex &v, int m)
{
  if(m<0)   return pow(1/v,-m);
  if(m==0)  return complex(1,0);
  if(m==1)  return v;
  if(m==2)  return v*v;
  if(m%2)   return v*sqr(pow(v,m/2));
            return   sqr(pow(v,m/2));
}

complex pow(const complex  &a, numero x)
{
  if(x<0)   return pow(1/a,-x);
  if(x==0)  return complex(1,0);
            return euler(pow(mod(a),x),x*arg(a));
}

complex pow(const complex  &a, const complex &b)
{
  numero r=mod(a),fi=arg(a);
  return euler(pow(r,b.re)*exp(-fi*b.im),fi*b.re+b.im*log(r));
}

complex sqrt(const complex &v)  // the returned real part is >=0
{
  if(!v.im)  {if(v.re<0)  return complex(0,sqrt(-v.re));
              else        return complex(sqrt(v.re), 0);}

  numero a;

  if(v.re>0) {
    a=sqrt((v.re+sqrt(v.re*v.re+v.im*v.im))/2);
    return complex(a,v.im/(2*a));
  } else {
    a=sqrt((-v.re+sqrt(v.re*v.re+v.im*v.im))/2);
    if(v.im>=0)  return complex(v.im/(2*a),a);
                 return complex(-v.im/(2*a),-a);
} }

numero arg(complex a)
{
  if(!a.re && !a.im)  return 0;

  return atan2(a.im,a.re);
}

complex log(const complex  &v)
{
  complex b=0;
  if(!v.re && !v.im)  on_error("log","log(0)\n");
  numero va=mod(v);
  if(va<1e-20)  return llog(v);

  b.re=log(va);  b.im=atan2(v.im,v.re);
  while(b.im>pi)  b.im-=2*pi;
  while(b.im<-pi) b.im+=2*pi;

  return b;
}

complex llog(const complex  &v)   // for very small or very large arguments
{
  numero magnitude=fabs(v.re)+fabs(v.im);
  complex b=0, w=v/magnitude;
  if(!v.re && !v.im)  on_error("llog","log(0)\n");

  b.re=log(mod(w))+log(magnitude);  b.im=atan2(w.im,w.re);
  while(b.im>pi)  b.im-=2*pi;
  while(b.im<-pi) b.im+=2*pi;

  return b;
}

complex atan(const complex &v)
{
  complex i_c(0,1);
  complex sq=sqrt(1-v*v);
  complex x=v+i_c*sq;
  numero  r=mod(x);
  numero  ar=arg(x);
  complex val=complex(ar, -log(r));

  return  val;
}

// --------------------------------------------------------------------------

complex i_c(0,1);

complex i_l(int l)
{
  l=l%4;  if(l<0) l+=4;

  if(l==0)  return complex( 1, 0);
  if(l==1)  return complex( 0, 1);
  if(l==2)  return complex(-1, 0);
  if(l==3)  return complex( 0,-1);

  return 0;
}

complex read_complex(FILE *f)
{
  numero x=read_numero(f);
  return complex(x,read_numero(f));
}

// **************************************************************************
// *** complex arrays manipulations                                       ***
// **************************************************************************

complex *copy(complex *x, int n)
{
  complex *val;  val=new complex [n];
  int i;

  for(i=0; i<n; i++)  val[i]=x[i];

  return val;
}

#endif  // ******************************************************************
