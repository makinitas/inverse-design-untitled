#ifndef jga_bessel        // ************************************************
#define jga_bessel 1      // ***  jga/bessel.h                            ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Bessel functions (include jga/complex for complex functions)      ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  FUNCTION         VALUE     PARAMETERS                INCLUDE      ***
// ***                                                                    ***
// ***  besselJ(n,x)     J_n(x),   int  n,  real x                        ***
// ***  besselJp(n,x)    J'_n(x),  int  n,  real x                        ***
// ***                                                                    ***
// ***  besselI(n,x)     I_n(x),   int  n,  real/complex x                ***
// ***  besselK(n,x)     K_n(x),   int  n,  real/complex x                ***
// ***                                                                    ***
// ***  besselIp(n,x)    I'_n(x),  int  n,  real/complex x                ***
// ***  besselKp(n,x)    K'_n(x),  int  n,  real/complex x                ***
// ***                                                                    ***
// ***  besselj(l,x)     j_l(x),   int  l,  real x                        ***
// ***  bessely(l,x)     y_l(x),   int  l,  real x                        ***
// ***  besselh(l,x)     h+_l(x),  int  l,  real x                        ***
// ***                                                                    ***
// ***  besselj(l,x)     j_l(x),   int  l,  complex x        jga/qfint    ***
// ***  bessely(l,x)     y_l(x),   int  l,  complex x                     ***
// ***  besselh(l,x)     h+_l(x),  int  l,  complex x        jga/qfint    ***
// ***  besselilh(l,x)   i^l h+_l(x), int l, complex x       jga          ***
// ***  besselilhp(l,x)  i^l h'+_l(x), int l, complex x      jga          ***
// ***                                                                    ***
// ***  besseljp(l,x)    j'_l(x),  int  l,  complex x        jga/qfint    ***
// ***  besselyp(l,x)    y'_l(x),  int  l,  complex x                     ***
// ***  besselhp(l,x)    h'+_l(x), int  l,  complex x        jga/qfint    ***
// ***                                                                    ***
// ***  besselh1(l,x)    h1_l(x),  int  l,  complex x        jga/qfint    ***
// ***  besselh2(l,x)    h2_l(x),  int  l,  complex x        jga/qfint    ***
// ***  besselh1p(l,x)   h1'_l(x), int  l,  complex x        jga/qfint    ***
// ***  besselh2p(l,x)   h2'_l(x), int  l,  complex x        jga/qfint    ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  ROUTINE                PARAMETERS and VALUES ON OUTPUT            ***
// ***                                                                    ***
// ***  besselJY(n,x,J,Y,Jp,Yp)   input:  real n, real x                  ***
// ***                            output: J=J_n(x),  Jp=J'_n(x)           ***
// ***                                    Y=Y_n(x),  Yp=Y'_n(x)           ***
// ***                                                                    ***
// ***  besselIK(n,x,I,K,Ip,Kp)   input:  real n, real x                  ***
// ***                            output: I=I_n(x),  Ip=I'_n(x)           ***
// ***                                    K=K_n(x),  Kp=K'_n(x)           ***
// ***                                                                    ***
// ***  besseljy(l,x,j,y)         input:  int l,  real x                  ***
// ***                            output: j=j_l(x),  y=y_l(x)             ***
// ***                                                                    ***
// ***  besseljyp(l,x,jp,yp)      input:  int l,  real x                  ***
// ***                            output: jp=j'_l(x),  yp=y'_l(x)         ***
// ***                                                                    ***
// ***  besseljh(l,x,j,jp,h,hp)   input:  int l,  real x                  ***
// ***                            output: j=j_l(x), jp=j'_l(x),           ***
// ***                                    h=h_l^(+)(x), hp=h'_l^(+)(x)    ***
// ***                                                                    ***
// **************************************************************************


// **************************************************************************
// *** J_n(x), Y_n(x), J'_n(x), Y'_n(x),  real n,  real x                 ***
// **************************************************************************

numero chebev(numero a, numero b, numero c[], int m, numero x)
{
  numero d=0,dd=0,sv,y,y2;  int j;

  if((x-a)*(x-b)>0)  on_warning("chebev", "x not in range");
  y=(2*x-a-b)/(b-a);  y2=2*y;
  for(j=m; j>=1; j--) {sv=d;  d=y2*d-dd+c[j];  dd=sv;}

  return y*d-dd+c[0]/2;
}

// --------------------------------------------------------------------------

void beschb (numero x, numero &g1, numero &g2, numero &gamp, numero &gamm)
{
  numero xx,c1[7],c2[8];

  c1[0] = -1.142022680371172;  c1[1] = 0.006516511267076;
  c1[2] = 0.000308709017308;   c1[3] = -3.470626964e-6;
  c1[4] = 6.943764e-9;         c1[5] = 3.678e-11;
  c1[6] = -1.36e-13;
  c2[0] = 1.843740587300906;   c2[1] = -0.076852840844786;
  c2[2] = 0.001271927136655;   c2[3] = -4.971736704e-6;
  c2[4] = -3.3126120e-8;       c2[5] = 2.42310e-10;
  c2[6] = -1.7e-13;            c2[7] = -1.0e-15;

  xx=8*x*x-1;
  g1=chebev(-1,1,c1,6,xx);     g2=chebev(-1,1,c2,7,xx);
  gamp=g2-x*g1;                gamm=g2+x*g1;
}

// --------------------------------------------------------------------------

int besselJY(numero xnu, numero x, numero &rj,  numero &ry,
             numero &rjp, numero &ryp)
{
  int i,isign=1,l,nl, maxit=100000, salir;
  numero xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
	 a,b,br,bi,c,cr,ci,d=0,del,del1,den,di,dlr,dli,
       	 dr,e,f,fct,fct2,fct3,ff,gam,gam1,gam2,gammi,gampl,h,
      	 p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,
       	 rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2;

  if(x<=0)  on_error("besselJY", "wrong argument x =", x, "");
  if(xnu<0)  on_error("besselJY", "wrong order n =", xnu, "");
  if(x<xmin)    nl=int(xnu+0.5);
  else {nl=int(xnu-x+1.5); nl=0>nl?0:nl;}
  xmu=xnu-nl;  xmu2=xmu*xmu;  xi=1/x;  xi2=2*xi;  w=xi2/pi;  h=xnu*xi;

  if(h<fpmin) h=fpmin;     c=h;  b=xi2*xnu;

  for(i=1, salir=1; i<=maxit && salir; i++) {
    b+=xi2;  d=b-d;             if(fabs(d)<fpmin) d=fpmin;
    c=b-1/c;                    if(fabs(c)<fpmin) c=fpmin;
    d=1/d;  del=c*d;  h=del*h;  if(d<0) isign=-isign;
    if(fabs(del-1)<eps) salir=0; else salir=1;
  }
  if(salir) on_warning("besselJY", "x too large", x, "");

  rjl1=rjl=isign*fpmin;   rjp1=rjpl=h*rjl;  fct=xnu*xi;

  for(l=nl; l; l--) {rjtemp=fct*rjl+rjpl;
                     fct-=xi;  rjpl=fct*rjtemp-rjl;  rjl=rjtemp;}

  if(rjl==0.0) rjl=eps;
  f=rjpl/rjl;

  if(x<xmin) {
    x2=x/2;  pimu=pi*xmu;
    if(fabs(pimu)<eps)  fct=1;  else fct=pimu/sin(pimu);
    d=-log(x2);  e=xmu*d;
    if(fabs(e)<eps)     fct2=1; else fct2=sinh(e)/e;

    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=2/pi*fct*(gam1*cosh(e)+gam2*fct2*d);
    e=exp(e);   p=e/(gampl*pi);
    q=1/(e*pi*gammi);
    pimu2=pimu/2;
    if(fabs(pimu2)<eps) fct3=1; else fct3=sin(pimu2)/pimu2;
    r=pi*pimu2*fct3*fct3;
    c=1;  d=-x2*x2;  sum=ff+r*q;  sum1=p;

    for(i=1, salir=1; i<=maxit && salir; i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);   c*=d/i;
      p/=(i-xmu);   q/=(i+xmu);   del=c*(ff+r*q);
      sum+=del;   del1=c*p-i*del;
      sum1+=del1;
      if(fabs(del)<(1+fabs(sum))*eps) salir=0; else salir=1;
    }

    if(salir)  on_warning("besselJY", "series failed to converge");

    rymu=-sum;  ry1=-sum1*xi2;  rymup=xmu*xi*rymu-ry1;
    rjmu=w/(rymup-f*rymu);

  } else {

    a=0.25-xmu2;  p=-xi/2;  q=1;  br=2*x;  bi=2;
    fct=a*xi/(p*p+q*q);  cr=br+q*fct;  ci=bi+p*fct;
    den=br*br+bi*bi;  dr=br/den;  di=-bi/den;  dlr=cr*dr-ci*di;
    dli=cr*di+ci*dr;  temp=p*dlr-q*dli;  q=p*dli+q*dlr;
    p=temp;
    for(i=2, salir=1; i<=maxit && salir; i++) {
      a+=2*(i-1);  bi+=2;  dr=a*dr+br;  di=a*di+bi;
      if(fabs(dr)+fabs(di)<fpmin) dr=fpmin;
      fct=a/(cr*cr+ci*ci);  cr=br+cr*fct;  ci=bi-ci*fct;
      if(fabs(cr)+fabs(ci)<fpmin) cr=fpmin;
      den=dr*dr+di*di;  dr=dr/den;  di=-di/den;
      dlr=cr*dr-ci*di;  dli=cr*di+ci*dr;  temp=p*dlr-q*dli;
      q=p*dli+q*dlr;  p=temp;
      if(fabs(dlr-1)+fabs(dli)<eps) salir=0; else salir=1;
    }
    if(salir)  on_warning("besselJY", "cf2 failed");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    if(rjl<0) rjmu=-rjmu;
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }

  fct=rjmu/rjl;
  rj=rjl1*fct;
  rjp=rjp1*fct;
  for(i=1; i<=nl; i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;  ry1=rytemp;
  }

  ry=rymu;  ryp=xnu*xi*rymu-ry1;

  return 0;
}

// **************************************************************************
// *** I_n(x), K_n(x), I'_n(x), K'_n(x),  real n,  real x                 ***
// **************************************************************************

void besselIK(numero xnu, numero x, numero &ri, numero &rk, numero &rip,
              numero &rkp)
{
  int maxit=100000;
  int i,l,nl;
  numero xmin=2.0, eps=1.0e-10, fpmin=1.0e-30,
         a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
         gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
         ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
  
  if(x<=0)  on_error("besselIK", "wrong argument x =", x, "");
  if(xnu<0)  on_error("besselIK", "wrong order n =", xnu, "");
  nl=(int)(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  if (h < fpmin) h=fpmin;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=1;i<=maxit;i++) {
    b += xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) < eps) break;
  }
  if (i > maxit)  on_warning("besselIK", "x too large: try asymptotic expan.");
  ril=fpmin;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl;l>=1;l--) {
    ritemp=fact*ril+ripl;
    fact -= xi;
    ripl=fact*ritemp+ril;
    ril=ritemp;
  }
  f=ripl/ril;
  if (x < xmin) {
    x2=0.5*x;
    pimu=pi*xmu;
    fact = (fabs(pimu) < eps ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < eps ? 1.0 : sinh(e)/e);
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for (i=1;i<=maxit;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*ff;
      sum += del;
      del1=c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*eps) break;
    }
    if (i > maxit) on_warning("besselIK", "series failed to converge");
    rkmu=sum;
    rk1=sum1*xi2;
  } else {
    b=2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=2;i<=maxit;i++) {
      a -= 2*(i-1);
      c = -a*c/i;
      qnew=(q1-b*q2)/a;
      q1=q2;
      q2=qnew;
      q += c*qnew;
      b += 2.0;
      d=1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels=q*delh;
      s += dels;
      if (fabs(dels/s) < eps) break;
    }
    if (i > maxit)  on_warning("besselIK", "cf2 failed");
    h=a1*h;
    rkmu=sqrt(pi/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu-rk1;
  rimu=xi/(f*rkmu-rkmup);
  ri=(rimu*ril1)/ril;
  rip=(rimu*rip1)/ril;
  for (i=1;i<=nl;i++) {
    rktemp=(xmu+i)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  rk=rkmu;
  rkp=xnu*xi*rkmu-rk1;
}

// **************************************************************************
// *** J_n(x),  int n,  real x                                            ***
// **************************************************************************

numero besselJ0(numero x)
{
  numero  ax,xx,z,y,ans1,ans2;

  if(x==0)  return 1;

  if(x<8 && x>-8)  {
      y = sqr(x);
      ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      return ans1/ans2;
  }

  ax=x<0?-x:x; z = 8.0/ax; y = sqr(z); xx = ax-0.785398164;
  ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
     +y*(-0.2073370639e-5+y*0.2093887211e-6)));
  ans2 = -0.1562499995e-1+y*(0.1430488765e-3
     +y*(-0.6911147651e-5+y*(0.7621095161e-6
     -y*0.934945152e-7)));
  return  sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
}

// --------------------------------------------------------------------------

numero besselJ1(numero x)
{
  numero  ax,xx,z,y,ans1,ans2;

  if (x==0)  return 0;
  if (x<8 && x>-8)  {
    y = sqr(x);
    ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
       +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74
       +y*(99447.43394+y*(376.9991397+y*1.0))));
    return  ans1/ans2;
  }

  ax = x<0?-x:x; z = 8.0/ax; y = sqr(z); xx = ax-2.356194491;
  ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
     +y*(0.2457520174e-5+y*(-0.240337019e-6))));
  ans2 = 0.04687499995+y*(-0.2002690873e-3
     +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
  return  sqrt(0.636619772/ax)*(cos(xx)*ans1
          -z*sin(xx)*ans2)*(x>0?1:-1);
}

// --------------------------------------------------------------------------

numero besselJJ(int n, numero x)
{
  if(n<0) {if(n%2) return -besselJJ(-n,x); else return besselJJ(-n,x);}
  if(ABS(x)<1e-12) {if(n==0) return 1; else return 0;}
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return J;
}

numero besselYY(int n, numero x)
{
  if(n<0) {if(n%2) return -besselYY(-n,x); else return besselYY(-n,x);}
  if(ABS(x)<1e-12) {if(n==0) return 1; else return 0;}
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return Y;
}

numero besselJp(int n, numero x)
{
  if(n<0) {if(n%2) return -besselJp(-n,x); else return besselJp(-n,x);}
  if(ABS(x)<1e-12) {if(n==1)  return 0.5;  else
                   {if(n==-1) return -0.5; else return 0;}}
  numero J,Y,Jp,Yp;  besselJY(n,x,J,Y,Jp,Yp);
  return Jp;
}

numero besselJ(int n, numero x)
{
  int iacc=40, j,jsum,m;
  numero  bigno=1.0e10, bigni=1.0e-10,
          bj,bjm,bjp,sum,tox,ans=0;

  if(ABS(x)<1e-12) {if(n==0) return 1; else return 0;}
  if(n==0)  return  besselJ0(x);
  if(n==1)  return  besselJ1(x);
  if(n<0)   {if(n%2)  return -besselJ(-n,x);  else  return besselJ(-n,x);}

  tox = 2.0/x;
  if(x>n)  {
    bjm=besselJ0(x);
    bj=besselJ1(x);
    for(j=1; j<n; j++) {
      bjp=j*tox*bj-bjm;
      bjm=bj;
      bj=bjp;
    }
    return  bj;
  }

  m= 2*((n+ int(floor(sqrt(1.0*(iacc*n)))) ) / 2);
  jsum=0;  sum=bjp=0;  bj=1;
  for(j=m; j>=1; j--) {
    bjm=j*tox*bj-bjp;
    bjp=bj;
    bj=bjm;
    if(bj>bigno || bj<-bigno) {
      bj=bj*bigni;
      bjp=bjp*bigni;
      ans=ans*bigni;
      sum=sum*bigni;
    }
    if(jsum)  sum = sum+bj;
    jsum = 1-jsum;
    if(j==n)  ans=bjp;
  }
  sum = 2.0*sum-bj;
  return  ans/sum;
}


// **************************************************************************
// *** j_l(x), y_l(x),  int l,  real x                                    ***
// **************************************************************************

void besseljy(int l, numero x, numero &j, numero &y)
{
  numero jp,yp;
  numero factor=sqrt(pi/(2*x));

  besselJY(l+0.5,x,j,y,jp,yp);

  j*=factor;  y*=factor;
}

// --------------------------------------------------------------------------

numero besselj(int l, numero x)
{
  if(ABS(x)<1e-10)  {if(l==0) return 1.0;
                     if(l==1) return x/3;
                     if(l==2) return x*x/15;  return 0;}
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return jl;
}

// --------------------------------------------------------------------------

numero bessely(int l, numero x)
{
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return yl;
}

// --------------------------------------------------------------------------

void besseljyp(int l, numero x, numero &jp, numero &yp)
{
  numero j,y;
  numero factor=sqrt(pi/(2*x));

  besselJY(l+0.5,x,j,y,jp,yp);

  jp=factor*(jp-j/(2*x));  yp=factor*(yp-y/(2*x));
}


// **************************************************************************
// *** I_n(x), K_n(x), I'_n(x), K'_n(x),  int n,  real x                  ***
// **************************************************************************

numero besselI0(numero x)
{
  numero  ax,y;

  if(-3.75<x && x<3.75) {
    y=x/3.75;  y=y*y;
    return  1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*
                  (0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    ax=(x<0)?-x:x;  y=3.75/ax;
    return  (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
              +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2))))))));
  }
}

// --------------------------------------------------------------------------

numero besselI1(numero x)
{
  numero  ax,y,ans;

  if(-3.75<x && x<3.75) {
    y=x/3.75;  y=y*y;
    return  x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
               +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    ax=(x<0)?-x:x;  y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                  +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    return  (exp(ax)/sqrt(ax))*ans;
  }
}

// --------------------------------------------------------------------------

numero besselK0(numero x)
{
  numero  y;

  if(x<=2) {
    y=x*x/4;
    return  -log(x/2.0)*besselI0(x)+(-0.57721566+y*(0.42278420
               +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
               +y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y=2/x;
    return  (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
               +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
               +y*(-0.251540e-2+y*0.53208e-3))))));
  }
}

// --------------------------------------------------------------------------

numero besselK1(numero x)
{
  numero  y;

  if(x<=2) {
      y=x*x/4;
      return  log(x/2.0)*besselI1(x)+(1.0/x)*(1.0+y*(0.15443144
                  +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
                  +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
      y=2/x;
      return  (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
                +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
                 +y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
}

// --------------------------------------------------------------------------

numero besselI(int n, numero x)
{
  if(ABS(x)<1e-20)  {if(n==0) return 1.0;  return 0.0;}
  int     iacc=40, j,m;
  numero  bigno=1.0e10,
          bigni=1.0e-10,
          bi=1, bim, bip=0, tox=2/x, ans=0;

  if(n==0)  return  besselI0(x);
  if(n==1)  return  besselI1(x);
  if(n<0)   n=-n;

  m=2*(n+int(floor(sqrt(iacc*n))));

  for(j=m; j>=1; j--) {
    bim=bip+j*tox*bi;    bip=bi;    bi=bim;
    if(bi<-bigno || bigno<bi) {
         ans=ans*bigni;
         bi=bi*bigni;
         bip=bip*bigni;
    }
    if(j==n)  ans=bip;
  }

  return  ans*besselI0(x)/bi;
}

// --------------------------------------------------------------------------

numero besselK(int n, numero x)
{
  numero  tox,bkp,bkm,bk;
  int  j;

  if(n==0)  return  besselK0(x);
  if(n==1)  return  besselK1(x);
  if(n<0)   n=-n;

  tox=2/x;
  bkm=besselK0(x);
  bk=besselK1(x);

  for(j=1; j<n; j++) {
    bkp=bkm+j*tox*bk;    bkm=bk;    bk=bkp;
  }

  return  bk;
}

// --------------------------------------------------------------------------

numero besselIp(int n, numero x)
{
  if(n==0)  return besselI(1,x);
  return  (besselI(n-1,x)+besselI(n+1,x))/2;
}

// --------------------------------------------------------------------------

numero besselKp(int n, numero x)
{
  if(n==0)  return -besselK(1,x);
  return -(besselK(n-1,x)+besselK(n+1,x))/2;
}


// **************************************************************************
// *** complex functions                                                  ***
// **************************************************************************

#ifdef jga_complex

// **************************************************************************
// *** I_n(x), K_n(x), I'_n(x), K'_n(x),  int n,  complex x               ***
// **************************************************************************

int besselIKexp_flag=1;

complex besselI0(complex x)
{
  complex  y,eax,deax;  if(real(x)<0) on_error("besselI0", "real(x)<0");
  if(besselIKexp_flag) {eax=exp(x); deax=1;} else {eax=1; deax=exp(-x);}

  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
                +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))))*deax;
  } else {
    y=3.75/x;
    return  (eax/sqrt(x))*(0.39894228+y*(0.1328592e-1
              +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2))))))));
} }

// --------------------------------------------------------------------------

complex besselI1(complex x)
{
  complex  y,ans,eax,deax;  if(real(x)<0) on_error("besselI1", "real(x)<0");
  if(besselIKexp_flag) {eax=exp(x); deax=1;} else {eax=1; deax=exp(-x);}

  if(mod(x)<3.75) {
    y=x/3.75;  y=y*y;
    return  (x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
                   +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))))*deax;
  } else {
    y=3.75/x;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
                  +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    return  (eax/sqrt(x))*ans;
} }

// --------------------------------------------------------------------------

complex besselK0(complex x)
{
  complex  y,eax,deax;
  if(besselIKexp_flag) {eax=exp(-x); deax=1;} else {eax=1; deax=exp(x);}

  if(mod(x)<=2) {
    y=x*x/4;
    return  (-log(x/2.0)*besselI0(x)*deax+(-0.57721566+y*(0.42278420
               +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
               +y*(0.10750e-3+y*0.74e-5)))))))*deax;
  } else {
    y=2/x;
    return  eax/sqrt(x)*(1.25331414+y*(-0.7832358e-1
               +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
               +y*(-0.251540e-2+y*0.53208e-3))))));
} }

// --------------------------------------------------------------------------

complex besselK1(complex x)
{
  complex  y,eax,deax;
  if(besselIKexp_flag) {eax=exp(-x); deax=1;} else {eax=1; deax=exp(x);}

  if(mod(x)<=2) {
      y=x*x/4;
      return  (log(x/2.0)*besselI1(x)*deax+(1.0/x)*(1.0+y*(0.15443144
                  +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
                  +y*(-0.110404e-2+y*(-0.4686e-4))))))))*deax;
  } else {
    y=2/x;
    return  eax/sqrt(x)*(1.25331414+y*(0.23498619
              +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
              +y*(0.325614e-2+y*(-0.68245e-3)))))));
} }

// --------------------------------------------------------------------------

complex besselI(int n, complex x)
{
  if(mod(x)<1e-20)  {if(n==0) return 1.0;  return 0.0;}
  int     iacc=40, j,m;
  numero   bigno=1.0e10, bigni=1.0e-10;
  complex   bi(1,0), bim, bip(0,0), tox, ans=0;   tox=2/x;

  if(n<0)   n=-n;
  if(n==0)  return  besselI0(x);
  if(n==1)  return  besselI1(x);

  m=2*(n+int(floor(sqrt(iacc*n))));

  for(j=m; j>=1; j--) {
    bim=bip+j*tox*bi;    bip=bi;    bi=bim;
    if(bigno<mod(bi)) {
         ans=ans*bigni;
         bi=bi*bigni;
         bip=bip*bigni;
    }
    if(j==n)  ans=bip;
  }

  return  ans*besselI0(x)/bi;
}

// --------------------------------------------------------------------------

complex besselK(int n, complex x)
{
  complex  tox,bkp,bkm,bk;
  int  j;

  if(n<0)   n=-n;
  if(n==0)  return  besselK0(x);
  if(n==1)  return  besselK1(x);

  tox=2/x;
  bkm=besselK0(x);
  bk=besselK1(x);

  for(j=1; j<n; j++) {
    bkp=bkm+j*tox*bk;    bkm=bk;    bk=bkp;
  }

  return  bk;
}

// --------------------------------------------------------------------------

complex besselIp(int n, complex x)
{
  if(n==0)  return besselI(1,x);
  return  (besselI(n-1,x)+besselI(n+1,x))/2;
}

// --------------------------------------------------------------------------

complex besselKp(int n, complex x)
{
  if(n==0)  return -besselK(1,x);
  return -(besselK(n-1,x)+besselK(n+1,x))/2;
}


// **************************************************************************
// *** j_l(x), y_l(x), j'_l(x), y'_l(x),                                  ***
// *** h+_l(x), h+'_l(x)                    int n,  real/complex x        ***
// **************************************************************************

complex besselh(int l, numero x)
{
  numero jl,yl;
  besseljy(l,x,jl,yl);
  return i_c*jl-yl;
}

// --------------------------------------------------------------------------

complex besselhp(int l, numero x)
{
  numero j,y,jp,yp;
  numero factor=sqrt(pi/(2*x));
  besselJY(l+0.5,x,j,y,jp,yp);
  return factor*((i_c*jp-yp)-(i_c*j-y)/(2*x));
}

// --------------------------------------------------------------------------

complex ak(int k,int n)  // n=n+1/2 -- auxiliary funct. in bessely|j|h
{
  // if(k<=n-1/2) return factorial(n+k)/(pow(2,k)*factorial(k)*factorial(n-k));
  if(k<=n-1/2) return sqr(((fact2[n+k]/pow(2,k/2.0))/fact2[k])/fact2[n-k]);
  else         return 0;
}

complex bessely(int l, complex z)
{
/*
  int k,kmax; complex bb, bb1=0, bb2=0;
  kmax=int(l/2);
  for(k=0; k<=kmax; k++) { bb1+=pow(-1,k)*ak(2*k,l)/pow(z,2*k+1); }
  kmax=int((l-1)/2);
  for(k=0; k<=kmax; k++) { bb2+=pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+2); }

  bb=-cos(z-l*pi/2)*bb1+sin(z-l*pi/2)*bb2;
  return bb;
*/
  // -- old implementation
  if(l==0)  return -cos(z)/z;
  if(l==1)  return -(cos(z)/z+sin(z))/z;
  //            return (2*l-1)/z*bessely(l-1,z)-bessely(l-2,z);
  int n; complex yl_2, yl_1=bessely(0,z), yl=bessely(1,z); // this is faster!
  for(n=2; n<=l; n++) {yl_2=yl_1;  yl_1=yl;  yl=(n+n-1)/z*yl_1-yl_2;}
  return yl;
}

// --------------------------------------------------------------------------

complex besselyp(int l, complex z)
{
/*
  int k,kmax; complex bb, bb1=0, bb2=0, bb3=0, bb4=0;
  kmax=int(l/2);
  for(k=0; k<=kmax; k++) {
    bb1+=-pow(-1,k)*ak(2*k,l)/pow(z,2*k+1);
    bb2+= pow(-1,k)*ak(2*k,l)/pow(z,2*k+2)*(1+2*k);
  }
  kmax=int((l-1)/2);
  for(k=0; k<=kmax; k++) {
    bb3+= pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+2);
    bb4+= pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+3)*(2+2*k);
  }
  bb=cos(z-l*pi/2)*(bb2+bb3)+sin(z-l*pi/2)*(bb1+bb4);
  return bb;
*/
  // -- old implementation
  if(l==0) return -bessely(1,z);
  else     return  bessely(l-1,z)-((l+1)/z)*bessely(l,z);
}

// --------------------------------------------------------------------------

complex besselilh(int l, complex z)
{                                        // i^l h^(+)_l(z): series expansion
  complex val=1, fct=1, fct1=i_c/(2*z);  // Eq.(11.152) of Arfken
  int s;

  init_fact();

  for(s=1; s<=l; s++) {
    fct=fct*fct1/s;
    val=val+fct*fact[l+s]/fact[l-s];
  }

  return  exp(i_c*z)/z * val;
}

// --------------------------------------------------------------------------

complex besselilhp(int l, complex z)
{                                        // i^l h'^(+)_l(z): series expansion

  if(l>0) return i_c*besselilh(l-1,z)-(l+1)/z*besselilh(l,z);
          return i_c*besselilh(l+1,z)+    l/z*besselilh(l,z);
}

// --------------------------------------------------------------------------

#ifdef jga_qfint

complex besselj_int_z;
int     besselj_int_l;

// --------------------------------------------------------------------------

complex besselj_int(numero th)
{
  return cos(besselj_int_z*cos(th)) * pow(sin(th),2*besselj_int_l+1.0);
}

// --------------------------------------------------------------------------

complex besselj(int l, complex z)
{
  /*
  int k,kmax; complex bb, bb1=0, bb2=0;
  kmax=int(l/2);
  for(k=0; k<=kmax; k++) bb1+=pow(-1,k)*ak(2*k,l)/pow(z,2*k+1);
  kmax=int((l-1)/2);
  for(k=0; k<=kmax; k++) bb2+=pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+2);
  bb=sin(z-l*pi/2)*bb1+cos(z-l*pi/2)*bb2;
  return bb;
  */
  // -- old implementation
  if(l==0)  {if(!real(z) && !imag(z))  return complex(1,0);
             else                      return sin(z)/z;}
  if(l==1)  {if(!real(z) && !imag(z))  return complex(0,0);
             else                      return (sin(z)/z-cos(z))/z;}
  complex coef(1,0);   int i;
  besselj_int_z=z;
  besselj_int_l=l;
  for(i=1; i<=l; i++)  coef=coef*z/(2*i);
  return coef * qfint(0,pi/2,besselj_int,0.00001,2048);
}

// --------------------------------------------------------------------------

complex besseljp(int l, complex z)
{
  /*
  int k,kmax; complex bb, bb1=0, bb2=0, bb3=0, bb4=0;
  kmax=int(l/2);
  for(k=0; k<=kmax; k++) {
    bb1+= pow(-1,k)*ak(2*k,l)/pow(z,2*k+1);
    bb2+=-pow(-1,k)*ak(2*k,l)/pow(z,2*k+2)*(1+2*k);
  }
  kmax=int((l-1)/2);
  for(k=0; k<=kmax; k++) {
    bb3+=-pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+2);
    bb4+=-pow(-1,k)*ak(2*k+1,l)/pow(z,2*k+3)*(2+2*k);
  }
  bb=sin(z-l*pi/2)*(bb2+bb3)+cos(z-l*pi/2)*(bb1+bb4);
  return bb;
  */
  // -- old implementation
  if(l==0) return -besselj(1,z);
  else     return  besselj(l-1,z)-((l+1)/z)*besselj(l,z);
}

// --------------------------------------------------------------------------

complex besselh(int l, complex z)
{
  /*
  int k; complex bb=0;
  for(k=0; k<=l; k++) bb+=pow(i_c,k-l-1)*ak(k,l)/pow(z,k+1);
  bb=i_c*exp(i_c*z)*bb;  // i_c added due to Messiah's definition (h=i*j-y)   
  return bb;
  */
  // -- old implementation
  if(ABS(imag(z))<infinitesimal)  return besselh(l, real(z));
  return i_c*besselj(l,z) -bessely(l,z);
}

// --------------------------------------------------------------------------

complex besselhp (int l, complex z)
{
  /*
  int k; complex bb=0;
  for(k=0; k<=l; k++) bb+=pow(i_c,k-l-1)*ak(k,l)*(1+k-i_c*z)/pow(z,k+2);
  bb=-i_c*exp(i_c*z)*bb;  // i_c added due to Messiah's definition (h=i*j-y)
  return bb;
  */
  // -- old implementation
  if(ABS(imag(z))<infinitesimal)  return besselhp(l, real(z));
  return i_c*besseljp(l,z)-besselyp(l,z);
}

// --------------------------------------------------------------------------

complex besselh1 (int l, complex z)  {return besselj(l,z) +i_c*bessely(l,z);}
complex besselh1p(int l, complex z)  {return besseljp(l,z)+i_c*besselyp(l,z);}
complex besselh2 (int l, complex z)  {return besselj(l,z) -i_c*bessely(l,z);}
complex besselh2p(int l, complex z)  {return besseljp(l,z)-i_c*besselyp(l,z);}

// --------------------------------------------------------------------------

void besseljh(int l, numero x, numero &j, numero &jp, complex &h, complex &hp)
{
  numero y,yp, x2=2*x;
  numero factor=sqrt(pi/x2);
  besselJY(l+0.5,x,j,y,jp,yp);

  j=factor*j;         y=factor*y;         h=i_c*j-y;      // mind this order
  jp=factor*jp-j/x2;  yp=factor*yp-y/x2;  hp=i_c*jp-yp;
}

// --------------------------------------------------------------------------

void besseljh(int l, complex z, complex &j, complex &jp,
                                complex &h, complex &hp)
{
  complex y,yp, hh,hhp, b=complex(0,imag(z));
  numero  rj,rjp,rjj,rjjp, x=real(z);

  if(ABS(imag(z))<1e-2*ABS(real(z)) && ABS(imag(z))<1e-2) {
    besseljh(l,   real(z), rj, rjp, h, hp);
    besseljh(l+1, real(z), rjj,rjjp,hh,hhp);
    j=rj+b*rjp;  h=h+b*hp;
    jp=rjp+b*(l/x*(rjp-rj/x)-rjjp);
    hp= hp+b*(l/x*( hp- h/x)- hhp);
  } else {
    on_warning("besseljh", "complex z");
    j=besselj(l,z);  jp=besseljp(l,z);
    y=bessely(l,z);  yp=besselyp(l,z);
    h=i_c*j-y;       hp=i_c*jp-yp;
} }

// --------------------------------------------------------------------------

#endif  // qfint

#endif  // complex

#endif  // ******************************************************************
