#include "jga/bessel.h"

// all lengths in units of lambda;
// particle electrostatic polarizability

int M,lmax;  // lmax=N in the presentation
numero *bj=NULL,*gj, *grad;

numero alphal(int l)  // alpha_l
{
  int j;  numero val=0;
  for(j=0; j<M; j++) val+=bj[j]*besselJ(2*l,2*ABS(gj[j]));
  return val;
}

numero RR(void)       // function to be minimized
{
  int l;  numero val=0;
  for(l=-lmax; l<=lmax; l++) val+=sqr(alphal(l)-1);
  return val;
}

void normalize(void)
{
  int j,l;  numero val=0;
  //  for(l=-lmax; l<=lmax; l++) val+=sqr(alphal(l));  val=sqrt(val);
  for(j=0; j<M; j++) val+=sqr(bj[j]);  val=sqrt(val);
  for(j=0; j<M; j++) bj[j]=bj[j]/val;
}

void init(int M_, int lmax_)
{
  M=M_;  lmax=lmax_;
  if(bj!=NULL) {delete[] bj;  delete[] gj;  delete[] grad;}
  bj=new numero[M];  gj=new numero[M];  grad=new numero[2*M];
  int j;  for(j=0; j<M; j++) {    // random initialization
    bj[j]=rand()/(0.0+RAND_MAX);  // with values from 0 to 1
    gj[j]=rand()/(0.0+RAND_MAX);
  }
  normalize();
}

int main(int argc, char **argv)
{
  if(argc<6) on_syntax_error("a.cpp", "a.out M lmax h hh niter\n");
  int    M    =atoi(argv[1]);
  int    lmax =atoi(argv[2]);
  numero h    =atof(argv[3]);  // small step for integration
  numero hh   =atof(argv[4]);  // small step to move along the gradient
  int    niter=atoi(argv[5]);
  numero max,val,tmp;  int n,j,l;

  init(M,lmax);
  //--- steepest descent
  for(n=0; n<niter; n++) {
    max=-1;
    for(j=0; j<M; j++) {  // gradient
      bj[j]+=h;  tmp=RR();                    bj[j]-=h;
      bj[j]-=h;  grad[j+0]=(tmp-RR())/(2*h);  bj[j]+=h;  // derivative of bj
      gj[j]+=h;  tmp=RR();                    gj[j]-=h;
      gj[j]-=h;  grad[j+M]=(tmp-RR())/(2*h);  gj[j]+=h;  // derivative of gj
      val=sqrt(sqr(grad[j])+sqr(grad[j+M]));
      if(val>max) max=val;
    }
    if(max<1e-10) max=1e-10;
    for(j=0; j<M; j++) {       // normalization of the gradient
      bj[j]-=grad[j+0]*hh/max; // and displacement along the gradient
      gj[j]-=grad[j+M]*hh/max;
    }
    if(n==niter-1) {
      printf("*** l, alpha_l -->\n");
      for(l=-lmax; l<=lmax; l++) printf("%d %g\n", l, alphal(l));
      printf("*** j, b_j, g_j -->\n");
      for(j=0; j<M; j++) printf("%d %g %g\n", j, bj[j], gj[j]);
    } else printf("%4d %g\n", n, RR());
  }

  return 0;
}
