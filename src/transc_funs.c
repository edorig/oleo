#include <math.h>
#include <stdio.h>

#include "funcdef.h"
          #include "sysdef.h"
          #include "global.h"
          #include "cell.h"
          #include "eval.h"
          #include "errors.h"
          
          struct value
            {
              int type;
              union vals x;
            };
          
          #define Float	x.c_d
          #define String	x.c_s
          #define Int	x.c_l
          #define Value	x.c_i
          #define Rng	x.c_r
/* Prototypes for our DIY functions. It would be better to use CEPHES, cernlib, NSWC or SLATEC mathematical libraries. But according to repology.org, there is no CEPHES or NSWC package, and SLATEC and CERNLIB are available for only a handful of Linux distros. There does not seem to be an automatic way to detect these mathematical libraries with GNU configure, so for the moment it will be DIY. */ 
double bessi0(double);
double bessk0(double);
double lambertW(double);
double chisquare(double,int);
double student(double,int); 
double kolmogorov(double);
double ellk(double);
double ellec(double);
double ibeta(double,double,double);
double fisherQ(double,double,double);

/* Modified Bessel function I_0(x) see Abramowitz and Stegun Chapter 9 */
double bessi0(double x)
{double y,z;

 if (fabs(x)<3.75) {
   z=(x/3.75)*(x/3.75);
   y=(((((0.0045813*z+0.0360768)*z+0.2659732)*z+1.2067492)*z+3.08899424)*z+3.5156229)*z+1.0;
 } else {
   z=3.75/fabs(x);
   y=(((((((0.00392337*z-0.01647633)*z+0.02635537)*z-0.02057706)*z+0.00916281)*z-0.00157565)*z+0.00225319)*z+0.01328592)*z+0.39894228;
   y*=exp(fabs(x))/sqrt(fabs(x));
 };

 return(y);
}

/* Modified Bessel function K_0(x) (also called MacDonald function) see Abramowitz and Stegun Chapter 9 */

double bessk0(double x) 
{double y,z;
 int i;
 double gamma,u,h;

 gamma=0.5772156649;
 if (fabs(x)<2.0){
   z=(x/2.0)*(x/2.0);
   u=1.0;
   h=0.;
   y=-(log(x/2.0)+gamma)*bessi0(x);
   for (i=1;i<=20;i++){
     u*=z/((float)i*(float)i);
     h+=1.0/((float)i);
     y+=u*h;
   };
 } else {
   z=2.0/fabs(x);
   y=(((((0.00053208*z-0.00251540)*z+0.00587872)*z-0.01062446)*z+0.02189568)*z-0.07832358)*z+1.25331414; 
  y/=exp(fabs(x))*sqrt(x);
 };

 return(y);
}

/* Computes the complete elliptic integral of the second kind E(m)  using approximation 17.3.36 of 
   Abramowitz and Stegun p. 591. Accuracy is 2e-8. Tested against the table 17.1  p. 603  */


double ellec(double m)
{
  double m1;
  double ec; 
  const double a[4]={0.4432514153,
		     0.05260601220,
		     0.04757383546,
		     0.01736506451
  }; 
  const double b[4]={0.24998368310,
		     0.09200180037,
		     0.04069697526,
		     0.00526449639
  };
  if (m>1.0) return (NAN); 
  if (m>=0.0) {
    m1=1.0-m;
  } else {
    m1=1.0/(1.0-m); 
      };
  
  ec=(((a[4]*m1+a[3])*m1+a[2])*m1+a[1])*m1+1.0-(((b[3]*m1+b[2])*m1+b[1])*m1+b[0])*m1*log(m1);
  
  if (m>=0.0) {
    return(ec);
  } else {
    return (sqrt(1.0-m)*ec);
  }
}

/* Computes the complete elliptic integral of the first kind using approximation 17.3.34 of 
   Abramowitz and Stegun p. 591. Accuracy is 2e-8. Tested against the table 17.1  p. 603. m<=1.0  */ 
double ellk(double m)
{
  double m1;
  double k; 
  const double a[5]={1.38629436112,
		     0.09666344259,
		     0.03590092383,
		     0.03742563713,
		     0.01451196212
  }; 
  const double b[5]={0.5,
		     0.12498593597,
		     0.06880248576,
		     0.03328355346,
		     0.00441787012
  };
  if (m>1.0) return(NAN); 
  if (m>=0.0) {
    m1=1.0-m;
  } else {
    m1=1.0/(1.0-m);
  };
  
  k=(((a[4]*m1+a[3])*m1+a[2])*m1+a[1])*m1+a[0]-((((b[4]*m1+b[3])*m1+b[2])*m1+b[1])*m1+b[0])*log(m1);

  if (m>=0.0) {
  return(k);
  } else {
    return (k/sqrt(1.0-m));
  }
}


/* Calculates the principal branch of the Lambert W function, i.e. the unique solution of W(x) exp(W(x))=x with W(x)>=-1 */
double lambertW(double x)
{
  int i;
  double u=0.0;
  /* Lambert W is undefined for x<exp(-1) */ 
  if (x< -exp(-1.0)) return (NAN);
  fprintf(stderr,"In lambertW()\n"); 
  if (x<600.) {
    /* A few steps to get close to the solution faster than with Newton */ 
    for (i=1;i<=5;i++) {
      u=u+(sqrt(1.0+(2.0+u)*x*exp(-u))-1.0-u)/(2.0+u);
    }
    while (fabs(u*exp(u)-x)>1e-13*fabs(x)) {
      /* Newton method */ 
      u=u-(u-x*exp(-u))/(1+u); 
    }
    return(u); 
  } else {
    u=log(x);
    i=0;
    while ((fabs(u*exp(u)-x)>1e-13*fabs(x)) && (i<40)) {
      u=log(x)-log(u);
      i++; 
    }
    return(u); 
  } 
}



/* Calculates Q(chi2|n) cumulative distribution function of the chi squared 
   distribution, defined in Abramowitz and Stegun (A&S) Eq. (26.4.2) p. 940. 
Reproduces table 26.7 p. 978 of A&S */ 
double chisquare(double chi2,int n)
{
int r; 
 double k,q,z,y,w,chi;

 chi=sqrt(chi2); 
 
if (n%2) {
  /* Odd n case: apply A&S (26.4.4) p. 941 */ 
    y=chi;
    k=1.0;
    w=0.0;
    z=exp(-0.5*chi*chi)/sqrt(2.0*M_PI);
    q=(1.0-erf(chi/sqrt(2.0)));   
    for (r=1;r<=((n-1)/2);r++) {
      w+=y/k;  
    y*=(chi*chi);
    k*=(double)(2*r+1);
    }
    return (q+2.0*z*w);
} else {
  /* Even n case apply A&S (26.4.5) */ 
    z=exp(-0.5*chi*chi);
    k=1.0; 
    w=y=z; 
    for (r=1;r<=(n/2-1);r++) { 
      k*=(double)(2*r);
    y*=(chi*chi);
    w+=y/k; 
    }
    return (w);
    }
}

/* Computes the cumulative distribution function of Student's t-distribution 
defined in Eq. (26.7.1) of Abramowitz and Stegun (A&S) or Eq. (14.4.2) of Ventsel, Elena Sergeevna. Probability theory (first steps). Mir Publishers, 1982. 
Tested against table 4 of the latter reference. */ 
double student(double t,int n)
{
  double theta;
  double y,z,s,w;
  int i; 
  theta = atan(t/sqrt((double)n));
  y=2.0*theta/M_PI;
  if (n==1) return (y);
  z=pow(cos(theta),2);
  /* Odd case Eq. (26.7.3) of A&S */ 
  if (n%2) { 
    w=1.0;
    s=1.0; 
    for (i=1;i<=((n-3)/2);i++)
      {
	w*=(((double)(2*i)/(double)(2*i+1))*z);
	s+=w; 
      }
    return (y+s*sin(2.0*theta)/M_PI);
  } else {
    /* Even case Eq. (26.7.4) of A&S */
    w=1.0;
    s=1.0;
    for (i=1;i<=(n/2-1);i++) {      
      w*=((double)(2*i-1)/(double)(2*i))*z;
      s+=w;
    }
    return(s*sin(theta)); 
  }
    
}

/* Algorithm for Kolmogorov distribution: Paul van Mulbregt arXiv:1803.00426 */
double kolmogorov(double x)
{
  double q,s;
  int i,k;
  s=0; 
  if (x>=0.8) {
    q=exp(-2.0*x*x);
    k=1; 
    for (i=1;i<=20;i++) {
      s+=(double)k*pow(q,i*i);
      k=-k;
    }
    return (2*s); 
  } else {
    q=exp(-M_PI*M_PI/(8.0*x*x));
    for (i=1;i<=41;i+=2) {
      s+=pow(q,i*i);
    }
return(1.0-(sqrt(2*M_PI)*s/x)); 
  }

}
/* Incomplete beta function See Abramowitz and Stegun p. 263 Eq. (6.6.1) */ 
double ibeta(double a,double b,double x)
{
  int i;
  double s,y;
  /* Eq. (6.6.8) and chapter 15 for series expansion */ 
  if (x<=0.5) {
    s=1.0/a;
    y=1.0;
    for (i=1;i<=30;i++) {
      y*=(1.0-b/(double)i)*x;
      s+=y/((double)i+a);
    }
    return (s*pow(x,a));
  } else {
    /* Eqs. (6.6.2) and (6.6.3) */  
    s=1.0/b;
    y=1.0; 
    for (i=1;i<=30;i++) 
      {
	y*=(1.0-a/(double)i)*(1.0-x);
	s+=y/((double)i+b);
      }
    return (exp(lgamma(a)+lgamma(b)-lgamma(a+b))-s*pow(1.0-x,b));
  }
	    
}

/* Ckecked against the results in table 26.9 p. 986 of A&S */ 
double fisherQ(double nu1,double nu2,double F)
{
  double x,z;
  x=nu2/(nu2+nu1*F);
  z=exp(lgamma(0.5*nu1)+lgamma(0.5*nu2)-lgamma(0.5*(nu1+nu2)));
  return (ibeta(0.5*nu2,0.5*nu1,x)/z); 
}

/* Functions from the C Library in SVr4,  4.3BSD,  POSIX.1-2001,
   and POSIX.1-2008 */ 

void do_cbrt(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=cbrt(arg0);
}

void do_cosh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=cosh(arg0);
}
void do_sinh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=sinh(arg0);
}

void do_tanh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=tanh(arg0);
}

void do_acosh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=acosh(arg0);
  p->type= TYP_FLT;
}

void do_asinh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=asinh(arg0);
  // p->Type= TYP_FLT;
}

void do_atanh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=atanh(arg0);
  // p->Type= TYP_FLT;
}

void do_lambertw(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=lambertW(arg0); 
}
void do_erf(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=erf(arg0);
  // p->Type= TYP_FLT;
}

void do_chi2(p)
     struct value *p; 
{
  double arg0=p[0].Float;
  int arg1=p[1].Int; 
  p->Float=chisquare(arg0,arg1);
  p->type=TYP_FLT;
}

void do_student(p)
     struct value *p; 
{
  double arg0=p[0].Float;
  int arg1=p[1].Int; 
  p->Float=student(arg0,arg1);
  p->type=TYP_FLT; 
}

void do_kolmogorov(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=kolmogorov(arg0);
  p->type=TYP_FLT; 
}

void do_gamma(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=tgamma(arg0);
  p->type= TYP_FLT;
}

void do_beta(p)
     struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  p->Float=exp(lgamma(arg0)+lgamma(arg1)-lgamma(arg0+arg1));
	       p->type=TYP_FLT;
}

void do_ibeta(p)
     struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  double arg2=p[2].Float;
  p->Float=ibeta(arg0,arg1,arg2);
  p->type=TYP_FLT;
} 

void do_fisherQ(p)
  struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  double arg2=p[2].Float;
  p->Float=fisherQ(arg0,arg1,arg2);
  p->type=TYP_FLT;
} 

void do_binomial(p)
     struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  p->Float=exp(lgamma(arg0+1.0)-lgamma(arg1+1.0)-lgamma(arg0-arg1+1.0));
  p->type=TYP_FLT; 
    }

void do_ellk(p)
     struct value *p;
{double arg0=p[0].Float;
  p->Float=ellk(arg0);
  p->type=TYP_FLT;
}

void do_ellec(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=ellec(arg0);
  p->type=TYP_FLT;
}

void do_j0(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=j0(x); 
}

void do_j1(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=j1(x); 
}

void do_y0(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=y0(x); 
}

void do_y1(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=y1(x); 
}

void do_jn(p)
     struct value *p;
{
  double x=p[1].Float;
  int n=p[0].Int;
  /*  fprintf(stderr,"j_n(%d,%f) called.\n",n,x); */ 
  p->Float=jn(n,x);
  p->type=TYP_FLT; 
}

void do_yn(p)
     struct value *p;
{
  double x=p[1].Float;
  int n=p[0].Int;
  /*  fprintf(stderr,"y_n(%d,%f) called.\n",n,x);  */ 
  p->Float=yn(n,x);
  p->type=TYP_FLT; 
}

/* Modified Bessel functions */

void do_i0(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=bessi0(x); 
}
void do_k0(p)
     struct value *p;
{
  double x=p[0].Float;
  p->Float=bessk0(x); 
}


struct function transc_funs[]=
  {
    {C_FN1,X_A1,"F",do_cbrt,"cbrt"},
    {C_FN1,X_A1,"F",do_cosh,"cosh"},
    {C_FN1,X_A1,"F",do_sinh,"sinh"},
    {C_FN1,X_A1,"F",do_tanh,"tanh"},
    {C_FN1,X_A1,"F",do_acosh,"acosh"},
    {C_FN1,X_A1,"F",do_asinh,"asinh"},
    {C_FN1,X_A1,"F",do_atanh,"atanh"},
    {C_FN1,X_A1,"F",do_lambertw,"lambw"}, 
    {C_FN1,X_A1,"F",do_erf,"erf"},
    {C_FN2,X_A2,"FI",do_chi2,"chi2Q"},
    {C_FN2,X_A2,"FI",do_student,"studentA"},
    {C_FN1,X_A1,"F",do_kolmogorov,"kolmQ"},
    {C_FN1,X_A1,"F",do_gamma,"gamma"},
    {C_FN2,X_A2,"FF",do_beta,"beta"},
    {C_FN2,X_A2,"FF",do_binomial,"binomial"},
    {C_FN3,X_A3,"FFF",do_ibeta,"ibeta"},
    {C_FN3,X_A3,"FFF",do_fisherQ,"fisherQ"},
    {C_FN1,X_A1,"F",do_j0,"besj0"},
    {C_FN1,X_A1,"F",do_j1,"besj1"},
    {C_FN1,X_A1,"F",do_y0,"besy0"},
    {C_FN1,X_A1,"F",do_y1,"besy1"},
    {C_FN2,X_A2,"IF",do_jn,"besjn"},
    {C_FN2,X_A2,"IF",do_yn,"besyn"},
    {C_FN1,X_A1,"F",do_i0,"besi0"},
    {C_FN1,X_A1,"F",do_k0,"besk0"},
    {C_FN1,X_A1,"F",do_ellk,"ellK"},
    {C_FN1,X_A1,"F",do_ellec,"ellE"},
    {0, 0, "", 0, 0},
          };

int init_transc_function_count(void) 
{
        return sizeof(transc_funs) / sizeof(struct function) - 1;
}
