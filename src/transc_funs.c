#include <math.h>

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
/* Prototypes for modified Bessel functions */ 
double bessi0(double);
double bessk0(double);

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

void do_erf(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=erf(arg0);
  // p->Type= TYP_FLT;
}

void do_gamma(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=tgamma(arg0);
  // p->Type= TYP_FLT;
}

void do_beta(p)
     struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  p->Float=exp(lgamma(arg0)+lgamma(arg1)-lgamma(arg0+arg1));
	       p->type=TYP_FLT;
}

void do_binomial(p)
     struct value *p;
{double arg0=p[0].Float;
  double arg1=p[1].Float;
  p->Float=exp(lgamma(arg0+1.0)-lgamma(arg1+1.0)-lgamma(arg0-arg1+1.0));
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
    {C_FN1,X_A1,"F",do_erf,"erf"},
    {C_FN1,X_A1,"F",do_gamma,"gamma"},
    {C_FN2,X_A2,"FF",do_beta,"beta"},
    {C_FN2,X_A2,"FF",do_binomial,"binomial"}, 
    {C_FN1,X_A1,"F",do_j0,"besj0"},
    {C_FN1,X_A1,"F",do_j1,"besj1"},
    {C_FN1,X_A1,"F",do_y0,"besy0"},
    {C_FN1,X_A1,"F",do_y1,"besy1"},
    {C_FN2,X_A2,"IF",do_jn,"besjn"},
    {C_FN2,X_A2,"IF",do_yn,"besyn"},
    {C_FN1,X_A1,"F",do_i0,"besi0"},
    {C_FN1,X_A1,"F",do_k0,"besk0"},
    {0, 0, "", 0, 0},
          };

int init_transc_function_count(void) 
{
        return sizeof(transc_funs) / sizeof(struct function) - 1;
}
