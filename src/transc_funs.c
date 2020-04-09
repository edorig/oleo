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
    {C_FN1,X_A1,"F",do_j0,"besj0"},
    {C_FN1,X_A1,"F",do_j1,"besj1"},
    {C_FN1,X_A1,"F",do_y0,"besy0"},
    {C_FN1,X_A1,"F",do_y1,"besy1"},
    {C_FN2,X_A2,"IF",do_jn,"besjn"},
    {C_FN2,X_A2,"IF",do_yn,"besyn"},
    {0, 0, "", 0, 0},
          };

int init_transc_function_count(void) 
{
        return sizeof(transc_funs) / sizeof(struct function) - 1;
}
