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
  p->Type= TYP_FLT;
}

void do_sinh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=sinh(arg0);
  p->Type= TYP_FLT;
}

void do_tanh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=tanh(arg0);
  p->Type= TYP_FLT;
}
void do_cosh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=cosh(arg0);
  p->Type= TYP_FLT;
}

void do_sinh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=sinh(arg0);
  p->Type= TYP_FLT;
}

void do_tanh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=tanh(arg0);
  p->Type= TYP_FLT;
}

void do_acosh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=acosh(arg0);
  p->Type= TYP_FLT;
}

void do_asinh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=asinh(arg0);
  p->Type= TYP_FLT;
}

void do_atanh(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=atanh(arg0);
  p->Type= TYP_FLT;
}

void do_erf(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=erf(arg0);
  p->Type= TYP_FLT;
}

void do_gamma(p)
     struct value *p;
{
  double arg0=p[0].Float;
  p->Float=gamma(arg0);
  p->Type= TYP_FLT;
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
    {0, 0, "", 0, 0},
          };
