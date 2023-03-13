#include <math.h>
#ifndef MAXDOUBLE
#define MAXDOUBLE 1e300
#endif
//Замечание: Возвращение значения MAXDOUBLE может быть признаком ошибки
	
/*************************************************************************
Вычисление функций Г(х) и обратной к ней
*************************************************************************/
// точность - 10-11 значащих цифр
double gamma(double x)
{
 if(x<=0)return MAXDOUBLE;
 if(x>171)return MAXDOUBLE;
 double f=1.,z;
  while(x<7.0)
  {
   f*=x;
   x+=1.0;
  }
 z=1.0/(x*x);
 return 2.506628274631*exp(x*(log(x)-1.)+(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z+0.083333333333333)/x)/(f*sqrt(x));
}
double inv_gamma(double x)
{
 double f=1.,z;
  while(x<7.0)
  {
   f*=x;
   x+=1.0;
  }
 z=1.0/(x*x);
 return f*sqrt(x)*exp(-(x*(log(x)-1.)+(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z+0.083333333333333)/x))/2.506628274631;
}
/*************************************************************************
Вычисление полигамма-функции d^n/dx^n(digamma(x))
*************************************************************************/
// точность: 10-11 значащих цифр для n=1, 11 значащих цифр для n=2,3
double polygamma(int n,double x)
{
    if((n<1)||(x<=0))return MAXDOUBLE;
	double f=0.,t,x2;
    int sfct=1,i;
    for(i=1;i<n;i++)sfct*=-i;
    t= n<4 ? 6.:9.*n;
	if(x<t){ do{f+=pow(x,-n-1);x++;}while(x<t);
    f*=n;}
    x2=x*x;
    return sfct*(f+(((((( (2.08767569878681e-08-5.28419013868749e-10*
    		(n+10.)*(n+11.)/x2)*(n+8.)*(n+9.)/x2-0.000000826719576719577)*
            (n+6.)*(n+7.)/x2+0.0000330687830687831)*(n+4.)*(n+5.)/x2-
            0.00138888888888889)*(n+2.)*(n+3.)/x2+0.0833333333333333)*
            (n+1.)/x+0.5)*n/x+1.)/pow(x,n));
}
//trigamma, pentagamma
double trigamma(double x){return polygamma(1,x);}
double pentagamma(double x){return polygamma(3,x);}
/****************************************************************************
Неполная Г-функция
****************************************************************************/
double x2aux(double x,int upper,double v,double vfac,int* ifail)
{
 const  double acc=1e-15, big=1.3e154;
 double _x2aux;
 if(x==0.0)
 {
   _x2aux= upper ? 1.0 : 0.0; *ifail=0;
 }
 else
 {
  double gin,factor,a,b,term,rn,dif,pn[6];
  int i,undone;
  factor=v*log(x)-x-vfac;
  if(x>1.0&&x>=v)
  {
   a=2.0-v;b=a+x;
   term=pn[0]=1.0;pn[1]=x;pn[2]=x+1.0;pn[3]=x*b;
   gin=pn[2]/pn[3];undone=1;
   b+=2.0;
   while(undone)
   {
    pn[4]=b*pn[2]-a*term*pn[0];pn[5]=b*pn[3]-a*term*pn[1];
    if(pn[5]!=0)
    {
     rn=pn[4]/pn[5];dif=fabs(gin-rn);
     if(dif<=acc&&dif<=acc*rn)undone=0;else gin=rn;
    }
    if(undone)
    {
     for(i=0;i<4;i++)pn[i]=pn[i+2];
     if(fabs(pn[4])>=big)for(i=0;i<4;i++)pn[i]=pn[i]/big;
     a+=1.0;term+=1.0;
    }
    b+=2.0;
   }
   gin=exp(factor+log(gin));
   *ifail= !upper&&gin>0.9 ? 5 : 0;
   _x2aux= upper ? gin : 1.0-gin;
  }
  else
  {
   gin=term=1.0;
   rn=v+1.0;
   do
   {
    term*=x/rn;gin+=term;
    rn+=1.0;
   }while(term>acc);
   gin=exp(factor)*gin/v;
   *ifail= upper&&gin>0.9 ? 5 : 0;
   _x2aux= upper ? 1.0-gin : gin;
  }
 }
 return(_x2aux);
}
//неполная гамма-функция гамма-малое г(а,х)
double igamma(double a,double x)
{
 if((a<=0)||(x<0))return MAXDOUBLE;
 int ifail;double gam=gamma(a);
 return x2aux(x,0,a,log(gam),&ifail)*gam;
}
/*************************************************************************
Вычисление бета-функции
*************************************************************************/
double beta(double x,double y)
{
	double g1=gamma(x),g2=gamma(y);
    if(g1==MAXDOUBLE||g2==MAXDOUBLE)return MAXDOUBLE;
    return g1*g2*inv_gamma(x+y);
}
/***********************************************************************************************
Вычисление интегральной показательной функции E1(x) через аппроксимацию многочленами
***********************************************************************************************/
//0<x<1 - точность 2е-7, x>=1 - точность 5е-5
double E1(double x)
{
 double a[]={-0.57721566,0.99999193,-0.24991055,0.05519968,-0.00976004,0.00107857},
 b[]={2.334733,0.250621,3.330657,1.681534};
 if(x>0&&x<1)return -log(x)+a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x+a[5]*x*x*x*x*x;
 else if(x>=1)return exp(-x)*(x*x+b[0]*x+b[1])/(x*x+b[2]*x+b[3])/x;
 else return MAXDOUBLE;
}