#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <cmath>
#include <math.h>

#include "Integ.h"

/*
gaulegÖ÷ÒªÊÇ¸ø³öÇø¼äµÄ»®·ÖºÍÈ¨ÖØ
 
*/

#define EPS 3.0e-15
void gauleg(double x1, double x2, double x[], double w[], int n)  
//w[]Ó¦¸ÃÔõÃ´¸ø¶¨£¿x[]ºÍw[]¶¼ÊÇÍâ²¿Ëù¸ø³öµÄ¿ÕÊ¸Á¿
//Return arrays x[1..n] and w[1..n] of length n, containing the absciaasa and weights of the Gauss-Legendre n-point quadrature formula
{
	int m, j, i;
	double z1, z, xm, xl, pp, p3, p2, p1;

	m=(n + 1)/2;
	xm=0.5*(x2 + x1);
	xl=0.5*(x2 - x1);
	for (i=1;i<=m;i++)
	{
		z=cos(3.141592653589793238462643*(i-0.25)/(n+0.5));
		do 
		{
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++)
			{
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		}
		while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}
#undef EPS




/*method--method of integration; a, b-integration limits*/
/*n_steps--number of panels, func(t, x)-integrand at given x*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Simpson·¨ 

double Simpson_integ (int n_steps, double a, double b, double (*func)( double t, double x, double xp, int xpp), double x, double xp, int xpp)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
 
    dt=(b-a)/n_steps;
 
    for(n=0;n<n_steps;n++)
    { 
      tk=a+n*dt;
      f1=(*func)(tk,x, xp, xpp);
      f2=(*func)(tk+dt/2., x, xp, xpp);
      f3=(*func)(tk+dt, x, xp, xpp);
      sum+=dt*(f1+4*f2+f3)/6.; 
    }
  
   return sum;
}

double Simpson_integ1 (int n_steps, double a, double b, double (*func)( double t, double x), double x)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
  
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,x );
        f2=(*func)(tk+dt/2., x);
        f3=(*func)(tk+dt, x);
        sum+=dt*(f1+4*f2+f3)/6;
    }
  
   return sum;
}

double Simpson_integ2 (int n_steps, double a, double b, double (*func)( double t))
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk );
        f2=(*func)(tk+dt/2);
        f3=(*func)(tk+dt);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ3 (int n_steps, double a, double b, double (*func)( double t, int i, double x[], double Nmb[], double v0, double beta ), int i, double x[], double Nmb[], double v0, double beta)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,i,x,Nmb,v0,beta);
        f2=(*func)(tk+dt/2,i,x,Nmb,v0,beta);
        f3=(*func)(tk+dt,i,x,Nmb,v0,beta);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ4 (int n_steps, double a, double b, double (*func)( double t, int i, int num, double x[], double Nmb[], double beta ),int i, int num, double x[], double Nmb[], double beta)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,i,num,x,Nmb,beta);
        f2=(*func)(tk+dt/2,i,num,x,Nmb,beta);
        f3=(*func)(tk+dt,i,num,x,Nmb,beta);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ5 (int n_steps, double a, double b, double (*func)( double t, int i, int k, double x[], double beta ),int i, int k, double x[], double beta)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,i,k,x,beta);
        f2=(*func)(tk+dt/2,i,k,x,beta);
        f3=(*func)(tk+dt,i,k,x,beta);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ6 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, double beta), double x, double xp, double beta)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,x,xp,beta);
        f2=(*func)(tk+dt/2,x,xp,beta);
        f3=(*func)(tk+dt,x,xp,beta);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ7 (int n_steps, double a, double b, double (*func)( double t, int i, int nn, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R), int i, int nn, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,i,nn,x,w,Nmb,Beta,dBeta,injectionmode,R);
        f2=(*func)(tk+dt/2,i,nn,x,w,Nmb,Beta,dBeta,injectionmode,R);
        f3=(*func)(tk+dt,i,nn,x,w,Nmb,Beta,dBeta,injectionmode,R);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ8 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, double xpp, double xppp), double x, double xp, double xpp, double xppp)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,x,xp,xpp,xppp);
        f2=(*func)(tk+dt/2,x,xp,xpp,xppp);
        f3=(*func)(tk+dt,x,xp,xpp,xppp);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

double Simpson_integ9 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, int j, double xpp[], double xppp), double x, double xp, int j, double xpp[], double xppp)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
   
    dt=(b-a)/n_steps;
    for(n=0;n<n_steps;n++)
    {
        tk=a+n*dt;
        f1=(*func)(tk,x,xp,j,xpp,xppp);
        f2=(*func)(tk+dt/2,x,xp,j,xpp,xppp);
        f3=(*func)(tk+dt,x,xp,j,xpp,xppp);
        sum+=dt*(f1+4*f2+f3)/6;/*simpson f-la for Ik*/
    }
  
    return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Ö´ÐÐChebyshevÇó»ý·¨

double Chebyshev_integ (double a, double b, double (*func)( double t, double x, double xp, double Rp), double x, double xp, double Rp)   

  { 
	  int m,i,j;
      double h,d,p,ep,g,aa,bb,s,tt;
      
	  double eps= 0.000001;
      static double t[5]={-0.8324975, -0.3745414, 0.0, 0.3745414, 0.8324975};
      
      m=1;
      h=b-a; d=fabs(0.001*h);
      
      p=1.0e+35; ep=1.0+eps;
      
      while ((ep>=eps)&&(fabs(h)>d))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              
              for (j=0;j<=4;j++)
              { 
				  tt=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+(*func)(tt,x,xp,Rp);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      
      return g;
  }

double Chebyshev_integ1 (double a, double b, double (*func)( double t, double x), double x)   
  { 
	  int m,i,j;
      double h,d,p,ep,g,aa,bb,s,tt;
      double eps= 0.000001;
      static double t[5]={-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
      m=1;
      h=b-a; d=fabs(0.001*h);
      
      p=1.0e+35; ep=1.0+eps;
      
      while ((ep>=eps)&&(fabs(h)>d))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              
              for (j=0;j<=4;j++)
              { 
				  tt=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+(*func)(tt,x);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      
      return g;
  }
  
  double Chebyshev_integ2 (double a, double b, double (*func)(double t))   
  { 
	  int m,i,j;
      double h,d,p,ep,g,aa,bb,s,tt;
      double eps= 0.000001;
      static double t[5]={-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
      m=1;
      h=b-a; d=fabs(0.001*h);
      
      p=1.0e+35; ep=1.0+eps;
      
      while ((ep>=eps)&&(fabs(h)>d))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              
              for (j=0;j<=4;j++)
              { 
				  tt=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+(*func)(tt);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      
      return g;
  }

  double Chebyshev_integ3 (double a, double b, double (*func)(double t, double x, double xp, double beta, double R), double x, double xp, double beta, double R)
  { 
    int m,i,j;
      double h,d,p,ep,g,aa,bb,s,tt;
      double eps= 0.000001;
      static double t[5]={-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
      m=1;
      h=b-a; d=fabs(0.001*h);
      
      p=1.0e+35; ep=1.0+eps;
      
      while ((ep>=eps)&&(fabs(h)>d))
      { 
      g=0.0;
          for (i=1;i<=m;i++)
          { 
        aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              
              for (j=0;j<=4;j++)
              { 
          tt=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+(*func)(tt,x,xp,beta,R);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      
      return g;
  }

  double Chebyshev_integ4 (double a, double b, double (*func)( double t, double x, double xp, int ii, double xpp[], double R), double x, double xp, int ii, double xpp[], double R)
  { 
    int m,i,j;
      double h,d,p,ep,g,aa,bb,s,tt;
      double eps= 0.000001;
      static double t[5]={-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
      m=1;
      h=b-a; d=fabs(0.001*h);
      
      p=1.0e+35; ep=1.0+eps;
      
      while ((ep>=eps)&&(fabs(h)>d))
      { 
      g=0.0;
          for (i=1;i<=m;i++)
          { 
        aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              
              for (j=0;j<=4;j++)
              { 
          tt=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+(*func)(tt,x,xp,ii,xpp,R);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      
      return g;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  //Ö´ÐÐ±ä²½³¤ÌÝÐÎÇó»ý·¨
  double integration1 (double a, double b, double (*func)(double t, double xx),double xx)       
  { 
	  int n,k;
      double fa,fb,h,t1,p,s,x,t;
      double eps= 0.000001;
       
      fa=func (a,xx);  fb=func (b,xx);
      n=1; h=b-a;
      t1=h*(fa+fb)/2.0;
      p=eps+1.0;
      while (p>=eps)
      { 
		  s=0.0;
          for (k=0;k<=n-1;k++)
          { 
			  x=a+(k+0.5)*h;
              s=s+func (x,xx);
          }
          t=(t1+h*s)/2.0;
          p=fabs(t1-t);
          t1=t; n=n+n; h=h/2.0;
      }
      return t;
  }
  
  double integration2 (double a, double b, double (*func)(double t))       
  { 
	  int n,k;
      double fa,fb,h,t1,p,s,x,t;
      double eps= 0.000001;
       
      fa=func (a);  fb=func (b);
      n=1; h=b-a;
      t1=h*(fa+fb)/2.0;
      p=eps+1.0;
      while (p>=eps)
      { 
		  s=0.0;
          for (k=0;k<=n-1;k++)
          { 
			  x=a+(k+0.5)*h;
              s=s+func (x);
          }
          t=(t1+h*s)/2.0;
          p=fabs(t1-t);
          t1=t; n=n+n; h=h/2.0;
      }
      return t;
  }

  double integration3 (double a, double b, double (*func)(double t, double xx, double xxx),double xx, double xxx)
  { 
    int n,k;
      double fa,fb,h,t1,p,s,x,t;
      double eps= 0.000001;
       
      fa=func (a,xx,xxx);  fb=func (b,xx,xxx);
      n=1; h=b-a;
      t1=h*(fa+fb)/2.0;
      p=eps+1.0;
      while (p>=eps)
      { 
      s=0.0;
          for (k=0;k<=n-1;k++)
          { 
        x=a+(k+0.5)*h;
              s=s+func (x,xx,xxx);
          }
          t=(t1+h*s)/2.0;
          p=fabs(t1-t);
          t1=t; n=n+n; h=h/2.0;
      }
      return t;
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  //RombergÇó»ý·¨
  
  double Romberg1 (double a, double b, double (*func)(double t, double xx),double xx)       
  { 
	  int m,n,i,k;
      double y[10],h,ep,p,x,s,q;
      double eps= 0.00000001;
      
      h=b-a;
      y[0]=h*(func (a,xx)+func (b,xx))/2.0;
      
	  m=1; n=1; ep=eps+1.0;
      
	  while ((ep>=eps)&&(m<=9))
      { 
		  p=0.0;
          for (i=0;i<=n-1;i++)
          { 
			  x=a+(i+0.5)*h;
              p=p+func (x,xx);
          }
          p=(y[0]+h*p)/2.0;
          s=1.0;
          for (k=1;k<=m;k++)
          { 
			  s=4.0*s;
              q=(s*p-y[k-1])/(s-1.0);
              y[k-1]=p; p=q;
          }
          ep=fabs(q-y[m-1]);
          m=m+1; y[m-1]=q; n=n+n; h=h/2.0;
      }
      return q;
  }
  
  double Romberg2 (double a, double b, double (*func)(double t))       
  { 
	  int m,n,i,k;
      double y[10],h,ep,p,x,s,q;
      double eps= 0.00000001;
      
      h=b-a;
      y[0]=h*(func (a)+func (b))/2.0;
      m=1; n=1; ep=eps+1.0;
      while ((ep>=eps)&&(m<=9))
      { 
		  p=0.0;
          for (i=0;i<=n-1;i++)
          { 
			  x=a+(i+0.5)*h;
              p=p+func (x);
          }
          p=(y[0]+h*p)/2.0;
          s=1.0;
          for (k=1;k<=m;k++)
          { 
			  s=4.0*s;
              q=(s*p-y[k-1])/(s-1.0);
              y[k-1]=p; p=q;
          }
          ep=fabs(q-y[m-1]);
          m=m+1; y[m-1]=q; n=n+n; h=h/2.0;
      }
      return q;
  }
  


