#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <cmath>
#include <math.h>

#include "Integ.h"

/*
gauleg主要是给出区间的划分和权重
 
*/

#define EPS 3.0e-15
void gauleg(double x1, double x2, double x[], double w[], int n)  
//w[]应该怎么给定？x[]和w[]都是外部所给出的空矢量
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
//Simpson法 

double Simpson_integ (int n_steps, double a, double b, double (*func)( double t, double x, double xp), double x, double xp)
{
    int n;
    double sum=0;
    double tk=0., dt;
    double f1, f2, f3;
 
    dt=(b-a)/n_steps;
 
    for(n=0;n<n_steps;n++)
    { 
      tk=a+n*dt;
      f1=(*func)(tk,x, xp );
      f2=(*func)(tk+dt/2., x, xp);
      f3=(*func)(tk+dt, x, xp);
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


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//执行Chebyshev求积法

double Chebyshev_integ (double a, double b, double (*func)( double t, double x, double xp), double x, double xp)   

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
                  s=s+(*func)(tt,x,xp);
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
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  //执行变步长梯形求积法
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


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  //Romberg求积法
  
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
  


