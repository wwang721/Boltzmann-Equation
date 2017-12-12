#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
using namespace std;

#include "nrutil.h"
#include "Boltzmann.h"

static double tc=th*tm/(th+tm);

void s1_Vector(double x[], double s1[], int n) 
//s1
{
	for(int i=1; i<=n; i++)
	{
		s1[i]=0.0;
		for(double v=0.1; v<N; v+=0.1)
			s1[i]+=(c1*exp(-x[i]*Mu/(v*tc))+c2*exp(-(distance-x[i])*Mu/(v*tc)))*0.1;
		s1[i]+=(tc/th)*Nm_0;	
	}	//积分

	FILE *fp;
	if((fp=fopen("s1.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1;i<=n;i++)
		fprintf(fp,"%lf\n",s1[i]);
	fclose(fp);	//输出
}

void input_fak(double x[], double **fak, int n)
//fak
{
	cout<<"<fak>\n";
	int p(1),q=n*n/1000;
	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=n; j++)
		{
			fak[i][j]=0.0;
			for(double v=0.1;v<N;v+=0.1)
				fak[i][j]+=(1/(v*tm))*exp(-fabs(x[i]-x[j])*Mu/(v*tc))*0.1*distance*Mu/n;
			if(((i-1)*n+j)==p*q)
			{	
				printf("已完成：%4.1f%%\n",(float)p/10);
				p++;
			}	
		}
	}	//积分

	FILE *fp;
	if((fp=fopen("fak.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1;i<=n;i++)
	{
		for(int j=1; j<=n; j++)
		{
			fprintf(fp,"%lf",fak[i][j]);
			if(j==n)
				fprintf(fp,"\n");
			else
				fprintf(fp,"\t");
		}
	}
	fclose(fp);	//以矩阵形式输出
}

#define EPS 3.0e-15
void gauleg(double x1, double x2, double x[], double w[], int n)  
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

#define NRANSI
#define TINY 1.0e-20
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) 
	{
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) 
	{
		for (i=1;i<j;i++) 
		{
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) 
		{
			sum=a[i][j];
			for (k=1;k<j;k++)
			sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) 
			{
				big=dum;
				imax=i;
			}
		}
		if (j != imax) 
		{
			for (k=1;k<=n;k++) 
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) 
		{
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY 
#undef NRANSI

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) 
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
		for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) 
	{
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i]; 
	}
}

void fred2(int n, double t[], double w[], double **omk, double g[], double f[])
{
	int i,j,*indx;
	double d;

	indx=ivector(1,n);
	
	for (i=1;i<=n;i++) 
	{ 
		for (j=1;j<=n;j++) 
			omk[i][j]=(double)(i == j) - omk[i][j];
		f[i]=g[i];
	}
	
	ludcmp(omk,n,indx,&d); 
	lubksb(omk,n,indx,f);
	
	free_ivector(indx,1,n);
}

