#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>

#include "nrutil.h"
#include "Integ.h" 
#include "Boltzmann.h"

using namespace std;

static double tc=th*tm/(th+tm);

double s1_integrand(double v, double x)
//s1被积函数
{
	return (c1*exp(-x*Mu/(v*tc))+c2*exp(-(distance-x)*Mu/(v*tc)));
}

void s1_vector(double x[], double s1[], int n) 
//s1
{
	FILE *fp;
	if((fp=fopen("s1.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1; i<=n; i++)
	{
		s1[i] = Chebyshev_integ1(Low, Up, s1_integrand, x[i]);
		s1[i]+=(tc/th)*Nm_0;
		s1[i]=s1[i]/intvx;
		fprintf(fp,"%lf\n",s1[i]);	
	}	
	fclose(fp);	//积分并输出
}

double fak_integrand(double v, double x, double xp)
//fak被积函数
{
	return ((1/(v*tm))*exp(-fabs(x-xp)*Mu/(v*tc)));
}

void fak_matrix(double x[], double w[], double **fak, int n)
//fak
{
	cout<<"<fak>\n";
	FILE *fp;
	if((fp=fopen("fak.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	int p(1),q=n*n/100;
	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=n; j++)
		{
			fak[i][j] = Chebyshev_integ(Low, Up, fak_integrand, x[i], x[j])*w[j]*Mu/intvx;   
			fprintf(fp,"%lf",fak[i][j]);
			if(j==n)
				fprintf(fp,"\n");
			else
				fprintf(fp,"\t");	
			if( ((i-1)*n+j) == p*q)
			{
				cout<<"已完成:"<<p<<"%\n";
				p++;						
			}							//积分并输出
		}
	}
	fclose(fp);	//以矩阵形式输出
}

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

void fred2(int n, double t[], double **omk, double g[], double f[])
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

