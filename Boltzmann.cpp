#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "string.h"

#include "nrutil.h"
#include "Integ.h" 
#include "Boltzmann.h"

using namespace std;

double s1_integrand(double v, double x, double v0)
//s1被积函数
{
	return (c1*exp(-(v-v0)*(v-v0)/(2*sigma*sigma))*exp(-x*Mu/(v*th))+c2*exp(-(distance-x)*Mu/(v*th)));
}

void s1_vector(double x[], double s1[], int n, double v0) 
//s1
{
	cout<<"<s1>"<<endl;
	FILE *fp;
	if((fp=fopen("s1.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1,p=1; i<=n; i++)
	{
		s1[i] = Simpson_integ(1e7,Low,Up,s1_integrand,x[i],v0);
		s1[i]+= Nm_0;
		s1[i] = s1[i]/intvx;
		if(i==p*(n/100))
		{
			cout<<"已完成:"<<p<<"%"<<endl;
			p++;
		}
		fprintf(fp,"%lf\n",s1[i]);	
	}	
	fclose(fp);	//积分并输出
	cout<<"<s1> done..."<<endl;
}

double fak_integrand1(double v, double x, double xp, double beta)
//fak被积函数1
{
	return (-(beta/(v*th))*exp(-fabs(x-xp)*Mu/(v*th)));
}

double fak_integrand2(double v, double beta)
//fak被积函数2
{
	return (2*beta);
}

double fak_integrand3(double v, double x, double beta)
//fak被积函数3
{
	return (-beta*exp(-x*Mu/(v*th)));
}

double fak_integrand4(double v, double x, double beta)
//fak被积函数4 
{
	return (-beta*exp(-(distance-x)*Mu/(v*th)));
}

void fak_matrix(double x[], double w[], double **fak, int n, double beta, int k)
//fak
{
	cout<<"<fak>\n";
	FILE *fp;
	char *name,*order;
	name=new char;
	order=new char;
	strcpy(name,"fak-");
	sprintf(order,"%d",k+1);
	name=strcat(name,order);
	name=strcat(name,".txt");
	if((fp=fopen(name,"w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	int p(1),q=n*n/100;
	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=n; j++)
		{
			fak[i][j] = Chebyshev_integ3(Low, Up, fak_integrand1, x[i], x[j], beta)*w[j]*Mu;   
			if(i==j)
				fak[i][j]+=Chebyshev_integ1(Low, Up, fak_integrand2, beta);
			if(j==1)
				fak[i][j]+=Chebyshev_integ(Low, Up, fak_integrand3, x[i], beta);
			if(j==n)
				fak[i][j]+=Chebyshev_integ(Low, Up, fak_integrand4, x[i], beta);
			fak[i][j]=fak[i][j]/intvx;

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
	cout<<"<fak> done..."<<endl;
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

double vxNm1(double v, int i, double x[], double Nmb[], double v0, double beta)
//Nm(x,vx>0)
{
	double sum;
	sum = c1*exp(-(v-v0)*(v-v0)/(2*sigma*sigma))*exp(-x[i]*Mu/(v*th));
	for(int k=1;k<i;k++)
		sum	+= beta*(Nmb[k+1]-Nmb[k])*exp(-(x[i]-x[k])*Mu/(v*th));
	return v*sum;
}

double vxNm2(double v, int i, int n, double x[], double Nmb[], double beta)
//Nm(x,vx<0)
{
	double sum;
	sum = c2*exp((distance-x[i])*Mu/(v*th));
	for(int k=i+1;k<n;k++)
	{
		sum -= beta*(Nmb[k+1]-Nmb[k])*exp(-(x[i]-x[k])*Mu/(v*th));
		
	}
	
	return v*sum;
}

void jm_vector(double x[], double jm[], double Nmb[], int n, double v0, double beta, int k)
//jm
{
	cout<<"<jm>"<<endl;
	FILE *fp;
	char *name,*order;
	name=new char;
	order=new char;
	strcpy(name,"jm-");
	sprintf(order,"%d",k+1);
	name=strcat(name,order);
	name=strcat(name,".txt");
	if((fp=fopen(name,"w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1,p=1; i<=n; i++)
	{
		jm[i] = Simpson_integ3(1e6,Low,Up,vxNm1,i,x,Nmb,v0,beta);		//vx>0积分
		jm[i] += Simpson_integ4(1e7,-1e6,-0.1,vxNm2,i,n,x,Nmb,beta);		//vx<0积分
		jm[i] = jm[i]/intvx;
		cout<<jm[i]<<endl;

		if(i==p*(n/100))
		{
			cout<<"已完成:"<<p<<"%"<<endl;
			p++;
		}	
		fprintf(fp,"%lf\n",jm[i]);	
	}	
	fclose(fp);			//积分并输出
	cout<<"<jm> done..."<<endl;
}
