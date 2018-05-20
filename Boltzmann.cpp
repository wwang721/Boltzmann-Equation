#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "nrutil.h"
#include "Integ.h" 
#include "Boltzmann.h"

using namespace std;

double cq=Pi*2;

double p=0.1;       	// tc/th 

double a=1.0/(1.0-p);	// tm/tc
double b=a-1.0;			// tm/th
  
double L=vm*tc;

double sgn(double x)
//阶跃函数
{
	if(x>0)
		return 1.0;
	else if(x==0)
		return 0.0;
	else
		return -1.0;
}

double nm0(double aq)
{
	return (1.0/(exp(30.1*aq*aq/(1.23*1.23*Temperature))-1));
}

double nm0p = Simpson_integ2(1e6,Low, Up, nm0);

double g_injection(double aq, double R, int injectionmode)
{
	double A=100.0;
	double Im=100.0;
	double aq0=Pi/2.0;
	double aq1=Pi/4.0;
	A = R*A;
	switch(injectionmode)
	{
		case 0: return A; break;
		case 1: return A+Im; break;
		case 2: return A+Im*exp(-(aq-aq0)*(aq-aq0)); break;
		case 3: return A+Im*exp(-(aq-aq0)*(aq-aq0))+Im*exp(-(aq-aq1)*(aq-aq1)); break;
		default: return 0; break;
	}
}

double S_integrand(double aq, double x, double R, int injectionmode)
//s被积函数
{

	return (g_injection(aq,R,injectionmode)*exp(-x/(L*sin(aq)))+(p*(nm0(aq)/nm0p)+(1.0-p))*(R-1.0)*exp(-x/(L*sin(aq))));
}

void S_vector(double x[], double S[], int n, double R, int injectionmode) 
//s
{
	cout<<"<S>"<<endl;
	char path[100];
	char fileName[100];
	
	sprintf(path, "./Gx");
    if(access(path,0)==-1)
        {
           printf("%s is set up\n",path);
           mkdir(path,0777);
        } 
    
	sprintf(path, "./Gx/Gx_mode(%d)", injectionmode);
    if(access(path,0)==-1)
        {
           printf("%s is set up\n",path);
           mkdir(path,0777);
        } 
        
	sprintf(fileName, "./Gx/Gx_mode(%d)/gx-mode(%d)R(%g).txt", injectionmode, injectionmode, R); 
	
	ifstream fin_gx(fileName);
	if(!(fin_gx))
	{	
		printf("the %s is set up.\n",fileName);
		
		ofstream fout_gx(fileName);

		for(int i=1,p=1; i<=n; i++)
		{
			S[i] = Simpson_integ(1e5,Low,Up,S_integrand,x[i],R,injectionmode);

			S[i] = S[i]/cq;

			if(i==p*(n*0.1))
			{
				cout<<"已完成:"<<p*10<<"%"<<endl;
				p++;
			}
			fout_gx<<S[i]<<endl;	
		}	
		fout_gx.close();
	}
	else
	{
		for(int i=1; i<=n; i++)
			fin_gx>>S[i];
		fin_gx.close();
	}
	cout<<"<S> done...\n"<<endl;
}

double fak_integrand_relax(double aq, double x, double xp, double R)
{
	if(x!=xp)
		return (1.0/(L*a*sin(aq)))*(R*exp(-(x+xp)/(L*sin(aq)))+exp(-fabs(x-xp)/(L*sin(aq))));
	else
		return (1.0/(L*a*sin(aq)))*(R*exp(-(x+xp)/(L*sin(aq))));
}

double fak_integrand_Beta(double aq, double x, double xp, double Beta, double R)
{
	double A,B;
	A = exp(-(x+xp)/(L*sin(aq)));
	B = exp(-fabs(x-xp)/(L*sin(aq)));

	if(x!=xp)
		return Beta*(R*A+B)/(L*sin(aq));
	else
		return Beta*R*A/(L*sin(aq));
}

double fak_integrand_dBeta(double aq, double x, double xp, int j, double dBeta[], double R)
{
	double A,B;
	A = exp(-(x+xp)/(L*sin(aq)));
	B = exp(-fabs(x-xp)/(L*sin(aq)));

	return (B*sgn(x-xp)-R*A)*dBeta[j];
}

double fak_integrand1(double aq, double Beta)
//fak被积函数2
{
	return (-2*Beta);
}

double fak_integrand2(double aq, double x, double Beta, double R)
//fak被积函数3
{
	return (1-R)*Beta*exp(-x/(L*sin(aq)));
}

void fak_relax(double x[], double w[], double **fak, int n, double R)
{   
    char fileName[100];
    sprintf(fileName, "./Fak_Relax/Fak-R(%g).txt", R);
	ofstream fout (fileName); 

	for (int i=1; i<=n; i++)
	{                       
		for (int j=1; j<=n; j++)
	   {
 	   		fak[i][j] = 0.0;
		 	fak[i][j] = Chebyshev_integ(Low, Up, fak_integrand_relax, x[i], x[j], R)*w[j];

			fak[i][j] = fak[i][j]/cq;

		  fout<<fak[i][j]<<"  ";
	   }
	   if(i%20==0) printf("%g%% of ak_F_Relax is ok!\n",100.0*i/n);
	}
	printf("\n");
    
    fout.close ();
}


void fak_Beta(double x[], double w[], double **fak, int n, double Beta, double R)
{   
    char fileName[100];
    sprintf(fileName, "./Fak_Beta/Fak_Beta_R(%g)/Fak-R(%g)Beta(%g).txt", R, R, Beta);
	ofstream fout (fileName); 

	for (int i=1; i<=n; i++)
	{                       
		for (int j=1; j<=n; j++)
	   {
 	   		fak[i][j] = 0.0;
		 	fak[i][j] = Chebyshev_integ3(Low, Up, fak_integrand_Beta, x[i], x[j], Beta, R)*w[j];
			
			if(i == j)
				fak[i][j] += Chebyshev_integ1(Low, Up, fak_integrand1, Beta);
			if(j == 1)
				fak[i][j] += Chebyshev_integ(Low, Up, fak_integrand2, x[i], Beta, R);

			fak[i][j] = fak[i][j]/cq;
			
			fout<<fak[i][j]<<"  ";
	   }
	   if(i%20==0) printf("%g%% of ak_F_Beta is ok!\n",100.0*i/n);
	}
	printf("\n");
    
    fout.close ();
}

void fak_dBeta(double x[], double w[], double **fak, int n, double dBeta[], double R)
{   
    char fileName[100];
    sprintf(fileName, "./Fak_dBeta/Fak_dBeta_R(%g)/Fak-R(%g)dBeta(%g).txt", R, R, dBeta[1]);
	ofstream fout (fileName); 

	for (int i=1; i<=n; i++)
	{                       
		for (int j=1; j<=n; j++)
	   {
 	   		fak[i][j] = 0.0;
		 	fak[i][j] = Chebyshev_integ4(Low, Up, fak_integrand_dBeta, x[i], x[j], j, dBeta, R)*w[j];

			fak[i][j] = fak[i][j]/cq;
			
			fout<<fak[i][j]<<"  ";
	   }
	   if(i%20==0) printf("%g%% of ak_F_dBeta is ok!\n",100.0*i/n);
	}
	printf("\n");
    
    fout.close ();
}

void fak_matrix(double x[], double w[], double **fak, int n, double Beta, double dBeta[], double R)
//fak
{   
	cout<<"<fak>"<<endl;

    char path[100];
    
    char fileName[100];
    char fileName_Relax[100];
    char fileName_Beta[100];
//  char fileName_dBeta[100];
    
    sprintf(path, "./Fak");
    if(access(path,0)==-1)
     {
         printf("new folder is set up\n");
         mkdir(path,0777);
      }  

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	sprintf(fileName_Relax, "./Fak_Relax/Fak-R(%g).txt", R);	
    ifstream fin_omk_Relax (fileName_Relax);
    if (!(fin_omk_Relax))
	{
		printf("The file %s is not created, being creating...\n",fileName_Relax); 
      
        sprintf(path, "./Fak_Relax"); 
        if(access(path,0)==-1)
           mkdir(path,0777);
        
        fak_relax(x, w, fak, n, R);
	}
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	sprintf(fileName_Beta, "./Fak_Beta/Fak_Beta_R(%g)/Fak-R(%g)Beta(%g).txt", R, R, Beta);	
    ifstream fin_omk_Beta (fileName_Beta);
    if (!(fin_omk_Beta))
	{
		printf("The file %s is not created, being creating...\n",fileName_Beta); 
      
        sprintf(path, "./Fak_Beta");
        if(access(path,0)==-1)
        	mkdir(path,0777); 
        
        sprintf(path, "./Fak_Beta/Fak_Beta_R(%g)", R); 
        if(access(path,0)==-1)
        	mkdir(path,0777);
        
        fak_Beta(x, w, fak, n, Beta, R);
	}
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
/*	sprintf(fileName_dBeta, "./Fak_dBeta/Fak_dBeta_R(%g)/Fak-R(%g)dBeta(%g).txt", R, R, dBeta[1]);
    ifstream fin_omk_dBeta (fileName_dBeta);
    if (!(fin_omk_dBeta))
	{
		printf("The file %s is not created, being creating...\n",fileName_dBeta); 
      
        sprintf(path, "./Fak_dBeta");
        if(access(path,0)==-1)
        	mkdir(path,0777);
    
        sprintf(path, "./Fak_dBeta/Fak_dBeta_R(%g)", R);
		if(access(path,0)==-1) 
       		mkdir(path,0777);
         
        fak_dBeta(x, w, fak, n, dBeta, R);
	}
*/
	fin_omk_Relax.close ();
	fin_omk_Beta.close ();
//	fin_omk_dBeta.close ();

	ifstream finp_omk_Relax (fileName_Relax);
	ifstream finp_omk_Beta (fileName_Beta);
//	ifstream finp_omk_dBeta (fileName_dBeta);

	double **fak_Relax;
    double **fak_Beta;
//  double **fak_dBeta;
   	fak_Relax = dmatrix(1,n,1,n);
    fak_Beta = dmatrix(1,n,1,n);
//  fak_dBeta = dmatrix(1,n,1,n);

    sprintf(fileName, "./Fak/Fak.txt");
	ofstream fout (fileName);
	for (int i=1; i<=n; i++)
	{                       
	   for (int j=1; j<=n; j++)
	   {
	   	  finp_omk_Relax >>fak_Relax[i][j];
	   	  finp_omk_Beta >>fak_Beta[i][j];
//	   	  finp_omk_dBeta >>fak_dBeta[i][j];
	   	  
	   	  fak[i][j]=0.0;
	   	  fak[i][j]=fak_Relax[i][j]+fak_Beta[i][j]/*+fak_dBeta[i][j]*/;

		  fout<<x[i]<<"  "<<x[j]<<"  "<<fak[i][j]<<endl;     
	   }
	}
	fout.close ();
	
	free_dmatrix(fak_Relax,1,n,1,n);
	free_dmatrix(fak_Beta,1,n,1,n);
//	free_dmatrix(fak_dBeta,1,n,1,n);

	finp_omk_Relax.close ();
	finp_omk_Beta.close ();
//	finp_omk_dBeta.close ();	
	
	
	cout<<"<fak> done...\n"<<endl;
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
		char nrerror_text[]="singular matrix in routine ludcmp";
		if (big == 0.0) nrerror(nrerror_text);
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

void Nmb_vector(double x[], double S[], double **fak, double Nmb[], int n, double Beta, double dBeta[], double R, int injectionmode,int k)
//Nmb(x)
{
	cout<<"<Nmb>"<<endl;	
		
	char path[100];
	char fileName[100];

	sprintf(path,"./Fx");
	if(access(path,0)==-1)
    {
        printf("Fx is set up\n");
    	mkdir(path,0777);
	}

	sprintf(path,"./Fx/Fx_injectionmode%d",injectionmode);
	if(access(path,0)==-1)
	{
        printf("Fx/Fx_injectionmode%d is set up\n",injectionmode);
    	mkdir(path,0777);
	}

	sprintf(path,"./Fx/Fx_injectionmode%d/Fx_R%g",injectionmode,R);
	if(access(path,0)==-1)
	{
        printf("Fx/Fx_injectionmode%d/Fx_R%g is set up\n",injectionmode,R);
    	mkdir(path,0777);
	}

	sprintf(path, "./Fx/Fx_injectionmode%d/Fx_R%g/Fx_Beta%g",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("Fx/Fx_injectionmode%d/Fx_R%g/Fx_Beta%g is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}
		
	sprintf(fileName, "./Fx/Fx_injectionmode%d/Fx_R%g/Fx_Beta%g/%d-Fx-mode(%d)R(%g)Beta(%g)dBeta(%g).txt", injectionmode,R,Beta,k,injectionmode,R,Beta,dBeta[1]);
	ifstream fin_fx(fileName);
	if(!fin_fx)
	{
		fred2(n,x,fak,S,Nmb);
		printf("%s is set up.\n",fileName);
		ofstream fout_fx(fileName);
		for(int i=1; i<=n; i++)
			fout_fx<<Nmb[i]<<endl;	
		fout_fx.close();
	}
	else
	{
		for(int i=1; i<=n; i++)
			fin_fx>>Nmb[i];	
	}
	fin_fx.close();
	
	cout<<"<Nmb> done...\n"<<endl;
}

double Nm(double aq, int i, int n, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R)
//Nm(x,q)
{
	double sum(0.0);
	if(aq>=0)
	{	
		sum = g_injection(aq,R,injectionmode)*exp(-x[i]/(L*sin(aq)));
		sum += (p*(nm0(aq)/nm0p)+(1.0-p))*(1-(1-R)*exp(-x[i]/(L*sin(aq)))-R*exp(-(distance+x[i])/(L*sin(aq))));
		sum += (1-R)*Beta*Nmb[1]*exp(-x[i]/(L*sin(aq)))-Beta*Nmb[i]+R*Beta*Nmb[n]*exp(-(distance+x[i])/(L*sin(aq)));

		double sump=0.0;
		for(int k=1; k<=n; k++)
		{
			sump += R*((1/(L*a*sin(aq)))+(Beta/(L*sin(aq))))*exp(-(x[k]+x[i])/(L*sin(aq)))*Nmb[k]*w[k];
		}
		sum += sump;

		for(int k=1;k<i;k++)
		{
			sum	+= (Beta/(L*sin(aq))+1/(L*a*sin(aq)))*Nmb[i]*exp(-(x[i]-x[k])/(L*sin(aq)))*w[k];
		}
	}
	else
	{
		sum = c2*exp((distance-x[i])/(L*sin(aq)))+(p*(nm0(aq)/nm0p)+(1.0-p))*(1-exp((distance-x[i])/(L*sin(aq))));
		sum += Beta*Nmb[n]*exp(-(x[i]-distance)/(L*sin(aq)))-Beta*Nmb[i];
		for(int k=i+1;k<=n;k++)	
			sum -= (Beta/(L*sin(aq))+1/(L*a*sin(aq)))*Nmb[i]*exp(-(x[i]-x[k])/(L*sin(aq)))*w[k];
	}
	return sum;
}

void Nm_3D(double x[], double w[], double Nmb[], int n, double Beta, double dBeta[], int injectionmode, double R, int k)
{
	cout<<"<Nm3D>"<<endl;

	char path[100];
    
    char fileName[100];
    
    sprintf(path, "./Nm3D");
    if(access(path,0)==-1)
    {
        printf("new folder is set up\n");
        mkdir(path,0777);
    } 

	sprintf(path,"./Nm3D/Nm3D_injectionmode%d",injectionmode);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d is set up\n",injectionmode);
    	mkdir(path,0777);
	}

	sprintf(path,"./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g",injectionmode,R);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g is set up\n",injectionmode,R);
    	mkdir(path,0777);
	}

	sprintf(path, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}

	sprintf(path, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Positive",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Positive is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}

	sprintf(path, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Negtive",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Negtive is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}

	sprintf(path, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Approx",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Approx is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}

	sprintf(fileName, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Positive/%d-Nm3D-mode(%d)R(%g)Beta(%g)dBeta(%g).txt",injectionmode,R,Beta,k,injectionmode,R,Beta,dBeta[1]);
	ifstream fin_Nm_Positive(fileName);
	if(!fin_Nm_Positive)
	{
		ofstream fout_Nm(fileName);
		for(int i=1;i<=n;i++)
			for(int j=1;j<=100;j++)
			{
				double aq=0.0+j*0.0314;
				fout_Nm<<x[i]<<"  "<<aq<<"  "<<Nm(aq,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R)<<endl;
			}
		fout_Nm.close();	
	}
	fin_Nm_Positive.close();

	sprintf(fileName, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Negtive/%d-Nm3D-mode(%d)R(%g)Beta(%g)dBeta(%g).txt",injectionmode,R,Beta,k,injectionmode,R,Beta,dBeta[1]);
	ifstream fin_Nm_Negtive(fileName);
	if(!fin_Nm_Negtive)
	{
		ofstream fout_Nm(fileName);
		for(int i=1;i<=n;i++)
			for(int j=1;j<=99;j++)
			{
				double aq=-3.14+j*0.0314;
				fout_Nm<<x[i]<<"  "<<aq<<"  "<<Nm(aq,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R)<<endl;
			}
		fout_Nm.close();	
	}
	fin_Nm_Negtive.close();

	sprintf(fileName, "./Nm3D/Nm3D_injectionmode%d/Nm3D_R%g/Nm3D_Beta%g/Approx/%d-Nm3D-mode(%d)R(%g)Beta(%g)dBeta(%g).txt",injectionmode,R,Beta,k,injectionmode,R,Beta,dBeta[1]);
	ifstream fin_Nm_Approx(fileName);
	if(!fin_Nm_Approx)
	{
		ofstream fout_Nm(fileName);
		for(int i=1;i<=n;i++)
			for(int j=0;j<=100;j++)
			{
				double aq=-0.05+j*0.001;
				if(aq != 0.0)
					fout_Nm<<x[i]<<"  "<<aq<<"  "<<Nm(aq,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R)<<endl;
			}
		fout_Nm.close();	
	}
	fin_Nm_Approx.close();

	cout<<"<Nm3D> done...\n"<<endl;
}

double aqNm(double aq, int i, int n, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R)
//aq*Nm
{
	return sin(aq)*Nm(aq,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R);
}

void jm_vector(double x[], double w[], double jm[], double Nmb[], int n, double Beta, double dBeta[], int injectionmode, double R, int k)
//jm
{
	cout<<"<jm>"<<endl;
	
	char path[100];
	char fileName[100];

	sprintf(path,"./jm");
	if(access(path,0)==-1)
	{
		cout<<"New folder is set up."<<endl;
		mkdir(path,0777);
	}

	sprintf(path,"./jm/jm_injectionmode%d",injectionmode);
	if(access(path,0)==-1)
	{
        printf("jm/jm_injectionmode%d is set up\n",injectionmode);
    	mkdir(path,0777);
	}

	sprintf(path,"./jm/jm_injectionmode%d/jm_R%g",injectionmode,R);
	if(access(path,0)==-1)
	{
        printf("jm/jm_injectionmode%d/jm_R%g is set up\n",injectionmode,R);
    	mkdir(path,0777);
	}

	sprintf(path, "./jm/jm_injectionmode%d/jm_R%g/jm_Beta%g",injectionmode,R,Beta);
	if(access(path,0)==-1)
	{
        printf("jm/jm_injectionmode%d/jm_R%g/jm_Beta%g is set up\n",injectionmode,R,Beta);
    	mkdir(path,0777);
	}
	
	sprintf(fileName, "./jm/jm_injectionmode%d/jm_R%g/jm_Beta%g/%d-jm-mode(%d)R(%g)Beta(%g)dBeta(%g).txt", injectionmode,R,Beta,k,injectionmode,R,Beta,dBeta[1]);
	ifstream fin_jm(fileName);
	if(!fin_jm)
	{
		ofstream fout_jm(fileName);
		for(int i=1,p=1; i<=n; i++)
		{
			jm[i] = Simpson_integ7(1e3,-1.0*Up,-1.0*Low,aqNm,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R);

			jm[i] = +Simpson_integ7(1e3,Low,Up,aqNm,i,n,x,w,Nmb,Beta,dBeta,injectionmode,R);		

			fout_jm<<jm[i]<<endl;
			if(i==p*(n/10))
			{
				cout<<"已完成:"<<p*10<<"%"<<endl;
				p++;
			}		
		}
		fout_jm.close();	
	}
	fin_jm.close();
	cout<<"<jm> done...\n"<<endl;
}
