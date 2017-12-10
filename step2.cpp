#include<iostream>
#include<time.h>
#include<math.h>
#include<stdio.h>
#include"data.h"
#include"nrutil.h"
using namespace std;
void fred2(int n, float a, float b, float t[], float f[], float w[], float (*g)(float), float (*ak)(float, float))
{
	void gauleg(float x1, float x2, float x[], float w[], int n);
	void lubksb(float **a, int n, int *indx, float b[]);
	void ludcmp(float **a, int n, int *indx, float *d);
	int i,j,*indx;
	float d,**omk;
	indx=ivector(1,n);
	omk=matrix(1,n,1,n);
	gauleg(a,b,t,w,n); 
	for (i=1;i<=n;i++) 
	{ 
		for (j=1;j<=n;j++) 
			omk[i][j]=(float)(i == j)-(*ak)(t[i],t[j])*w[j];
		f[i]=(*g)(t[i]);
	}
	ludcmp(omk,n,indx,&d); 
	lubksb(omk,n,indx,f);
	free_matrix(omk,1,n,1,n);
	free_ivector(indx,1,n);
}
int main()
{
	time_t start,end;
	start=time(NULL);	//简单计时
	
	int i,j;
	float v;
	static float s[N],x[N],fak[N][N];
	
	for(i=0;i<N;i++)
		x[i]=(i+1)*d/N;	//x

//s1 
	for(i=0;i<N;i++)
	{
		for(v=0.1;v<n;v+=0.1)
			s[i]+=(c1*exp(-x[i]/(v*tc))+c2*exp(-(d-x[i])/(v*tc)))*0.1;
		s[i]+=(tc/th)*Nm_0;	
		printf("%3d. %4.1fμm %14f\n",i+1,x[i]/Mu,s[i]);
	}	

//fak
	for(i=0;i<2;i++)
		for(j=0;j<N;j++)
		{
			for(v=0.1;v<n;v+=0.1)
				fak[i][j]+=(1/(v*tm))*exp(-fabs(x[i]-x[j])/(v*tc))*0.1*d/N;
			printf("已完成:%.2f%%\n",float(i*N+j+1)/(N*N)*100);
		}
/*
//输出数据
	FILE *fp;
	if((fp=fopen("data.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			fprintf(fp,"%f",fak[i][j]);
			if(j==N-1)
				fprintf(fp,"\n");
			else
				fprintf(fp,"\t");
		}
	}
	fclose(fp);
*/
	end=time(NULL);
	cout<<"\n用时:"<<end-start<<"s\n";	//计时结束并输出所用时间
	return 0;
}
