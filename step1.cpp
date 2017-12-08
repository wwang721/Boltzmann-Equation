#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
using namespace std;
void creatfile(const char *name,const char *s,double *x,double *n)
{
//定义生成数据文件的函数
	FILE *fp;
	if((fp=fopen(name,s))==NULL)
	{
		cout<<"Can't open file1\n";
		exit(1);
	}
	for(int i=0;i<99;i++)
		fprintf(fp,"%.1lf\t%lf\n",x[i]*1e6,n[i]);
	fclose(fp);
}

int main()
{
	time_t start,end;
	start=time(NULL);	//简单计时

 	int i,j,imax(0),imin(0);
	double a,b(0),v,c1,c2,nm[99],x[99],jm[99],k;
	const double d=4e-5,t=1e-11;
	a=1e6;

	for(v=0;v<1e6;v+=0.1)
		b+=exp(-d/(v*t))*0.1;		//经过验算步长取1和取0.01时b值相差0.01，取1以节省时间
	c1=a/(a*a-b*b);					//实际上c1还要乘系数n0
	c2=b/(b*b-a*a);					//实际上还要乘系数n0
	
	for(i=0;i<99;i++)
		x[i]=0.4*(i+1)*1e-6;			//x 

//计算nm(χ)的表达式
	cout<<"nm(x):\n";
	for(i=0;i<99;i++)
	{
		for(v=0;v<1e6;v+=0.1)
			nm[i]+=(c1*exp(-x[i]/(v*t))+c2*exp((x[i]-d)/v*t))*0.1;//nm(x),实际上还要乘系数n0
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]*1e6,nm[i]);
	}

//计算jm(χ)的表达式	
	cout<<"\njm(x):\n";
	for(i=0;i<99;i++)
	{
		for(v=0;v<1e6;v+=0.1)
			jm[i]+=v*(c1*exp(-x[i]/(v*t))+c2*exp((x[i]-d)/v*t))*0.1;//jm(x),实际上还要乘系数n0
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]*1e6,jm[i]);
	}

//根据论文图像中数据计算得到的jm'(χ)
	double lambda,C,jm_1[99];
	lambda=9.4*1e-6;
	C=-4.2*1e-4;
	cout<<"\njm'(x):\n";	
	for(i=0;i<99;i++)
	{
		jm_1[i]=C*exp(x[i]/lambda)/(lambda*(1-exp(2*x[i]/lambda))); //jm'(x),实际上还要乘系数n0/K
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]*1e6,jm_1[i]);
	}

//计算比例系数K
	cout<<"\n比例系数K:\n";
	double K[99];
	for(i=0;i<99;i++)
	{		
		K[i]=jm_1[i]/jm[i];			
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]*1e6,K[i]);	
	}										//分别计算各个点的比例系数
	for(i=0;i<99;i++)
	{
		if(K[i]>K[imax])
			imax=i;	        
		if(K[i]<K[imin])
			imin=i;
	}	       
	int p=(int)round((K[imax]-K[imin])*1e5);	
	double pn[p];
	for(i=0;i<p;i++)
		pn[i]=0;
	for(i=0;i<99;i++)
		for(j=0;j<p;j++)
			if(K[i]>=(K[imin]+j*1e-5) && K[i]<(K[imin]+(j+1)*1e-5))
				pn[j]++;
	for(i=1;i<p;i++)
		if(pn[imax]<pn[i])	
			imax=i;
	k=(double)round((K[imin]+imax*1e-5)*1e5)*1e-5;
	cout<<"最符合的比例系数K:"<<k<<endl;  //取最匹配的比例系数k，精度为1e-5

//生成数据文件
 	creatfile("nmp.txt","w",x,nm);
 	creatfile("jmp.txt","w",x,jm);
	cout<<"\n比例系数调整过的jm'(x):\n";
	for(i=0;i<99;i++)
	{
		jm_1[i]=jm_1[i]/k;		
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]*1e6,jm_1[i]);
	}
	cout<<"最符合的比例系数K:"<<k<<endl;  //取最匹配的比例系数k，精度为1e-5
	creatfile("jm_1p.txt","w",x,jm_1); //输出已经用比例系数调整过的jm'(x)

	end=time(NULL);
	cout<<"\n用时:"<<end-start<<"s\n";		//计时结束并输出

	system("pause");	
	return 0;
}
