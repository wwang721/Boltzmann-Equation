#include <iostream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "nrutil.h"
#include "Integ.h"
#include "Boltzmann.h"

#define num 100	//取点数

using namespace std;

int main()
{
	time_t start,end;
	start=time(NULL);					//简单计时
	cout<<"程序运行中..."<<endl;

	//定义及初始化
	double *x,*w,*s1,*Nmb,*jm;
	double **fak,**Nm1,**Nm2;
	double beta[4]={0.1,0.004,0.007,0.01};
	double v0=5e5;
	
	x=dvector(1,num);
	w=dvector(1,num);
	s1=dvector(1,num);
	Nmb=dvector(1,num);
	jm=dvector(1,num);
	fak=dmatrix(1,num,1,num);
	Nm1=dmatrix(1,num,1,1e7);
	Nm2=dmatrix(1,num,1,1e7);

	//x取点并输出
	gauleg(0,distance,x,w,num);
	FILE *fp1;
	if((fp1 = fopen("x.txt","w")) == NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i = 1; i <= num; i++)
		fprintf(fp1,"%lf\n",x[i]);
	fclose(fp1);
	cout<<"<x> done..."<<endl;			

	s1_vector(x,s1,num,v0); 			//s1
	
	FILE *fp4;
	char *name,*order;
	name=new char;
	order=new char;
	
	for(int k=0;k<4;k++)
	{
		cout<<"正在进行第"<<k+1<<"重循环，共4重循环..."<<endl;

		fak_matrix(x,w,fak,num,beta[k],k); 			//fak
	
		//求出Nmb并输出
		fred2(num,x,fak,s1,Nmb);	
		cout<<"<Nmb> done..."<<endl;
		strcpy(name,"Nmb-");
		sprintf(order,"%d",k+1);
		name=strcat(name,order);
		name=strcat(name,".txt");
		if((fp4=fopen(name,"w"))==NULL)
		{
			cout<<"Can't open file.\n";			
			exit(1);
		}
		for(int i=1;i<=num;i++)
			fprintf(fp4,"%lf\n",Nmb[i]);
		fclose(fp4);					

		jm_vector(x,jm,Nmb,num,v0,beta[k],k);		//jm(x)
	}
	//释放变量空间
	free_dvector(x,1,num);
	free_dvector(w,1,num);
	free_dvector(s1,1,num);
	free_dvector(Nmb,1,num);
	free_dvector(jm,1,num);
    free_dmatrix(fak,1,num,1,num);
    free_dmatrix(Nm1,1,num,1,1e7);
    free_dmatrix(Nm2,1,num,1,1e7);

    //计时结束并输出所用时间
	end=time(NULL);
	cout<<"程序结束！"<<endl;
	int T,second,minute,hour;
	T=end-start;
	second=T%60;
	minute=T/60;
	hour=minute/60;
	minute=minute%60; 
	cout<<"用时:"<<hour<<" h "<<minute<<" min "<<second<<" s.\n";	
	
	return 0;
}
