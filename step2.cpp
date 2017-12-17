#include<iostream>
#include<time.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
using namespace std;

#include"nrutil.h"
#include"Boltzmann.h"

#define num 100	//取点数

int main()
{
	time_t start,end;
	start=time(NULL);	//简单计时
	cout<<"程序运行中..."<<endl;

	//定义及初始化
	double *x,*w,*s1,*Nmb;
	double **fak;
	
	x=dvector(1,num);
	w=dvector(1,num);
	s1=dvector(1,num);
	Nmb=dvector(1,num);
	fak=dmatrix(1,num,1,num);

	gauleg(0, distance, x, w, num);	
	FILE *fp1;
	if((fp1 = fopen("x.txt","w")) == NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i = 1; i <= num; i++)
		fprintf(fp1,"%lf\n",x[i]);
	fclose(fp1);
	cout<<"<x> done..."<<endl;	//x取点并输出

	s1_Vector(x,s1,num); //s1
	cout<<"<s1> done..."<<endl;

	input_fak(x,w,fak,num); //fak
	cout<<"<fak> done..."<<endl;

	fred2(num,x,w,fak,s1,Nmb);	
	cout<<"<Nmb> done..."<<endl;
	FILE *fp4;
	if((fp4=fopen("Nmb.txt","w"))==NULL)
	{
		cout<<"Can't open file.\n";
		exit(1);
	}
	for(int i=1;i<=num;i++)
		fprintf(fp4,"%lf\n",Nmb[i]);
	fclose(fp4);	//求出Nmb并输出

	//释放变量空间
	free_dvector(x,1,num);
	free_dvector(w,1,num);
	free_dvector(s1,1,num);
	free_dvector(Nmb,1,num);
    free_dmatrix(fak,1,num,1,num);

	end=time(NULL);
	cout<<"程序结束！"<<endl;
	int T,second,minute,hour;
	T=end-start;
	second=T%60;
	minute=T/60;
	hour=minute/60;
	minute=minute%60;
	cout<<"用时:"<<hour<<"h "<<minute<<"min "<<second<<"s.\n";	//计时结束并输出所用时间
	
	return 0;
}
