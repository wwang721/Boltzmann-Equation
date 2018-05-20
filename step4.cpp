#include <iostream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "nrutil.h"
#include "Integ.h"
#include "Boltzmann.h"

#define num 500	//取点数
extern double cq;

using namespace std;

int main()
{
	time_t start,end;
	start=time(NULL);					//简单计时
	cout<<"程序运行中..."<<endl;

	//定义及初始化
	double *x,*w,*S,*Nmb,*jm,*dBeta;
	double **fak;
	double Beta[33];
	double R;
	int injectionmode;

	x=dvector(1,num);
	w=dvector(1,num);
	
	dBeta=dvector(1,num);

	for(int i=0; i<=22; i++)
	{
		if(i<=12)
			Beta[i] = -1.2+0.1*i;
		else
			Beta[i] =0.5*(i-12);

		if(fabs(Beta[i]) <= 0.000001)
			Beta[i] = 0.0;
	}
	
	for(int i=1; i<=num; i++)
		dBeta[i]=0.0;

	//x取点并输出
	gauleg(0,distance,x,w,num);
	ifstream fin_x("x.txt");
	if(!fin_x)
	{
		ofstream fout_x("x.txt");
		for(int i=1; i<=num; i++)
			fout_x<<x[i]<<endl;
		fout_x.close();
	}
	fin_x.close();
	cout<<"<x> done...\n"<<endl;			

	//循环
	for(injectionmode=0; injectionmode<=3; injectionmode++)
	{
		for(int j=0; j<=5; j++)
		{
			R = 0.2*j;

			S=dvector(1,num);

			//S_vector(x,S,num,R,injectionmode);

			for(int k=0; k<=22; k+=2)
			{
				int counter=0;
					
				fak=dmatrix(1,num,1,num);
				Nmb=dvector(1,num);
				jm=dvector(1,num);
					
				//fak_matrix(x,w,fak,num,Beta[k],dBeta,R); 			//fak
			
				Nmb_vector(x,S,fak,Nmb,num,Beta[k],dBeta,R,injectionmode,counter);    //求出Nmb并输出
				
				Nm_3D(x,w,Nmb,num,Beta[k],dBeta,injectionmode,R,counter);    //Nm3D

				jm_vector(x,w,jm,Nmb,num,Beta[k],dBeta,injectionmode,R,counter);  //jm
			
				free_dmatrix(fak,1,num,1,num);
				free_dvector(Nmb,1,num);
				free_dvector(jm,1,num);
				
			}
			free_dvector(S,1,num);
		}
	}
	//释放变量空间
	free_dvector(x,1,num);
	free_dvector(w,1,num);
	
	
	free_dvector(dBeta,1,num);
    
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
