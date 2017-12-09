#include<iostream>
#include<time.h>
#include<math.h>
#include"data.h"
using namespace std;
int main()
{
	time_t start,end;
	start=time(NULL);	//简单计时
	
	int i;
	double v,s[N],x[N];
	
	for(i=0;i<N;i++)
		x[i]=0.4*(i+1)*Mu;	//x
//s1 
	for(i=0;i<N;i++)
	{
		for(v=0;v<n;v+=0.1)
			s[i]+=(c1*exp(-x[i]/(v*tc))+c2*exp(-(d-x[i])/(v*tc)))*0.1;
		s[i]+=(tc/th)*Nm_0;	
		printf("%3d. %4.1lfμm %14lf\n",i+1,x[i]/Mu,s[i]);
	}		

	end=time(NULL);
	cout<<"\n用时:"<<end-start<<"s\n";	//计时结束并输出所用时间
	return 0;
}
