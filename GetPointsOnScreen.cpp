//This program is aimed to get the relative x and y coordinate of the points on the computer screen.
#include<iostream>
#include"windows.h"
#include"conio.h"
using namespace std;
int main()
{
	FILE *fp;
	int num(1),xmax,ymax,x0,y0,xm,ym;
	double xf,yf;
	POINT  pt;   
	fp=fopen("Point.txt","w");
	if(fp==NULL)
	{
		cout<<"can't open file\n";
		exit(1);
	}
	cout<<"1.请输入x轴最大值和y轴最大值"<<endl;
	cin>>xmax>>ymax;
	cout<<"2.确定原点位置，按空格键记录原点坐标;\n";
	getch();
	GetCursorPos(&pt); 
	x0=pt.x;
	y0=pt.y;
	cout<<"Finished...\n3.确定x轴最大值点位置，按空格键记录该点坐标;\n";
	getch();
	GetCursorPos(&pt); 
	xm=pt.x;
	cout<<"Finished...\n4.确定y轴最大值点位置，按空格键记录该点坐标;\n";
	getch();
	GetCursorPos(&pt); 
	ym=pt.y;
	cout<<"Finished...\n按任意键继续:"<<endl;
	getch();
	system("cls");
	cout<<"请按空格记录点坐标:"<<endl;
	while(1)
	{
		getch();
		GetCursorPos(&pt); 
		xf=(double)(pt.x-x0)*xmax/(xm-x0);
		yf=(double)(pt.y-y0)*ymax/(ym-y0);
		printf("%d.(%.2lf,%.2lf)\n",num,xf,yf);
		fprintf(fp,"%.2lf\t%.2lf\n",xf,yf);
		num++;
	}
	fclose(fp);
	system("pause");
	return 0;
}
