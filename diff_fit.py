import os
import numpy as np
import matplotlib.pyplot as plt 	#python绘图 	
from scipy import optimize		#最小二乘法拟合

def func1(x, p):
	C,lam = p
	return np.log10(C*np.exp(x/lam)/(lam*(1-np.exp(2*x/lam))))

def func2(x, p):
	'''
	数据拟合所用的函数：A*exp(-x/lambda)
	'''
	C,lam = p;
	return np.log10(C*np.exp(-x/lam))

def residuals(p, y, x, func):
	'''
	实验得到y和拟合函数之间的差，p为拟合需要找到的所有系数
	'''
	return y-func(x,p)



fitmode=2	#选择拟合函数

R=0
mode = 2
#------------------------------------------------------------------------
#读取T和Beta值
c5=1
fp=open("../beta_{}.txt".format(c5))
lines=fp.readlines()
Ts=[float(line.strip().split()[0]) for line in lines]
Betas=[float(line.strip().split()[1]) for line in lines]
fp.close()

#------------------------------------------------------------------------

Betas=[Beta+1.0 for Beta in Betas]

#Ts=[Ts[i] for i in range(12)]
#Betas=[Betas[i] for i in range(12)]
print('fx')
ds=[0.2,1,2,2.5,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]


fp=open("data1.txt","w")
yppp=[]
for n in range(len(Ts)):
	os.chdir("./Temperature{:g}".format(Ts[n]))
	yppp.append([])
	for d in ds:
		os.chdir("./Distance{:g}".format(d))

		#-----------------------------------------------------------------------------
		'''导入x数据'''

		f1 = open(r'./x.txt','r')
		lines = f1.readlines()
		x = [float(line.strip()) for line in lines]	#string经过strip或split处理后强制转化为float组成list
		f1.close()
		#-----------------------------------------------------------------------------

		N=range(0,999)
		xp = [(((x[Np+1]-x[Np])/2)+x[Np]) for Np in N]
		#=============================================================
		#手动求导
		
		Beta=Betas[n]
		f2 = open(r'./Fx/Fx_injectionmode{}/Fx_R{}/Fx_Beta{:g}/0-Fx-mode({})R({})Beta({:g})dBeta(0).txt'.format(mode,R,Beta,mode,R,Beta),'r')
		lines = f2.readlines()
		y = [float(line.strip()) for line in lines]
		#y = [Y_E/y[0] for Y_E in y]
		yp = [(y[Np+1]-y[Np])/(x[Np]-x[Np+1]) for Np in N]
		#yp = [Y_E/yp[0] for Y_E in yp]
		#yp = [ypp*1200/yp[0] for ypp in yp]
		f2.close()
		#==============================================================
		
		yppp[n].append(yp[960])	
		os.chdir("../")
	os.chdir("../")


for m in range(len(ds)):		
	for n in range(len(Ts)):
		fp.write("{:.2f}\t".format(yppp[n][m]))
	fp.write("\n")


		

fp.close()
#----------------------------------------------------------
#----------------------------------------------------------

#----------------------------------------------------------
#----------------------------------------------------------
#拟合求出lambda并作图
lams=[]
C=[]
fp=open("data1.txt","r")
lines=fp.readlines()

for num in range(len(Ts)):
	yy=[np.log10(float(line.strip().split()[num])) for line in lines]

	x=np.array(ds)
	y=np.array(yy)  

	if fitmode == 1:
		plsq = optimize.leastsq(residuals, [-400,10], args=(y,x,func1))
	elif fitmode == 2:
		plsq = optimize.leastsq(residuals, [400,10], args=(y,x,func2))

	c,lam = plsq[0]

	print("T:{}K:\t{}\t{}".format(Ts[num],c,lam))
	lams.append(lam)
	C.append(c)


plt.figure(figsize=(8,6), dpi=300)

plt.plot(Ts, lams, "o-", color="r", linewidth=2)

plt.xlabel("T")
plt.ylabel("$\lambda$")
plt.grid()
plt.savefig("fxlambda_fitmode{}.png".format(fitmode))

#----------------------------------------------------------------
#----------------------------------------------------------------
#读入数据和拟合曲线作图
x=np.arange(0.01,30,0.01)
N=range(0,2999)
plt.figure(figsize=(8,6), dpi=300)
y=np.arange(0.01,30,0.01)

fp=open("data1.txt","r")
lines=fp.readlines()

color="rygmckb"
point="sohdvp*D4321><^HD|_"

#Tindex=range(12)
Tindex=[0,1,3,6,8,9,10,11,13,16,21]

for n in range(len(Tindex)):
	num=Tindex[n]
	lam=lams[num]
	c=C[num]
	for Np in N:
		if fitmode == 1:
			y[Np]=(c*np.exp(x[Np]/lam)/(lam*(1-np.exp(2*x[Np]/lam))))
		elif fitmode == 2:
			y[Np]=c*np.exp(-x[Np]/lam)
	
	plt.plot(x, y, color=color[n%7], linewidth=1.5)
	
	yy=[float(line.strip().split()[num]) for line in lines]
	plt.semilogy(ds, yy, point[n], color=color[n%7], label="T={:g}T".format(Ts[num]), linewidth=1.5)


plt.xlim(0,30)
#plt.ylim(0,1200)
plt.xlabel("x/μm")
plt.ylabel("fx")
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0)
plt.grid()
plt.savefig("./fx_fitmode{}.png".format(fitmode), bbox_inches="tight")








