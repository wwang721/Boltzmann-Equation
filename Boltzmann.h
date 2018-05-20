#define Pi 3.1415926535

//可改变的常量
#define Temperature 300.0 //Temperature

#define Up 3.1415926535 //积分上限
#define Low 1.0e-4	//积分下限

#define tc 1.0e-9

#define vm 5.5e9   //单位：μm/s

#define c2 0.0

#define distance 80.0 //单位：μm

double sgn(double x);	//阶跃函数

double nm0(double aq);
double g_injection(double aq, double R, int injectionmode);

double S_integrand(double aq, double x, double R, int injectionmode);
double fak_integrand_relax(double aq, double x, double xp, double R);
double fak_integrand_Beta(double aq, double x, double xp, double Beta, double R);
double fak_integrand_dBeta(double aq, double x, double xp, int j, double dBeta[], double R);
double fak_integrand1(double aq, double Beta);
double fak_integrand2(double aq, double x, double Beta, double R);


void S_vector(double x[], double S[], int n, double R, int injectionmode);  
void fak_relax(double x[], double w[], double **fak, int n, double R);
void fak_Beta(double x[], double w[], double **fak, int n, double Beta, double R);
void fak_dBeta(double x[], double w[], double **fak, int n, double dBeta[], double R);
void fak_matrix(double x[], double w[], double **fak, int n, double Beta, double dBeta[], double R);

void Nmb_vector(double x[], double S[], double **fak, double Nmb[], int n, double Beta, double dBeta[], double R, int injectionmode,int k);

double Nm(double aq, int i, int n, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R);
void Nm_3D(double x[], double w[], double Nmb[], int n, double Beta, double dBeta[], int injectionmode, double R, int k);

double aqNm(double aq, int i, int n, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R);
void jm_vector(double x[], double w[], double jm[], double Nmb[], int n, double Beta, double dBeta[], int injectionmode, double R, int k);
//jm
//以上是自己定义的

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void fred2(int n, double t[], double **omk, double g[], double f[]);
