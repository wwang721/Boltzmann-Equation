//可改变的常量
#define Up 1e6	//积分上限
#define Low 0	//积分下限

#define th 1e-11
#define tm 1e-11

#define Nm_0 0

#define c1 1e9
#define c2 0

#define sigma 1 //正态分布方差

#define distance 40 //单位：μm

#define Mu 1e-6 //微米数量级

#define intvx 2e6 //对Vx的积分

double s1_integrand(double v, double x, double v0);
double fak_integrand1(double v, double x, double xp, double beta);
double fak_integrand2(double v, double beta);
double fak_integrand3(double v, double x, double beta);
double fak_integrand4(double v, double x, double beta);

void s1_vector(double x[], double s1[], int n, double v0); 
void fak_matrix(double x[], double w[], double **fak, int n, double beta, int k);
double vxNm1(double v, int i, double x[], double Nmb[], double v0, double beta);
double vxNm2(double v, int i, int n, double x[], double Nmb[], double beta);
void jm_vector(double x[], double jm[], double Nmb[], int n, double v0, double beta, int k);
//以上是自己定义的

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void fred2(int n, double t[], double **omk, double g[], double f[]);
