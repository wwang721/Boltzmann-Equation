//可改变的常量
#define Up 1e6	//积分上限
#define Low 0	//积分下限

#define th 1e-11
#define tm 1e-11

#define Nm_0 0

#define c1 2000
#define c2 0

#define distance 40 //单位：μm

#define Mu 1e-6 //微米数量级

#define intvx 2e9 //对Vx的积分

double s1_integrand(double v, double x);
double fak_integrand(double v, double x, double xp);

void s1_vector(double x[], double s1[], int n); 
void fak_matrix(double x[], double w[], double **fak, int n);

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void fred2(int n, double t[], double **omk, double g[], double f[]);
