#define N 1e6	//积分上限

#define th 1e-11
#define tm 1e-11

#define Nm_0 0

#define c1 -1
#define c2 0

#define distance 40 //单位：μm

#define Mu 1e-6 //微米数量级

void s1_Vector(double x[],double s1[],int n); 
void input_fak(double x[], double **fak, int n);

void gauleg(double x1, double x2, double x[], double w[], int n);  

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void fred2(int n, double t[], double w[], double **omk, double g[], double f[]);
