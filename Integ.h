/*method--method of integration; a, b-integration limits*/
/*n_steps--number of panels, func(t, x)-integrand at given x*/

void gauleg(double x1, double x2, double x[], double w[], int n);  

//
double Simpson_integ (int n_steps, double a, double b, double (*func)( double t, double x, double xp, int xpp), double x, double xp, int xpp);
double Simpson_integ1 (int n_steps, double a, double b, double (*func)( double t, double x), double x);
double Simpson_integ2 (int n_steps, double a, double b, double (*func)( double t));
double Simpson_integ3 (int n_steps, double a, double b, double (*func)( double t, int i, double x[], double Nmb[], double v0, double beta ), int i, double x[], double Nmb[], double v0, double beta);
double Simpson_integ4 (int n_steps, double a, double b, double (*func)( double t, int i, int num, double x[], double Nmb[], double beta ),int i, int num, double x[], double Nmb[], double beta);
double Simpson_integ5 (int n_steps, double a, double b, double (*func)( double t, int i, int k, double x[], double beta ),int i, int k, double x[], double beta);
double Simpson_integ6 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, double beta), double x, double xp, double beta);
double Simpson_integ7 (int n_steps, double a, double b, double (*func)( double t, int i, int nn, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R),int i, int nn, double x[], double w[], double Nmb[], double Beta, double dBeta[], int injectionmode, double R);
double Simpson_integ8 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, double xpp, double xppp), double x, double xp, double xpp, double xppp);
double Simpson_integ9 (int n_steps, double a, double b, double (*func)( double t, double x, double xp, int j, double xpp[], double xppp), double x, double xp, int j, double xpp[], double xppp);
//
double Chebyshev_integ (double a, double b, double (*func)( double t, double x, double xp, double Rp), double x, double xp, double Rp);
double Chebyshev_integ1 (double a, double b, double (*func)( double t, double x), double x);
double Chebyshev_integ2 (double a, double b, double (*func)( double t));
double Chebyshev_integ3 (double a, double b, double (*func)( double t, double x, double xp, double beta, double R), double x, double xp, double beta, double R);
double Chebyshev_integ4 (double a, double b, double (*func)( double t, double x, double xp, int ii, double xpp[], double R), double x, double xp, int ii, double xpp[], double R);
//
double integration1 (double a, double b, double (*func)(double t, double x), double xx); 
double integration2 (double a, double b, double (*func)(double t));
double integration3 (double a, double b, double (*func)(double t, double xx, double xxx),double xx, double xxx);
//Ä¿Ç°ÕâÒ»»ý·ÖÊÇÕýÈ·µÄ£¡ 

double Romberg1 (double a, double b, double (*func)(double t, double xx),double xx);
double Romberg2 (double a, double b, double (*func)(double t));



  
