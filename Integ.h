/*method--method of integration; a, b-integration limits*/
/*n_steps--number of panels, func(t, x)-integrand at given x*/

void gauleg(double x1, double x2, double x[], double w[], int n);  

//
double Simpson_integ (int n_steps, double a, double b, double (*func)( double t, double x, double xp), double x, double xp);
double Simpson_integ1 (int n_steps, double a, double b, double (*func)( double t, double x), double x);
double Simpson_integ2 (int n_steps, double a, double b, double (*func)( double t));
double Simpson_integ3 (int n_steps, double a, double b, double (*func)( double t, int i, double x[], double Nmb[], double v0, double beta ), int i, double x[], double Nmb[], double v0, double beta);
double Simpson_integ4 (int n_steps, double a, double b, double (*func)( double t, int i, int num, double x[], double Nmb[], double beta ),int i, int num, double x[], double Nmb[], double beta);

//
double Chebyshev_integ (double a, double b, double (*func)( double t, double x, double xp), double x, double xp);
double Chebyshev_integ1 (double a, double b, double (*func)( double t, double x), double x);
double Chebyshev_integ2 (double a, double b, double (*func)(double t));
double Chebyshev_integ3 (double a, double b, double (*func)(double t, double x, double xp, double xpp), double x, double xp, double xpp);

//
double integration1 (double a, double b, double (*func)(double t, double x), double xx); 
double integration2 (double a, double b, double (*func)(double t));
double integration3 (double a, double b, double (*func)(double t, double xx, double xxx),double xx, double xxx);
//目前这一积分是正确的！ 

double Romberg1 (double a, double b, double (*func)(double t, double xx),double xx);
double Romberg2 (double a, double b, double (*func)(double t));



  




  
