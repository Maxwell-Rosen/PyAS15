#ifndef FIRSTORDER
#define FIRSTORDER

// Define a function pointer type for the derivative function
typedef double (*derivative_func)(double, double);

void push(double* xf, double* tf, double* B, double xi, double ti, double dt, derivative_func derivative);

void calculateNodePositions(double* xh, double* th, double* B, double xi, double ti,  derivative_func derivative, double dt, const double* hp, int len);

void calculateNewPosition(double* xn, double*tn, double* B, double xi, double ti, derivative_func derivative, double dt);

void calculateDerivatives(double* F, derivative_func derivative, double* x, double* t, int len);

void calculateB(double* B, double* G);

void calculateGFromF(double* G, double* F);

#endif