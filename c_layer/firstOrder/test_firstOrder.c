#include "firstOrder.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Define a function pointer type for the derivative function
static const double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

typedef double (*derivative_func)(double, double);

double t_plus_t5(double y, double t) {
    return t + pow(t, 5);
}

void test_individual() {
    double xi = 0.0;
    double ti = 0.0;
    double dt = 1;
    double B[7] = {0,0,0,0,0,0,0};
    double xf, tf;
    int len = 8;
    double* xh = (double*)malloc(len * sizeof(double));
    double* th = (double*)malloc(len * sizeof(double));
    double* F  = (double*)malloc(len * sizeof(double));
    double* G  = (double*)malloc(len * sizeof(double));

    for (int i = 0; i < 12; i++) {
        calculateNodePositions(xh, th, B, xi, ti, t_plus_t5, dt, h, len);
        calculateDerivatives(F, t_plus_t5, xh, th, len);
        calculateGFromF(G, F);
        calculateB(B, G);
    }

    // Check that B is correct
    double B_check[7] = {1,0,0,0,1,0,0};
    for (int i = 0; i < 7; i++) {printf("error[%i]=%1.4e \n",i,B_check[i]-B[i]);}
}

void test_push() {
    double xi = 0.0;
    double ti = 0.0;
    double dt = 1;
    double B[7] = {0,0,0,0,0,0,0};
    double xf, tf;

    push(&xf, &tf, B, xi, ti, dt, t_plus_t5);

    // Check that B is correct
    double B_check[7] = {1,0,0,0,1,0,0};
    for (int i = 0; i < 7; i++) {printf("error[%i]=%1.4e \n",i,B_check[i]-B[i]);}
}

int main(){
    test_individual();
    test_push();
    return 0;
}