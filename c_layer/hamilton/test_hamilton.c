#include "hamilton.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Define a function pointer type for the derivative function
static const double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

typedef double (*derivative_func)(double, double);

double dxds(double x, double y){
    return x;
}

double dyds(double x, double y){
    return -y;
}

void test_parts(){
    double xi = 0.0;
    double yi = 1.0;
    double ds = 1;
    double Bx[7] = {0,0,0,0,0,0,0};
    double By[7] = {0,0,0,0,0,0,0};
    double xf, yf;

    int len_hF = 8;
    int len_BG = 7;
    double* xh = (double*)malloc(len_hF * sizeof(double));
    double* yh = (double*)malloc(len_hF * sizeof(double));
    double* Fx  = (double*)malloc(len_hF * sizeof(double));
    double* Fy  = (double*)malloc(len_hF * sizeof(double));
    double* Gx  = (double*)malloc(len_BG * sizeof(double));
    double* Gy  = (double*)malloc(len_BG * sizeof(double));
    for (int i = 0; i < 1; i++) {
        calculateNodePositions(xh, yh, Bx, By, xi, yi, dyds, dxds, ds, h, len_hF);
        for (int i = 0; i < 7; i++) {printf("xh[%i]=%1.4e \n",i,xh[i]);}
        for (int i = 0; i < 7; i++) {printf("yh[%i]=%1.4e \n",i,yh[i]);}
        calculateDerivatives(Fx, Fy, dyds, dxds, xh, yh, len_hF);
        for (int i = 0; i < 7; i++) {printf("Fx[%i]=%1.4e \n",i,Fx[i]);}
        for (int i = 0; i < 7; i++) {printf("Fy[%i]=%1.4e \n",i,Fy[i]);}

        calculateGFromF(Gx, Fx);
        calculateGFromF(Gy, Fy);
        for (int i = 0; i < 7; i++) {printf("Gx[%i]=%1.4e \n",i,Gx[i]);}
        for (int i = 0; i < 7; i++) {printf("Gy[%i]=%1.4e \n",i,Gy[i]);}
        calculateB(Bx, Gx);
        calculateB(By, Gy);
        for (int i = 0; i < 7; i++) {printf("Bx[%i]=%1.4e \n",i,Bx[i]);}
        for (int i = 0; i < 7; i++) {printf("By[%i]=%1.4e \n",i,By[i]);}
    }
    double hf = 1.0;
    calculateNodePositions(&xf, &yf, Bx, By, xi, yi, dyds, dxds, ds, &hf, 1);
    printf("xf=%1.4e \n",xf);
    printf("yf=%1.4e \n",yf);
    free(xh);
    free(yh);
    free(Fx);
    free(Fy);
    free(Gx);
    free(Gy);

    
}

void test_push() {
    double xi = 0.0;
    double yi = 1.0;
    double ds = 1;
    double Bx[7] = {0,0,0,0,0,0,0};
    double By[7] = {0,0,0,0,0,0,0};
    double xf, yf;

    push(&xf, &yf, Bx, By, xi, yi, ds, dxds, dyds);

        // Check that B is correct
    for (int i = 0; i < 7; i++) {printf("Bx[%i]=%1.4e \n",i,Bx[i]);}
    for (int i = 0; i < 7; i++) {printf("By[%i]=%1.4e \n",i,By[i]);}
}

int main(){
    test_parts();
    if (0) {
    test_push();
}
    return 0;
}