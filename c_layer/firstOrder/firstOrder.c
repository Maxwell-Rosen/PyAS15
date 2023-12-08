// Author: Maxwell Rosen 12/7/2023
// Description: This file contains the C implementation of the first order
//              Gauss Radau integrator. So far I'm just translating it
//              directly from the python implementation. It will be a standalone
//              C file.
#include "firstOrder.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Gauss Radau spacings
static const double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};
// Other constants
static const double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};

// Define a function pointer type for the derivative function

void push(double* xf, double* tf, double* B, double xi, double ti, double dt, derivative_func derivative) {
    int len_hF = 8;
    int len_BG = 7;
    double* xh = (double*)malloc(len_hF * sizeof(double));
    double* th = (double*)malloc(len_hF * sizeof(double));
    double* F  = (double*)malloc(len_hF * sizeof(double));
    double* G  = (double*)malloc(len_BG * sizeof(double));

    for (int i = 0; i < 12; i++) {
        calculateNodePositions(xh, th, B, xi, ti, derivative, dt, h, len_hF);
        calculateDerivatives(F, derivative, xh, th, len_hF);
        calculateGFromF(G, F);
        calculateB(B, G);
    }

    calculateNewPosition(xf, tf, B, xi, ti, derivative, dt);

    free(xh);
    free(th);
    free(F);
    free(G);
}


void calculateNodePositions(double* xh, double* th, double* B, double xi, double ti,
 derivative_func derivative, double dt, const double* hp, int len) {
    double F1 = derivative(xi, ti);
    for (int n = 0; n < len; n++) {
        xh[n] = xi + hp[n]*dt*(F1   + hp[n]*1/2*\
                             ( B[0] + hp[n]*2/3*\
                             ( B[1] + hp[n]*3/4*\
                             ( B[2] + hp[n]*4/5*\
                             ( B[3] + hp[n]*5/6*\
                             ( B[4] + hp[n]*6/7*\
                             ( B[5] + hp[n]*7/8*\
                               B[6])))))));
        th[n] = ti + hp[n]*dt;
    }
}

void calculateNewPosition(double* xn, double*tn, double* B, double xi, double ti, derivative_func derivative, double dt) {
    double hp = 1.0;
    calculateNodePositions(xn, tn, B, xi, ti, derivative, dt, &hp, 1);
}

void calculateDerivatives(double* F, derivative_func derivative, double* x, double* t, int len) {
    for (int i = 0; i < len; i++) {
        F[i] = derivative(x[i], t[i]);
    }
}

void calculateB(double* B, double* G) {
    B[0] = G[0] + c[0]*G[1] + c[1]*G[2] + c[3]*G[3] + c[6]*G[4] + c[10]*G[5] + c[15]*G[6];
    B[1] =             G[1] + c[2]*G[2] + c[4]*G[3] + c[7]*G[4] + c[11]*G[5] + c[16]*G[6];
    B[2] =                         G[2] + c[5]*G[3] + c[8]*G[4] + c[12]*G[5] + c[17]*G[6];
    B[3] =                                     G[3] + c[9]*G[4] + c[13]*G[5] + c[18]*G[6];
    B[4] =                                                 G[4] + c[14]*G[5] + c[19]*G[6];
    B[5] =                                                              G[5] + c[20]*G[6];
    B[6] =                                                                           G[6];
}

void calculateGFromF(double* G, double* F) {
    G[0] =       (F[1] - F[0])/rr[0];
    G[1] =      ((F[2] - F[0])/rr[1]  - G[0])/rr[2];
    G[2] =     (((F[3] - F[0])/rr[3]  - G[0])/rr[4]  - G[1])/rr[5];
    G[3] =    ((((F[4] - F[0])/rr[6]  - G[0])/rr[7]  - G[1])/rr[8]  - G[2])/rr[9];
    G[4] =   (((((F[5] - F[0])/rr[10] - G[0])/rr[11] - G[1])/rr[12] - G[2])/rr[13] - G[3])/rr[14];
    G[5] =  ((((((F[6] - F[0])/rr[15] - G[0])/rr[16] - G[1])/rr[17] - G[2])/rr[18] - G[3])/rr[19] - G[4])/rr[20];
    G[6] = (((((((F[7] - F[0])/rr[21] - G[0])/rr[22] - G[1])/rr[23] - G[2])/rr[24] - G[3])/rr[25] - G[4])/rr[26] - G[5])/rr[27];
}