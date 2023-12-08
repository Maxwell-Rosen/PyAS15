import PyAS15.firstOrder2D as ias15
import numpy as np


def test_1():
    poly = lambda x,t: t+t**5
    xi = 0
    ti = 0
    dt = 1
    B = np.zeros(7)
    for i in range(1):
        xh, th = ias15.calculateNodePositions(B, xi, ti, poly, dt)
        print("xh = ", xh)
        print("th = ", th)
        F = ias15.calculateDerivatives(poly, xh, th)
        print("F = ", F)
        G = ias15.calculateGFromF(F)
        print("G = ", G)
        B = ias15.calculateB(G)
        print("B = ", B)
    B_correct = [1,0,0,0,1,0,0]
    np.testing.assert_allclose(B, B_correct, atol=1e-11)

    #Where is the 0th order term?
    #The 0th order term is accounted for in the pusher by hp[n]*dt*F1
    #Why is it limited to 1e-12?
    # Print(B) returns [ 1.00000000e+00  4.47402872e-14 -3.66173813e-13  1.31312788e-12 1.00000000e+00  1.94918152e-12 -6.31849588e-13]
    # Larger number of itterations does not reduce the size of the negligable terms.

def test_2():
    poly = lambda x,t: 1
    xi = 0
    ti = 0
    dt = 1
    B = np.zeros(7)
    for i in range(12):
        xh, th = ias15.calculateNodePositions(B, xi, ti, poly, dt)
        F = ias15.calculateDerivatives(poly, xh, th)
        G = ias15.calculateGFromF(F)
        B = ias15.calculateB(G)
    B_correct = [0,0,0,0,0,0,0]
    np.testing.assert_allclose(B, B_correct, atol=1e-11)

test_1()
# test_2()