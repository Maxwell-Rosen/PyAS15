import ias15PDE as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: x
dxds = lambda x,y: -y

x0 = 1
y0 = 0

dt = 1
Bx = np.zeros(7)
By = np.zeros(7)
for i in range(12):
    xh, yh = ias15.calculateNodePositions(Bx, By, x0, y0, dxds, dyds, dt)
    Fx, Fy = ias15.calculateDerivatives(dxds, dyds, xh, yh)
    Gx, Gy = ias15.calculateGFromF_xy(Fx, Fy)
    Bx, By = ias15.calculateB_xy(Gx, Gy)
    