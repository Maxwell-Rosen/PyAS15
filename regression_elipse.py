import ias15PDE as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: x
dxds = lambda x,y: -y

x0 = 1
y0 = 0

ds = .001
s_final = np.pi*2/100
nSteps = int(s_final/ds)
s = np.linspace(0,s_final,nSteps+1, endpoint=True)
x = np.zeros(nSteps+1)
y = np.zeros(nSteps+1)
x[0] = x0
y[0] = y0
Bxi = np.zeros(7)
Byi = np.zeros(7)
for i in range(1,nSteps+1):
    x[i], y[i], Bxf, Byf = ias15.push(dyds, dxds, x[i-1], y[i-1], ds, Bxi, Byi)
    Bxi = Bxf
    Byi = Byf

radius  = np.sqrt(x**2 + y**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
ax1.set_title("Elliptical Orbit")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
#Set aspect ratio to be square
ax1.set_aspect('equal')
ax1.grid()

ax2.plot(s,radius-1,'r.')
ax2.set_title("Radius - 1")
ax2.set_xlabel("s")
ax2.set_ylabel("radius - 1")
ax2.grid()
plt.savefig("figures/elliptical_orbit.png")
plt.show()
