import ias15PDE as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: x
dxds = lambda x,y: -y

xi = 2
yi = 0
end = np.pi

x, y, s = ias15.adaptiveTrace(dyds, dxds, xi, yi, end)

radius  = np.sqrt(x**2 + y**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
# ax1.set_xlim(-1.5,1.5)
# ax1.set_ylim(-1.5,1.5)
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
plt.savefig("figures/elliptical_orbit_adaptive.png")
plt.show()
