import PyAS15.hamilton as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: 2*x
dxds = lambda x,y: 2*y

#Top Right
xi = .1
yi = .1
end = 1

x, y, s = ias15.adaptiveTrace(dyds, dxds, xi, yi, end)

error  = np.sqrt(x**2 - y**2) - np.sqrt(xi**2 - yi**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
ax1.plot(xi,yi,'g.')
ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
ax1.set_title("Hyperboloic Orbit")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
#Set aspect ratio to be square
ax1.set_aspect('equal')
ax1.grid()

ax2.plot(s,error,'r.')
ax2.set_title("x^2 - y^2 - initial")
ax2.set_xlabel("s")
ax2.set_ylabel("x^2 - y^2 - initial")
ax2.grid()
# plt.savefig("figures/hyperbolic_orbit_adaptive.png")
plt.show()

#Bottom Left
xi = -.1
yi = -.1
end = 1

x, y, s = ias15.adaptiveTrace(dyds, dxds, xi, yi, end)

error  = np.sqrt(x**2 - y**2) - np.sqrt(xi**2 - yi**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
ax1.plot(xi,yi,'g.')
ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
ax1.set_title("Hyperboloic Orbit")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
#Set aspect ratio to be square
ax1.set_aspect('equal')
ax1.grid()

ax2.plot(s,error,'r.')
ax2.set_title("x^2 - y^2 - initial")
ax2.set_xlabel("s")
ax2.set_ylabel("x^2 - y^2 - initial")
ax2.grid()
# plt.savefig("figures/hyperbolic_orbit_adaptive.png")
plt.show()

# Bottom right
xi = 1
yi = -1
end = 1

x, y, s = ias15.adaptiveTrace(dyds, dxds, xi, yi, end)

error  = np.sqrt(x**2 - y**2) - np.sqrt(xi**2 - yi**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
ax1.plot(xi,yi,'g.')
ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
ax1.set_title("Hyperboloic Orbit")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
#Set aspect ratio to be square
ax1.set_aspect('equal')
ax1.grid()

ax2.plot(s,error,'r.')
ax2.set_title("x^2 - y^2 - initial")
ax2.set_xlabel("s")
ax2.set_ylabel("x^2 - y^2 - initial")
ax2.grid()
# plt.savefig("figures/hyperbolic_orbit_adaptive.png")
plt.show()

# Top left
xi = -1
yi = 1
end = 1

x, y, s = ias15.adaptiveTrace(dyds, dxds, xi, yi, end)

error  = np.sqrt(x**2 - y**2) - np.sqrt(xi**2 - yi**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x,y,'r.')
ax1.plot(xi,yi,'g.')
ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
ax1.set_title("Hyperboloic Orbit")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
#Set aspect ratio to be square
ax1.set_aspect('equal')
ax1.grid()

ax2.plot(s,error,'r.')
ax2.set_title("x^2 - y^2 - initial")
ax2.set_xlabel("s")
ax2.set_ylabel("x^2 - y^2 - initial")
ax2.grid()
# plt.savefig("figures/hyperbolic_orbit_adaptive.png")
plt.show()