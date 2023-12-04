import ias15PDE as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: 2*x
dxds = lambda x,y: 2*y

end = 1
x1, y1, s1 = ias15.trace(dyds, dxds, .75, -1, 1, 100, end)
x2, y2, s2 = ias15.trace(dyds, dxds, .25, -1, 1, 100, end)

x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, .25, -1, end)
x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, .75, -1, end)

plt.plot(x1,y1,'.r')
plt.plot(x2,y2,'.r')
plt.plot(x3,y3,'.g')
plt.plot(x4,y4,'.g')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Hyperbolic Orbit Comparing Uniform Steps to Adaptive Steps")
plt.grid()
plt.legend(["Uniform .75","Uniform .25","Adaptive .25","Adaptive .75"])
plt.show()