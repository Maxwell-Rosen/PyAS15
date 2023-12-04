import ias15PDE as ias15
import numpy as np
import matplotlib.pyplot as plt

dyds = lambda x,y: 2*x
dxds = lambda x,y: 2*y

#Top Right
end = 1


x1, y1, s1 = ias15.adaptiveTrace(dyds, dxds, .1, .1, end)
x2, y2, s2 = ias15.adaptiveTrace(dyds, dxds, -.1, -.1, end)
x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, 1, -1, end)
x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, -1, 1, end)

plt.plot(x1,y1,'.r')
plt.plot(x2,y2,'.r')
plt.plot(x3,y3,'.r')
plt.plot(x4,y4,'.r')

# x1, y1, s1 = ias15.adaptiveTrace(dyds, dxds, -.75, -1, end)
# x2, y2, s2 = ias15.adaptiveTrace(dyds, dxds, -.25, -1, end)
x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, .25, -1, end)
x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, .75, -1, end)

# plt.plot(x1,y1,'.g')
# plt.plot(x2,y2,'.g')
plt.plot(x3,y3,'.g')
plt.plot(x4,y4,'.g')

x1, y1, s1 = ias15.adaptiveTrace(dyds, dxds, -.75, 1, end)
x2, y2, s2 = ias15.adaptiveTrace(dyds, dxds, -.25, 1, end)
# x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, .25, 1, end)
# x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, .75, 1, end)

plt.plot(x1,y1,'.b')
plt.plot(x2,y2,'.b')
# plt.plot(x3,y3,'.b')
# plt.plot(x4,y4,'.b')

x1, y1, s1 = ias15.adaptiveTrace(dyds, dxds, 1,-.75,end)
x2, y2, s2 = ias15.adaptiveTrace(dyds, dxds, 1,-.25, end)
# x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, 1, .25, end)
# x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, 1, .75, end)

plt.plot(x1,y1,'.m')
plt.plot(x2,y2,'.m')
# plt.plot(x3,y3,'.m')
# plt.plot(x4,y4,'.m')

# x1, y1, s1 = ias15.adaptiveTrace(dyds, dxds, -1,-.75,end)
# x2, y2, s2 = ias15.adaptiveTrace(dyds, dxds, -1,-.25, end)
x3, y3, s3 = ias15.adaptiveTrace(dyds, dxds, -1, .25, end)
x4, y4, s4 = ias15.adaptiveTrace(dyds, dxds, -1, .75, end)

# plt.plot(x1,y1,'.k')
# plt.plot(x2,y2,'.k')
plt.plot(x3,y3,'.k')
plt.plot(x4,y4,'.k')

plt.show()