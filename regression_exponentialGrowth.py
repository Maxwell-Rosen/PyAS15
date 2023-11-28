import ias15FirstOrder1D as ias15
import rk4 as rk4
import euler as euler
import numpy as np
import matplotlib.pyplot as plt

def exponentialGrowth(y):
    return k*y

global k
k = 1

x_i = 1
dt = .3
tFinal = 5
nSteps = int(tFinal/dt)
t = np.linspace(0,tFinal,nSteps)
x = np.zeros(nSteps)
x_rk4 = np.zeros(nSteps)
x_euler = np.zeros(nSteps)
x[0] = x_i
x_rk4[0] = x_i
x_euler[0] = x_i
Bi = np.zeros(7)
for i in range(1,nSteps):
    x[i], Bf = ias15.push(exponentialGrowth, x[i-1], dt, B=Bi)
    x_rk4[i] = rk4.push(exponentialGrowth, x_rk4[i-1], dt)
    x_euler[i] = euler.push(exponentialGrowth, x_euler[i-1], dt)
    Bi = Bf
    print("Bi", Bi, "step", i)

print("Error", np.abs(x[-1]-np.exp(k*tFinal)))

#Make 2 side by side figures
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(t,x,'r.')
ax1.plot(t,np.exp(k*t),'b-')
ax1.plot(t,x_rk4,'gx')
ax1.plot(t,x_euler,'bo')

ax2.plot(t,np.abs(x-np.exp(k*t)),'r.')
ax2.plot(t,np.abs(x_rk4-np.exp(k*t)),'gx')
plt.show()

