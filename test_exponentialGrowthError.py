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
dt = np.array([5,1,.5,.1,.05,.01,.005,.001,])
tFinal = 10
error_ias15 = np.zeros(len(dt))
error_rk4 = np.zeros(len(dt))
error_euler = np.zeros(len(dt))

for idt in range(len(dt)):
    dtloop = dt[idt]
    nSteps = int(tFinal/dtloop)
    t = np.linspace(0,tFinal,nSteps)
    x_ias15 = np.zeros(nSteps)
    x_rk4   = np.zeros(nSteps)
    x_euler = np.zeros(nSteps)
    x_ias15[0] = x_i
    x_rk4[0]   = x_i
    x_euler[0] = x_i
    Bi = np.zeros(7)
    for i in range(1,nSteps):
        x_ias15[i], Bf = ias15.push(exponentialGrowth, x_ias15[i-1], dtloop, B=Bi)
        x_rk4[i]         = rk4.push(exponentialGrowth, x_rk4[i-1]  , dtloop)
        x_euler[i]     = euler.push(exponentialGrowth, x_euler[i-1], dtloop)
        Bi = Bf
    # plt.plot(t,x_ias15,'r.')
    # plt.plot(t,x_rk4,'gx')
    # plt.plot(t,np.exp(k*t),'b-')
    # plt.show()
    error_ias15[idt] = np.abs(x_ias15[-1]-np.exp(k*tFinal)) 
    error_rk4[idt]   = np.abs(  x_rk4[-1]-np.exp(k*tFinal))
    error_euler[idt] = np.abs(x_euler[-1]-np.exp(k*tFinal))

plt.loglog(dt,error_ias15,'r.')
plt.loglog(dt,error_rk4,'gx')
plt.loglog(dt,error_euler,'bo')
plt.xlabel('dt')
plt.ylabel('Error')
plt.legend(['IAS15','RK4','Euler'])
plt.title('Error vs dt over t ='+str(tFinal)+' k='+str(k)+' at end point')
plt.savefig('figures/test_exponentialGrowthError.png')
plt.show()

# print("Error", np.abs(x[-1]-np.exp(k*tFinal)))

# #Make 2 side by side figures
# fig = plt.figure()
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)

# ax1.plot(t,x,'r.')
# ax1.plot(t,np.exp(k*t),'b-')
# ax1.plot(t,x_rk4,'gx')

# ax2.plot(t,np.abs(x-np.exp(k*t)),'r.')
# ax2.plot(t,np.abs(x_rk4-np.exp(k*t)),'gx')
# plt.show()

