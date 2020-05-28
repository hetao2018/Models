import numpy as np
import sfa_interface as sfa
import matplotlib.pyplot as plt
import time
import os

start_time = time.time()

nn = 15
E1, E2, w, nc = 0.08, 0.046, 0.057, nn

def fuhao(theta):
    if theta < 0:
        return -1
    else:
        return 1

#定义基频光电场：S方向
def Efs_field_sin2(t, theta):
    duration = 2. * nc * np.pi / w
    result = np.zeros(len(t))
    idx = np.nonzero((t >= -duration/2) & (t <= duration/2))
    result[idx] = (E1*np.sqrt(1/(1+np.tan(theta*np.pi/180)**2))*
          np.sin(w*(t[idx]-duration/2)/(2.*nc))**2*np.cos(w*t[idx]))
    return result

#定义基频光电场：P方向
def Efp_field_sin2(t,theta):
    duration = 2. * nc * np.pi / w
    result = np.zeros(len(t))
    idx = np.nonzero((t >= -duration/2) & (t <= duration/2))
    result[idx] = (E1*np.sqrt(np.tan(theta*np.pi/180)**2/(1+np.tan(theta*np.pi/180)**2))*
          np.sin(w*(t[idx]-duration/2)/(2.*nc))**2*fuhao(theta)*np.sin(w*t[idx]))
    return result

#定义倍频光电场矢势，倍频光是圆偏场
def Eds_field_sin2(t, DT):
    duration = 2. * nc * np.pi / w
    result = np.zeros(len(t))
    idx = np.nonzero((t >=- duration/2+DT*41.34) & (t <= duration/2+DT*41.34))
    result[idx] = (E2*np.cos(np.pi/4)* np.sin(w*(t[idx]-duration/2-DT*41.34)/(2.*nc))**2*
          np.cos(2*w*(t[idx]-DT*41.34)))
    return result

def Edp_field_sin2(t, DT):
    duration = 2. * nc * np.pi / w
    result = np.zeros(len(t))
    idx = np.nonzero((t >=- duration/2+DT*41.34) & (t <= duration/2+DT*41.34))
    result[idx] = (E2*np.sin(np.pi/4)* np.sin(w*(t[idx]-duration/2-DT*41.34)/(2.*nc))**2*
          np.sin(2*w*(t[idx]-DT*41.34)))
    return result

time = np.linspace(-900, 900, 1801)
nt = len(time)
dt = time[1] - time[0]

ES = Efs_field_sin2(time, 45) + Eds_field_sin2(time, 0)
EP = Efp_field_sin2(time, 45) - Edp_field_sin2(time, 0)

AS, AP =np.zeros(nt), np.zeros(nt)
alphaS, alphaP =  np.zeros(nt), np.zeros(nt)

AS[0] = ES[0]
for i in range(nt-1):
    AS[i+1] = AS[i] + 0.5 * (ES[i] + ES[i+1]) * dt
        
alphaS[0] = AS[0]
for i in range(nt-1):
    alphaS[i+1] = alphaS[i] + 0.5*(AS[i] + AS[i+1])*dt

AP[0] = EP[0]
for i in range(nt-1):
    AP[i+1] = AP[i] + 0.5 * (EP[i] + EP[i+1]) * dt
        
alphaP[0] = AP[0]
for i in range(nt-1):
    alphaP[i+1] = alphaP[i] + 0.5*(AP[i] + AP[i+1])*dt

dipole_d = sfa.calc_dipole_CC_d2 (time, ES, AS, alphaS, EP, AP, alphaP)
np.savetxt('counter_dipole_d.dat', np.transpose((time, dipole_d[0], dipole_d[1])))

#plt.plot(ES,EP)
#plt.plot(EP)

# dipole_d = sfa.calc_wt (time, ES, AS, alphaS, EP, AP, alphaP)
# os.rename('result.dat', 'Wt.dat')

#dipole_d = sfa.calc_Nt (time, ES, AS, alphaS, EP, AP, alphaP)

import time
end_time = time.time()
print('Operation time:{:.2f} min'.format((end_time-start_time)/60))
plt.show()













