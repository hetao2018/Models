import numpy as np
import ctypes as ct
import sfa_interface as sfa
import matplotlib.pyplot as plt
import os
from time import time as Time

time_start = Time()

nn = 15   #nc=10是勉强的，nc=8是不行的
E1, E2, w, nc = 0.08, 0.046, 0.057, nn

def fuhao(theta):
    if theta < 0:
        return -1
    else:
        return 1

#定义基频光电场：S方向
def Efs_field_sin2(t, theta=45):
    duration = 2. * nc * np.pi / w
    result = np.zeros(len(t))
    idx = np.nonzero((t >= -duration/2) & (t <= duration/2))
    result[idx] = (E1*np.sqrt(1/(1+np.tan(theta*np.pi/180)**2))*
          np.sin(w*(t[idx]-duration/2)/(2.*nc))**2*np.cos(w*t[idx]))
    return result

#定义基频光电场：P方向
def Efp_field_sin2(t,theta=45):
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


time = np.linspace(-1400., 1400., 2801)
nt = len(time)
dt = time[1] - time[0]

tau0 = np.linspace(-15., 15., 420)[140:280]

data_S, data_P = [], []
no = 0

for ta in tau0:
    ES = Efs_field_sin2(time) + Eds_field_sin2(time, ta)    
    EP = Efp_field_sin2(time) + Edp_field_sin2(time, ta)

    AS, AP, alphaS, alphaP = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
    
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


    # dipole_d = sfa.calc_wt (time, ES, AS, alphaS, EP, AP, alphaP)
    # os.rename('result.dat', 'Wt_{}.dat'.format(no))
    # no = no +1
    
    dipole_d = sfa.calc_dipole_CC_d2 (time, ES, AS, alphaS, EP, AP, alphaP)
    data_S.append(dipole_d[0])
    data_P.append(dipole_d[1])

np.savetxt('dipole_d_S2_term12.dat', np.transpose(data_S))
np.savetxt('dipole_d_P2_term12.dat', np.transpose(data_P))

from time import time as Time
end_time = Time()
operation  = (time_start -end_time)/60
print("operation time is {}".format(operation))















