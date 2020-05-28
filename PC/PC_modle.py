import numpy as np
import matplotlib.pyplot as plt
import time
import os
from scipy.fftpack import fft,fftshift,ifft


E1, E2, w, nc = 0.08, 0.046, 0.057, 15

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
    
def W(E):
    '''电离率 经验公式'''
    Ee = 0.0000001 + abs(E)             #  防止分母为零
    Ei, Eh = 15.576/27.2116, 0.58       ##Eh and ⁢𝐸i are the ionization potentials of hydrogen and the atom in question.
    W = 4 * (Ei / Eh) ** 2.5 * (1 / Ee) * np.exp(-2 / 3 * (Ei / Eh) ** 1.5 * (1 / Ee))
    return W

def Ne(E):
    '''等离子体浓度'''
    Ng = 2.4 * 10 ** 19 * (0.5292 * 10 ** (-8)) ** 3   #initial density of neutral particles
    Ne = np.zeros(len(time),float)
    Ne[0] = Ng * W(E[0])
    # Ne[1] = Ne[0] + (Ng-Ne[0])*W(E[1])
    # for i in range(len(time)-1):

    for i in range(len(time) - 1):
        Ne[i + 1] = Ne[i] + (Ng - Ne[i]) * W(E[i + 1]) * dt
    return Ne

def Je(E):
    dJ_dt = Ne(E) * E
    J = np.zeros(len(time),float)
    J[0] = Ne(E)[0] * E[0]
    for i in range(len(time)-1):
        J[i+1] = J[i] + dJ_dt[i] * dt
    return J

def Fourier(E):
    num = 5 * len(time)
    EE = np.pad(E, (num, num), 'constant')
    FFt = fft(EE)
    # return fftshift(FFt)
    return FFt
    
    
time = np.linspace(-836, 900, 1801)
nt = len(time)
dt = time[1] - time[0]

ES = Efs_field_sin2(time, 45) + Eds_field_sin2(time, 0)
EP = Efp_field_sin2(time, 45) - Edp_field_sin2(time, 0)

E = np.zeros(nt,dtype=complex)
E.real , E.imag = ES , EP

dj_dt_S = Ne(E) * ES
dj_dt_P = Ne(E) * EP

np.savetxt('PC_counter.dat',np.transpose((time,dj_dt_S,dj_dt_P)))

FFt_S = Fourier(dj_dt_S)
FFt_P = Fourier(dj_dt_P)
	
plt.plot(np.log10(abs(FFt_S)))
plt.plot(np.log10(abs(FFt_P)))
plt.xlim(0,8000)
plt.title("PC co-rotating")
plt.show()












