import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftshift,ifft


def E_FW_P(t):
    return (1/np.sqrt(1+ep_FW))*E1*np.exp(-2*np.log(2)*(t/sigma)**2)*np.cos(w*t)

def E_FW_S(t):
    return (1 / np.sqrt(1 + ep_FW)) * E1 * np.exp(-2 * np.log(2) * (t / sigma) ** 2) * np.sin(w * t)

def E_DW_P(t):
    return (1 / np.sqrt(1 + ep_DW)) * E2 * np.exp(-2 * np.log(2) * ((t-tau) / sigma) ** 2) * np.cos(2*w * (t-tau) + sita)

def E_DW_S(t):
    return (1 / np.sqrt(1 + ep_DW)) * E2 * np.exp(-2 * np.log(2) * ((t-tau) / sigma) ** 2) * np.sin(2*w * (t-tau) + sita)

def co_P(t):
    return E_DW_P(t) + E_FW_P(t)

def co_S(t):
    return E_DW_S(t) + E_FW_S(t)

def counter_P(t):
    return E_FW_P(t) + E_DW_P(t)

def counter_S(t):
    return E_FW_S(t) - E_DW_S(t)


def W(E):
    '''电离率 经验公式'''
    Ee = 0.0000001 + abs(E)             #  防止分母为零
    Ei, Eh = 15.6, 13.6       ##Eh and ⁢𝐸i are the ionization potentials of hydrogen and the atom in question.
    W = 4 * (Ei / Eh) ** 2.5 * (1 / Ee) * np.exp(-2 / 3 * (Ei / Eh) ** 1.5 * (1 / Ee))
    return W

def Ne(E):
    '''等离子体浓度'''
    Ng = 2.4 * 10 ** 19 * (0.5292 * 10 ** (-8)) ** 3   #initial density of neutral particles
    Ne = np.zeros(len(t),float)
    Ne[0] = Ng * W(E)[0]
    # for i in range(len(time)-1):

    for i in range(len(t) - 1):
        Ne[i + 1] = Ne[i] + (Ng - Ne[i]) * W(E[i + 1]) * (t[i+1] - t[i])
    return Ne

def Je(E):
    dJ_dt = Ne(E) * E
    J = np.zeros(len(t),float)
    J[0] = Ne(E)[0] * E[0]
    for i in range(len(t)-1):
        J[i+1] = J[i] + dJ_dt[i] * (t[i+1]-t[i])
    return J

def Fourier(E):
    num = 5 * len(t)
    EE = np.pad(E, (num, num), 'constant')
    FFt = fft(EE)
    # return fftshift(FFt)
    return FFt

def inFourier(E):
    return ifft(E)


'''参数设置'''
E1 = 0.0534      #基频光电场，对于光强10^14 W/cm2
k = 0.5           #倍频光与基频光光强比，
E2 = np.sqrt(k)*E1
w = 0.057             #800nm光频率
sigma = 35*41.34      #Gauss 半高宽
ep_FW = 1
ep_DW = 1               #epliticity of DW
sita = 0                #relative phase
tau = 0
def t():
    return np.linspace(-3000,3000,10000)
t = t()
dt = 1

co_rotating, counter_rotating = np.zeros(t.size, dtype=complex), np.zeros(t.size, dtype=complex)
co_rotating.real = co_P(t)
co_rotating.imag = co_S(t)
counter_rotating.real = counter_P(t)
counter_rotating.imag = counter_S(t)

dJ_dt_co_P = Ne(co_rotating) * co_P(t)
dJ_dt_co_S = Ne(co_rotating) * co_S(t)
dJ_dt_counter_P = Ne(counter_rotating) * counter_P(t)
dJ_dt_counter_S = Ne(counter_rotating) * counter_S(t)

FFt_dJ_dt_co_P = Fourier(dJ_dt_co_P)
FFt_dJ_dt_co_S = Fourier(dJ_dt_co_S)
FFt_dJ_dt_counter_P = Fourier(dJ_dt_counter_P)
FFt_dJ_dt_counter_S = Fourier(dJ_dt_counter_S)

plt.plot(FFt_dJ_dt_co_S)
plt.show()

