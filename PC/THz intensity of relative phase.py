import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft


E1, E2, w, sigma, nc = 0.08, 0.046, 0.057, 30*41.34, 15
tl, tr = -3000, 3000
time = np.linspace(tl, tr, 10001)
nt = len(time)

def E_800x(theta,t):
    return E1*np.cos(theta)*np.exp(-2*np.log(2)*(t/sigma)**2)*np.cos(w*t)
def E_800y(theta,t):
    return E1*np.sin(theta)*np.exp(-2*np.log(2)*(t/sigma)**2)*np.cos(w*t)
def E_400(t,tau):
    return E2*np.exp(-2*np.log(2)*((t-tau)/sigma)**2)*np.cos(2*w*(t-tau))


def W(E):
    Ei, Eh = 15.6, 13.6
    W = 4 * (Ei / Eh) ** 2.5 * (1 / abs(E)) * np.exp(-2 / 3 * (Ei / Eh) ** 1.5 * (1 / abs(E)))
    return W


def Ne(E,t):
    Ng = 2.4 * 10 ** 19 * (0.5292 * 10 ** (-8)) ** 3
    Ne = np.zeros(len(t))
    Ne[0] = Ng * W(E[0])
    # for i in range(len(time)-1):

    for i in range(len(t) - 1):
        Ne[i + 1] = Ne[i] + (Ng - Ne[i]) * W(E[i + 1]) * (t[1] - t[0])
    return Ne


def Fourier(dJ_dt,t):
    num = 10 * len(t)
    dJ_dt_N = np.pad(dJ_dt, (num, num), 'constant')
    Four_fft = fft(dJ_dt_N)
    return Four_fft



plt.figure(figsize=(11,6))

plt.ion()
space = np.linspace(0,2*np.pi,100)
for theta in space:
    E = E_800x(theta,time)+E_400(time,0)
    dJ_dt = Ne(E,time)*E
    plt.xlim(0,8)
    plt.ylim(0,7e-7)
    plt.plot(np.linspace(0,3000,210021),abs(Fourier(dJ_dt,time)))
    plt.pause(0.01)
    plt.clf()
plt.ioff()
# THz_new = abs(Fourier(dJ_dt,time))
# THz_new[560:209460] = 0




# plt.ion()
# space = np.linspace(0,1002*np.pi,100)
# for theta in space:
#     plt.ylim(-0.1,0.1)
#     plt.plot(time,E_800x(0,time))
#     plt.plot(time,E_400(time,theta))
#     plt.pause(0.01)
#     plt.clf()
# plt.ioff()
plt.show()


