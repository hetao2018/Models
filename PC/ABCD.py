import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft


def E_function(t, tau):
    sigma_800, sigma_THz = 30 * 41.34, 600 * 41.34
    E1, E2 = 0.0534, 0.046
    w_800 = 0.057
    w_THz = 0.9 * 1.52e-4  # 1.52e-4  1THz
    E_800 = E1 * np.exp(-(t - tau) ** 2 / sigma_800 ** 2) * np.cos(w_800 * (t - tau))
    E_THz = E2 * np.exp(-((t) / sigma_THz) ** 2) * np.cos(w_THz * (t - 3500))
    E = E_800 + E_THz
    return E


def E_THz(t):
    sigma_800, sigma_THz = 30 * 41.34, 600 * 41.34
    E1, E2 = 0.0534, 0.046
    w_800 = 0.057
    w_THz = 0.9 * 1.52e-4  # 1.52e-4  1THz
    return E2 * np.exp(-((t) / sigma_THz) ** 2) * np.cos(w_THz * (t - 3500))


def W(E):
    Ei, Eh = 15.6, 13.6
    W = 4 * (Ei / Eh) ** 2.5 * (1 / abs(E)) * np.exp(-2 / 3 * (Ei / Eh) ** 1.5 * (1 / abs(E)))
    return W


def Ne(E):
    Ng = 2.4 * 10 ** 19 * (0.5292 * 10 ** (-8)) ** 3
    Ne = np.zeros(len(time))
    Ne[0] = Ng * W(E[0])
    # for i in range(len(time)-1):

    for i in range(len(time) - 1):
        Ne[i + 1] = Ne[i] + (Ng - Ne[i]) * W(E[i + 1]) * (time[1] - time[0])
    return Ne


def Fourier(dJ_dt):
    num = 10 * len(time)
    dJ_dt_N = np.pad(dJ_dt, (num, num), 'constant')
    Four_fft = fft(dJ_dt_N)
    return Four_fft


tl, tr = -1000, 1000
time = np.linspace(tl * 41.34, tr * 41.34, 5000)
plt.figure(figsize=(18, 8))
taus = np.linspace(-1000 * 41.34, 1000 * 41.34, 40)
peak = []
count = 0
for tau in taus:
    E = E_function(time, tau)
    dt = time[1] - time[0]
    dJ_dt = Ne(E) * E
    J = np.cumsum(dJ_dt * dt)

    Four_fft = Fourier(dJ_dt)

    f = 2 * np.pi / (time[1] - time[0])
    Hz = np.linspace(0, f, 21 * len(time))

    for i in range(21 * len(time) - 1380):
        Four_fft[i + 690] = complex(0, 0)

    # print(Four_fft[1000])

    Four_ifft = ifft(Four_fft)

    # plt.figure(figsize=(18,8))
    plt.subplot(1, 2, 1)
    plt.plot(time / 41.34, E,color='r',linewidth=2,linestyle='-',label='THz and Prob')
    plt.legend(loc=1,prop={'size':25,'family':'Time New Roman'},frameon=False)

    plt.subplot(1, 2, 2)
    plt.xlim(-0.2, 0.5)
    plt.ylim(0, 5e-7)
    Four = Fourier(dJ_dt)
    plt.plot(Hz, abs(Four))
    peak.append((max(abs(Four)[28000:36000])) ** 2)
    # plt.plot(time/41.34,W(E))

    # plt.subplot(2,2,3)
    # plt.xlim(-50,60)
    # plt.plot(np.linspace(tl,tr,21*len(time)),Four_ifft)
    # #plt.plot(time/41.34,Ne(E))
    #
    # plt.subplot(2,2,4)
    # plt.plot(time/41.34, J, 'r')
    # plt.xlim(-100, 1200)
    # plt.axhline(y=0, lw=1,  c='black', ls='--', alpha=0.3)

    # plt.savefig('D:picture')

    # plt.savefig('D:\\picture\\text%d.png'% count)
    count += 1
    plt.pause(0.01)
    plt.clf()
plt.ioff()

def bias():
    pass

plt.figure(figsize=(12, 7))
ax = plt.gca()
max_peak = max(peak)
plt.plot(taus / 41.34, peak / max_peak, linestyle='-', linewidth=2, label='SH', color='r')
plt.legend(loc=1, prop={'size': 25}, frameon=False)
# plt.legend(loc=)
# plt.legend()
E_THz = (abs(E_THz(time))) ** 2
E_THz_max = max(E_THz)
plt.plot(time / 41.34, E_THz / E_THz_max, linestyle='--', linewidth=2, label='THz', color='g')
plt.legend(loc=1, prop={'size': 25}, frameon=False)
plt.tick_params(labelsize=20)

font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 30,
         }
plt.xlabel('Time',font2)
plt.ylabel('Intensity', font2)

ax = plt.gca()
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
# plt.legend()
plt.show()
