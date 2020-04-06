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

def low_harmonic():
    co_rotating,counter_rotating = np.zeros(t.size,dtype=complex),np.zeros(t.size,dtype=complex)
    co_rotating.real = co_P(t)
    co_rotating.imag = co_S(t)
    counter_rotating.real = counter_P(t)
    counter_rotating.imag = counter_S(t)

    dJ_dt_co_P = Ne(co_rotating) * co_P(t)
    dJ_dt_co_S = Ne(co_rotating) * co_S(t)
    dJ_dt_counter_P = Ne(counter_rotating) * counter_P(t)
    dJ_dt_counter_S = Ne(counter_rotating) * counter_S(t)

    '''
    # THz时域谱
    FFT_co_S = Fourier(dJ_dt_co_S)    
    for i in range(len(FFT_co_S)-70):
        FFT_co_S[i+35] = complex(0,0)
    plt.plot(ifft(FFT_co_S))
    plt.show()
    plt.xlim(4000,6000)
    plt.plot(abs(fftshift(Fourier(dJ_dt_co_S))))
    A,B = 0,30
    print(max(abs(Fourier(dJ_dt_co_S)[A:B])),max(abs(Fourier(dJ_dt_co_P)[A:B])),
          max(abs(Fourier(dJ_dt_counter_S)[A:B])),max(abs(Fourier(dJ_dt_counter_P)[A:B])))
    '''

    FFt_dJ_dt_co_P = Fourier(dJ_dt_co_P)
    FFt_dJ_dt_co_S = Fourier(dJ_dt_co_S)
    FFt_dJ_dt_counter_P = Fourier(dJ_dt_counter_P)
    FFt_dJ_dt_counter_S = Fourier(dJ_dt_counter_S)

    '''
    plt.figure(figsize=(24,14))
    plt.subplot(2,2,1)
    plt.plot(FFt_dJ_dt_co_P,'r',linestyle='-',label='co_P')
    plt.legend(prop={'size':12},loc='best')
    plt.subplot(2,2,2)
    plt.plot(FFt_dJ_dt_co_S,'r',linestyle='-',label='co_S')
    plt.legend(prop={'size':12},loc='best')
    plt.subplot(2,2,3)
    plt.plot(FFt_dJ_dt_counter_P,'b',linestyle='-',label='counter_P')
    plt.legend(prop={'size':12},loc='best')
    plt.subplot(2,2,4)
    plt.plot(FFt_dJ_dt_counter_P,'b',linestyle='-',label='counter_P')
    plt.legend(prop={'size':12},loc='best')
    plt.show()
    '''

    co_ , counter_ = {},{}
    co_['co_P_THz'] = max(abs(FFt_dJ_dt_co_P)[0:250])
    co_['co_P_three'] = max(abs(FFt_dJ_dt_co_P[1500:2100]))
    co_['co_P_four'] = max(abs(FFt_dJ_dt_co_P[2100:2700]))
    co_['co_P_five'] = max(abs(FFt_dJ_dt_co_P[2700:3300]))
    co_['co_P_six'] = max(abs(FFt_dJ_dt_co_P[3300:3900]))
    co_['co_S_THz'] = max(abs(FFt_dJ_dt_co_S)[0:250])
    co_['co_S_three'] = max(abs(FFt_dJ_dt_co_S[1500:2100]))
    co_['co_S_four'] = max(abs(FFt_dJ_dt_co_S[2100:2700]))
    co_['co_S_five'] = max(abs(FFt_dJ_dt_co_S[2700:3300]))
    co_['co_S_six'] = max(abs(FFt_dJ_dt_co_S[3300:3900]))
    counter_['counter_P_THz'] = max(abs(FFt_dJ_dt_counter_P)[0:250])
    counter_['counter_P_three'] = max(abs(FFt_dJ_dt_counter_P[2100:2700]))
    counter_['counter_P_four'] = max(abs(FFt_dJ_dt_counter_P[2700:3600]))
    counter_['counter_P_five'] = max(abs(FFt_dJ_dt_counter_P[3600:4500]))
    counter_['counter_P_six'] = max(abs(FFt_dJ_dt_counter_P[4500:5400]))
    counter_['counter_S_THz'] = max(abs(FFt_dJ_dt_counter_S)[0:250])
    counter_['counter_S_three'] = max(abs(FFt_dJ_dt_counter_S[2100:2700]))
    counter_['counter_S_four'] = max(abs(FFt_dJ_dt_counter_S[2700:3600]))
    counter_['counter_S_five'] = max(abs(FFt_dJ_dt_counter_S[3600:4500]))
    counter_['counter_S_six'] = max(abs(FFt_dJ_dt_counter_S[4500:5400]))
    co_['co_THz'] = np.sqrt(co_['co_P_THz']**2 + co_['co_S_THz']**2)
    co_['co_three'] = np.sqrt(co_['co_P_three']**2 + co_['co_S_three']**2)
    co_['co_four'] = np.sqrt(co_['co_P_four']**2 + co_['co_S_four']**2)
    co_['co_five'] = np.sqrt(co_['co_P_five']**2 + co_['co_S_five']**2)
    co_['co_six'] = np.sqrt(co_['co_P_six']**2 + co_['co_S_six']**2)
    counter_['counter_THz'] = np.sqrt(counter_['counter_P_THz']**2 + counter_['counter_S_THz']**2)
    counter_['counter_three'] = np.sqrt(counter_['counter_P_three']**2 + counter_['counter_S_three']**2)
    counter_['counter_four'] = np.sqrt(counter_['counter_P_four']**2 + counter_['counter_S_four']**2)
    counter_['counter_five'] = np.sqrt(counter_['counter_P_five']**2 + counter_['counter_S_five']**2)
    counter_['counter_six'] = np.sqrt(counter_['counter_P_six']**2 + counter_['counter_S_six']**2)
    return co_,counter_


'''参数设置'''
E1 = 0.0534      #基频光电场，对于光强10^14 W/cm2
# k = 0.5           #倍频光与基频光光强比，
# E2 = np.sqrt(k)*E1
w = 0.057             #800nm光频率
sigma = 35*41.34      #Gauss 半高宽
ep_FW = 1
ep_DW = 1               #epliticity of DW
sita = 0                #relative phase
def t():
    return np.linspace(-3000,3000,10000)
t = t()
dt = 1


THz_co_P,THz_co_S,THz_co,THz_counter_P,THz_counter_S,THz_counter = [],[],[],[],[],[]
three_co,three_co_P,three_co_S,three_counter,three_counter_P,three_counter_S   = [],[],[],[],[],[]
four_co,four_co_P,four_co_S,four_counter,four_counter_P,four_counter_S   = [],[],[],[],[],[]
five_co,five_co_P,five_co_S,five_counter,five_counter_P,five_counter_S   = [],[],[],[],[],[]
six_co,six_co_P,six_co_S,six_counter,six_counter_P,six_counter_S   = [],[],[],[],[],[]

N = 100
k_x = np.linspace(0.1,10,N)
# k_x = 0.5
# E2 = np.sqrt(k_x) * E1
# tau_ = np.linspace(0,1.333*41.34,N)
tau = 0
for k in k_x:
    E2 = np.sqrt(k) * E1
    # tau = k
    co_,counter_ = low_harmonic()
    THz_co.append(co_['co_THz'])
    THz_co_P.append(co_['co_P_THz'])
    THz_co_S.append(co_['co_S_THz'])
    THz_counter.append(counter_['counter_THz'])
    THz_counter_P.append(counter_['counter_P_THz'])
    THz_counter_S.append(counter_['counter_S_THz'])

    three_co.append(co_['co_three'])
    three_co_P.append(co_['co_P_three'])
    three_co_S.append(co_['co_S_three'])
    three_counter.append(counter_['counter_three'])
    three_counter_P.append(counter_['counter_P_three'])
    three_counter_S.append(counter_['counter_S_three'])

    four_co.append(co_['co_four'])
    four_co_P.append(co_['co_P_four'])
    four_co_S.append(co_['co_S_four'])
    four_counter.append(counter_['counter_four'])
    four_counter_P.append(counter_['counter_P_four'])
    four_counter_S.append(counter_['counter_S_four'])

    five_co.append(co_['co_five'])
    five_co_P.append(co_['co_P_five'])
    five_co_S.append(co_['co_S_five'])
    five_counter.append(counter_['counter_five'])
    five_counter_P.append(counter_['counter_P_five'])
    five_counter_S.append(counter_['counter_S_five'])

    six_co.append(co_['co_six'])
    six_co_P.append(co_['co_P_six'])
    six_co_S.append(co_['co_S_six'])
    six_counter.append(counter_['counter_six'])
    six_counter_P.append(counter_['counter_P_six'])
    six_counter_S.append(counter_['counter_S_six'])

'''
# fig,ax = plt.subplots()
# ax.set_yscale('log')
plt.figure(figsize=(12,7))
plt.plot(np.linspace(0,2*np.pi,N),THz_co_S,c='r',linestyle='-',label='THz_co_P',lw=2)
plt.plot(np.linspace(0,2*np.pi,N),three_co_S,c='b',linestyle='-',label='three_co_P',lw=2)
plt.plot(np.linspace(0,2*np.pi,N),four_co_S,c='g',linestyle='-',label='four_co_P',lw=2)
plt.plot(np.linspace(0,2*np.pi,N),five_co_S,c='y',linestyle='-',label='five_co_P',lw=2)
plt.plot(np.linspace(0,2*np.pi,N),six_co_S,c='k',linestyle='-',label='six_co_P',lw=2)
plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','270$\degree$','360$\degree$'))
plt.xlabel(r'Relative phase',size=15,**{'fontname':'serif'})
plt.ylabel(r'P-polarized co-rotating',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()
'''

# plt.figure(figsize=(12,7))
# plt.plot(k_x,THz_counter_S,c='r',linestyle='-',label='THz_counter_S',lw=2)
# plt.plot(k_x,three_counter_S,c='b',linestyle='-',label='three_counter_S',lw=2)
# plt.plot(k_x,four_counter_S,c='g',linestyle='-',label='four_counter_S',lw=2)
# plt.plot(k_x,five_counter_S,c='y',linestyle='-',label='five_counter_S',lw=2)
# plt.plot(k_x,six_counter_S,c='k',linestyle='-',label='six_counter_S',lw=2)
# # plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
# plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
# plt.ylabel(r'S_counter-rotating',size=15,**{'fontname':'serif'})
# plt.legend(loc='best',prop={'size':12})
# plt.show()

plt.figure(figsize=(12,7))
plt.plot(k_x,THz_counter,c='r',linestyle='-',label='THz_counter',lw=2)
plt.plot(k_x,THz_counter_P,c='b',linestyle='-',label='THz_counter_P',lw=2)
plt.plot(k_x,THz_counter_S,c='g',linestyle='-',label='THz_counter_S',lw=2)
plt.plot(k_x,THz_co,c='y',linestyle='-',label='THz_co',lw=2)
plt.plot(k_x,THz_co_P,c='k',linestyle='-',label='THz_co_P',lw=2)
plt.plot(k_x,THz_co_S,c='m',linestyle='-',label='THz_co_S',lw=2)
# plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
plt.ylabel(r'Intensity of THz',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()

plt.figure(figsize=(12,7))
plt.plot(k_x,three_counter,c='r',linestyle='-',label='three_counter',lw=2)
plt.plot(k_x,three_counter_P,c='b',linestyle='-',label='three_counter_P',lw=2)
plt.plot(k_x,three_counter_S,c='g',linestyle='-',label='three_counter_S',lw=2)
plt.plot(k_x,three_co,c='y',linestyle='-',label='three_co',lw=2)
plt.plot(k_x,three_co_P,c='k',linestyle='-',label='three_co_P',lw=2)
plt.plot(k_x,three_co_S,c='m',linestyle='-',label='three_co_S',lw=2)
# plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
plt.ylabel(r'Intensity of three',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()

plt.figure(figsize=(12,7))
plt.plot(k_x,four_counter,c='r',linestyle='-',label='four_counter',lw=2)
plt.plot(k_x,four_counter_P,c='b',linestyle='-',label='four_counter_P',lw=2)
plt.plot(k_x,four_counter_S,c='g',linestyle='-',label='four_counter_S',lw=2)
plt.plot(k_x,four_co,c='y',linestyle='-',label='four_co',lw=2)
plt.plot(k_x,four_co_P,c='k',linestyle='-',label='four_co_P',lw=2)
plt.plot(k_x,four_co_S,c='m',linestyle='-',label='four_co_S',lw=2)
# plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
plt.ylabel(r'Intensity of four',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()

plt.figure(figsize=(12,7))
plt.plot(k_x,five_counter,c='r',linestyle='-',label='five_counter',lw=2)
plt.plot(k_x,five_counter_P,c='b',linestyle='-',label='five_counter_P',lw=2)
plt.plot(k_x,five_counter_S,c='g',linestyle='-',label='five_counter_S',lw=2)
plt.plot(k_x,five_co,c='y',linestyle='-',label='five_co',lw=2)
plt.plot(k_x,five_co_P,c='k',linestyle='-',label='five_co_P',lw=2)
plt.plot(k_x,five_co_S,c='m',linestyle='-',label='five_co_S',lw=2)
# plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
plt.ylabel(r'Intensity of five',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()

plt.figure(figsize=(12,7))
plt.plot(k_x,six_counter,c='r',linestyle='-',label='six_counter',lw=2)
plt.plot(k_x,six_counter_P,c='b',linestyle='-',label='six_counter_P',lw=2)
plt.plot(k_x,six_counter_S,c='g',linestyle='-',label='six_counter_S',lw=2)
plt.plot(k_x,six_co,c='y',linestyle='-',label='six_co',lw=2)
plt.plot(k_x,six_co_P,c='k',linestyle='-',label='six_co_P',lw=2)
plt.plot(k_x,six_co_S,c='m',linestyle='-',label='six_co_S',lw=2)
# plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),('0$\degree$','90$\degree$','180$\degree$','360$\degree$''360$\degree$'))
plt.xlabel(r'Intensity ratio of $E_{2w}/E_w$',size=15,**{'fontname':'serif'})
plt.ylabel(r'Intensity of six',size=15,**{'fontname':'serif'})
plt.legend(loc='best',prop={'size':12})
plt.show()




'''
co_theta = np.angle(co_rotating)
co_abs = abs(co_rotating)
counter_theta = np.angle(counter_rotating)
counter_abs = abs(counter_rotating)
Size = 12
plt.figure(figsize=(10,5))
plt.subplot(121,polar=True)
plt.plot(co_theta,co_abs,linewidth=1.5,color="y",label="Co-rotating")
plt.rgrids(np.arange(0.01,1.1*max(co_abs),0.03),angle=180)
plt.legend(loc=(0.7,1),prop={'size': Size})
plt.subplot(122,polar=True)
plt.plot(counter_theta,counter_abs,linewidth=1.5,color="c",label="Counter-rotating")
plt.legend(loc=(0.7,1),prop={'size': Size})
plt.rgrids(np.arange(0.01,1.1*max(counter_abs),0.03),angle=180)
plt.tight_layout(pad=1.0, w_pad=4)
plt.show()
'''


