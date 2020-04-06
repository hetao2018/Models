import numpy as np
import matplotlib.pyplot as plt
from photocurrent_model import low_harmonic

'''参数设置'''
E1 = 0.0534      #基频光电场，对于光强10^14 W/cm2
k = 0.5           #倍频光与基频光光强比，
E2 = np.sqrt(k)*E1
w = 0.057             #800nm光频率
sigma = 35*41.34      #Gauss 半高宽
ep_FW = 1
ep_DW = 1               #epliticity of DW
sita = 0                #relative phase
def t():
    return np.linspace(-3000,3000,10000)
t = t()
dt = 1


def E_FW_P(t):
    return (1/np.sqrt(1+ep_FW))*E1*np.exp(-2*np.log(2)*(t/sigma)**2)*np.cos(w*t)

def E_FW_S(t):
    return (1 / np.sqrt(1 + ep_FW)) * E1 * np.exp(-2 * np.log(2) * (t / sigma) ** 2) * np.sin(w * t)

def E_DW_P(t):
    return (1 / np.sqrt(1 + ep_DW)) * E2 * np.exp(-2 * np.log(2) * (t / sigma) ** 2) * np.cos(2*w * t + sita)

def E_DW_S(t):
    return (1 / np.sqrt(1 + ep_DW)) * E2 * np.exp(-2 * np.log(2) * (t / sigma) ** 2) * np.sin(2*w * t + sita)

def co_P(t):
    return E_DW_P(t) + E_FW_P(t)

def co_S(t):
    return E_DW_S(t) + E_FW_S(t)

def counter_P(t):
    return E_FW_P(t) + E_DW_P(t)

def counter_S(t):
    return E_FW_S(t) - E_DW_S(t)


'''。。。。。矢势A（初速度0）。。。。。。'''
def co_A_P(t):
    co_A_p = np.zeros(t.size)
    co_p = co_P(t)
    for i in range(t.size):
        for j in range(t.size-i):
            co_A_p[i] = co_A_p[i] + co_p[i+j]*dt
    return co_A_p

def co_A_S(t):
    co_A_s = np.zeros(t.size)
    co_s = co_S(t)
    for i in range(t.size):
        for j in range(t.size - i):
            co_A_s[i] = co_A_s[i] + co_s[i + j] * dt
    return co_A_s

def counter_A_P(t):
    counter_A_p = np.zeros(t.size)
    counter_p = counter_P(t)
    for i in range(t.size):
        for j in range(t.size - i):
            counter_A_p[i] = counter_A_p[i] + counter_p[i + j] * dt
    return counter_A_p

def counter_A_S(t):
    counter_A_s = np.zeros(t.size)
    counter_s = counter_S(t)
    for i in range(t.size):
        for j in range(t.size - i):
            counter_A_s[i] = counter_A_s[i] + counter_s[i + j] * dt
    return counter_A_s




