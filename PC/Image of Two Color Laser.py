import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from Two_color_laser import *



def im_FW_DW():
    plt.figure(figsize=(12, 8))
    plt.subplot(4,2,1)
    plt.plot(t,E_FW_P(t),"r",label='FW_P',linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 2)
    plt.plot(t, E_FW_S(t), "r", label='FW_S', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 3)
    plt.plot(t, E_DW_P(t), "b", label='DW_P', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 4)
    plt.plot(t, E_DW_S(t), "b", label='DW_S', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 5)
    plt.plot(t, co_P(t), "m", label='co_P', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 6)
    plt.plot(t, co_S(t), "m", label='co_S', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 7)
    plt.plot(t, counter_P(t), "g", label='counter_P', linestyle='-')
    plt.legend(prop={'size': 12})
    plt.subplot(4, 2, 8)
    plt.plot(t, counter_S(t), "g", label='counter_S', linestyle='-')
    plt.legend(prop={'size': 12})
    matplotlib.rc('xtick', labelsize=12)
    matplotlib.rc('ytick', labelsize=12)
    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
    plt.show()
    return None

def im_Co_Counter():
    plt.figure(figsize=(12, 7))
    plt.subplot(121)
    plt.plot(co_P(t), co_S(t), "g", linestyle='-')
    plt.title("Co-rotating", size=12)
    plt.subplot(122)
    plt.plot(counter_P(t), counter_S(t), "purple", linestyle='-')
    plt.title("Counter-rotating", size=12)
    plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
    plt.show()
    return None

def im_3D_co():
    fig = plt.figure(figsize=(12, 7))
    ax = fig.gca(projection='3d')
    ax.plot(t, co_P(t), co_S(t))
    plt.show()
    return None

def im_3D_counter():
    fig = plt.figure(figsize=(12, 7))
    ax = fig.gca(projection='3d')
    ax.plot(t, counter_P(t), counter_S(t))
    plt.show()
    return None

def polar_coo():
    Co_dcf = np.zeros(t.size, dtype=complex)
    Counter_dcf = np.zeros(t.size, dtype=complex)
    Co_dcf.real = co_P(t)
    Co_dcf.imag = co_S(t)
    Counter_dcf.real = counter_P(t)
    Counter_dcf.imag = counter_S(t)
    co_theta = np.angle(Co_dcf)
    co_abs = abs(Co_dcf)
    counter_theta = np.angle(Counter_dcf)
    counter_abs = abs(Counter_dcf)
    plt.figure(figsize=(10, 5))
    plt.subplot(121, polar=True)
    plt.plot(co_theta, co_abs, linewidth=1.5, color="y", label="Co-rotating")
    plt.rgrids(np.arange(0.01, 1.1 * max(co_abs), 0.03), angle=180)
    plt.legend(loc=(0.7, 1), prop={'size': 12})
    plt.subplot(122, polar=True)
    plt.plot(counter_theta, counter_abs, linewidth=1.5, color="c", label="Counter-rotating")
    plt.legend(loc=(0.7, 1), prop={'size': 12})
    plt.rgrids(np.arange(0.01, 1.1 * max(counter_abs), 0.03), angle=180)
    plt.tight_layout(pad=1.0, w_pad=4)
    plt.show()
    return None

def polar_coo_A():
    co_A = np.zeros(t.size, dtype=complex)
    counter_A = np.zeros(t.size, dtype=complex)
    co_A.real = -co_A_P(t)
    co_A.imag = -co_A_S(t)
    counter_A.real = -counter_A_P(t)
    counter_A.imag = -counter_A_S(t)
    co_A_theta = np.angle(co_A)
    co_A_abs = abs(co_A)
    counter_A_theta = np.angle(counter_A)
    counter_A_abs = abs(counter_A)
    plt.figure(figsize=(10, 5))
    plt.subplot(121, polar=True)
    plt.plot(co_A_theta, co_A_abs, linewidth=1.5, color="y", label="Co-rotating")
    plt.rgrids(np.arange(0.1, 1.1 * max(co_A_abs), 0.3), angle=180)
    plt.legend(loc=(0.7, 1), prop={'size': 12})
    plt.subplot(122, polar=True)
    plt.plot(counter_A_theta, counter_A_abs, linewidth=1.5, color="c", label="Counter-rotating")
    plt.rgrids(np.arange(0.1, 1.1 * max(counter_A_abs), 0.3), angle=180)
    plt.legend(loc=(0.7, 1), prop={'size': 12})
    plt.tight_layout(pad=1.0, w_pad=4)
    plt.show()
    return None

im_FW_DW()
im_Co_Counter()
im_3D_co()
im_3D_counter()
# polar_coo_A()
plt.plot(t,counter_P(t)+counter_S(t))
plt.show()


