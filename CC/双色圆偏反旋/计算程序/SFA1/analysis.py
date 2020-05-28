import matplotlib.pyplot as plt
import numpy as np

def Second_derivative(Dx,Dy):
    Dx_1 = np.zeros(nt)
    Dy_1 = np.zeros(nt)
    Dx_1[0] = Dx[0]
    Dy_1[0] = Dy[0]
    for i in range(1,nt):
        Dx_1[i] = (Dx[i] - Dx[i-1])/dt
        Dy_1[i] = (Dx[i] - Dx[i-1])/dt
    Dx_2 = np.zeros(nt)
    Dy_2 = np.zeros(nt)
    Dx_2[0] = Dx_1[0]
    Dy_2[0] = Dy_1[0]
    for j in range(1,nt):
        Dx_2[j] = (Dx_1[j] - Dx_1[j-1])/dt
        Dy_2[j] = (Dy_1[j] - Dy_1[j-1])/dt
    return Dx_2,Dy_2

data = np.loadtxt("./dipole_d.dat")
time , D_s , D_p = data[:,0] , data[:,1] , data[:,2]

t0,t1 = -836,900
nt = 1801
dt = (t1-t0)/(nt-1)
E0 = 0.08
w0 = 0.057
Ip = 0.5724
Up = Ip + 3.17 * E0**2 /(4 * w0**2)
r = np.sqrt(Ip/(2*E0**2/(4*w0**2)))
print(r)

D_S , D_P = Second_derivative(D_s , D_p)

fre_s = np.fft.fft(D_s)
fre_p = np.fft.fft(D_p)

w = 2* np.pi* np.fft.fftfreq(len(time), dt)

plt.plot(np.log10(abs(fre_s)**2))
plt.plot(np.log10(abs(fre_p)**2))
plt.xlim(0,1000)
plt.axvline(w0*2)
plt.axvline(w0*3)

#plt.plot(time,D_s)
#plt.plot(time,D_p)
plt.show()




