import ctypes as ct
import numpy as np

libsfa = ct.CDLL('./libsfa.so')

def arr (data):
    return data.ctypes.data_as(ct.POINTER(ct.c_double))


def calc_dipole_CC_d2 (time, ES, AS, alphaS, EP, AP, alphaP):
    func = libsfa.calc_dipole_CC_d2
    func.argtypes = [ct.c_long, ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double),
                     ct.POINTER(ct.c_double)]
    func.restype = None
    nt = len(time)
    dipole_d_S_d2 = np.empty(nt)
    dipole_d_P_d2 = np.empty(nt)
    func(nt, arr(time), arr(ES), arr(AS), arr(alphaS), arr(EP), arr(AP), arr(alphaP), arr(dipole_d_S_d2), arr(dipole_d_P_d2))
    return dipole_d_S_d2, dipole_d_P_d2


# def calc_wt (time, ES, AS, alphaS, EP, AP, alphaP):
#     func = libsfa.calc_wt
#     func.argtypes = [ct.c_long, ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double)]
#     func.restype = None
#     nt = len(time)
#     dipole_d_d2 = np.empty(nt)
#     func(nt, arr(time), arr(ES), arr(AS), arr(alphaS), arr(EP), arr(AP), arr(alphaP), arr(dipole_d_d2))
#     return dipole_d_d2


# def calc_Nt (time, ES, AS, alphaS, EP, AP, alphaP):
#     func = libsfa.calc_Nt
#     func.argtypes = [ct.c_long, ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double),
#                      ct.POINTER(ct.c_double)]
#     func.restype = None
#     nt = len(time)
#     dipole_d_d2 = np.empty(nt)
#     func(nt, arr(time), arr(ES), arr(AS), arr(alphaS), arr(EP), arr(AP), arr(alphaP), arr(dipole_d_d2))
#     return dipole_d_d2
