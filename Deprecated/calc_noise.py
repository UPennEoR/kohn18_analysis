import numpy as np
from wedgie import cosmo_utils as cu
from astropy import units as u

def Tsys(f,Trcvr):
    Tsky = 180*u.K*np.power(f/(0.180*u.GHz),-2.55)
    return Tsky + Trcvr

def delta2_noise(k,Tsys,B,Opp,t_int,Nbl,Npol):
    delta2_noise = 150 * (k/(0.1*0.7/u.Mpc))**2.5 *\
                   (6*u.MHz/B)**0.5 * (Opp/(0.76*u.steradian))**1.25 *\
                   (Tsys/(500*u.K))**2 * (6*u.hour/t_int)**0.5 * Nbl**-0.5 * Npol**-1
    return delta2_noise*(u.mK**2)

def delta2_noise_2(k,B,Opp,T_sys,t_days,u_bl,Npol):
    noise = 2.48e4 * (k/(0.1*0.7/u.Mpc))**2.5 *\
            (6*u.MHz/B)**0.5 * (Opp/(0.76*u.steradian))**1.5 *\
            (T_sys/(500*u.K))**2 * (120*u.day/t_days) * (u_bl/20.) * (1/Npol)
    return noise*(u.mK**2)

bl_occ = [(14.6, 25), (25.2879, 19), (29.2, 15), (38.628, 22), (43.8, 9), (50.5759, 6), (52.641, 8), (58.4, 1)]

f_lo = 0.125*u.GHz
f_hi = 0.160*u.GHz
Tsys_lo = Tsys(f_lo,600*u.K)
Tsys_hi = Tsys(f_hi,300*u.K)
omega_eff = np.polyval(cu.HERA_BEAM_POLY,f_hi.to(u.GHz).value)*u.steradian
#t_int = 10 * (60*u.second/u.minute) * (60*u.minute/u.hour) * (10*u.hour) * 8
t_int = 8*u.day * (10*u.hour/u.day)
B = 10*u.MHz

noisearr = []
for (bl,occ) in bl_occ:
    kperp_hi = cu.uv2kpr(bl*u.m,f_hi)
    kpl = (0.2/0.7) * (1/u.Mpc)
    k_hi = np.sqrt(kperp_hi**2 + kpl**2)
    #d2k_noise_hi = delta2_noise(k_hi, Tsys_hi, 10*u.MHz,omega_eff*u.steradian,  t_int, occ, 2)
    #pk_noise_hi = (2*(np.pi**2)/(k_hi**3))*d2k_noise_hi
    u_bl = bl/1.5
    Npol = 2.
    d2k_noise_hi = delta2_noise_2(k_hi,B,omega_eff,Tsys_hi,t_int,u_bl,Npol)
    noisearr.append(d2k_noise_hi)

import IPython;IPython.embed()
