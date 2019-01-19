import numpy as np
from .constants import *
from .thermo import esat

K_t     = 2.5e-2   # conductivity of heat [J/(sKm)]
D_v     = 3.e-5    # diffusivity of water vapor [m2/s]
rho_w   = 1.e3     # Density water
rho_0   = 1.225    # SB06, p48
pirhow  = np.pi * rho_w / 6.
xc_min  = 4.2e-15  # Min mean mass of cloud droplet (DALES)
xc_max  = 2.6e-10  # Max mean mass of cloud droplet (DALES)
xr_min  = xc_max   # Min mean mass of precipitation drop (DALES) 
xr_max  = 5e-6     # Max mean mass of precipitation drop, (DALES)
ql_min  = 1.e-7    # Min cloud liquid water for which calculations are performed (DALES) 
qr_min  = 1.e-20   # Min rain liquid water for which calculations are performed (UCLA-LES) 

def calc_mur(Dr):
    return 10. * (1. + np.tanh(1200 * (Dr - 0.0014))) # SS08

def calc_drop_properties(qr, nr, rho):
    #nr_lim = np.maximum(nr, 1)                          # Lower limit of one drop (because qr>0, we should have one drop...)
    xr     = rho * qr / (nr+1e-12)                       # Mean mass of prec. drops (kg)
    xr     = np.minimum(np.maximum(xr, 1e-16), xr_max)   # Limit mass
    Dr     = (xr / pirhow)**(1./3.)                      # Mean diameter of prec. drops (m)

    return (xr,Dr)

def int2(a,b):
    return 0.5*a + 0.5*b

def int4(a,b,c,d):
    ci0  = -1./16.; ci1  =  9./16.;
    return ci0*a + ci1*b + ci1*c + ci0*d

def minmod(a,b):
    return np.sign(a) * max(0., min(np.abs(a), np.sign(a)*b))

class Micro:
    def __init__(self):
        self.scheme = 'KK00'

        self.max_cfl = 1e-3
        self.nc    = 70e6

        self.sw_auto = True
        self.sw_evap = True
        self.sw_accr = True
        self.sw_scbr = True
        self.sw_sedi = True

    def get_time_limit(self, cfl, dt):
        if(self.max_cfl > 0):
            return cfl / self.max_cfl * dt
        else:
            return 1e12

    def autoconversion(self, qr_tend, nr_tend, thl_tend, qt_tend, \
                       qr, nr, ql, rho, exn, kstart, kend):
        
        if(self.sw_auto): 
            x_star = 2.6e-10    # SB06, list of symbols, same as UCLA
            k_cc   = 9.44e+9    # UCLA-LES (Long, 1974), 4.44e9 in SB06, p48  
            nu_c   = 1          # SB06, Table 1; same as UCLA
            kccxs  = k_cc / (20 * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, 2) 

            for k in range(kstart, kend):
                if(ql[k] > ql_min):
                    xc           = rho[k] * ql[k] / self.nc           # Mean mass of cloud drops
                    tau          = max(0, 1 - ql[k] / (ql[k] + qr[k]))  # SB06, Eq 5
                    phi_au       = 400 * pow(tau, 0.7) * \
                                   pow(1. - pow(tau, 0.7), 3)           # SB06, Eq 6

                    # Tendencies
                    # Seifert and Beheng (SB06), eq 4: rho_0 is the surface density; in SB06 rho_0 / rho, which simply becomes rho_0 in the conversion to kg/kg
                    # Khairoutdinov and Kogan (KK00), eq 29 (which has nc in cm-3)
                    qr_t = self.auto_tuning_prefac * kccxs * pow(ql[k], 2) * pow(xc, 2) * (1. + phi_au / pow(1 - tau, 2)) * rho_0 if self.scheme == 'SB06' else self.auto_tuning_prefac * 1350 * pow(ql[k], 2.47) * pow(self.nc * 1e-6, self.auto_exponent_KK)
                    nr_t = qr_t * rho[k] / x_star

                    qr_tend[k]  += qr_t 
                    nr_tend[k]  += nr_t
                    thl_tend[k] += Lv / (cp * exn[k]) * qr_t 
                    qt_tend[k]  -= qr_t

    def evaporation(self, qr_tend, nr_tend, thl_tend, qt_tend, qr, nr, qt, qs, ql, T, rho, exner, kstart, kend):
        if(self.sw_evap):
            lambda_evap = 1.

            # Constants ventilation term:
            nu_a   = 1.4086e-5  # SB06, page 62
            a      = 9.65       # S08, p3618
            b      = 9.8        # S08, p3618
            c      = 600        # S08, p3618
            a_vent = 0.78       # SB06, p62
            b_vent = 0.308      # SB06, p62
            Nsc    = 0.71       # SB06, p63

            for k in range(kstart, kend):
                if(qr[k] > qr_min):
                    xr, Dr = calc_drop_properties(qr[k], nr[k], rho[k])

                    G      = 1. / ( (Rv * T[k] / (D_v * esat(T[k]))) + (Lv * Lv / (K_t * Rv * T[k] * T[k])) )
                    S      = (qt[k] - ql[k]) / qs[k] - 1 # Supersaturation (-)

                    # Ventilation:
                    #Vr     = a - b * np.exp(-c * Dr)     # Terminal velocity drop of diamter Dr [m s-1]
                    #Nre    = max(0, Vr * Dr / nu_a)          
                    #F0     = a_vent + b_vent * 0.71**(1./3.) * Nre**(1/2) 
                    F      = 1. #F0

                    qr_t   = 2. * np.pi * Dr * G * S * F * nr[k] / rho[k]
                    nr_t   = lambda_evap * qr_t * rho[k] / xr

                    # Accumulate tendencies
                    qr_tend[k]  += qr_t
                    nr_tend[k]  += nr_t
                    thl_tend[k] += Lv / (cp * exner[k]) * qr_t 
                    qt_tend[k]  -= qr_t

    def accretion(self, qr_tend, nr_tend, thl_tend, qt_tend, qr, nr, ql, rho, exner, kstart, kend):
        if(self.sw_accr):
            k_cr  = 5.25 # SB06, p49

            for k in range(kstart, kend):
                if(ql[k] > ql_min and qr[k] > qr_min):
                    tau    = 1 - ql[k] / (ql[k] + qr[k]) # SB06, Eq 5
                    phi_ac = pow(tau / (tau + 5e-5), 4)  # SB06, Eq 8

                    # accretion rate:
                    # Seifert and Behen (SB06): eq 7
                    # Khairoutdinov and Kogan (KK00): eq 33
                    qr_t   = k_cr * ql[k] *  qr[k] * phi_ac * pow(rho_0 / rho[k], 0.5) if self.scheme == 'SB06' else 67. * pow(ql[k] * qr[k], 1.15)

                    # Accumulate tendencies
                    qr_tend[k]  += qr_t 
                    qt_tend[k]  -= qr_t
                    thl_tend[k] += Lv / (cp * exner[k]) * qr_t 

    def selfcollection_breakup(self, qr_tend, nr_tend, thl_tend, qt_tend, qr, nr, rho, kstart, kend):
        if(self.sw_scbr):
            k_rr     = 7.12   # SB06, p49
            kappa_rr = 60.7   # SB06, p49 

            D_eq     = 0.9e-3 # SB06, list of symbols 
            k_br1    = 1.0e3  # SB06, p50, for 0.35e-3 <= Dr <= D_eq
            k_br2    = 2.3e3  # SB06, p50, for Dr > D_eq
   
            for k in range(kstart, kend):
                if(qr[k] > qr_min):
                    xr, Dr   = calc_drop_properties(qr[k], nr[k], rho[k])
                    mur      = calc_mur(Dr)
                    lambda_r = pow((mur+3)*(mur+2)*(mur+1), 1./3.) / Dr

                    # self-collection tendency
                    sc_t     = -k_rr * nr[k] * qr[k]*rho[k] * pow(1. + kappa_rr / lambda_r * pow(pirhow, 1./3.), -9) * pow(rho_0 / rho[k], 0.5)

                    br_t = 0
                    dDr  = Dr - D_eq
                    if(Dr > 0.35e-3):
                        if(Dr <= D_eq):
                            phi_br = k_br1 * dDr
                        elif(Dr > D_eq):
                            phi_br = 2. * np.exp(k_br2 * dDr) - 1. 
                          
                        # Breakup tendency 
                        br_t = -(phi_br + 1) * sc_t

                    nr_tend[k] += sc_t + br_t

    def sedimentation_velocity(self, w_qr, w_nr, qr, nr, rho, kstart, kend):
        w_max = 9.65    # Maximum sedimentation velocity
        w_min = 0.01    # Minimum sedimentation velocity
        a_R   = 9.65    # SB06, p51
        c_R   = 600     # SB06, p51
        Dv    = 25.e-6  # UCLA-LES
        b_R   = a_R * np.exp(c_R*Dv) # UCLA-LES

        for k in range(kstart, kend):
            if(qr[k] > qr_min):
                xr, Dr   = calc_drop_properties(qr[k], nr[k], rho[k])
                mur      = calc_mur(Dr)
                lambda_r = pow((mur+3)*(mur+2)*(mur+1), 1./3.) / Dr

                # SS08:
                w_qr[k]  = -max(w_min, min(w_max, a_R - b_R * pow(1. + c_R/lambda_r, -(mur+4))))
                w_nr[k]  = -max(w_min, min(w_max, a_R - b_R * pow(1. + c_R/lambda_r, -(mur+1))))

    def sedimentation_ss08(self, qr_tend, nr_tend, qr, nr, rho, dz, dzi, dzhi, dt, subdt, kstart, kend, kcells, rainrate=None):
        if(self.sw_sedi):
            # 1. Calculate sedimentation velocity
            w_qr = np.zeros(kcells)
            w_nr = np.zeros(kcells)
            self.sedimentation_velocity(w_qr, w_nr, qr, nr, rho, kstart, kend)
    
            # 2. Calculate CFL number based on interpolated velocity
            c_qr = np.zeros(kcells)
            c_nr = np.zeros(kcells)

            c_qr[kstart] = -w_qr[kstart] * dzi[kstart] * dt 
            c_nr[kstart] = -w_nr[kstart] * dzi[kstart] * dt 

            for k in range(kstart+1, kend):
                c_qr[k] = -0.25 * (w_qr[k+1] + 2*w_qr[k] + w_qr[k-1]) * dzi[k] * dt
                c_nr[k] = -0.25 * (w_nr[k+1] + 2*w_nr[k] + w_nr[k-1]) * dzi[k] * dt

            self.max_cfl = c_qr.max()

            # Calculate slopes
            qr_slope = np.zeros(kcells)
            nr_slope = np.zeros(kcells)

            for k in range(kstart, kend):
                qr_slope[k] = minmod(qr[k]-qr[k-1], qr[k+1]-qr[k])
                nr_slope[k] = minmod(nr[k]-nr[k-1], nr[k+1]-nr[k])

            qr_slope[kstart-1] = qr_slope[kstart]
            nr_slope[kstart-1] = nr_slope[kstart]

            # calculate flux and tendency
            qr_flux = np.zeros(kcells)
            nr_flux = np.zeros(kcells)

            for k in range(kend-1, kstart-1, -1):
                kk  = k
                tot = 0
                zz  = 0
                cc  = min(1., c_qr[k])
                while(cc > 0 and kk < kend):
                    tot += rho[kk] * (qr[kk] + qr_slope[kk]*(1.-cc)) * cc * dz[kk]
                    zz  += dz[kk]
                    kk  += 1
                    cc  = min(1., c_qr[kk] - zz * dzi[kk])

                lim = rho[k] * dz[k] * qr[k] - qr_flux[k+1] * dt + small
                #if(lim < tot):
                #    print('limiter qr!!')
                tot = min(tot, lim)
                qr_flux[k] = -tot / dt 

                kk  = k
                tot = 0
                zz  = 0
                cc  = min(1., c_nr[k])
                while(cc > 0 and kk < kend):
                    tot += rho[kk] * (nr[kk] + nr_slope[kk]*(1.-cc)) * cc * dz[kk]
                    zz  += dz[kk]
                    kk  += 1
                    cc  = min(1., c_nr[kk] - zz * dzi[kk])

                lim = rho[k] * dz[k] * nr[k] - nr_flux[k+1] * dt + small
                #if(lim < tot):
                #    print('limiter nr!!')
                tot = min(tot, lim)
                nr_flux[k] = -tot / dt 

                qr_tend[k] += -(qr_flux[k+1] - qr_flux[k]) * dzi[k] / rho[k] 
                nr_tend[k] += -(nr_flux[k+1] - nr_flux[k]) * dzi[k] / rho[k]

                if rainrate is not None: rainrate[k] = - qr_flux[k] / rho_w
