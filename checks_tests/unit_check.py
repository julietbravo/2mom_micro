# Unit checking microphysics scheme
# Output:
# --------------------------------------------
# > Base units: qr = dimensionless, Nr = 1 / meter ** 3
# > Units autoconversion: dqr/dt = 1 / second, dNr/dt = 1 / meter ** 3 / second
# > Units autoconversion: dqt/dt = 1 / second, dtheta/dt = kelvin / second
# > Units evaporation: dqr/dt = 1 / second, dNr/dt = 1 / meter ** 3 / second
# > Units evaporation: dqt/dt = 1 / second, dtheta/dt = kelvin / second
# --------------------------------------------

import numpy as np
import pint
ureg = pint.UnitRegistry()

# J  = kg * m**2 * s**-2
# pa = kg * m**-1 * s**-2 

# Constants
Rv  = 287.   * ureg['kg * m**2 * s**-2 * kg**-1 * K**-1']
Lv  = 2.45e6 * ureg['kg * m**2 * s**-2 * kg**-1']
cp  = 1004.  * ureg['kg * m**2 * s**-2 * kg**-1 * K**-1']

D_v = 3.e-5  * ureg['m**2 * s**-1']                                # Diffusivity of water vapor [m2/s]
K_t = 2.5e-2 * ureg['kg * m**2 * s**-2 * s**-1 * K**-1 * m**-1 ']  # conductivity of heat [J/(sKm)]

# Define base units for calculations
rho   = np.random.random() * ureg['kg * m**-3']         # Density
qt    = np.random.random() * ureg['kg * kg**-1']        # Total specfic humidity
qs    = np.random.random() * ureg['kg * kg**-1']        # Saturation specfic humidity
qr    = np.random.random() * ureg['kg * kg**-1']        # Rain liquid water
ql    = np.random.random() * ureg['kg * kg**-1']        # Cloud liquid water
Nc0   = np.random.random() * ureg['m**-3']              # Number density cloud droplets
Nr    = np.random.random() * ureg['m**-3']              # Number density rain
T     = np.random.random() * ureg['K']                  # Temperature
e     = np.random.random() * ureg['kg * m**-1 * s**-2'] # Vapor pressure
exner = np.random.random() * ureg['']                   # Exner

print('Base units: qr = {0.units}, Nr = {1.units}'.format(qr, Nr))

# ---------------------------
# Autoconversion
# ---------------------------
if(True):
    x_star = 2.6e-10 * ureg['kg']                       # Separating drop mass (SB06, list of symbols)
    k_cc   = 4.44e9  * ureg['m**3 * kg**-2 * s**-1']    # SB06, p48 
    nu_c   = 1       * ureg['']                         # SB06, Table 1
    kccxs  = k_cc / (20 * x_star) * (nu_c+2)*(nu_c+4)/((nu_c+1)**2) 

    xc     = rho * ql / Nc0
    tau    = 1 - (ql / (ql+qr))  
    phi_au = 400. * tau**0.7 * (1-tau**0.7)**3. 

    qr_t   = kccxs * ql**2. * xc**2. * (1. + phi_au / (1 - tau)**2.) * rho # SB06, eq 4
    nr_t   = qr_t * rho / x_star
    qt_t   = -qr_t
    th_t   = Lv / (cp * exner) * qr_t 

    print('Units autoconversion: dqr/dt = {0.units}, dNr/dt = {1.units}'.format(qr_t, nr_t))
    print('Units autoconversion: dqt/dt = {0.units}, dtheta/dt = {1.units}'.format(qt_t, th_t))

# ---------------------------
# Evaporation
# ---------------------------
if(True):
    pi      = 1. * ureg['']
    rho_w   = 1. * ureg['kg * m**-3']
    pirhow  = pi * rho_w / 6.
    nu_air  = 1.4086e-5 * ureg['m**2 * s**-1'] # SB06, page 62
    lambda_evap = 1. * ureg['']

    xr  = rho * qr / Nr               # Mean mass of prec. drops (kg)
    Dr  = (xr / pirhow)**(1./3.)      # Mean diameter of prec. drops (m)

    # Ventilation term rain drop (unitless)
    # ------------------------
    a   = 9.65 * ureg['m * s**-1']  # S08, p3618
    b   = 9.8  * ureg['m * s**-1']  # S08, p3618
    c   = 600  * ureg['m**-1']      # S08, p3618

    a_vent = 0.78    # SB06, p62
    b_vent = 0.308   # SB06, p62
    Nsc    = 0.71    # SB06, p63

    Vr   = a - b * np.exp(-c * Dr)     # Terminal velocity drop of diamter Dr [m s-1]
    Nre  = Vr * Dr / nu_air          
    Fv   = a_vent + b_vent * 0.71**(1./3.) * Nre**(1/2) 
    # ------------------------

    # Saturation
    S   = (qt - ql) / qs - 1 # Supersaturation (-)

    G1   = Rv * T / (D_v * e)
    G2   = Lv * Lv / (K_t * Rv * T * T)
    Glv  = 1. / (G1 + G2) 

    qr_t = 2. * np.pi * Dr * Glv * Fv * S * Nr / rho
    nr_t = lambda_evap * qr_t * rho / xr 
    qt_t = -qr_t
    th_t = Lv / (cp * exner) * qr_t

    print('Units evaporation: dqr/dt = {0.units}, dNr/dt = {1.units}'.format(qr_t, nr_t))
    print('Units evaporation: dqt/dt = {0.units}, dtheta/dt = {1.units}'.format(qt_t, th_t))
