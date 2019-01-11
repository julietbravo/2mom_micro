import numpy as np
import matplotlib.pyplot as pl
from itertools import cycle

# Import the microphysics model code
from src.model import Model
from src.tools import *
from src.thermo import calc_zLCL
from src.plot_util import plot_tendencies, plot_profiles

def run_column_model(run_name='run', thl=298, qt=16e-3, nc=60e6, sw_auto=True, sw_accr=True, sw_evap=True, sw_scbr=False, sw_sedi=True):
    case = Case(
              thl=thl,      # liquid water potential temperature, unit: K
              qt=qt,     # total water mixing ratio, unit: kg/kg
              nc=nc,      # cloud droplet number concentration, unit: 1/m3
              sw_auto=sw_auto, # autoconverion: on
              sw_accr=sw_accr, # accretion: on
              sw_evap=sw_evap, # rain drop evaporation: on
              sw_sedi=True, # rain drop sedimentation: on
              )
    return Model(case, run_name)


class Case:
    def __init__(self, thl=298, qt=16e-3, nc=60e6, sw_auto=True, sw_accr=True, sw_evap=True, sw_scbr=False, sw_sedi=True):
        # Initial vertical profiles:
        self.zi1     = calc_zLCL(thl, qt, 101325)       # Mixed-layer depth (m)
        self.thl     = thl       # Mixed-layer liquid water potential temperature (K)
        self.qt      = qt     # Mixed-layer total specific humidity (kg kg-1)

        self.dthldz  = 2.3e-3    # Cloud layer liquid water potential temperature lapse rate (K m-1)
        self.dqtdz   = -2.8e-6   # Cloud layer specific humidity lapse rate (kg kg-1 m-1)
        self.zi2     = 2200      # Cloud top inversion height (m)

        # Large scale tendencies:
        self.dthldt  = 0.        # Large scale warming/cooling tendency (K s-1)
        self.dqtdt   = 0.        # Large scale moistening/drying tendency (kg kg-1 s-1)

        # Microphysics settings:
        self.scheme = 'KK00'
        self.nc      = nc        # Cloud droplet number
        self.sw_auto = sw_auto   # Enable/disable   autoconversion
        self.sw_evap = sw_evap     #    "     "       evaporation
        self.sw_accr = sw_accr   #    "     "       accretion
        self.sw_scbr = sw_scbr     #    "     "       self-collection and breakup
        self.sw_sedi = sw_sedi    #    "     "       sedimentation

        # Time:
        self.ttot    = 3600      # Total integration time (s)

        # Vertical grid:
        self.zsize   = 2600      # Vertical extent domain
        self.ktot    = 100       # Number of vertical grid levels

        # Output / statistics
        self.dt_stat = 90        # Time interval of the statistics output

if (__name__ == "__main__"):
    pl.close('all')
   
    if("case" not in locals()):     
        # Create default case:
        case = Case()
        
        # Edit (if necessary) and run model
        case.nc = 50e6
        r1 = Model(case, 'nc=50e6')

        #case.nc = 60e6
        #r2 = Model(case, 'nc=60e6')
        
        case.nc = 100e6
        r3 = Model(case, 'nc=100e6')

    # Predefined functions to plot vertical profiles and tendencies
    plot_tendencies([r1,r3], plot_every=10)
    plot_profiles([r1,r3], plot_every=10)

    # Model output is available in run.stat; e.g. r1.stat.qr, r1.stat.z, etc.
    pl.figure()
    pl.subplot(121)
    pl.title('qr (g kg-1), ncc=50e6')
    pl.pcolormesh(r1.stat.var['t'], r1.stat.z, r1.stat.qr.transpose()*1e3, vmin=0, vmax=0.5, cmap=pl.cm.gist_earth_r)
    pl.xlim(0,r1.stat.var['t'].max())
    pl.ylim(0,r1.stat.z.max())
    pl.xlabel('t (s)')
    pl.ylabel('z (m)')
    pl.colorbar()

    pl.subplot(122)
    pl.title('qr (g kg-1), ncc=100e6')
    pl.pcolormesh(r3.stat.var['t'], r3.stat.z, r3.stat.qr.transpose()*1e3, vmin=0, vmax=0.5, cmap=pl.cm.gist_earth_r)
    pl.xlim(0,r3.stat.var['t'].max())
    pl.ylim(0,r3.stat.z.max())
    pl.xlabel('t (s)')
    pl.ylabel('z (m)')
    pl.colorbar()
    pl.savefig('qr.png')
