import numpy as np
import matplotlib.pyplot as pl
from itertools import cycle

# Import the microphysics model code
from src.model import Model
from src.tools import *
from src.plot_util import plot_tendencies, plot_profiles



class Case:
    def __init__(self, nc=60e6, sw_auto=True, sw_accr=True):
        # Initial vertical profiles:
        self.zi1     = 500       # Mixed-layer depth (m)
        self.thl     = 298       # Mixed-layer liquid water potential temperature (K)
        self.qt      = 16e-3     # Mixed-layer total specific humidity (kg kg-1)

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
        self.sw_evap = True      #    "     "       evaporation
        self.sw_accr = sw_accr   #    "     "       accretion
        self.sw_scbr = False     #    "     "       self-collection and breakup
        self.sw_sedi = True      #    "     "       sedimentation

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
