import numpy as np
import matplotlib.pyplot as pl
from itertools import cycle

# Import the microphysics model code
from src.model import Model
from src.tools import *

def plot_tendencies(cases=[], plot_every=1):
    """
    Function to plot the tendencies of qr and nr 
    """

    # Define unique colors based on the number of output times
    cc = pl.cm.jet(np.linspace(0,1,cases[0].stat.t.size))

    # Define line types
    lines = ['-', '--', '-.', ':']
    linetypes = cycle(lines)

    # Vertical extent plot
    ymax = cases[0].grid.zsize

    fig = pl.figure(figsize=(16,10)) #
    #                   L    B    R    T    ws  hs
    fig.subplots_adjust(0.06,0.07,0.97,0.97,0.11,0.21)
   
    # qr-tendencies -------------------
    ax=pl.subplot(2,4,1)
    remove_top_right_axis()
    pl.title('Autoconversion q_r', loc='left')
    linetypes = cycle(lines)
    times = []
    names = []
    for ic,case in enumerate(cases):
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            p,=pl.plot(case.stat.au_qr[t,:]*1e6*3600, case.stat.z, color=cc[t], linestyle=lt)

            if (ic == 0): times.append(p)
            if (t  == 0): names.append(p)
    
    legend1 = pl.legend(times, ['t={0:.0f} s'.format(cases[0].stat.t[t]) for t in range(0,cases[0].stat.t.size,plot_every)], frameon=False, loc=4)
    ax.add_artist(legend1)
    legend2 = pl.legend(names, ['{}'.format(case.name) for case in cases], frameon=False, loc=2)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    pl.ylabel('z (m)')
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,2)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Evaporation q_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.ev_qr[t,:]*1e6*3600, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,3)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Accretion q_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.ac_qr[t,:]*1e6*3600, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
    
    ax=pl.subplot(2,4,4)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Sedimentation q_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.se_qr[t,:]*1e6*3600, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    # nr-tendencies -------------------
    ax=pl.subplot(2,4,5)
    remove_top_right_axis()
    pl.title('Autoconversion n_r', loc='left')
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.au_nr[t,:]*3600/1000., case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m3-1 h-1)')
    pl.ylabel('z (m)')
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,6)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Evaporation n_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.ev_nr[t,:]*3600/1000., case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,7)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Self-collection and breakup n_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.sb_nr[t,:]*3600/1000., case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
    
    ax=pl.subplot(2,4,8)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Sedimentation n_r', loc='left')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.se_nr[t,:]*3600/1000., case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 hour-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    pl.savefig('tendencies_qr_nr.pdf') 

def plot_profiles(cases=[], plot_every=1):
    """
    Function to plot all vertical profiles 
    """

    # Define unique colors based on the number of output times
    cc = pl.cm.jet(np.linspace(0,1,cases[0].stat.t.size))

    # Define line types
    lines = ['-', '--', '-.', ':']
    linetypes = cycle(lines)

    # Vertical extent plot
    ymax = cases[0].grid.zsize

    fig = pl.figure(figsize=(16,6))
    #                   L    B    R    T    ws  hs
    fig.subplots_adjust(0.06,0.11,0.98,0.96,0.13,0.21)
   
    ax = pl.subplot(161)
    linetypes = cycle(lines)
    times = []
    names = []
    for ic,case in enumerate(cases):
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            p,=pl.plot(case.stat.thl[t,:]-273.15, case.stat.z, color=cc[t], linestyle=lt)

            if (ic == 0): times.append(p)
            if (t  == 0): names.append(p)
    
    legend1 = pl.legend(times, ['t={0:.0f} s'.format(cases[0].stat.t[t]) for t in range(0,cases[0].stat.t.size,plot_every)], frameon=False, loc=4)
    ax.add_artist(legend1)
    legend2 = pl.legend(names, ['{}'.format(case.name) for case in cases], frameon=False, loc=2)

    pl.xlabel('thl (C)')
    pl.ylabel('z (m)')
    remove_top_right_axis()
    pl.ylim(0,ymax)

    ax=pl.subplot(162)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.qt[t,:]*1000, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('qt (g kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(163)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.ql[t,:]*1000, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('ql (g kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(164)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.qr[t,:]*1e6, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('qr (mg kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(165)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            pl.plot(case.stat.nr[t,:]/1e3, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('nr (dm-3)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
 
    ax=pl.subplot(166)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.t.size, plot_every):
            if (np.any(case.stat.qr[t,:]>0)):
                # Calculate mean mass/drop
                m = (case.stat.rho[t,:]*case.stat.qr[t,:])/(case.stat.nr[t,:]+1e-12)
                # Calculate mean diameter
                pirhow = np.pi * 1e3 / 6.
                D = (m / pirhow)**(1./3.)
                # Mask small drops:
                D = np.ma.masked_less(D, 1e-5)

                pl.semilogx(D*1e3, case.stat.z, color=cc[t], linestyle=lt)
    pl.xlabel('mean diameter (mm)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    pl.savefig('profiles.pdf')

class Case:
    def __init__(self):
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
        self.nccn    = 60e6      # Cloud droplet number
        self.sw_auto = True      # Enable/disable   autoconversion
        self.sw_evap = True      #    "     "       evaporation
        self.sw_accr = True      #    "     "       accretion
        self.sw_scbr = True      #    "     "       self-collection and breakup
        self.sw_sedi = True      #    "     "       sedimentation

        # Time:
        self.ttot    = 3600      # Total integration time (s)

        # Vertical grid:
        self.zsize   = 2600      # Vertical extent domain
        self.ktot    = 100       # Number of vertical grid levels

        # Output / statistics
        self.dt_stat = 90       # Time interval of the statistics output

if (__name__ == "__main__"):
    pl.close('all')
   
    if("case" not in locals()):     
        # Create default case:
        case = Case()
        
        # Edit (if necessary) and run model
        case.nccn = 50e6
        r1 = Model(case, 'nccn=50e6')

        #case.nccn = 60e6
        #r2 = Model(case, 'nccn=60e6')
        
        case.nccn = 100e6
        r3 = Model(case, 'nccn=100e6')

    # Predefined functions to plot vertical profiles and tendencies
    plot_tendencies([r1,r3], plot_every=10)
    plot_profiles([r1,r3], plot_every=10)

    # Model output is available in run.stat; e.g. r1.stat.qr, r1.stat.z, etc.
    pl.figure()
    pl.subplot(121)
    pl.title('qr (g kg-1), ncc=50e6', loc='left')
    pl.pcolormesh(r1.stat.t, r1.stat.z, r1.stat.qr.transpose()*1e3, vmin=0, vmax=0.5, cmap=pl.cm.gist_earth_r)
    pl.xlim(0,r1.stat.t.max())
    pl.ylim(0,r1.stat.z.max())
    pl.xlabel('t (s)')
    pl.ylabel('z (m)')
    pl.colorbar()

    pl.subplot(122)
    pl.title('qr (g kg-1), ncc=100e6', loc='left')
    pl.pcolormesh(r3.stat.t, r3.stat.z, r3.stat.qr.transpose()*1e3, vmin=0, vmax=0.5, cmap=pl.cm.gist_earth_r)
    pl.xlim(0,r3.stat.t.max())
    pl.ylim(0,r3.stat.z.max())
    pl.xlabel('t (s)')
    pl.ylabel('z (m)')
    pl.colorbar()
