import numpy as np
import pylab as pl
import sys
import matplotlib.pyplot as plt

from pylab import *
close('all')

from grid          import Grid
from field         import Fields
from thermo        import Thermo
from timerk3       import Time
from statistics    import Statistics
from micro_kernels import Micro

class Model:
    def __init__(self):
        pass

def init_case(model):
    zi     = 1400
    thl    = 290
    dthldz = 6e-3
    qt     = 7e-3
    dqtdz  = -3.e-6
    labda  = 4000

    for k in range(model.grid.kstart, model.grid.kend):
        if(model.grid.z[k] < zi):
            model.fields.thl.data[k] = thl
            model.fields.qt. data[k] = 7e-3 
        elif(model.grid.z[k] < model.grid.zsize-500):    
            model.fields.thl.data[k] = thl + (model.grid.z[k] - zi) * dthldz 
            model.fields.qt .data[k] = qt  * np.exp((-(model.grid.z[k] - zi))/labda) 
        else:
            model.fields.thl.data[k] = thl + (model.grid.z[k] - zi) * dthldz + (model.grid.z[k] - (model.grid.zsize-500)) * 3 * dthldz
            model.fields.qt .data[k] = qt  * np.exp((-(model.grid.z[k] - zi))/labda) 

    model.fields.thl.data[model.grid.kstart-1] = model.fields.thl.data[model.grid.kstart]
    model.fields.qt .data[model.grid.kstart-1] = model.fields.qt .data[model.grid.kstart]
    model.fields.thl.data[model.grid.kend]     = model.fields.thl.data[model.grid.kend-1]
    model.fields.qt .data[model.grid.kend]     = model.fields.qt .data[model.grid.kend-1]

if(__name__ == "__main__"):
    # Empty class to wrap all components
    model = Model()

    # Create grid, and initialize arrays
    grid = Grid(zsize=4000, ktot=100)
    model.grid = grid
    grid.create()

    # Create fields (thl, qt, qr, nr)
    fields = Fields(model)
    model.fields = fields
    init_case(model)
    
    # Create thermodynamics
    thermo = Thermo(model)
    model.thermo = thermo

    # Time settings (some are defined in time.py):
    time = Time(model)
    model.time = time
    time.total_time = 7200
    time.dt_max = 60
    time.cfl_max = 2

    # Statistics
    stat = Statistics(model)
    model.stat = stat
    stat.sampletime = 600
   
    # Microphysics kernels
    micro = Micro()
    model.micro = micro
    micro.sw_auto = True
    micro.sw_evap = True
    micro.sw_accr = True
    micro.sw_scbr = True
    micro.sw_sedi = True
    
    # Time loop
    while(not time.finished()):
        # Set the time step based on stat, sedimentation, etc.
        time.set_time_step()

        # Slowly decrease qt to disolve cloud
        fields.qt.data[:] -= 1.e-3 / 3600. * model.time.dt

        # Calculate cloud liquid water, pressure, etc.
        thermo.execute()

        # Save statistics before calculating tendencies and integrating the fields
        stat.execute()

        # Loop over time steps
        for rk in range(time.nsubstep):
            T   = fields.thl.data * fields.exn.data
            dt  = time.dt  # Full time step
            sdt = time.get_sub_dt() # RK3 sub time step

            fields.clip_at_zero()

            micro.autoconversion(fields.qr.tend, fields.nr.tend, fields.thl.tend, fields.qt.tend, \
                                 fields.qr.data, fields.nr.data, fields.ql.data, fields.rho.data, \
                                 fields.exn.data, grid.kstart, grid.kend)

            micro.evaporation(fields.qr.tend, fields.nr.tend, fields.thl.tend, fields.qt.tend, \
                              fields.qr.data, fields.nr.data, fields.qt.data, fields.qs.data, \
                              fields.ql.data, T, fields.rho.data, fields.exn.data, grid.kstart, grid.kend)

            micro.accretion(fields.qr.tend, fields.nr.tend, fields.thl.tend, fields.qt.tend, \
                            fields.qr.data, fields.nr.data, fields.ql.data, \
                            fields.rho.data, fields.exn.data, grid.kstart, grid.kend)

            micro.selfcollection_breakup(fields.qr.tend, fields.nr.tend, fields.thl.tend, fields.qt.tend, \
                                         fields.qr.data, fields.nr.data, fields.rho.data, \
                                         grid.kstart, grid.kend)

            micro.sedimentation_ss08(fields.qr.tend, fields.nr.tend, \
                                     fields.qr.data, fields.nr.data, \
                                     fields.rho.data, grid.dz, grid.dzi, grid.dzhi, dt, sdt, \
                                     grid.kstart, grid.kend, grid.kcells)

            # Integrate!
            time.rk3()

        time.time += time.dt
        time.iter += 1

    stat.execute()
    stat.finish()

    # make colors...
    cc = plt.get_cmap(cm.jet)(np.linspace(0, 1, stat.t.size))
 
    if(True):
        figure()
        
        subplot(131)
        for t in range(stat.t.size):
            plot(stat.ql[t,:]*1e3, stat.z, color=cc[t])
        xlabel('ql [g/kg]')

        subplot(132)
        for t in range(stat.t.size):
            plot(stat.qr[t,:]*1e6, stat.z, color=cc[t])
        xlabel('qr [mg/kg]')

        subplot(133)
        for t in range(stat.t.size):
            plot(stat.nr[t,:], stat.z, color=cc[t])
        xlabel('Nr [1/m3]')
  
    if(True):
        figure()
        subplot(131)
        title('ql', loc='left')
        pcolormesh(stat.t/3600., stat.z, np.transpose(stat.ql), cmap=cm.ocean_r)
        colorbar()

        subplot(132)
        title('qr', loc='left')
        pcolormesh(stat.t/3600., stat.z, np.transpose(stat.qr), cmap=cm.ocean_r)
        colorbar()

        subplot(133)
        title('Nr', loc='left')
        pcolormesh(stat.t/3600., stat.z, np.transpose(stat.nr), cmap=cm.ocean_r)
        colorbar()

    if(True):
        figure()
        
        subplot(121)
        for t in range(stat.t.size):
            plot(stat.w_qr[t,:], stat.z, color=cc[t])
        xlabel('w_qr [m/s]')

        subplot(122)
        for t in range(stat.t.size):
            plot(stat.w_nr[t,:], stat.z, color=cc[t])
        xlabel('w_nr [m/s]')

    if(True):
        fig = figure(figsize=(16,10)) #
        #                   L    B    R    T    ws  hs
        fig.subplots_adjust(0.07,0.14,0.97,0.94,0.11,0.3)

        # qr-tendencies -------------------
        ax=subplot(2,5,1)
        title('autoconversion qr', loc='left')
        for t in range(stat.t.size):
            plot(stat.au_qr[t,:]*1e6*3600, stat.z, color=cc[t])
        #xlabel('dqr/dt [mg/kg/hour]')
        legend(frameon=False)

        ax=subplot(2,5,2)
        ax.tick_params(labelleft='off')  
        title('evaporation qr', loc='left')
        for t in range(stat.t.size):
            plot(stat.ev_qr[t,:]*1e6*3600, stat.z, color=cc[t])
        xlabel('dqr/dt [mg/kg/hour]')
        legend(frameon=False)

        ax=subplot(2,5,3)
        ax.tick_params(labelleft='off')  
        title('accretion qr', loc='left')
        for t in range(stat.t.size):
            plot(stat.ac_qr[t,:]*1e6*3600, stat.z, color=cc[t])
        xlabel('dqr/dt [mg/kg/hour]')
        legend(frameon=False)

        ax=subplot(2,5,5)
        ax.tick_params(labelleft='off')  
        title('sedimentation qr', loc='left')
        for t in range(stat.t.size):
            plot(stat.se_qr[t,:]*1e6*3600, stat.z, color=cc[t])
        xlabel('dqr/dt [mg/kg/hour]')
        legend(frameon=False)

        # nr-tendencies -------------------
        ax=subplot(2,5,6)
        title('autoconversion nr', loc='left')
        for t in range(stat.t.size):
            plot(stat.au_nr[t,:]*3600/1000., stat.z, color=cc[t])
        xlabel('dnr/dt [1000/m3/hour]')
        legend(frameon=False)

        ax=subplot(2,5,7)
        ax.tick_params(labelleft='off')  
        title('evaporation nr', loc='left')
        for t in range(stat.t.size):
            plot(stat.ev_nr[t,:]*3600/1000., stat.z, color=cc[t])
        xlabel('dnr/dt [1000/m3/hour]')
        legend(frameon=False)

        ax=subplot(2,5,9)
        ax.tick_params(labelleft='off')  
        title('self-c/breakup qr', loc='left')
        for t in range(stat.t.size):
            plot(stat.sb_nr[t,:]*3600/1000., stat.z, color=cc[t])
        xlabel('dnr/dt [1000/m3/hour]')
        legend(frameon=False)

        ax=subplot(2,5,10)
        ax.tick_params(labelleft='off')  
        title('sedimentation nr', loc='left')
        for t in range(stat.t.size):
            plot(stat.se_nr[t,:]*3600/1000., stat.z, color=cc[t])
        xlabel('dnr/dt [1000/m3/hour]')
        legend(frameon=False)
