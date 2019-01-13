import numpy as np
import matplotlib.pyplot as plt
from .tools import remove_top_right_axis
from itertools import cycle

def plot_tendencies(cases=[], plot_every=1):
    """
    Function to plot the tendencies of qr and nr 
    """

    # Define unique colors based on the number of output times
    cc = pl.cm.jet(np.linspace(0,1,cases[0].stat.var['t'].size))

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
    pl.title('Autoconversion q_r')
    linetypes = cycle(lines)
    times = []
    names = []
    for ic,case in enumerate(cases):
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            p,=pl.plot(case.stat.var['au_qr'][t,:]*1e6*3600, case.stat.var['z'], color=cc[t], linestyle=lt)

            if (ic == 0): times.append(p)
            if (t  == 0): names.append(p)
    
    legend1 = pl.legend(times, ['t={0:.0f} s'.format(cases[0].stat.var['t'][t]) for t in range(0,cases[0].stat.var['t'].size,plot_every)], loc=4)
    ax.add_artist(legend1)
#    legend2 = pl.legend(names, ['{}'.format(case.name) for case in cases], loc=2)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    pl.ylabel('z (m)')
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,2)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Evaporation q_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['ev_qr'][t,:]*1e6*3600, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,3)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Accretion q_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['ac_qr'][t,:]*1e6*3600, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
    
    ax=pl.subplot(2,4,4)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Sedimentation q_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['se_qr'][t,:]*1e6*3600, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dqr/dt (mg kg-1 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    # nr-tendencies -------------------
    ax=pl.subplot(2,4,5)
    remove_top_right_axis()
    pl.title('Autoconversion n_r')
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['au_nr'][t,:]*3600/1000., case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m3-1 h-1)')
    pl.ylabel('z (m)')
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,6)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Evaporation n_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['ev_nr'][t,:]*3600/1000., case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(2,4,7)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Self-collection and breakup n_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['sb_nr'][t,:]*3600/1000., case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 h-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
    
    ax=pl.subplot(2,4,8)
    remove_top_right_axis()
    linetypes = cycle(lines)
    pl.title('Sedimentation n_r')
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['se_nr'][t,:]*3600/1000., case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('dnr/dt (1000 m-3 hour-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    pl.savefig('tendencies_qr_nr.pdf') 

def plot_profiles(cases=[], plot_every=1):
    """
    Function to plot all vertical profiles 
    """

    # Define unique colors based on the number of output times
    cc = pl.cm.jet(np.linspace(0,1,cases[0].stat.var['t'].size))

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
        for t in range(0, case.stat.var['t'].size, plot_every):
            p,=pl.plot(case.stat.var['thl'][t,:]-273.15, case.stat.var['z'], color=cc[t], linestyle=lt)

            if (ic == 0): times.append(p)
            if (t  == 0): names.append(p)
    
    legend1 = pl.legend(times, ['t={0:.0f} s'.format(cases[0].stat.var['t'][t]) for t in range(0,cases[0].stat.var['t'].size,plot_every)], loc=4)
    ax.add_artist(legend1)
#    legend2 = pl.legend(names, ['{}'.format(case.name) for case in cases],loc=2)

    pl.xlabel('thl / C')
    pl.ylabel('z / m')
    remove_top_right_axis()
    pl.ylim(0,ymax)

    ax=pl.subplot(162)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['qt'][t,:]*1000, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('qt / (g kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(163)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['ql'][t,:]*1000, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('ql (g kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(164)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['qr'][t,:]*1e6, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('qr (mg kg-1)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    ax=pl.subplot(165)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            pl.plot(case.stat.var['nr'][t,:]/1e3, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('nr (dm-3)')
    ax.set_yticks([])
    pl.ylim(0,ymax)
     
    ax=pl.subplot(166)
    remove_top_right_axis()
    linetypes = cycle(lines)
    for case in cases:
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
            if (np.any(case.stat.var['qr'][t,:]>0)):
                # Calculate mean mass/drop
                m = (case.stat.var['rho'][t,:]*case.stat.var['qr'][t,:])/(case.stat.var['nr'][t,:]+1e-12)
                # Calculate mean diameter
                pirhow = np.pi * 1e3 / 6.
                D = (m / pirhow)**(1./3.)
                # Mask small drops:
                D = np.ma.masked_less(D, 1e-5)

                pl.semilogx(D*1e3, case.stat.var['z'], color=cc[t], linestyle=lt)
    pl.xlabel('mean diameter (mm)')
    ax.set_yticks([])
    pl.ylim(0,ymax)

    pl.savefig('profiles.pdf')



def plot_profile_evolution(cases, varname, plot_every=10):
    specs = {
            'thl': ['liquid water potential temperature / K', 1],
            'th': ['potential temperature / K', 1],
            'qt': ['total water mixing ratio / (g/kg)', 1e3],
            'qv': ['water vapor mixing ratio / (g/kg)', 1e3],
            'qs': ['saturation water vapor mixing ratio / (g/kg)', 1e3],
            'au_qr': ['autoconversion / (mg/kg/h)', 1e6*3600],
            'ac_qr': ['accretion / (mg/kg/h)', 1e6*3600],
            'ev_qr': ['rain evaporation / (mg/kg/h)', 1e6*3600],
            'ql': ['cloud droplet water mixing ratio / (g/kg)', 1e3],
            'qr': ['rain drop water mixing ratio / (mg/kg)', 1e6],
            }
    assert varname in specs.keys()
    xlabel, fac = specs[varname]

    # Define unique colors based on the number of output times
    cc = plt.cm.jet(np.linspace(0,1,cases[0].stat.var['t'].size))

    # Define line types
    lines = ['-', '--', '-.', ':']
    linetypes = cycle(lines)

    # Vertical extent plot
    ymax = cases[0].grid.zsize
   
    plt.figure(figsize=(3,5))
    remove_top_right_axis()
    
    times = []
    names = []
    for ic,case in enumerate(cases):
        lt = next(linetypes)
        for t in range(0, case.stat.var['t'].size, plot_every):
    #        pl.plot(case.stat.se_qr[t,:]*1e6*3600, case.stat.z, color=cc[t], linestyle=lt)
    #pl.xlabel('dqr/dt (mg kg-1 h-1)')
            p, = plt.plot(case.stat.var[varname][t,:]*fac, case.stat.var['z'], color=cc[t], linestyle=lt)

            if (ic == 0): times.append(p)
            if (t  == 0): names.append(p)

    plt.legend(times, ['t={0:.0f} min'.format(cases[0].stat.var['t'][t]/60.) for t in range(0,cases[0].stat.var['t'].size,plot_every)], bbox_to_anchor=(1.05, 1))
    plt.xlabel(xlabel)
    plt.ylabel('height / m')
    #plt.set_yticks([])
    plt.ylim(0,ymax)
    plt.savefig("profile_evolution_"+varname+".pdf", bbox_inches='tight')
    plt.show()

def plot_profile_comparison(case, varnames):
    specs = {
            'rainrate': ['rain rate', '(mm/h)', 3600*1e3],
            'thl': ['liquid water potential temperature', 'K', 1],
            'th': ['potential temperature', 'K', 1],
            'qt': ['total water mixing ratio', '(g/kg)', 1e3],
            'qv': ['water vapor mixing ratio', '(g/kg)', 1e3],
            'qs': ['saturation water vapor mixing ratio', '(g/kg)', 1e3],
            'au_qr': ['autoconversion', '(mg/kg/h)', 1e6*3600],
            'ac_qr': ['accretion', '(mg/kg/h)', 1e6*3600],
            'ev_qr': ['rain evaporation', '(mg/kg/h)', 1e6*3600],
            'ql': ['cloud droplet water mixing ratio', '(g/kg)', 1e3],
            'qr': ['rain drop water mixing ratio', '(mg/kg)', 1e6],
            }

    assert len(set([specs[varname][2] for varname in varnames])) == 1
    # Vertical extent plot
    ymax = case.grid.zsize
   
    plt.figure(figsize=(3,5))
    remove_top_right_axis()
    
    for ic,varname in enumerate(varnames):
        assert varname in specs.keys()
        xlabel, unit, fac = specs[varname]
        p, = plt.plot(case.stat.var[varname][0,:]*fac, case.stat.var['z'], label=xlabel)

    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.xlabel('variable / '+unit) 
    plt.ylabel('height / m')
    #plt.set_yticks([])
    plt.ylim(0,ymax)
    varnamestring = ''
    for varname in varnames:
        varnamestring += varname
    plt.savefig("initial_conditions_"+varnamestring+".pdf", bbox_inches='tight')
    plt.show()


def plot_initial_profiles(case):
    plot_profile_comparison(case, ['qt', 'qv', 'ql'])
    plot_profile_comparison(case, ['thl', 'th'])

def plot_timeseries(cases, varname):
    specs = {
            'cb_rain_rate': ['cloud-base rain rate / (mm/h)', 3600*1e3],
            'surf_rain_rate': ['surface rain rate / (mm/h)', 3600*1e3],
            'au_qr': ['maximum autoconversion rate / (mg/kg/h)', 1e6*3600],
            'ac_qr': ['maximum accretion rate / (mg/kg/h)', 1e6*3600],
            'ev_qr': ['maximum rain evaporation / (mg/kg/h)', 1e6*3600],
            'thl': ['maximum liquid water potential temperature / K', 1],
            'th': ['maximum potential temperature / K', 1],
            'qt': ['maximum total water mixing ratio / (g/kg)', 1e3],
            'qv': ['maximum water vapor mixing ratio / (g/kg)', 1e3],
            'qs': ['maximum saturation water vapor mixing ratio / (g/kg)', 1e3],
            'ql': ['maximum cloud droplet water mixing ratio / (g/kg)', 1e3],
            'qr': ['maximum rain drop water mixing ratio / (mg/kg)', 1e6],
            'processes': ['maximum rate / (mg/kg/h)', 1e6*3600],
            }
    assert varname in specs.keys()
    ylabel, fac = specs[varname]

    plt.figure(figsize=(5,3))
    remove_top_right_axis()
    
    for ic,case in enumerate(cases):
       
        if 'rain_rate' in varname:
            # Find cloud base
            cb_index = np.nonzero(case.stat.var['ql'])[0]
            
            rainrate = np.array([case.stat.var['rainrate'][t, np.nonzero(case.stat.var['ql'][t])[0][0]]*fac for t in range(len(case.stat.var['t']))]) if varname == 'cb_rain_rate' else case.stat.var['rainrate'][:,0]*fac
            p, = plt.plot(case.stat.var['t'] / 60., rainrate, label=case.name + ', tot = %.2f mm' % np.sum(rainrate * case.stat.sampletime/3600.) )
        elif varname == 'processes':
            p, = plt.plot(case.stat.var['t'] / 60., np.amax(case.stat.var['ac_qr'], axis=1)*fac, label=case.name)
            p, = plt.plot(case.stat.var['t'] / 60., np.amax(case.stat.var['au_qr'], axis=1)*fac, label=case.name)
        else:
            p, = plt.plot(case.stat.var['t'] / 60., np.amax(np.abs(case.stat.var[varname]), axis=1)*fac, label=case.name)

    plt.xlabel('time / min')
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.savefig("timeseries_"+varname+".pdf", bbox_inches='tight')
    plt.show()
