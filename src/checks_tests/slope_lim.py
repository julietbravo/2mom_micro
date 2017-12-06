import numpy as np
import matplotlib.pylab as pl

pl.close('all')

def phi_minmod(R):
    return max(0, min(1,R))

def phi_vanleer(R):
    return R + np.abs(R) / (1 + np.abs(R))

def phi_superbee(R):
    return np.min([max(1,R), 2, 2*R])


def interpol(vp, v, vm, kind='minmod'):
    # ratio of forward to backward difference:
    ri = (vp - v) / ((v - vm)+1e-16)

    # total slope:
    slope = 0.5 * (vp - vm)

    if (kind == 'minmod'):
        phi = phi_minmod
    if (kind == 'vanleer'):
        phi = phi_vanleer
    if (kind == 'superbee'):
        phi = phi_superbee
    if (kind == 'vanAlbada'):
        phi = phi_vanAlbada

    vph = v + 0.5 * phi(ri) * slope
    vmh = v - 0.5 * phi(1/ri) * slope
    
    pl.figure()
    pl.title(kind)
    pl.plot([vm, v, vp], [0,1,2], '-x', label='original')
    pl.scatter([vmh, vph], [0.5, 1.5], s=30, color='r', label='interpol')
    pl.legend(frameon=False)


interpol(0, 1, 2, 'minmod')
interpol(0, 1, 2, 'vanleer')
interpol(0, 1, 2, 'superbee')

