import numpy as np

def slope_mihh(a,b):
    # minmod() function
    return np.sign(a) * max(0., min(np.abs(a), np.sign(a)*b))

def slope_palm(qplus, q, qmin):
    d_mean = 0.5 * (qplus - qmin)
    d_min  = q - np.min([qplus, q, qmin])
    d_max  = np.max([qplus, q, qmin]) - q
    slope  = np.sign(d_mean) * np.min([2*d_min, 2*d_max, np.abs(d_mean)])

    return slope

def slope_ucla(qplus, q, qmin):
    dnp    = qplus - q
    dnm    = q - qmin
    sk     = 0.5 * (dnp + dnm)
    mini   = np.min([qmin, q, qplus])
    maxi   = np.max([qmin, q, qplus])
    slope  = 0.5 * np.sign(sk) * np.min([np.abs(sk), 2*(q-mini), 2*(maxi-q)]) 

    return slope

qp = 2.
q  = 1.
qm = 0.5
print('d/dz > 0:\n----')
print('PALM', slope_palm(qp, q, qm))
print('UCLA', slope_ucla(qp, q, qm))
print('MiHH', slope_mihh(q-qm, qp-q))

qp = 0.5
q  = 1.
qm = 2.
print('d/dz < 0:\n----')
print('PALM', slope_palm(qp, q, qm))
print('UCLA', slope_ucla(qp, q, qm))
print('MiHH', slope_mihh(q-qm, qp-q))


# UCLA-LES
# sk = 0.5 * (dn(k-1) + dn(k))
# mini = min(np(k-1),np(k),np(k+1))
# maxi = max(np(k-1),np(k),np(k+1))
# nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np(k)-mini), &
#      &                                     2.*(maxi-np(k)))

# PALM:
# d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
# d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
# d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)
# 
# qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  &
#                                            2.0_wp * d_max,  &
#                                            ABS( d_mean ) )32
