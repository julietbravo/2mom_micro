import numpy as np
from constants import *

def exner(p):
    return (p/p0)**(Rd/cp);

def esat(T):
    if(np.size(T) > 1):
        x = T - T0
        x[x<-80] = -80
    else:
        x = np.max((-80.,T - T0));
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));

def qsat(p, T):
    return ep*esat(T)/(p-(1-ep)*esat(T));

def sat_adjust(thl, qt, p, exn):
    niter    = 0
    nitermax = 30
    tnr_old  = 1.e9
    tl       = thl * exn
    tnr      = tl
    while (abs(tnr-tnr_old)/tnr_old> 1e-5 and niter < nitermax):
        niter  += 1
        tnr_old = tnr;
        qs      = qsat(p, tnr);
        tnr     = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(pow(Lv,2)*qs)/ (Rv*cp*pow(tnr,2)));
    ql = np.max((0.,qt - qs))
    return (ql, qs)

class Thermo:
    def __init__(self, model):
        self.model = model

        self.pbot  = 1e5 

    def execute(self):
        # Some shortcuts
        model  = self.model
        kstart = model.grid.kstart
        kend   = model.grid.kend
        thl    = model.fields.thl.data
        qt     = model.fields.qt.data

        ql     = model.fields.ql .data
        qs     = model.fields.qs .data
        p      = model.fields.p  .data
        exn    = model.fields.exn.data
        rho    = model.fields.rho.data
        thv    = model.fields.thv.data

        # Some tmp arrays 
        exh   = np.zeros(model.grid.kcells)
        prefh = np.zeros(model.grid.kcells)
        rhoh  = np.zeros(model.grid.kcells)
        thvh  = np.zeros(model.grid.kcells)
        exh   = np.zeros(model.grid.kcells)

        thlsurf       = 0.5*(thl[kstart-1] + thl[kstart]);
        qtsurf        = 0.5*(qt [kstart-1] + qt [kstart]);
        exh[kstart]   = exner(self.pbot);
        qltmp, qstmp  = sat_adjust(thlsurf, qtsurf, self.pbot, exh[kstart]); 
        thvh[kstart]  = (thlsurf + Lv*qltmp/(cp*exh[kstart])) * (1. - (1. - Rv/Rd)*qtsurf - Rv/Rd*qltmp);
        prefh[kstart] = self.pbot;
        rhoh[kstart]  = self.pbot / (Rd * exh[kstart] * thvh[kstart]);
    
        # First full grid level pressure
        p[kstart]  = pow((pow(self.pbot, rdcp) - grav * pow(p0, rdcp) * model.grid.z[kstart] / (cp * thvh[kstart])),(1./rdcp)); 
    
        for k in range(kstart+1, kend+1):
            # 1. Calculate values at full level below zh[k] 
            exn[k-1]          = exner(p[k-1]);
            ql[k-1], qs[k-1]  = sat_adjust(thl[k-1], qt[k-1], p[k-1], exn[k-1]); 
            thv[k-1]          = (thl[k-1] + Lv*ql[k-1]/(cp*exn[k-1])) * (1. - (1. - Rv/Rd)*qt[k-1] - Rv/Rd*ql[k-1]); 
            rho[k-1]          = p[k-1] / (Rd * exn[k-1] * thv[k-1]);
    
            # 2. Calculate half level pressure at zh[k] using values at z[k-1]
            prefh[k] = pow((pow(prefh[k-1],rdcp) - grav * pow(p0,rdcp) * model.grid.dz[k-1] / (cp * thv[k-1])),(1./rdcp));
    
            # 3. Interpolate conserved variables to half level
            thli     = 0.5 * (thl[k-1] + thl[k]);
            qti      = 0.5 * (qt[k-1]  + qt[k] );
    
            # 4. Calculate half level values
            exh[k]   = exner(prefh[k]);
            qli, qsi = sat_adjust(thli, qti, prefh[k], exh[k]);
            thvh[k]  = (thli + Lv*qli/(cp*exh[k])) * (1. - (1. - Rv/Rd)*qti - Rv/Rd*qli); 
            rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]); 
    
            # 5. Calculate full level pressure at z[k]
            p[k]  = pow((pow(p[k-1],rdcp) - grav * pow(p0,rdcp) * model.grid.dzh[k] / (cp * thvh[k])),(1./rdcp)); 
    
        # Fill bottom and top full level ghost cells 
        p[kstart-1]   = 2.*prefh[kstart] - p  [kstart];
        p[kend]       = 2.*prefh[kend]   - p  [kend-1];
        exn[kstart-1] = 2.*exh[kstart]   - exn[kstart];
        exn[kend]     = 2.*exh[kend]     - exn[kend-1];

