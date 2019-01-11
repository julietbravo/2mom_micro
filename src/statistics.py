import numpy as np

class Statistics:
    def __init__(self, model):
        self.model  = model
        self.grid   = model.grid
        self.fields = model.fields
        self.time   = model.time

        self.sampletime = 120

        self.var = {}
        ## Empty lists to store data or numpy arrays
        self.var['t']     = np.zeros( 0                 )
        self.var['thl']   = np.zeros((0, self.grid.ktot))
        self.var['th']   = np.zeros((0, self.grid.ktot))
        self.var['qt']    = np.zeros((0, self.grid.ktot))
        self.var['qs']    = np.zeros((0, self.grid.ktot))
        self.var['qv']    = np.zeros((0, self.grid.ktot))
        self.var['ql']    = np.zeros((0, self.grid.ktot))
        self.var['qr']    = np.zeros((0, self.grid.ktot))
        self.var['nr']    = np.zeros((0, self.grid.ktot))
        
        self.var['rho']   = np.zeros((0, self.grid.ktot))

        # Sedimentation velocities
        self.var['w_qr']  = np.zeros((0, self.grid.ktot))
        self.var['w_nr']  = np.zeros((0, self.grid.ktot))

        # Autoconversion
        self.var['au_qr'] = np.zeros((0, self.grid.ktot))
        self.var['au_nr'] = np.zeros((0, self.grid.ktot))
        self.var['au_qt'] = np.zeros((0, self.grid.ktot))
        self.var['au_th'] = np.zeros((0, self.grid.ktot))

        # Evaporation
        self.var['ev_qr'] = np.zeros((0, self.grid.ktot))
        self.var['ev_nr'] = np.zeros((0, self.grid.ktot))
        self.var['ev_qt'] = np.zeros((0, self.grid.ktot))
        self.var['ev_th'] = np.zeros((0, self.grid.ktot))

        # Accretion
        self.var['ac_qr'] = np.zeros((0, self.grid.ktot))
        self.var['ac_qt'] = np.zeros((0, self.grid.ktot))
        self.var['ac_th'] = np.zeros((0, self.grid.ktot))

        # Self-collection and breaukup
        self.var['sb_nr'] = np.zeros((0, self.grid.ktot))

        # Sedimentation
        self.var['se_qr'] = np.zeros((0, self.grid.ktot))
        self.var['se_nr'] = np.zeros((0, self.grid.ktot))
        self.var['rainrate'] = np.zeros((0, self.grid.ktot))

        # grid without ghost cells
        self.var['z']   = self.grid.z[self.grid.kstart:self.grid.kend]

        # Tmp fields for tendencies
        self.var['tmp1'] = np.zeros(self.grid.kcells)
        self.var['tmp2'] = np.zeros(self.grid.kcells)
        self.var['tmp3'] = np.zeros(self.grid.kcells)
        self.var['tmp4'] = np.zeros(self.grid.kcells)
        
    def get_time_limit(self):
        return self.sampletime - np.mod(self.time.time, self.sampletime)

    def execute(self, verbose=False):
        if(np.mod(self.time.time, self.sampletime) == 0):
            if verbose:
                print('saving stats for t=%7.1f'%self.time.time)

            kstart   = self.grid.kstart
            kend     = self.grid.kend
            self.var['t']   = np.append(self.var['t'], self.time.time) 

            # Mean fields --------------------------------------------------------
            self.var['thl'] = np.append(self.var['thl'], self.fields.thl.data[kstart:kend]) 
            self.var['th'] = np.append(self.var['th'], self.fields.th.data[kstart:kend]) 
            self.var['qt']  = np.append(self.var['qt'],  self.fields.qt. data[kstart:kend]) 
            self.var['qv']  = np.append(self.var['qv'],  self.fields.qv. data[kstart:kend]) 
            self.var['qs']  = np.append(self.var['qs'],  self.fields.qs. data[kstart:kend]) 
            self.var['ql']  = np.append(self.var['ql'],  self.fields.ql. data[kstart:kend]) 
            self.var['qr']  = np.append(self.var['qr'],  self.fields.qr. data[kstart:kend]) 
            self.var['nr']  = np.append(self.var['nr'],  self.fields.nr. data[kstart:kend]) 

            self.var['rho'] = np.append(self.var['rho'], self.fields.rho.data[kstart:kend]) 

            # Sedimentation velocity ---------------------------------------------
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0

            self.model.micro.sedimentation_velocity(self.var['tmp1'], self.var['tmp2'], \
                                                    self.fields.qr.data, self.fields.nr.data, 
                                                    self.fields.rho.data, self.grid.kstart, self.grid.kend)

            self.var['w_qr'] = np.append(self.var['w_qr'], self.var['tmp1'][kstart:kend])
            self.var['w_nr'] = np.append(self.var['w_nr'], self.var['tmp2'][kstart:kend])

            # Autoconversion -----------------------------------------------------
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0; self.var['tmp3'][:] = 0; self.var['tmp4'][:] = 0

            self.model.micro.autoconversion(self.var['tmp1'], self.var['tmp2'], self.var['tmp3'], self.var['tmp4'], \
                                            self.fields.qr.data, self.fields.nr.data, self.fields.ql.data, self.fields.rho.data, \
                                            self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.var['au_qr'] = np.append(self.var['au_qr'], self.var['tmp1'][kstart:kend])
            self.var['au_nr'] = np.append(self.var['au_nr'], self.var['tmp2'][kstart:kend])
            self.var['au_th'] = np.append(self.var['au_th'], self.var['tmp3'][kstart:kend])
            self.var['au_qt'] = np.append(self.var['au_qt'], self.var['tmp4'][kstart:kend])

            # Evaporation --------------------------------------------------------
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0; self.var['tmp3'][:] = 0; self.var['tmp4'][:] = 0

            self.model.micro.evaporation(self.var['tmp1'], self.var['tmp2'], self.var['tmp3'], self.var['tmp4'], \
                                         self.fields.qr.data, self.fields.nr.data, self.fields.qt.data, self.fields.qs.data, \
                                         self.fields.ql.data, self.fields.thl.data*self.fields.exn.data, self.fields.rho.data, \
                                         self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.var['ev_qr'] = np.append(self.var['ev_qr'], self.var['tmp1'][kstart:kend])
            self.var['ev_nr'] = np.append(self.var['ev_nr'], self.var['tmp2'][kstart:kend])
            self.var['ev_th'] = np.append(self.var['ev_th'], self.var['tmp3'][kstart:kend])
            self.var['ev_qt'] = np.append(self.var['ev_qt'], self.var['tmp4'][kstart:kend])

            # Accretion ----------------------------------------------------------
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0; self.var['tmp3'][:] = 0; self.var['tmp4'][:] = 0

            self.model.micro.accretion(self.var['tmp1'], self.var['tmp2'], self.var['tmp3'], self.var['tmp4'], \
                                       self.fields.qr.data, self.fields.nr.data, self.fields.ql.data, \
                                       self.fields.rho.data, self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.var['ac_qr'] = np.append(self.var['ac_qr'], self.var['tmp1'][kstart:kend])
            self.var['ac_th'] = np.append(self.var['ac_th'], self.var['tmp3'][kstart:kend])
            self.var['ac_qt'] = np.append(self.var['ac_qt'], self.var['tmp4'][kstart:kend])

            # Self-collection and breakup ----------------------------------------
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0; self.var['tmp3'][:] = 0; self.var['tmp4'][:] = 0

            self.model.micro.selfcollection_breakup(self.var['tmp1'], self.var['tmp2'], self.var['tmp3'], self.var['tmp4'], \
                                                self.fields.qr.data, self.fields.nr.data, \
                                                self.fields.rho.data, self.grid.kstart, self.grid.kend)
            
            self.var['sb_nr'] = np.append(self.var['sb_nr'], self.var['tmp2'][kstart:kend])

            # Sedimentation ------------------------------------------------------ 
            self.var['tmp1'][:] = 0; self.var['tmp2'][:] = 0; self.var['tmp3'][:] = 0; self.var['tmp4'][:] = 0

            self.model.micro.sedimentation_ss08(self.var['tmp1'], self.var['tmp2'], \
                                                self.fields.qr.data, self.fields.nr.data, \
                                                self.fields.rho.data, self.grid.dz, self.grid.dzi, self.grid.dzhi, self.time.dt, \
                                                self.time.dt, self.grid.kstart, self.grid.kend, self.grid.kcells, rainrate=self.var['tmp3'])
            
            self.var['se_qr'] = np.append(self.var['se_qr'], self.var['tmp1'][kstart:kend])
            self.var['se_nr'] = np.append(self.var['se_nr'], self.var['tmp2'][kstart:kend])
            self.var['rainrate'] = np.append(self.var['rainrate'], self.var['tmp3'][kstart:kend])

    def finish(self):
        nt = self.var['t'].size
        nz = self.grid.ktot

        self.var['thl']   = self.var['thl'].reshape([nt, nz])
        self.var['th']   = self.var['th'].reshape([nt, nz])
        self.var['qt']    = self.var['qt'].reshape([nt, nz])
        self.var['qv']    = self.var['qv'].reshape([nt, nz])
        self.var['qs']    = self.var['qs'].reshape([nt, nz])
        self.var['ql']    = self.var['ql'].reshape([nt, nz])
        self.var['qr']    = self.var['qr'].reshape([nt, nz])
        self.var['nr']    = self.var['nr'].reshape([nt, nz])

        self.var['rho']   = self.var['rho'].reshape([nt, nz])

        self.var['w_qr']  = self.var['w_qr'].reshape([nt, nz])
        self.var['w_nr']  = self.var['w_nr'].reshape([nt, nz])

        self.var['au_qr'] = self.var['au_qr'].reshape([nt, nz])
        self.var['au_nr'] = self.var['au_nr'].reshape([nt, nz])
        self.var['au_th'] = self.var['au_th'].reshape([nt, nz])
        self.var['au_qt'] = self.var['au_qt'].reshape([nt, nz])

        self.var['ev_qr'] = self.var['ev_qr'].reshape([nt, nz])
        self.var['ev_nr'] = self.var['ev_nr'].reshape([nt, nz])
        self.var['ev_th'] = self.var['ev_th'].reshape([nt, nz])
        self.var['ev_qt'] = self.var['ev_qt'].reshape([nt, nz])

        self.var['ac_qr'] = self.var['ac_qr'].reshape([nt, nz])
        self.var['ac_th'] = self.var['ac_th'].reshape([nt, nz])
        self.var['ac_qt'] = self.var['ac_qt'].reshape([nt, nz])

        self.var['sb_nr'] = self.var['sb_nr'].reshape([nt, nz])

        self.var['se_qr'] = self.var['se_qr'].reshape([nt, nz])
        self.var['se_nr'] = self.var['se_nr'].reshape([nt, nz])
        self.var['rainrate'] = self.var['rainrate'].reshape([nt, nz])
