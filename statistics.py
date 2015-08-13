import numpy as np

class Statistics:
    def __init__(self, model):
        self.model = model
        self.grid   = model.grid
        self.fields = model.fields
        self.time   = model.time

        self.sampletime = 120

        ## Empty lists to store data or numpy arrays
        self.t     = np.zeros( 0                 )
        self.thl   = np.zeros((0, self.grid.ktot))
        self.qt    = np.zeros((0, self.grid.ktot))
        self.ql    = np.zeros((0, self.grid.ktot))
        self.qr    = np.zeros((0, self.grid.ktot))
        self.nr    = np.zeros((0, self.grid.ktot))

        # Sedimentation velocities
        self.w_qr  = np.zeros((0, self.grid.ktot))
        self.w_nr  = np.zeros((0, self.grid.ktot))

        # Autoconversion
        self.au_qr = np.zeros((0, self.grid.ktot))
        self.au_nr = np.zeros((0, self.grid.ktot))
        self.au_qt = np.zeros((0, self.grid.ktot))
        self.au_th = np.zeros((0, self.grid.ktot))

        # Evaporation
        self.ev_qr = np.zeros((0, self.grid.ktot))
        self.ev_nr = np.zeros((0, self.grid.ktot))
        self.ev_qt = np.zeros((0, self.grid.ktot))
        self.ev_th = np.zeros((0, self.grid.ktot))

        # Accretion
        self.ac_qr = np.zeros((0, self.grid.ktot))
        self.ac_qt = np.zeros((0, self.grid.ktot))
        self.ac_th = np.zeros((0, self.grid.ktot))

        # Self-collection and breaukup
        self.sb_nr = np.zeros((0, self.grid.ktot))

        # Sedimentation
        self.se_qr = np.zeros((0, self.grid.ktot))
        self.se_nr = np.zeros((0, self.grid.ktot))

        # grid without ghost cells
        self.z   = self.grid.z[self.grid.kstart:self.grid.kend]

        # Tmp fields for tendencies
        self.tmp1 = np.zeros(self.grid.kcells)
        self.tmp2 = np.zeros(self.grid.kcells)
        self.tmp3 = np.zeros(self.grid.kcells)
        self.tmp4 = np.zeros(self.grid.kcells)
        
    def get_time_limit(self):
        return self.sampletime - np.mod(self.time.time, self.sampletime)

    def execute(self):
        if(np.mod(self.time.time, self.sampletime) == 0):
            print('saving stats for t=%7.1f'%self.time.time)

            kstart   = self.grid.kstart
            kend     = self.grid.kend
            self.t   = np.append(self.t, self.time.time) 

            # Mean fields --------------------------------------------------------
            self.thl = np.append(self.thl, self.fields.thl.data[kstart:kend]) 
            self.qt  = np.append(self.qt,  self.fields.qt. data[kstart:kend]) 
            self.ql  = np.append(self.ql,  self.fields.ql. data[kstart:kend]) 
            self.qr  = np.append(self.qr,  self.fields.qr. data[kstart:kend]) 
            self.nr  = np.append(self.nr,  self.fields.nr. data[kstart:kend]) 

            # Sedimentation velocity ---------------------------------------------
            self.tmp1[:] = 0; self.tmp2[:] = 0

            self.model.micro.sedimentation_velocity(self.tmp1, self.tmp2, \
                                                    self.fields.qr.data, self.fields.nr.data, 
                                                    self.fields.rho.data, self.grid.kstart, self.grid.kend)

            self.w_qr = np.append(self.w_qr, self.tmp1[kstart:kend])
            self.w_nr = np.append(self.w_nr, self.tmp2[kstart:kend])

            # Autoconversion -----------------------------------------------------
            self.tmp1[:] = 0; self.tmp2[:] = 0; self.tmp3[:] = 0; self.tmp4[:] = 0

            self.model.micro.autoconversion(self.tmp1, self.tmp2, self.tmp3, self.tmp4, \
                                            self.fields.qr.data, self.fields.nr.data, self.fields.ql.data, self.fields.rho.data, \
                                            self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.au_qr = np.append(self.au_qr, self.tmp1[kstart:kend])
            self.au_nr = np.append(self.au_nr, self.tmp2[kstart:kend])
            self.au_th = np.append(self.au_th, self.tmp3[kstart:kend])
            self.au_qt = np.append(self.au_qt, self.tmp4[kstart:kend])

            # Evaporation --------------------------------------------------------
            self.tmp1[:] = 0; self.tmp2[:] = 0; self.tmp3[:] = 0; self.tmp4[:] = 0

            self.model.micro.evaporation(self.tmp1, self.tmp2, self.tmp3, self.tmp4, \
                                         self.fields.qr.data, self.fields.nr.data, self.fields.qt.data, self.fields.qs.data, \
                                         self.fields.ql.data, self.fields.thl.data*self.fields.exn.data, self.fields.rho.data, \
                                         self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.ev_qr = np.append(self.ev_qr, self.tmp1[kstart:kend])
            self.ev_nr = np.append(self.ev_nr, self.tmp2[kstart:kend])
            self.ev_th = np.append(self.ev_th, self.tmp3[kstart:kend])
            self.ev_qt = np.append(self.ev_qt, self.tmp4[kstart:kend])

            # Accretion ----------------------------------------------------------
            self.tmp1[:] = 0; self.tmp2[:] = 0; self.tmp3[:] = 0; self.tmp4[:] = 0

            self.model.micro.accretion(self.tmp1, self.tmp2, self.tmp3, self.tmp4, \
                                       self.fields.qr.data, self.fields.nr.data, self.fields.ql.data, \
                                       self.fields.rho.data, self.fields.exn.data, self.grid.kstart, self.grid.kend)

            self.ac_qr = np.append(self.ac_qr, self.tmp1[kstart:kend])
            self.ac_th = np.append(self.ac_th, self.tmp3[kstart:kend])
            self.ac_qt = np.append(self.ac_qt, self.tmp4[kstart:kend])

            # Self-collection and breakup ----------------------------------------
            self.tmp1[:] = 0; self.tmp2[:] = 0; self.tmp3[:] = 0; self.tmp4[:] = 0

            self.model.micro.selfcollection_breakup(self.tmp1, self.tmp2, self.tmp3, self.tmp4, \
                                                self.fields.qr.data, self.fields.nr.data, \
                                                self.fields.rho.data, self.grid.kstart, self.grid.kend)
            
            self.sb_nr = np.append(self.sb_nr, self.tmp2[kstart:kend])

            # Sedimentation ------------------------------------------------------ 
            self.tmp1[:] = 0; self.tmp2[:] = 0; self.tmp3[:] = 0; self.tmp4[:] = 0

            self.model.micro.sedimentation_ss08(self.tmp1, self.tmp2, \
                                                self.fields.qr.data, self.fields.nr.data, \
                                                self.fields.rho.data, self.grid.dz, self.grid.dzi, self.grid.dzhi, self.time.dt, \
                                                self.time.dt, self.grid.kstart, self.grid.kend, self.grid.kcells)
            
            self.se_qr = np.append(self.se_qr, self.tmp1[kstart:kend])
            self.se_nr = np.append(self.se_nr, self.tmp2[kstart:kend])

    def finish(self):
        nt = self.t.size
        nz = self.grid.ktot

        self.thl   = self.thl.reshape([nt, nz])
        self.qt    = self.qt.reshape([nt, nz])
        self.ql    = self.ql.reshape([nt, nz])
        self.qr    = self.qr.reshape([nt, nz])
        self.nr    = self.nr.reshape([nt, nz])

        self.w_qr  = self.w_qr.reshape([nt, nz])
        self.w_nr  = self.w_nr.reshape([nt, nz])

        self.au_qr = self.au_qr.reshape([nt, nz])
        self.au_nr = self.au_nr.reshape([nt, nz])
        self.au_th = self.au_th.reshape([nt, nz])
        self.au_qt = self.au_qt.reshape([nt, nz])

        self.ev_qr = self.ev_qr.reshape([nt, nz])
        self.ev_nr = self.ev_nr.reshape([nt, nz])
        self.ev_th = self.ev_th.reshape([nt, nz])
        self.ev_qt = self.ev_qt.reshape([nt, nz])

        self.ac_qr = self.ac_qr.reshape([nt, nz])
        self.ac_th = self.ac_th.reshape([nt, nz])
        self.ac_qt = self.ac_qt.reshape([nt, nz])

        self.sb_nr = self.sb_nr.reshape([nt, nz])

        self.se_qr = self.se_qr.reshape([nt, nz])
        self.se_nr = self.se_nr.reshape([nt, nz])
