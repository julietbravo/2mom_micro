import numpy as np

class Time:
    def __init__(self, model):
        self.model = model

        self.total_time = 3600
        self.dt_min     = 1
        self.dt_max     = 30
        self.time       = 0
        self.dt         = self.dt_min
        self.cfl_max    = 2.0

        self.iter       = 0
        self.substep    = 0
        self.nsubstep   = 3

    def get_time_limit(self):
        return max(self.dt_min, min(self.total_time - self.time, self.dt_max))
    
    def set_time_step(self):
        dt1 = self.get_time_limit()
        dt2 = self.model.micro.get_time_limit(self.cfl_max, self.dt)
        dt3 = self.model.stat.get_time_limit()
        self.dt = np.array([dt1,dt2,dt3]).min()

        if(self.iter % 10 == 0):
            print('it=%4i [-], dt=%6.2f [s], cfl_max=%6.2f [-]'%(self.iter, self.dt, self.model.micro.max_cfl))

    def finished(self):
        return True if self.time >= self.total_time else False 

    def get_sub_dt(self):
        cB = ([1./3., 15./16., 8./15.])   
        return self.dt * cB[self.substep] 

    def rk3(self):
        cA = ([0., -5./9., -153./128.])
        cB = ([1./3., 15./16., 8./15.])   

        substepn = (self.substep+1)%3

        for pfield in self.model.fields.ap:
            pfield.data[:] += cB[self.substep] * self.dt * pfield.tend[:]
            pfield.tend[:] *= cA[substepn    ]

        self.substep = substepn 
