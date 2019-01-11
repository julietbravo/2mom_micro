import numpy as np

class Field_1d_d:
    def __init__(self, size):
        self.data = np.zeros(size)

class Field_1d_p:
    def __init__(self, size):
        self.data = np.zeros(size)
        self.tend = np.zeros(size)

class Fields:
    def __init__(self, model):
        self.ap = []
        self.init_prognostic("thl", model.grid.kcells)
        self.init_prognostic("qt",  model.grid.kcells)
        self.init_prognostic("qr",  model.grid.kcells)
        self.init_prognostic("nr",  model.grid.kcells)

        self.ad =[]
        self.init_diagnostic("ql",  model.grid.kcells)  
        self.init_diagnostic("qv",  model.grid.kcells)  
        self.init_diagnostic("qs",  model.grid.kcells)  
        self.init_diagnostic("p" ,  model.grid.kcells)   
        self.init_diagnostic("exn", model.grid.kcells) 
        self.init_diagnostic("rho", model.grid.kcells) 
        self.init_diagnostic("thv", model.grid.kcells) 
        self.init_diagnostic("th", model.grid.kcells) 

    def init_prognostic(self, name, kcells):
        setattr(self, name, Field_1d_p(kcells)) 
        self.ap.append(getattr(self, name))

    def init_diagnostic(self, name, kcells):
        setattr(self, name, Field_1d_d(kcells)) 
        self.ad.append(getattr(self, name))

    def clip_at_zero(self):
        for f in self.ap:
            f.data[f.data<0] = 0.

    def set_bc(self):
        for f in self.ap:
            f.data[0] = f.data[1]
            f.data[-1] = f.data[-2]
