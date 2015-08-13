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
        self.model = model

        self.ap = []
        self.init_prognostic("thl")
        self.init_prognostic("qt")
        self.init_prognostic("qr")
        self.init_prognostic("nr")

        self.ad =[]
        self.init_diagnostic("ql")  
        self.init_diagnostic("qs")  
        self.init_diagnostic("p")   
        self.init_diagnostic("exn") 
        self.init_diagnostic("rho") 
        self.init_diagnostic("thv") 

    def init_prognostic(self, name):
        setattr(self, name, Field_1d_p(self.model.grid.kcells)) 
        self.ap.append(getattr(self, name))

    def init_diagnostic(self, name):
        setattr(self, name, Field_1d_d(self.model.grid.kcells)) 
        self.ad.append(getattr(self, name))

    def clip_at_zero(self):
        for f in self.ap:
            f.data[f.data<0] = 0.
