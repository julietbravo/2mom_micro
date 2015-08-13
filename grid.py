import numpy as np

class Grid:
    def __init__(self, zsize, ktot):
        self.zsize  = zsize
        self.ktot   = ktot
        
    def create(self):
        self.kcells = self.ktot + 2
        self.kstart = 1
        self.kend   = self.ktot + 1 

        dz0 = self.zsize / self.ktot

        # Arrays 
        self.z    = np.zeros(self.kcells)
        self.zh   = np.zeros(self.kcells)
        self.dz   = np.zeros(self.kcells)
        self.dzh  = np.zeros(self.kcells)

        # Create uniform grid
        self.z[self.kstart:self.kend] = np.linspace(0.5*dz0, self.zsize-0.5*dz0, self.ktot)
        self.z[0] = -self.z[1]
        self.z[-1] = self.z[-2] + dz0

        # Calculate half levels
        for k in range(self.kstart+1, self.kend):
            self.zh[k] = 0.5 * (self.z[k-1] + self.z[k])
        self.zh[self.kstart] = 0.;
        self.zh[self.kend]   = self.zsize;

        # Calculate grid spacing
        for k in range(self.kstart, self.kcells):
            self.dzh[k] = self.z[k] - self.z[k-1]
        self.dzh[self.kstart-1] = self.dzh[self.kstart+1]

        for k in range(self.kstart, self.kcells-1):
            self.dz[k] = self.zh[k+1] - self.zh[k]
        self.dz[self.kstart-1] = self.dz[self.kstart]
        self.dz[self.kend] = self.dz[self.kend-1]

        self.dzi  = 1./self.dz
        self.dzhi = 1./self.dzh
