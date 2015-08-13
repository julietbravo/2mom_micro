kappa = 0.4        # von Karman constant
grav  = 9.81       # Gravitational acceleration [m s-2]
Rd    = 287.04     # Gas constant for dry air [J K-1 kg-1] 
Rv    = 461.5      # Gas constant for water vapor [J K-1 kg-1]
cp    = 1005       # Specific heat of air at constant pressure [J kg-1 K-1]
Lv    = 2.5e6      # Latent heat of condensation or vaporization [J kg-1]
T0    = 273.15     # Freezing / melting temperature [K]
p0    = 1.e5       # Reference pressure [pa]

# Short cuts
ep    = Rd/Rv
rdcp  = Rd/cp
small = 1e-16

# Coefficients esat
c0    = 0.6105851e+03 
c1    = 0.4440316e+02 
c2    = 0.1430341e+01 
c3    = 0.2641412e-01 
c4    = 0.2995057e-03 
c5    = 0.2031998e-05 
c6    = 0.6936113e-08 
c7    = 0.2564861e-11 
c8    = -.3704404e-13 
