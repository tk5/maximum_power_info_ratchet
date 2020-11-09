#!/usr/bin/env python3

from numpy import array, sqrt, pi, exp, zeros
from scipy.special import erf
from scipy.constants import g as g_earth, Boltzmann as kB



# ============= DEFINE PROPERTIES ==============

# define bead properties
diameters = array([0.5, 1.5, 3.0, 5.0])*1e-6
r = diameters*0.5 #radius
mass = array([6.54e-17, 1.73e-15, 1.41e-14, 6.54e-14])
gamma = array([4.20e-9, 1.25e-8, 2.52e-8, 4.19e-8])

# define medium properties
eta = 8.9e-4 # dynamic viscosity of water
rho_f = 1e3 # density of water
kT = kB*(293.0) # assume temperature of 20 degrees celcius
beta = kT**(-1)

# ============= COMPUTE QUANTITIES ==============

# define S^{2} = tau_{f}/tau_{r} and Q^{2} = tau_{v}/tau_{r}
Sval = 1.0
Qval = 1.0

# compute critical k values for Q = Qval and S = Sval
kf = ((Sval**2)*gamma*eta)/((r**2)*rho_f)
kv = ((Qval**2)*(gamma**2))/mass

# compute delta values associated with kf and kv
delta_s = mass*g_earth*sqrt(beta/kf)
delta_q = mass*g_earth*sqrt(beta/kv)

# compute the fluid and velocity relaxation time scales
tau_f = (r**2)*rho_f/eta
tau_v = mass/gamma

# ============= DISPLAY QUANTITIES ==============

# print diameters in units of um
print("d:", diameters/1e-6)

# print time scales in units of us
print("tau_f: ", tau_f/1e-6)
print("tau_v: ", tau_v/1e-6)

# print the critical values of kappa in units of pN/um
print("k_f:", kf/1e-6)

# print velocity and power associated with critical kappa
print("v_s:", v_s := sqrt(2.0/pi)*exp(-0.5*(delta_s**2))*(1.0/sqrt(beta*kf)) / ((gamma/kf)*(1.0+erf(delta_s/sqrt(2.0)))) / 1e-6)
print("P_s:", pows_s := sqrt(2.0/pi)*delta_s*exp(-0.5*(delta_s**2)) / ((gamma/kf)*(1.0+erf(delta_s/sqrt(2.0)))))

# print("k_v:", kv/1e-6)
# print("v_q:", v_q := sqrt(2.0/pi)*exp(-0.5*(delta_q**2))*(1.0/sqrt(beta*kv)) / ((gamma/kv)*(1.0+erf(delta_q/sqrt(2.0)))) / 1e-6)
# print("P_q:", pows_q := sqrt(2.0/pi)*delta_q*exp(-0.5*(delta_q**2)) / ((gamma/kv)*(1.0+erf(delta_q/sqrt(2.0)))))
