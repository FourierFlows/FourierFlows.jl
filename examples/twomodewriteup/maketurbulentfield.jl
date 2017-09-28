include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot,
      JLD

import FourierFlows.TwoDTurb,
       FourierFlows.TwoModeBoussinesq

# Rocha Wagner Young parameters
nx   = 512                                 # Resolution
Lx   = 2*pi*200e3                          # Domain extent
Ue   = 5e-2
ke   = 20*pi/Lx
f0   = 1e-4

# Numerical parameters
dt   = 2e-2 * 2*pi/f0
nnu  = 8                                   # Hyperviscous order
nu0 = nu1 = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)  # Hyperviscosity
name = @sprintf("macroturb_niw")


#Z0 = TwoDTurb.makematureturb(nx, Lx; nu=nu0, nnu=nnu,
#  k0=ke*Lx/(2*pi), E0=0.5*Ue^2, tf=200/(Ue*ke), plots=true)

Z0 = TwoDTurb.makematureturb(nx, Lx; nu=nu0, nnu=nnu,
  k0=256, qf=0.1*f0, q0=0.12*f0, tf = 1000/f0, plots=true)

fig, axs = subplots()
imshow(Z0)
show()

save("testturb.jld", "Z", Z0)
