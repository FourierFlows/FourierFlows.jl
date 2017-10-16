include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot,
      JLD

import FourierFlows.TwoDTurb,
       FourierFlows.TwoModeBoussinesq

# Rocha Wagner Young parameters
#nx   = 1024                                 # Resolution
#nnu  = 8                                    # Hyperviscous order

Lx   = 2*pi*200e3                           # Domain extent
Ue   = 5e-2
ke   = 20*pi/Lx
f0   = 1e-4

for nx in [2048, 4096]
  for nnu in [2, 4, 8]

    # Numerical parameters
    dt   = 2e-2 * 2*pi/f0
    nu0  = nu1 = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)  # Hyperviscosity
    name = @sprintf("macroturb_niw")

    qf   = 0.1*f0                               # Final vorticity 

    k0 = 256 #max(nx/4, 256)

    Z0 = TwoDTurb.makematureturb(nx, Lx; nu=nu0, nnu=nnu,
      k0=k0, qf=qf, q0=1.2*qf, tf=20/qf, plots=true)

    maxz = maximum(Z0) 

    savename = @sprintf("twodturb_nx%04d_nnu%d.jld", nx, nnu)
    plotname = @sprintf("twodturb_nx%04d_nnu%d.png", nx, nnu)
    titlname = @sprintf("twodturb: \$n=%d\$, \$n_{\\nu}=%d\$, max(Z)=%.2e", 
      nx, nnu, maxz)

    fig, axs = subplots()
    imshow(Z0)
    title(titlname)
    savefig(plotname, dpi=240)

    save(savename, "Z", Z0)

  end
end
