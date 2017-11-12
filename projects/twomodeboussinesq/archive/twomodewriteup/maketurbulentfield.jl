include("../../src/fourierflows.jl")

using  FourierFlows, PyPlot, JLD2
import FourierFlows.TwoDTurb


nx, Lx, f = 512,   2*pi*1600e3, 1e-4
nnu, dt   = 4,     0.02*2pi/f
qf, k0    = 0.1*f, 192

nu = 2e-2/(dt*(0.65*pi*nx/Lx)^nnu) 

Z = TwoDTurb.makematureturb(nx, Lx; nu=nu, nnu=nnu,
  k0=k0, qf=qf, q0=1.2*qf, tf=80/qf, plots=true)

savename = @sprintf("twodturb_nx%04d_nnu%d.jld2", nx, nnu)
plotname = @sprintf("twodturb_nx%04d_nnu%d.png", nx, nnu)
titlname = @sprintf("twodturb: \$n=%d\$, \$n_{\\nu}=%d\$, max(Z)=%.2e", 
  nx, nnu, maximum(Z))

fig, axs = subplots()
imshow(Z)
title(titlname)
savefig(plotname, dpi=240)

@save savename Z
