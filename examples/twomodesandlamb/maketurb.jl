#include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")
#using FourierFlows,
#      PyPlot

import FourierFlows.TwoDTurb

substeps = 100

function matureturb(nx, Lx, dt, nu, nnu, f; kpeak=16, Ro=0.1, Ro0=1.2*Ro)

  #nu = 5e-2/(dt*(0.65*pi*nx/Lx)^nnu)

  g  = TwoDTurb.Grid(nx, Lx)
  pr = TwoDTurb.Params(nu, nnu)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LC)

  q0 = FourierFlows.peaked_isotropic_spectrum(nx, kpeak; maxval=Ro0*f)
  TwoDTurb.set_q!(vs, g, q0)
  maxRo = maximum(abs.(vs.q)/f)

  fig, axs = subplots()

  @printf("Stepping until flow has max(Ro) = %.2f...\n", Ro)
  while maxRo > Ro && ts.step < 10000
    stepforward!(vs, substeps, ts, eq, pr, g)
    TwoDTurb.updatevars!(vs, g)
    maxRo = maximum(abs.(vs.q)/f)

    @printf("step: %d, t*f/Ro0: %.2e, max Ro: %.3f, CFL: %.3f\n", 
      ts.step, vs.t*f/Ro0, maxRo, maximum([vs.U; vs.V])*ts.dt/g.dx)

    imshow(vs.q)
    pause(0.01)
  end


  return vs.q
end

#Z0 = matureturb(nx, Lx, 2*pi*Ro/(5*f0), 4, f0; kpeak=256, Ro=Ro)
