#__precompile__()

include("../../src/fourierflows.jl")

using PyPlot, PyCall, NullableArrays
using FourierFlows 

import FourierFlows, 
       FourierFlows.TracerAdvDiff

import FourierFlows: Cy2, Diagnostic, unpack


@pyimport warnings
warnings.filterwarnings("ignore")




@pyimport numpy.ma as ma
PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)




""" Initialize a tracer problem with constant diffusivity. """
function ConstDiffTracerProblem(g::AbstractGrid, kap::Real, 
  u::Function, w::Function, dt::Real)

  vs = TracerAdvDiff.Vars(g)
  pr = TracerAdvDiff.ConstDiffParams(kap, u, w)
  eq = TracerAdvDiff.Equation(pr, g)
  ts = FourierFlows.ETDRK4TimeStepper(dt, eq.LC)

  Problem(g, vs, pr, eq, ts)
end


function ConstDiffTracerProblem(g::AbstractGrid, eta::Real, kap::Real,
  u::Function, w::Function, dt::Real)

  vs = TracerAdvDiff.Vars(g)
  pr = TracerAdvDiff.ConstDiffParams(eta, kap, u, w)
  eq = TracerAdvDiff.Equation(pr, g)
  ts = FourierFlows.ETDRK4TimeStepper(dt, eq.LC)

  Problem(g, vs, pr, eq, ts)
end







""" Make a nice looking tracer plot. """
function makeplot!(axs, prob, tfinal, yvar, mask)

  vs, pr, g = unpack(prob)

  
  cplot = NullableArray(vs.c, mask)

  t_theory = linspace(0.0, tfinal, 100)

  yvar_theory = 2.0*kap*t_theory + yvar.data[1]

  axs[1][:cla]()
  axs[2][:cla]()

  axes(axs[1])
  pcolormesh(g.X, g.Y, PyObject(cplot))
  xlabel(L"$x$"); ylabel(L"$z$")

  axes(axs[2])
  plot(t_theory, yvar_theory;    color="C0", linestyle="-")
  plot(yvar.time, yvar.data; color="C1", linestyle="--")
  xlabel(L"$t$")
  ylabel(L"C_2 = \int (y-\bar y)^2 c \, \mathrm{d}y \, \mathrm{d}x")
    
  xlim(0.0, tfinal)

  # Formatting
  axs[1][:xaxis][:set_ticks_position]("top")
  axs[1][:xaxis][:set_label_position]("top")

  nothing
end




function writemessage(ep, kap, nsteps, prob, yvar)

  yvar_theory = 2.0*kap*yvar.time[end] + yvar.data[1]

  @printf(
  "(ep: %.1f) fraction complete: %.3f, variance/variance(flat): %.4f\n",
  ep, prob.step/nsteps, yvar.value/yvar_theory)

  nothing
end



function wavy_barotropic_tracer_run(ep, kap; nx=256, dcx=0.05, dcy=0.05,
  xi=0.0, yi=-0.5, Lx=2*pi, Ly=1.0, savefreq=40, U=1.0, CFL=0.2)
       
  # Names for plotting and saving
  plotname = @sprintf("./plots/wavychannel_ep%d_nx%04d", Int(ep*10), nx)
  filename = @sprintf("wavychannel_ep%03d_nx%04d.jld2", Int(ep*100), nx)

  # Parameters
  k1       = 2*pi/Lx      
  del      = Lx/nx
  umax     = U/(Ly*(1-ep))
  tfinal   = Ly*Lx/U 

  roughdt  = CFL*del/umax
  roughnsteps = tfinal/roughdt

  substeps = ceil(Int, roughnsteps/savefreq) # Number of snapshots and diags
  nsteps   = substeps*savefreq
  dt       = tfinal/nsteps


  # Initial condition
  ci(x, y) = (exp( -(x-xi)^2.0/(2.0*dcx^2.0) - (y-yi)^2.0/(2.0*dcy^2) )
      / (2.0*pi*dcx*dcy))


  # Sinsuidal topography and barotropic velocity
  h(x) = Ly*(1.0 - ep*sin(k1*x))           # Topography
  hx(x) = -Ly*k1*ep*cos(k1*x)              # Topographic x-gradient
  ub(x, y, t) = U/h(x)                     # Barotropic x-velocity
  vb(x, y, t) = y*ub(x, y, t)*hx(x)/h(x)   # Topography-following z-velocity


  # Initialize problem
  g = FourierFlows.TwoDGrid(nx, Lx, nx, Ly+ep; y0=-Ly-ep)
  prob = ConstDiffTracerProblem(g, kap, ub, vb, CFL*del/umax)
  TracerAdvDiff.set_c!(prob, ci)


  # Diagnostics
  calc_yvar(prob) = Cy2(prob.vars.c, prob.grid)
  yvar = Diagnostic(calc_yvar, prob)
  diags = [yvar]


  # Output
  get_c(prob) = prob.vars.c
  c_out = Output("c", get_c, prob, filename)
  outputs = [c_out]
  saveproblem!(prob, filename)
  saveoutput!(outputs)


  # Plotting
  fig, axs = subplots(nrows=2)
  mask = Array{Bool}(g.nx, g.ny)
  for j=1:g.ny, i=1:g.nx
    mask[i, j] = g.y[j] < -h(g.x[i]) ? true : false
  end


  # Integrate
  while prob.step < nsteps

    stepforward!(prob; nsteps=substeps)

    TracerAdvDiff.updatevars!(prob)
    increment!(diags)
    saveoutput!(outputs)

    writemessage(ep, kap, nsteps, prob, yvar)

    makeplot!(axs, prob, tfinal, yvar, mask)
    savefig(@sprintf("%s_%09d.png", plotname, prob.step), dpi=240)
  end

  return prob, yvar

end
