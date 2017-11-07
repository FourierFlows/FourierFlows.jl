__precompile__()

include("../../../../src/fourierflows.jl")


module SLSPlotTools

using JLD2

import FourierFlows.TwoModeBoussinesq



function load_boussinesq_ivp(filename)
  file = jldopen(filename, "r")

  ν0, nν0 = file["params"]["ν0"], file["params"]["nν0"]
  ν1, nν1 = file["params"]["ν1"], file["params"]["nν1"]

   n = file["grid"]["nx"]
   L = file["grid"]["Lx"]
   f = file["params"]["f"]
   N = file["params"]["N"]
   m = file["params"]["m"]
  dt = file["timestepper"]["dt"]

  steps = parse.(Int, keys(file["timeseries/t"]))

  prob = TwoModeBoussinesq.InitialValueProblem(
    nx=n, Lx=L, ν0=ν0, nν0=nν0, ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt)

  close(file)

  prob, steps
end





function get_snapshot(filename, step)

  file = jldopen(filename, "r")

     t = file["timeseries/t/$step"]
  solr = file["timeseries/solr/$step"]
  solc = file["timeseries/solc/$step"]

  close(file)

  t, solr, solc
end


# End module
end
