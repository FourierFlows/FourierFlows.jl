using JLD2, HDF5
using PyPlot, PyCall, NullableArrays
@pyimport numpy.ma as ma

# Plotting
PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)

# For integrals
ddxy(g) = g.dx*g.dy
ddx(g)  = g.dx
ddy(g)  = g.dy



# Moment-calculating functions
M0(c, g)     = ddxy(g)*sum(c)
Mxn(c, g, n) = ddxy(g)*sum(g.X.^n.*c)
Myn(c, g, n) = ddxy(g)*sum(g.Y.^n.*c) 


# Cumulants
Cx1(c, g) = ddxy(g)*sum(g.X.*c) / M0(c, g)
Cy1(c, g) = ddxy(g)*sum(g.Y.*c) / M0(c, g)

Cx2(c, g) = ddxy(g)*sum((g.X-Cx1(c, g)).^2.0.*c) / M0(c, g)
Cy2(c, g) = ddxy(g)*sum((g.Y-Cy1(c, g)).^2.0.*c) / M0(c, g)

Cx3(c, g) = ddxy(g)*sum((g.X-Cx1(c, g)).^3.0.*c) / M0(c, g)
Cy3(c, g) = ddxy(g)*sum((g.Y-Cy1(c, g)).^3.0.*c) / M0(c, g)

intx(II, g) = ddx(g)*sum(II)
inty(II, g) = ddy(g)*sum(II)

mxn(c, g, n) = ddx(g)*sum(g.Y.^n.*c, 1)
myn(c, g, n) = ddy(g)*sum(g.Y.^n.*c, 2)

cx1(c, g) = ddx(g)*sum(g.X.*c) ./ mxn(c, g, 0)
cy1(c, g) = ddy(g)*sum(g.Y.*c) ./ myn(c, g, 0)

delyc1(c, g) = g.Y - broadcast(*, cy1(c, g), ones(g.nx, g.ny))
cy2(c, g) = ddy(g) * sum(delyc1(c, g).*c, 2)./myn(c, g, 0)




""" Initialize a tracer problem with constant diffusivity. """
function TracerProblem(g::AbstractGrid, kap::Real, u::Function, w::Function,
  dt::Real)

  vs = TracerAdvDiff.Vars(g)
  pr = TracerAdvDiff.ConstDiffParams(kap, u, w)
  eq = TracerAdvDiff.Equation(pr, g)
  ts = FourierFlows.ETDRK4TimeStepper(dt, eq.LC)

  Problem(g, vs, pr, eq, ts)
end




abstract type AbstractDiagnostic end




""" A diagnostic type associated with FourierFlows.Problem types """
type Diagnostic{T} <: AbstractDiagnostic
  calc::Function
  prob::Problem
  num::Int
  data::Array{T, 1}
  time::Array{Float64, 1}
  value::T
  count::Int
end

""" Constructor for the ProblemDiagnostic type. """
function Diagnostic(calc::Function, prob::FourierFlows.Problem; num=1)
  value = calc(prob)
  T = typeof(value)

  data    = Array{T}(num)
  time    = Array{Float64}(num)

  data[1] = value
  time[1] = prob.vars.t

  Diagnostic{T}(calc, prob, num, data, time, value, 1)
end




function increment!(diag::Diagnostic)
  if diag.count < diag.num
    diag.data[diag.count+1] = diag.calc(diag.prob)
    diag.time[diag.count+1] = diag.prob.vars.t
  else
    push!(diag.data, diag.calc(diag.prob))
    push!(diag.time, diag.prob.vars.t)
  end
  diag.count += 1
  diag.value = diag.data[diag.count]
  nothing
end


function increment!(diags::AbstractArray)
  for d in diags
    increment!(d)
  end
  nothing
end




""" Output type for FourierFlows problems. """
type Output
  name::String
  calc::Function
  prob::Problem
  filename::String
end


""" Save output to file. """
function saveoutput!(out::Output)
  step = out.prob.step
  groupname = "timeseries"
  name = out.name

  jldopen(out.filename, "a+") do file
    file["$groupname/$name/$step"] = out.calc(out.prob)
    file["$groupname/t/$step"]     = out.prob.t
  end

  nothing
end


""" Save an array of outputs to file. """
function saveoutput!(outs::AbstractArray)

  step = outs[1].prob.step
  groupname = "timeseries"

  jldopen(outs[1].filename, "a+") do file
    file["$groupname/t/$step"] = outs[1].prob.t # save timestamp
    for out in outs # save output data 
      name = out.name
      file["$groupname/$name/$step"] = out.calc(out.prob)
    end
  end

  nothing
end






type Model
  name::String
  prob::AbstractProblem
  grid::AbstractGrid
  vars::AbstractVars
  params::AbstractParams
  eqn::AbstractEquation
  ts::AbstractTimeStepper
  t::Real
  step::Int
  diags::Array{AbstractDiagnostic, 1}
  outputs::Array{Output, 1}
  diagfreq::Int
  outputfreq::Int
end

""" Model Constructor. """
function Model(name::String, prob::Problem)
  diags = Array{AbstractDiagnostic, 1}(0)
  outputs = Array{Output, 1}(0)
  Model(name, prob, prob.grid, prob.vars, prob.params, prob.eqn, prob.ts,
    prob.t, prob.step, diags, outputs, 0, 0)
end

""" Model Constructor. """
function Model(name::String, prob::Problem, diags, outputs; 
  diagfreq=0, outputfreq=0)
  Model(name, prob, prob.grid, prob.vars, prob.params, prob.eqn, prob.ts,
    prob.t, prob.step, diags, outputs, diagfreq, outputfreq)
end
 
function Model(name::String, prob::Problem; 
  diags=nothing, outputs=nothing, diagfreq=0, outputfreq=0)

  if !(typeof(diags) <: AbstractArray);  diags   = [diags];   end
  if !(typeof(output) <: AbstractArray); outputs = [outputs]; end

  Model(name, prob, prob.grid, prob.vars, prob.params, prob.eqn, prob.ts,
    prob.t, prob.step, diags, outputs, diagfreq, outputfreq)
end
  


function groupsize(group::JLD2.Group)
  try
    value = length(group.unwritten_links) + length(group.written_links)
  catch
    value = length(group.unwritten_links)
  end
  return value
end
