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
type ProblemDiagnostic{T} <: AbstractDiagnostic
  calc::Function
  freq::Int64
  prob::Problem
  num::Int
  data::Array{T, 1}
  time::Array{Float64, 1}
  value::T
  count::Int
end


""" Constructor for the ProblemDiagnostic type. """
function ProblemDiagnostic(calc::Function, freq::Int, 
  prob::FourierFlows.Problem; num=1)
  value = calc(prob)
  T = typeof(value)

  data    = Array{T}(num)
  time    = Array{Float64}(num)

  data[1] = value
  time[1] = prob.vars.t

  ProblemDiagnostic{T}(calc, freq, prob, num, data, time, value, 1)
end




""" A Diagnostic type for FourierFlows problems with Vars, Params, Timestepper
and Grid explicltly associated. """
type Diagnostic{T} <: AbstractDiagnostic
  calc::Function
  freq::Int64
  num::Int
  count::Int
  data::Array{T, 1}
  time::Array{Float64, 1}
  value::T
  vs::AbstractVars
  pr::AbstractParams
  ts::AbstractTimeStepper
  g::AbstractGrid
end

""" Constructor for the Diagnostic type. """
function Diagnostic(calc::Function, freq, vs, pr, ts, g; num=1)
  value = calc(vs, pr, g)
  T = typeof(value)

  data    = Array{T}(num)
  time    = Array{Float64}(num)

  data[1] = value
  time[1] = vs.t

  Diagnostic{T}(calc, freq, num, 1, data, time, value, vs, pr, ts, g)
end





""" Calculate a diagnostic and store the result in Diagnostic.data. """
function increment!(diag::Diagnostic)
  if diag.count < diag.num
    diag.data[diag.count+1] = diag.calc(diag.vs, diag.pr, diag.g)
    diag.time[diag.count+1] = diag.vs.t
  else
    push!(diag.data, diag.calc(diag.vs, diag.pr, diag.g))
    push!(diag.time, diag.vs.t)
  end
  diag.count += 1
  diag.value = diag.data[diag.count]
  nothing
end


function increment!(diag::ProblemDiagnostic)
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




type Output
  name::String
  freq::Int
  calc::Function
  prob::Problem
end




type Model
  grid::AbstractGrid
  vars::AbstractVars
  params::AbstractParams
  eqn::AbstractEquation
  ts::AbstractTimeStepper
  diags::Array{AbstractDiagnostic, 1}
  outputs::Array{Any, 1}
end
