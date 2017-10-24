using PyPlot, PyCall, NullableArrays
@pyimport numpy.ma as ma

# Plotting
PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)

# For integrals
iixy(g) = g.dx*g.dy/(g.nx*g.ny)

# Moment-calculating functions
Mxn(c, g, n) = iixy(g)*sum(g.X.^n.*c)
Myn(c, g, n) = iixy(g)*sum(g.Y.^n.*c) 

# Cumulants
Cy1(c, g) = iixy(g)*sum(g.Y.*c) / Myn(c, g, 0)
Cy2(c, g) = iixy(g)*sum((g.Y-Cy1(c, g)).^2.0.*c) / Myn(c, g, 0)


abstract type AbstractDiagnostic end

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

function Diagnostic(calc::Function, freq, vs, pr, ts, g; num=1)

  value = calc(vs, pr, g)
  T = typeof(value)

  data    = Array{T}(num)
  time    = Array{Float64}(num)

  data[1] = value
  time[1] = vs.t

  Diagnostic{T}(calc, freq, num, 1, data, time, value, vs, pr, ts, g)
end


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

function increment!(diags::AbstractArray)
  for d in diags
    increment!(d)
  end
  nothing
end
