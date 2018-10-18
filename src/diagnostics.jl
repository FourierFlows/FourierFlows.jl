""" A diagnostic type associated with FourierFlows.Problem types """
mutable struct Diagnostic{T,N} <: AbstractDiagnostic
  calc::Function
  prob::Problem
  data::Array{T,1}
  time::Array{Float64,1}
  steps::Array{Int,1}
  freq::Int
  i::Int
end

""" Constructor for the Diagnostic type. """
function Diagnostic(calc::Function, prob::FourierFlows.Problem; freq=1, nsteps=1, N=ceil(Int, (nsteps+1)/freq))

  firstvalue = calc(prob)

   data = Array{typeof(firstvalue)}(undef, N)
  time = Array{Float64}(undef, N)
  steps = Array{Int}(undef, N)

  data[1] = firstvalue
  time[1] = prob.clock.t
  step[1] = prob.clock.step
  count = DiagStatus(1)

  Diagnostic{T,N}(calc, prob, data, time, steps, freq, 1)
end

"""
    update!(diag)

Update `diag`.
"""
function update!(d::Diagnostic, i)
   d.data[i] = d.calc(prob)
  d.time[i] = d.prob.clock.t
  d.steps[i] = d.prob.clock.step
  nothing
end

"""
    increment!(diag)

Increment the Diagnostic diag, or an array of Diagnostics diags.
"""
function increment!(d::AbstractDiagnostic)
  d.i += 1
  update!(d, d.i)
  nothing
end

"""
e.g. `plot(energy["time"], energy["data"])`.
"""
getindex(d, s::AbstractString) = getindex(getfield(d, Symbol(s)), 1:d.i)

# d() = current value
# e = Diagnostic(energy, ...)
# e() returns energy at the current time.
# sensible?
(d::Diagnostic)() = d.calc(d.prob)

