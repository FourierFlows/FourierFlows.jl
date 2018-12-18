""" 
    Diagnostic(calc, prob; freq=1, nsteps=100, N=floor(Int, (nsteps+1)/freq))

Construct a diagnostic which stores the result of `calc(prob)` with frequency `freq`
"""
mutable struct Diagnostic{T,N} <: AbstractDiagnostic
  calc::Function
  prob::Problem
  data::Vector{T}
  t::Vector{Float64}
  steps::Vector{Int}
  freq::Int
  i::Int
end

function Diagnostic(calc, prob; freq=1, nsteps=100, ndata=ceil(Int, (nsteps+1)/freq), N=ndata)
  firstvalue = calc(prob)
  T = typeof(firstvalue)

  data = Vector{T}(undef, N)
  t = Vector{Float64}(undef, N)
  steps = Vector{Int}(undef, N)

  data[1] = firstvalue
  t[1] = prob.clock.t
  steps[1] = prob.clock.step

  Diagnostic{T,N}(calc, prob, data, t, steps, freq, 1)
end

"""
    extend!(diag::AbstractDiagnostic, n)

Extend the `data`, `time`, and `steps` vectors of `diag` by `n`.
"""
function extend!(d, n)
  resize!(d.data, length(d.steps)+n)
  resize!(d.t, length(d.steps)+n)
  resize!(d.steps, length(d.steps)+n)
  nothing
end

extend!(diag::Diagnostic{T,N}) where {T,N} = extend!(diag, N)

"""
    update!(diag)

Update `diag`.
"""
function update!(d, i)
  i <= length(d.steps) || extend!(d)
   d.data[i] = d() #d.calc(prob)
      d.t[i] = d.prob.clock.t
  d.steps[i] = d.prob.clock.step
         d.i = i
  nothing
end

"""
    increment!(diag)

Increment the Diagnostic diag, or an array of Diagnostics diags.
"""
function increment!(d)
  d.prob.clock.step % d.freq == 0 && update!(d, d.i+1)
  nothing
end

function increment!(diags::AbstractVector)
  for d in diags
    increment!(d)
  end
  nothing
end

"""
e.g. `plot(energydiag[:t], energydiag[:data])`.
"""
getindex(d::Diagnostic, s::Union{AbstractString,Symbol}) = getfield(d, Symbol(s))[1:d.i]
getindex(d::Diagnostic, idx...) = getindex(d.data, idx...)

lastindex(d::Diagnostic) = d.i

# d() = current value
# e = Diagnostic(energy, ...)
# e() returns energy at the current time.
# sensible?
(d::Diagnostic)() = d.calc(d.prob)

