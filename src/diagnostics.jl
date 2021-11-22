"""
    mutable struct Diagnostic{T, N} <: AbstractDiagnostic

A diagnostic that includes `N` elements of type `T`.

$(TYPEDFIELDS)
"""
mutable struct Diagnostic{T, N} <: AbstractDiagnostic
    "function that returns the diagnostic via `calc(prob)`"
   calc :: Function
    "the relevant problem for this diagnostic"
   prob :: Problem
    "vector where the diagnostic time-series is saved"
   data :: Vector{T}
    "vector with the times for which the diagnostic was saved"
      t :: Vector{Float64}
    "vector with the problem's `step` for which the diagnostic was saved"
  steps :: Vector{Int}
    "integer denoting how often (every how many `problem.step`s) to save the diagnostic"
   freq :: Int
    "integer denoting how many times the diagnostic.data was updated"
      i :: Int
end

"""
    Diagnostic(calc, prob; freq=1, nsteps=100, ndata=ceil(Int, (nsteps+1)/freq))

Construct a diagnostic that stores the result of `calc(prob)` with frequency `freq`.

Keywords
========
  - `freq`: Diagnostic is saved every `freq` steps.
  - `nsteps`: The total number of steps in problem.
  - `ndata`: The number of diagnostics to be saved.
"""
function Diagnostic(calc, prob; freq=1, nsteps=100, ndata=ceil(Int, (nsteps+1)/freq))
  firstvalue = calc(prob)
  T = typeof(firstvalue)

  data = Vector{T}(undef, ndata)
  t = Vector{Float64}(undef, ndata)
  steps = Vector{Int}(undef, ndata)

  data[1] = firstvalue
  t[1] = prob.clock.t
  steps[1] = prob.clock.step

  return Diagnostic{T,ndata}(calc, prob, data, t, steps, freq, 1)
end

"""
    extend!(diag, n)

Extend the `data`, `time`, and `steps` vectors of the diagnostic `diag` by `n`.
"""
function extend!(diag, n)
  resize!(diag.data, length(diag.steps)+n)
  resize!(diag.t, length(diag.steps)+n)
  resize!(diag.steps, length(diag.steps)+n)

  return nothing
end

"""
    extend!(diag::Diagnostic{T,N})

Double the extend of the `data`, `time`, and `steps` vectors of the diagnostic `diag`.
"""
extend!(diag::Diagnostic{T,N}) where {T,N} = extend!(diag, N)

"""
    update!(diag)

Update `diag`.
"""
function update!(diag, i)
  i <= length(diag.steps) || extend!(diag)
  
   diag.data[i] = diag() #diag.calc(prob)
      diag.t[i] = diag.prob.clock.t
  diag.steps[i] = diag.prob.clock.step
         diag.i = i
         
  return nothing
end

"""
    increment!(diag)

Increment the Diagnostic `diag`, or an array of Diagnostics `diags`.
"""
function increment!(diag)
  diag.prob.clock.step % diag.freq == 0 && update!(diag, diag.i+1)
  
  return nothing
end

function increment!(diags::AbstractVector)
  for diag in diags
    increment!(diag)
  end
  
  return nothing
end

"""
e.g. `plot(energydiag[:t], energydiag[:data])`.
"""
getindex(d::Diagnostic, s::Union{AbstractString,Symbol}) = getfield(d, Symbol(s))[1:d.i]
getindex(d::Diagnostic, idx...) = getindex(d.data, idx...)

# d() = current value
# e = Diagnostic(energy, ...)
# e() returns energy at the current time.
# sensible?
(d::Diagnostic)() = d.calc(d.prob)

show(io::IO, d::Diagnostic{T, N}) where {T, N} =
     print(io, "Diagnostic\n",
               "  ├─── calc: ", d.calc, '\n',
               "  ├─── prob: ", summary(d.prob), '\n', 
               "  ├─── data: ", summary(d.data), '\n', 
               "  ├────── t: ", summary(d.t), '\n', 
               "  ├── steps: ", summary(d.steps), '\n', 
               "  ├─── freq: ", d.freq, '\n', 
               "  └────── i: ", d.i)
           