__precompile__()

export AbstractDiagnostic, Diagnostic
export increment!


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
