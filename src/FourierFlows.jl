__precompile__()

module FourierFlows

export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper,
       AbstractProblem

export Problem




# Abstract supertypes
abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractEquation end
abstract type AbstractTimeStepper end
abstract type AbstractProblem end


"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct Equation{dims} <: AbstractEquation
  LC::Array{Complex{Float64},dims} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualEquation{dimsc,dimsr} <: AbstractEquation
  LCc::Array{Complex{Float64},dimsc}
  LCr::Array{Complex{Float64},dimsr}
  calcN!::Function
end




# Problem type and associated functions
mutable struct Problem <: AbstractProblem
  grid::AbstractGrid
  vars::AbstractVars
  params::AbstractParams
  eqn::AbstractEquation
  ts::AbstractTimeStepper
  t::Real
  step::Int
end

function Problem(g::AbstractGrid, v::AbstractVars, p::AbstractParams,
  eq::AbstractEquation, ts::AbstractTimeStepper)
  Problem(g, v, p, eq, ts, 0.0, 0)
end



function unpack(prob::AbstractProblem)
  prob.vars, prob.params, prob.grid
end




# Include base functionality
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("timesteppers.jl")
include("utils.jl")


# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/traceradvdiff.jl")
include("physics/tracerpatcheqn.jl")

end
