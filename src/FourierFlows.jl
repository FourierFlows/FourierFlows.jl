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




# Problem type and associated functions
type Problem <: AbstractProblem
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
include("plotting.jl")


# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/niwqg.jl")
include("physics/traceradvdiff.jl")
include("physics/tracerpatcheqn.jl")

end
