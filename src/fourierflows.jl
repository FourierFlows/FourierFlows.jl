__precompile__()

module FourierFlows


export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper

abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractEquation end
abstract type AbstractTimeStepper end

type Problem
  g::AbstractGrid
  v::AbstractVars
  p::AbstractParams
  eq::AbstractEquation
  ts::AbstractTimeStepper
end




include("domains.jl")
include("timesteppers.jl")
include("utils.jl")

include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/niwqg.jl")

end
