__precompile__()

module FourierFlows


export AbstractVars,
       AbstractParams,
       AbstractGrid,
       AbstractEquation,
       AbstractTimeStepper

abstract type AbstractVars end
abstract type AbstractParams end
abstract type AbstractGrid end
abstract type AbstractEquation end
abstract type AbstractTimeStepper end




include("domains.jl")
include("timesteppers.jl")
include("utils.jl")

include("physics/twodturb.jl")
include("physics/barotropicqg.jl")

end
