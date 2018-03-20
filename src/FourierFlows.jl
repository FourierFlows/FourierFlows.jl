__precompile__()

module FourierFlows

export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper,
       AbstractProblem

# Abstract supertypes
abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractTimeStepper end
abstract type AbstractEquation end
abstract type AbstractState end
abstract type AbstractProblem end

"""
    cxeltype(a)
Returns Complex{eltype(a)} if eltype(a) <: Real; eltype(a) otherwise.
"""
cxeltype(a) = eltype(a) <: Real ? Complex{eltype(a)} : eltype(a)


# Include base functionality
include("problemstate.jl")
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("utils.jl")
include("timesteppers.jl")

# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/traceradvdiff.jl")
include("physics/tracerpatcheqn.jl")
include("physics/kuramotosivashinsky.jl")

end # module
