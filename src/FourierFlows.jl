module FourierFlows

export 
  OneDGrid, 
  TwoDGrid, 
  dealias!,

  unpack,

  Diagnostic,
  resize!, 
  update!, 
  increment!,
    
  Output, 
  saveoutput, 
  saveproblem, 
  groupsize, 
  savediagnostic,

  @createarrays, 

  ForwardEulerTimeStepper, 
  FilteredForwardEulerTimeStepper,
  RK4TimeStepper, 
  FilteredRK4TimeStepper,
  DualRK4TimeStepper, 
  DualFilteredRK4TimeStepper,
  ETDRK4TimeStepper, 
  FilteredETDRK4TimeStepper,
  DualETDRK4TimeStepper, 
  DualFilteredETDRK4TimeStepper,
  AB3TimeStepper, 
  FilteredAB3TimeStepper,
  stepforward!

using 
  Requires, 
  FFTW, 
  Statistics,
  JLD2,
  Interpolations

using SpecialFunctions: besselj

import Base: resize!, getindex, setindex!, push!, append!, fieldnames
import LinearAlgebra: mul!, ldiv!

abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractTimeStepper end
abstract type AbstractEquation end
abstract type AbstractState end
abstract type AbstractProblem end

abstract type AbstractTwoDGrid <: AbstractGrid end
abstract type AbstractOneDGrid <: AbstractGrid end

abstract type AbstractDiagnostic end

# ------------------
# Base functionality
# ------------------

include("problemstate.jl")
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("utils.jl")
include("timesteppers.jl")


# -------
# Physics modules
# -------

include("traceradvdiff.jl")
include("kuramotosivashinsky.jl")

# ----------------------
# CUDA/GPU functionality
# ----------------------

@require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
  using CuArrays
  include("cuda/cuutils.jl")
  include("cuda/cuproblemstate.jl")
  include("cuda/cudomains.jl")
  include("cuda/cutimesteppers.jl")
end

end # module
