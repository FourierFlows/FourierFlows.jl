module FourierFlows

export
  # Helper variables and macros for determining if machine is CUDA-enabled.
  @has_cuda,
  
  Device,
  CPU,
  GPU,
  ArrayType,

  cxtype,
  fltype,
  innereltype,

  AbstractVars,
  AbstractParams,

  AbstractGrid,
  ZeroDGrid,
  OneDGrid,
  TwoDGrid,
  dealias!,
  gridpoints,

  Diagnostic,
  resize!,
  update!,
  increment!,

  Output,
  saveoutput,
  saveproblem,
  groupsize,
  savediagnostic,

  @zeros,
  @devzeros,
  @createarrays,
  @superzeros,
  devzeros,
  superzeros,
  supersize,

  TimeStepper,
  ForwardEulerTimeStepper,
  FilteredForwardEulerTimeStepper,
  RK4TimeStepper,
  FilteredRK4TimeStepper,
  ETDRK4TimeStepper,
  FilteredETDRK4TimeStepper,
  AB3TimeStepper,
  FilteredAB3TimeStepper,
  stepforward!

using
  FFTW,
  JLD2,
  Statistics,
  Interpolations,
  CUDAapi,
  Requires

import Base: resize!, getindex, setindex!, lastindex, push!, append!

using Base: fieldnames

using LinearAlgebra: mul!, ldiv!

abstract type AbstractGrid{T, Ta} end
abstract type AbstractTimeStepper{T} end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractDiagnostic end

abstract type Device end
struct CPU <: Device end
struct GPU <: Device end

# The main show
include("problem.jl")
include("domains.jl")
include("utils.jl")
include("diagnostics.jl")
include("output.jl")
include("timesteppers.jl")

# Physics
include("diffusion.jl")

# Import CUDA utilities if cuda is detected.
if has_cuda()
  try
    using CuArrays # we have CUDA, so this should not fail
  catch ex
    # something is wrong with the user's set-up (or there's a bug in CuArrays)
    @warn "CUDA is installed, but CuArrays.jl fails to load" exception=(ex, catch_backtrace())
  end
end


macro has_cuda(ex)
  return has_cuda() ? :($(esc(ex))) : :(nothing)
end


function __init__()
  @require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("CuFourierFlows.jl")
end

end # module
