module FourierFlows

export
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
  @createarrays,
  @superzeros,
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
  Interpolations

import Base: resize!, getindex, setindex!, lastindex, push!, append!

using Base: fieldnames

using LinearAlgebra: mul!, ldiv!

abstract type AbstractGrid{T} end
abstract type AbstractTwoDGrid{T} <: AbstractGrid{T} end
abstract type AbstractOneDGrid{T} <: AbstractGrid{T} end
abstract type AbstractTimeStepper{T} end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractDiagnostic end

# The main show
include("problem.jl")
include("domains.jl")
include("utils.jl")
include("diagnostics.jl")
include("output.jl")
include("timesteppers.jl")

# Physics
include("diffusion.jl")

end # module
