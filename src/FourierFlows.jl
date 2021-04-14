"""
Main module for `FourierFlows.jl` -- an ecosystem for solving partial differential equations 
on periodic domains using Fourier-based pseudospectral methods."
 
# Exports
$(EXPORTS)
"""
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
  ThreeDGrid,
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

  AbstractTimeStepper,
  TimeStepper,
  ForwardEulerTimeStepper,
  FilteredForwardEulerTimeStepper,
  RK4TimeStepper,
  FilteredRK4TimeStepper,
  ETDRK4TimeStepper,
  FilteredETDRK4TimeStepper,
  AB3TimeStepper,
  FilteredAB3TimeStepper,
  stepforward!,
  step_until!

using
  FFTW,
  JLD2,
  Statistics,
  Interpolations,
  Reexport,
  DocStringExtensions

@reexport using FFTW: fft, ifft, rfft, irfft
@reexport using CUDA

import Base: resize!, getindex, setindex!, push!, append!, show, summary

using Base: fieldnames
using FFTW: fftfreq, rfftfreq

"Abstract supertype for grids."
abstract type AbstractGrid{T, A} end

"Abstract supertype for timesteppers."
abstract type AbstractTimeStepper{T} end

"Abstract supertype for parameters."
abstract type AbstractParams end

"Abstract supertype for variables."
abstract type AbstractVars end

"Abstract supertype for diagnostics."
abstract type AbstractDiagnostic end

"Abstract supertype for device."
abstract type Device end

"CPU device."
struct CPU <: Device end

"GPU device."
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

"""
    @has_cuda expr
A macro to compile and execute `expr` only if CUDA is installed and available.
"""
macro has_cuda(expr)
  return has_cuda() ? :($(esc(expr))) : :(nothing)
end

# CUDA functionality
@has_cuda include("CuFourierFlows.jl")

function __init__()
  threads = Threads.nthreads()
  if threads > 1
    @info "FourierFlows will use $threads threads"
    FFTW.set_num_threads(threads)
  end

  @has_cuda begin
    @debug "CUDA-enabled GPU(s) detected:"
    for (gpu, dev) in enumerate(CUDA.devices())
      @debug "$dev: $(CUDA.name(dev))"
    end

    CUDA.allowscalar(false)
  end
end


function show(io::IO, vars::AbstractVars)
  names = propertynames(vars)
  showstring = ""
  for name in names[1:end-1]
    field = getproperty(vars, name)
    showstring = string(showstring, "  ├───── variable: " * string(name) * " -> ", summary(field), "\n")
  end
  name = names[end]
  field = getproperty(vars, name)
  showstring = string(showstring, "  └───── variable: " * string(name) * " -> ", summary(field), "\n")
  
  return print(io, "Variables\n", showstring)
end

summary(::Function) = "Function"

function show(io::IO, params::AbstractParams)
  names = propertynames(params)
  showstring = ""
  for name in names[1:end-1]
    field = getproperty(params, name)
    showstring = string(showstring, "  ├───── parameter: " * string(name) * " -> ", summary(field), "\n")
  end
  name = names[end]
  field = getproperty(params, name)
  showstring = string(showstring, "  └───── parameter: " * string(name) * " -> ", summary(field), "\n")
  
  return print(io, "Parameters\n", showstring)
end

end # module
