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
  CUDA,
  DocStringExtensions

import Base: resize!, getindex, setindex!, push!, append!, show

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

if has_cuda()
  try
    include("CuFourierFlows.jl")
  catch ex
    # something is wrong with the user's set-up (or there's a bug in CuArrays)
    @warn "CUDA is installed, but CuArrays.jl fails to load" exception=(ex, catch_backtrace())
  end
end


"""
    @has_cuda expr
A macro to compile and execute `expr` only if CUDA is installed and available. Generally 
used to wrap expressions that can only be compiled if `CuArrays` and `CUDAnative` can be
loaded.
"""
macro has_cuda(expr)
  return has_cuda() ? :($(esc(expr))) : :(nothing)
end


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
    showstring = string(showstring, "  ├───── variable: " * string(name) * ", size: ", size(field), ", type: ", eltype(field), "\n")
  end
  name = names[end]
  field = getproperty(vars, name)
  showstring = string(showstring, "  └───── variable: " * string(name) * ", size: ", size(field), ", type: ", eltype(field), "\n")
  
  return print(io, "Variables\n", showstring)
end

function show(io::IO, params::AbstractParams)
  names = propertynames(params)
  showstring = ""
  for name in names[1:end-1]
    field = getproperty(params, name)
    showstring = string(showstring, "  ├───── parameter: " * string(name) * ", size: ", size(field), ", type: ", eltype(field), "\n")
  end
  name = names[end]
  field = getproperty(params, name)
  showstring = string(showstring, "  └───── parameter: " * string(name) * ", size: ", size(field), ", type: ", eltype(field), "\n")
  
  return print(io, "Parameters\n", showstring)
end

end # module
