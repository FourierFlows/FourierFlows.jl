module FourierFlows

export 
  AbstractProblem,
  AbstractVars,
  AbstractParams,

  AbstractGrid,
  ZeroDGrid, 
  OneDGrid, 
  TwoDGrid, 
  dealias!,
  gridpoints,

  State,
  DualState,
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
  @superzeros,

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
  Interpolations

import Base: resize!, getindex, setindex!, push!, append!, fieldnames
import LinearAlgebra: mul!, ldiv!

abstract type AbstractGrid{T} end
abstract type AbstractTwoDGrid{T} <: AbstractGrid{T} end
abstract type AbstractOneDGrid{T} <: AbstractGrid{T} end

abstract type AbstractTimeStepper{T} end

abstract type AbstractProblem end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractDiagnostic end

# Some utilities

"""
    innereltype(x)
    
Recursively determine the 'innermost' type in by the collection `x` (which may be, for example, 
a collection of a collection).
"""
function innereltype(x)
  T = eltype(x)
  T <: AbstractArray ? innereltype(T) : return T
end 

"""
    cxtype(T)

Returns `T` when `T` is `Complex`, or `Complex{T}` when `T` is `Real`.
"""
cxtype(::Type{T}) where T<:Number = T
cxtype(::Type{T}) where T<:Real = Complex{T}

"""
    fltype(T)

Returns `T` when `T<:AbstractFloat` or `Tf` when `T<:Complex{Tf}`.
"""
fltype(::Type{T})          where T<:AbstractFloat = T
fltype(::Type{Complex{T}}) where T<:AbstractFloat = T

cxeltype(x) = cxtype(innereltype(x))
fleltype(x) = fltype(innereltype(x))

"""
    superzeros(T, A)

Returns an array like `A`, but full of zeros. If `innereltype(A)` can be promoted to `T`, then
the innermost elements of the array will have type `T`.
"""
superzeros(T, A) = T(0)*A
superzeros(A) = superzeros(innereltype(A), A)

function superzeros(T, dims::Tuple)
  if eltype(tup) <: Tuple
    return [ superzeros(T, d) for d in dims ]
  else
    return zeros(T, dims)
  end
end

superzeros(dims::Tuple) = superzeros(Float64, dims)

"""
    @superzeros T a b c d...
    @superzeros a b c d...

Generate arrays `b, c, d...` with the super-dimensions of `a` and innereltype `T`. If `T` is not
supplied then b, c, d have the same innereltype as `a`.
"""
macro superzeros(T::Type, A, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = superzeros($(esc(T)), $(esc(A))); ) for var in vars])
  expr
end

macro superzeros(A, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = superzeros($(esc(A))); ) for var in vars])
  expr
end


# The meat and potatoes

include("problem.jl")
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("utils.jl")
include("timesteppers.jl")


# Physics

include("heatequation1d.jl")
#include("kuramotosivashinsky.jl")

end # module
