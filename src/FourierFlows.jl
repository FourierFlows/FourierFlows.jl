module FourierFlows

export 
  cxtype,
  fltype,
  innereltype,

  AbstractProblem,
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

import Base: resize!, getindex, setindex!, push!, append!
using Base: fieldnames
using LinearAlgebra: mul!, ldiv!

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
superzeros(T, A::AbstractArray) = T(0)*A

superzeros(A::AbstractArray) = superzeros(innereltype(A), A)
superzeros(T, dims::Tuple) = eltype(dims) <: Tuple ? [ superzeros(T, d) for d in dims ] : zeros(T, dims)
superzeros(dims::Tuple) = superzeros(Float64, dims) # default

"""
    @superzeros T a b c d...
    @superzeros T dims b c d...

Generate arrays `b, c, d...` with the super-dimensions of `a` and innereltype `T`.
"""
macro superzeros(T, ad, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = superzeros($(esc(T)), $(esc(ad))); ) for var in vars])
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

include("diffusion.jl")
#include("kuramotosivashinsky.jl")

end # module
