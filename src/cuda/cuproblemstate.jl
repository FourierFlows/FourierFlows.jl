export CuEquation, CuProblem, CuState, CuDualState

function CuProblem(g, v, p, eq, ts)
  st = CuState{eltype(eq.LC),ndims(eq.LC)}(0, 0, ts.dt, CuArray(zeros(eq.LC)))
  Problem(g, v, p, eq, ts, st)
end

mutable struct CuState{T,dim} <: AbstractState
  t::T
  step::Int
  dt::T
  sol::CuArray{T,dim}
end

function CuState(T::DataType, sz::Tuple, dt)
  sol = CuArray(zeros(T, sz))
  CuState(0, 0, dt, sol)
end

mutable struct CuDualState{T,dimc,dimr} <: AbstractState
  t::T
  step::Int
  dt::T
  solc::CuArray{T,dimc}
  solr::CuArray{T,dimr}
end

function CuDualState(T::DataType, sizec, sizer, dt)
  solc = CuArray(zeros(T, sizec))
  solr = CuArray(zeros(T, sizer))
  CuDualState(0, 0, dt, solc, solr)
end

# Equation Composite Type
"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct CuEquation{T,dim} <: AbstractEquation
  LC::CuArray{Complex{T},dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualCuEquation{T,dimc,dimr} <: AbstractEquation
  LCc::CuArray{Complex{T},dimc}
  LCr::CuArray{Complex{T},dimr}
  calcN!::Function
end

CuEquation(eq::AbstractEquation) = CuEquation(eq.LC, eq.calcN!)
