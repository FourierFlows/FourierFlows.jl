export CuEquation, CuProblem, CuState, CuDualState

function CuProblem(g, v, p, eq, ts)
  st = CuState(cxeltype(eq.LC), size(eq.LC), ts.dt)
  Problem(g, v, p, eq, ts, st)
end

mutable struct CuState{T,dim} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  sol::CuArray{T,dim}
end

CuState(T::DataType, sz::Tuple, dt) = CuState(0.0, 0, dt, CuArray(zeros(T, sz)))

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
  LC::CuArray{T,dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualCuEquation{T,dimc,dimr} <: AbstractEquation
  LCc::CuArray{Complex{T},dimc}
  LCr::CuArray{Complex{T},dimr}
  calcN!::Function
end

CuEquation(eq::AbstractEquation) = CuEquation(CuArray(eq.LC), eq.calcN!)
#CuEquation(eq::AbstractEquation) = CuEquation(eq.LC, eq.calcN!)
