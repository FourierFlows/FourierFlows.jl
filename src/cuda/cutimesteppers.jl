export CuForwardEulerTimeStepper, CuFilteredForwardEulerTimeStepper,
       CuRK4TimeStepper, CuFilteredRK4TimeStepper

struct CuForwardEulerTimeStepper{T,dim} <: AbstractForwardEulerTimeStepper
  dt::Float64
  N::CuArray{Complex{T},dim}    # Explicit linear and nonlinear terms
  CuForwardEulerTimeStepper(dt::T, N::CuArray{Complex{T},dim}) where {T,dim} = new{T,dim}(
    dt, CuArray(zeros(size(N))))
end

CuForwardEulerTimeStepper(dt, N::Array{Complex{T},dim}) where {T,dim} = CuForwardEulerTimeStepper(dt, CuArray(N))

struct CuFilteredForwardEulerTimeStepper{T,dim} <: AbstractFilteredForwardEulerTimeStepper
  dt::Float64
  N::CuArray{Complex{T},dim}    # Explicit linear and nonlinear terms
  filter::CuArray{T,dim}        # Filter for solution
end

function CuFilteredForwardEulerTimeStepper(dt, LC, g; filterkwargs...)
  N = CuArray(zeros(cxeltype(LC), size(LC)))
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  CuFilteredForwardEulerTimeStepper{ndims(N)}(dt, N, filter)
end

struct CuRK4TimeStepper{T,dim} <: AbstractRK4TimeStepper
  dt::Float64
  sol₁::CuArray{T,dim}
  RHS₁::CuArray{T,dim}
  RHS₂::CuArray{T,dim}
  RHS₃::CuArray{T,dim}
  RHS₄::CuArray{T,dim}
end

function CuRK4TimeStepper(dt, LC)
  @createcuarrays cxeltype(LC) size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  CuRK4TimeStepper{cxeltype(LC),ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end


struct CuFilteredRK4TimeStepper{T,dim} <: AbstractFilteredRK4TimeStepper
  dt::Float64
  sol₁::CuArray{T,dim}
  RHS₁::CuArray{T,dim}
  RHS₂::CuArray{T,dim}
  RHS₃::CuArray{T,dim}
  RHS₄::CuArray{T,dim}
  filter::CuArray{T,dim}    # Filter for solution
end

function CuFilteredRK4TimeStepper(dt, LC, g; filterkwargs...)
  @createcuarrays cxeltype(LC) size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  CuFilteredRK4TimeStepper{cxeltype(LC),ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄, filter)
end

ForwardEulerTimeStepper(dt, LC::CuArray{T,dim}) where {T,dim} = CuForwardEulerTimeStepper(dt, LC)
RK4TimeStepper(dt, LC::CuArray{T,dim}) where {T,dim} = CuRK4TimeStepper(dt, LC)

FilteredForwardEulerTimeStepper(
  dt, LC::CuArray{T,dim}, g; filterkwargs...) where {T,dim} = CuFilteredForwardEulerTimeStepper(dt, LC, g; 
                                                                                                filterkwargs...)
FilteredRK4TimeStepper(
  dt, LC::CuArray{T,dim}, g; filterkwargs...) where {T,dim} = CuFilteredRK4TimeStepper(dt, LC, g; filterkwargs...)
