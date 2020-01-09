struct Field{T<:AbstractFloat, Aphys<:AbstractArray, Atrans<:AbstractArray} <: AbstractField{T, Aphys, Atrans}
         grid_values :: Aphys
  coefficient_values :: Atrans
                grid :: AbstractGrid{T}
end

zerofield(g::OneDGrid{T, <:Array}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.x)), devzeros(CPU(), Complex{T}, size(g.kr)), g) : Field(devzeros(CPU(), Complex{T}, size(g.x)), devzeros(CPU(), Complex{T}, size(g.k)), g)

zerofield(g::Union{TwoDGrid{T, <:Array}, ThreeDGrid{T, <:Array}}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Krsq)), g) : Field(devzeros(CPU(), Complex{T}, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Ksq)), g)

function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  mul!(field.coefficient_values, field.grid.rfftplan, field.grid_values)
end

function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  ldiv!(field.grid_values, field.grid.rfftplan, field.coefficient_values)
end

function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  mul!(field.coefficient_values, field.grid.fftplan, field.grid_values)
end

function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  ldiv!(field.grid_values, field.grid.fftplan, field.coefficient_values)
end

function field_from_grid_valuesphysical(A::Aphys, g::AbstractGrid{T}) where {Aphys, T}
  F = Aphys<:AbstractArray{T} ? zerofield(g) :  zerofield(g; realvalued=false)
  F.grid_values .= A
  compute_coeff_space!(F)
  return F
end

function field_from_coefficient_values(A::Atrans, g::AbstractGrid{T}) where {Atrans, T}
  F = zerofield(g)
  F.coefficient_values .= A
  compute_grid_space!(F)
  return F
end
