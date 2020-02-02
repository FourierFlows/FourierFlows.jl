struct Field{T<:AbstractFloat, Aphys<:AbstractArray, Atrans<:AbstractArray} <: AbstractField{T, Aphys, Atrans}
         grid_values :: Aphys
  coefficient_values :: Atrans
                grid :: AbstractGrid{T}
end

"""
    zerofield(g)

Constructs a field that lives on grid g with zero values.
"""
zerofield(g::OneDGrid{T, <:Array}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.x)), devzeros(CPU(), Complex{T}, size(g.kr)), g) : Field(devzeros(CPU(), Complex{T}, size(g.x)), devzeros(CPU(), Complex{T}, size(g.k)), g)

zerofield(g::Union{TwoDGrid{T, <:Array}, ThreeDGrid{T, <:Array}}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Krsq)), g) : Field(devzeros(CPU(), Complex{T}, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Ksq)), g)

"""
    compute_coeff_space!(f)

Computes the coefficient values of a field f from its grid values.
"""
function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  mul!(field.coefficient_values, field.grid.rfftplan, field.grid_values)
end

function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  mul!(field.coefficient_values, field.grid.fftplan, field.grid_values)
end

"""
    compute_grid_space!(f)

Computes the grid values of a field f from its coefficient values.
"""
function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  ldiv!(field.grid_values, field.grid.rfftplan, deepcopy(field.coefficient_values))
end

function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  ldiv!(field.grid_values, field.grid.fftplan, field.coefficient_values)
end

"""
    field_from_grid_values(a, g)

Constructs a field that lives on grid g with grid values given by array a.
"""
function field_from_grid_values(a::Aphys, g::AbstractGrid{T}) where {T, Aphys}
  F = Aphys<:AbstractArray{T} ? zerofield(g) :  zerofield(g; realvalued=false)
  F.grid_values .= a
  compute_coeff_space!(F)
  return F
end

"""
    field_from_coefficient_values(ah, g)

Constructs a field that lives on grid g with coefficient values given by array ah.
"""
function field_from_coefficient_values(ah::Atrans, g::AbstractGrid{T}; realvalued=true) where {T, Atrans<:AbstractArray{Complex{T}}}
  F = zerofield(g; realvalued=realvalued)
  F.coefficient_values .= ah
  compute_grid_space!(F)
  return F
end
