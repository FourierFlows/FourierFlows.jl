struct Field{T<:AbstractFloat, Aphys<:AbstractArray, Atrans<:AbstractArray} <: AbstractField{T, Aphys, Atrans}
     physicalvalues :: Aphys
  transformedvalues :: Atrans
               grid :: AbstractGrid{T}
end

zerofield(g::OneDGrid{T, <:Array}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.x)), devzeros(CPU(), Complex{T}, size(g.kr)), g) : Field(devzeros(CPU(), Complex{T}, size(g.x)), devzeros(CPU(), Complex{T}, size(g.k)), g)

zerofield(g::Union{TwoDGrid{T, <:Array}, ThreeDGrid{T, <:Array}}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Krsq)), g) : Field(devzeros(CPU(), Complex{T}, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Ksq)), g)

function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  mul!(field.transformedvalues, field.grid.rfftplan, field.physicalvalues)
end

function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{T}, Atrans}
  ldiv!(field.physicalvalues, field.grid.rfftplan, field.transformedvalues)
end

function compute_coeff_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  mul!(field.transformedvalues, field.grid.fftplan, field.physicalvalues)
end

function compute_grid_space!(field::Field{T, Aphys, Atrans}) where {T, Aphys<:AbstractArray{Complex{T}}, Atrans}
  ldiv!(field.physicalvalues, field.grid.fftplan, field.transformedvalues)
end
