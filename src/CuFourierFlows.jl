using .CuArrays

# Discard `effort` argument for CuArrays
plan_flows_fft(a::CuArray, effort) = plan_fft(a)
plan_flows_rfft(a::CuArray, effort) = plan_rfft(a)

OneDGrid(dev::GPU, args...; kwargs...) = OneDGrid(args...; ArrayType=CuArray, kwargs...)
TwoDGrid(dev::GPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=CuArray, kwargs...)
ThreeDGrid(dev::GPU, args...; kwargs...) = ThreeDGrid(args...; ArrayType=CuArray, kwargs...)

function Base.zeros(::GPU, T, dims)
  a = CuArray{T}(undef, dims...)
  a .= 0
  return a
end

ArrayType(::GPU) = CuArray
ArrayType(::GPU, T, dim) = CuArray{T, dim}

supersize(a::CuArray) = size(a)

getetdcoeffs(dt, L::CuArray; kwargs...) = 
    (CuArray(ζ) for ζ in getetdcoeffs(dt, Array(L); kwargs...))

makefilter(K::CuArray; kwargs...) = CuArray(makefilter(Array(K); kwargs...))

function makefilter(g::AbstractGrid{Tg, <:CuArray}, T, sz; kwargs...) where Tg
  CuArray(ones(T, sz)) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
end

zerofield(g::OneDGrid{T, <:CuArray}; realvalued=true) where T = realvalued ? Field(devzeros(GPU(), T, size(g.x)), devzeros(GPU(), Complex{T}, size(g.kr)), g) : Field(devzeros(CPU(), Complex{T}, size(g.x)), devzeros(CPU(), Complex{T}, size(g.k)), g)

zerofield(g::Union{TwoDGrid{T, <:Array}, ThreeDGrid{T, <:Array}}; realvalued=true) where T = realvalued ? Field(devzeros(CPU(), T, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Krsq)), g) : Field(devzeros(CPU(), Complex{T}, size(g.Ksq)), devzeros(CPU(), Complex{T}, size(g.Ksq)), g)

