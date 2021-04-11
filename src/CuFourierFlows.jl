# Discard `effort` argument for CuArrays
plan_flows_fft(a::CuArray, args...; flags=nothing, kwargs...) = plan_fft(a, args...; kwargs...)
plan_flows_rfft(a::CuArray, args...; flags=nothing, kwargs...) = plan_rfft(a, args...; kwargs...)

OneDGrid(dev::GPU, args...; kwargs...) = OneDGrid(args...; ArrayType=CuArray, T=Float32, kwargs...)
TwoDGrid(dev::GPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=CuArray, T=Float32, kwargs...)
ThreeDGrid(dev::GPU, args...; kwargs...) = ThreeDGrid(args...; ArrayType=CuArray, T=Float32, kwargs...)

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
