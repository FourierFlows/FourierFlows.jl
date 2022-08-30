# Discard `effort` argument for CuArrays
plan_flows_fft(a::CuArray, args...; flags=nothing, kwargs...) = plan_fft(a, args...; kwargs...)
plan_flows_rfft(a::CuArray, args...; flags=nothing, kwargs...) = plan_rfft(a, args...; kwargs...)

Base.zeros(::GPU, T, dims) = CUDA.zeros(T, dims)

device_array(::GPU) = CuArray
device_array(::GPU, T, dim) = CuArray{T, dim}

supersize(a::CuArray) = size(a)

getetdcoeffs(dt, L::CuArray; kwargs...) = 
    (CuArray(ζ) for ζ in getetdcoeffs(dt, Array(L); kwargs...))

makefilter(K::CuArray; kwargs...) = CuArray(makefilter(Array(K); kwargs...))

makefilter(g::AbstractGrid{Tg, <:CuArray}, T, sz; kwargs...) where Tg =
    CuArray(ones(T, sz)) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
