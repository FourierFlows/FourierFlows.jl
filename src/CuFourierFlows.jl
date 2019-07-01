using .CuArrays

# Discard `effort` argument for CuArrays
plan_flows_fft(a::CuArray, effort) = plan_fft(a)
plan_flows_rfft(a::CuArray, effort) = plan_rfft(a)

OneDGrid(dev::GPU, args...; kwargs...) = OneDGrid(args...; ArrayType=CuArray, kwargs...)
TwoDGrid(dev::GPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=CuArray, kwargs...)
