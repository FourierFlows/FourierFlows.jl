export @cuconvertarrays, @createcuarrays

macro cuconvertarrays(vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = CuArray($(esc(var))); ) for var in vars])
  expr
end

macro createcuarrays(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = CuArray(zeros($(esc(T)), $(esc(dims)))); ) for var in vars])
  expr
end
