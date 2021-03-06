"""
    innereltype(x)

Recursively determine the 'innermost' type in by the collection `x` (which may be, for example,
a collection of a collection).
"""
function innereltype(x)
  T = eltype(x)
  T <: AbstractArray ? innereltype(T) : return T
end

"""
    cxtype(T)

Returns `T` when `T` is `Complex`, or `Complex{T}` when `T` is `Real`.
"""
cxtype(::Type{T}) where T<:Number = T
cxtype(::Type{T}) where T<:Real = Complex{T}

"""
    fltype(T)

Returns `T` when `T<:AbstractFloat` or `Tf` when `T<:Complex{Tf}`.
"""
fltype(::Type{T})          where T<:AbstractFloat = T
fltype(::Type{Complex{T}}) where T<:AbstractFloat = T
fltype(T::Tuple) = fltype(T[1])

"""
    superzeros(T, A)

Returns an array like `A`, but full of zeros. If `innereltype(A)` can be promoted to `T`, then
the innermost elements of the array will have type `T`.
"""
superzeros(T, A::AbstractArray) = T(0)*A
superzeros(A::AbstractArray) = superzeros(innereltype(A), A)
superzeros(T, dims::Tuple) = eltype(dims) <: Tuple ? [ superzeros(T, d) for d in dims ] : zeros(T, dims)
superzeros(dims::Tuple) = superzeros(Float64, dims) # default
superzeros(T::Tuple, dims::Tuple) = [ superzeros(T[i], dims[i]) for i=1:length(dims) ]

"""
    @superzeros T a b c d...
    @superzeros T dims b c d...

Generate arrays `b, c, d...` with the super-dimensions of `a` and innereltype `T`.
"""
macro superzeros(T, ad, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = superzeros($(esc(T)), $(esc(ad))); ) for var in vars])
  
  return expr
end
supersize(a) = Tuple([size(ai) for ai in a])
supersize(a::Array{T}) where T<:Number = size(a)

macro createarrays(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  
  return expr
end

"""
    @zeros T dims a b c...

Create arrays of all zeros with element type `T`, size `dims`, and global names
`a`, `b`, `c` (for example). An arbitrary number of arrays may be created.
"""
macro zeros(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  
  return expr
end

Base.zeros(::CPU, T, dims) = zeros(T, dims)

"""
    devzeros(dev, T, dims)

Returns an array like `A` of type `T`, but full of zeros.
"""
devzeros(dev, T, dims) = zeros(dev, T, dims)


"""
    @devzeros dev T dims a b c...

Create arrays of all zeros with element type `T`, size `dims`, and global names
`a`, `b`, `c` (for example) on device `dev`.
"""
macro devzeros(dev, T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = zeros($(esc(dev))(), $(esc(T)), $(esc(dims))); ) for var in vars])
  
  return expr
end


"""
    parsevalsum2(uh, grid)

Returns `∫ u² dxdy = Σ|uh|²` on the `grid`. More specifically, it returns
```math
\\int u(\\boldsymbol{x})^2 \\, \\mathrm{d}^2 \\boldsymbol{x} = \\sum_{\\boldsymbol{k}} |\\hat{u}_{\\boldsymbol{k}}|^2 L_x L_y
```
where ``\\hat{u}_{\\boldsymbol{k}} =`` `uh / grid.nx`. 
"""
function parsevalsum2(uh, grid::TwoDGrid)
  if size(uh, 1) == grid.nkr # uh is in conjugate symmetric form
    U = sum(abs2, uh[1, :])           # k=0 modes
    U += 2*sum(abs2, uh[2:end, :])    # sum k>0 modes twice
  else # count every mode once
    U = sum(abs2, uh)
  end

  norm = grid.Lx * grid.Ly / (grid.nx^2 * grid.ny^2) # normalization for dft

  return norm * U
end

function parsevalsum2(uh, grid::OneDGrid)
  if size(uh, 1) == grid.nkr                 # uh is conjugate symmetric
    U = sum(abs2, CUDA.@allowscalar uh[1])   # k=0 modes
    U += @views 2 * sum(abs2, uh[2:end])     # sum k>0 modes twice
  else # count every mode once
    U = sum(abs2, uh)
  end
  
  norm = grid.Lx / grid.nx^2 # normalization for dft
  
  return norm * U
end

"""
    parsevalsum(uh, grid)

Returns `real(Σ uh)` on the `grid`.
"""
function parsevalsum(uh, grid::TwoDGrid)
  if size(uh, 1) == grid.nkr    # uh is conjugate symmetric
    U = sum(uh[1, :])           # k=0 modes
    U += 2*sum(uh[2:end, :])    # sum k>0 modes twice
  else # count every mode once
    U = sum(uh)
  end

  norm = grid.Lx * grid.Ly / (grid.nx^2 * grid.ny^2) # weird normalization for dft
  
  return norm * real(U)
end

"""
    jacobianh(a, b, grid)

Returns the Fourier transform of the Jacobian of `a` and `b` on `grid`.
"""
function jacobianh(a, b, grid::TwoDGrid)
  if eltype(a) <: Real
    bh = rfft(b)
    bx = irfft(im * grid.kr .* bh, grid.nx)
    by = irfft(im * grid.l  .* bh, grid.nx)
    return im * grid.kr .* rfft(a .* by) .- im * grid.l .* rfft(a .* bx)
  else
    # J(a, b) = ∂(a * ∂b/∂y)/∂x - ∂(a * ∂b/∂x)/∂y
    bh = fft(b)
    bx = ifft(im * grid.k .* bh)
    by = ifft(im * grid.l .* bh)
    return im * grid.k .* fft(a .* by) .- im * grid.l .* fft(a .* bx)
  end
end

"""
    jacobian(a, b, grid)

Returns the Jacobian of `a` and `b` on `grid`.
"""
function jacobian(a, b, grid::TwoDGrid)
  if eltype(a) <: Real
    return irfft(jacobianh(a, b, grid), grid.nx)
  else
    return ifft(jacobianh(a, b, grid))
  end
end

"""
    radialspectrum(ah, grid; n=nothing, m=nothing, refinement=2)

Returns `aρ = ∫ ah(ρ, θ) ρ dρ dθ`, the radial spectrum of `ah` known on the
Cartesian wavenumber grid `(k, l)`.

`aρ` is found by intepolating `ah` onto a polar wavenumber grid (ρ, θ), and
then integrating over `θ` to find `aρ`. The default resolution (n, m) for the
polar wave number grid is `n=refinement * maximum(nk, nl),
m = refinement * maximum(nk, nl)`, where `refinement = 2` by default. If
`ah` is in conjugate symmetric form only the upper half plane in `θ` is
represented on the polar grid.
"""
function radialspectrum(ah, grid::TwoDGrid; n=nothing, m=nothing, refinement=2)

  n = n == nothing ? refinement * maximum([grid.nk, grid.nl]) : n
  m = m == nothing ? refinement * maximum([grid.nk, grid.nl]) : m

  # Calcualte shifted k and l
  lshift = range(-grid.nl/2+1, stop=grid.nl/2, length=grid.nl) * 2π/grid.Ly

  if size(ah)[1] == grid.nkr # conjugate symmetric form
    m = Int(m/2)                         # => half resolution in θ
    θ = range(-π/2, stop=π/2, length=m)  # θ-grid from k=0 to max(kr)
    ahshift = fftshift(ah, 2)            # shifted ah
    kshift = range(0, stop=grid.nkr-1, length=grid.nkr) * 2π/grid.Lx
  else # ordinary form
    θ = range(0, stop=2π, length=m)      # θ grid
    ahshift = fftshift(ah, [1, 2])       # shifted ah
    kshift = range(-grid.nk/2+1, stop=grid.nk/2, length=grid.nk) * 2π/grid.Lx
  end

  # Interpolator for ah
  itp = scale(interpolate(ahshift, BSpline(Linear())), kshift, lshift)

  # Get radial wavenumber vector
  ρmax = minimum([(grid.nk/2-1) * 2π/grid.Lx, (grid.nl/2-1) * 2π/grid.Ly])
  ρ = range(0, stop=ρmax, length=n)

  # Interpolate ah onto fine grid in (ρ,θ).
  ahρθ = zeros(eltype(ahshift), (n, m))

  for i₁=2:n, i₂=1:m # ignore zeroth mode; i₁≥2
    ahρθ[i₁, i₂] = itp(ρ[i₁] * cos(θ[i₂]), ρ[i₁] * sin(θ[i₂]))
  end

  # ahρ = ρ ∫ ah(ρ,θ) dθ  =>  Ah = ∫ ahρ dρ = ∫∫ ah dk dl
  dθ = θ[2]-θ[1]
  if size(ah)[1] == grid.nkr
    ahρ = 2ρ .* sum(ahρθ, dims=2) * dθ # multiply by 2 for conjugate symmetry
  else
    ahρ =  ρ .* sum(ahρθ, dims=2) * dθ
  end

  CUDA.@allowscalar ahρ[1] = ah[1, 1] # zeroth mode

  return ρ, ahρ
end

"""
    on_grid(func, grid)

Returns an array, of the ArrayType of the device `grid` lives on, that contains the values of
function `func` evaluated on the `grid`.
"""
on_grid(func, grid::OneDGrid) = CUDA.@allowscalar func.(grid.x)

function on_grid(func, grid::TwoDGrid)
  x, y = gridpoints(grid)
  return CUDA.@allowscalar func.(x, y)
end

function on_grid(func, grid::ThreeDGrid)
  x, y, z = gridpoints(grid)
  return CUDA.@allowscalar func.(x, y, z)
end

"""
    ArrayType(::Device)
    ArrayType(::Device, T, dim)

Returns the proper array type according to the Device chosen, i.e., `Array` for CPU and
`CuArray` for GPU.
"""
ArrayType(::CPU) = Array
ArrayType(::CPU, T, dim) = Array{T, dim}
