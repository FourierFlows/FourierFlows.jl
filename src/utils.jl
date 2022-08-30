"""
    innereltype(x)

Recursively determine the 'innermost' type in by the collection `x` (which may be, for example,
a collection of a collection).
"""
function innereltype(x)
  T = eltype(x)

  return T <: AbstractArray ? innereltype(T) : return T
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
superzeros(T, A::AbstractArray) = T(0) * A
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

Return `Σ |uh|²` on the `grid`, which is equal to the domain integral of `u`. More specifically, 
it returns
```math
\\sum_{𝐤} |û_{𝐤}|² L_x L_y = \\int u(𝐱)² \\, 𝖽x 𝖽y \\,,
```
where ``û_{𝐤} =`` `uh` ``/(`` `grid.nx` ``e^{- i 𝐤 ⋅ 𝐱₀})``, with ``𝐱₀`` the vector with components
the left-most position in each direction.
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

Return `real(Σ uh)` on the `grid`, i.e.
```math
ℜ [ \\sum_{𝐤} û_{𝐤} L_x L_y ] \\,,
```
where ``û_{𝐤} =`` `uh` ``/(`` `grid.nx` ``e^{- i 𝐤 ⋅ 𝐱₀})``, with ``𝐱₀`` the vector with components
the left-most position in each direction.
"""
function parsevalsum(uh, grid::TwoDGrid)
  if size(uh, 1) == grid.nkr    # uh is conjugate symmetric
    U = sum(uh[1, :])           # k=0 modes
    U += 2*sum(uh[2:end, :])    # sum k>0 modes twice
  else # count every mode once
    U = sum(uh)
  end

  norm = grid.Lx * grid.Ly / (grid.nx^2 * grid.ny^2) # normalization for dft

  return norm * real(U)
end

"""
    jacobianh(a, b, grid)

Return the Fourier transform of the Jacobian of `a` and `b` on `grid`.
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

Return the Jacobian of `a` and `b` on `grid`.
"""
function jacobian(a, b, grid::TwoDGrid)
  if eltype(a) <: Real
    return irfft(jacobianh(a, b, grid), grid.nx)
  else
    return ifft(jacobianh(a, b, grid))
  end
end

"""
    radialspectrum(fh, grid; n=nothing, m=nothing, refinement=2)

Return the radial spectrum of `fh`. `fh` lives on Cartesian wavenumber grid ``(k, l)``. To 
compute the radial spectrum, we first interpolate ``f̂(k, l)`` onto a radial wavenumber grid 
``(ρ, θ)``, where ``ρ² = k²+l²`` and ``θ = \\arctan(l/k)``. Note here that 
``f̂ =`` `fh` ``/(`` `grid.nx` ``e^{- i 𝐤 ⋅ 𝐱₀})``,  with ``𝐱₀`` the vector with components the 
left-most position in each direction. After interpolation, we integrate ``f̂``over angles ``θ`` 
to get `fρ`,

```math
  f̂_ρ = \\int f̂(ρ, θ) ρ 𝖽ρ 𝖽θ \\, .
```

The resolution `(n, m)` for the polar wavenumber grid is `n = refinement * maximum(nk, nl), 
m = refinement * maximum(nk, nl)`, where `refinement = 2` by default. If `fh` is in conjugate 
symmetric form then only the upper-half plane in ``θ`` is represented on the polar grid.
"""
function radialspectrum(fh, grid::TwoDGrid; n=nothing, m=nothing, refinement=2)

  n = n === nothing ? refinement * maximum([grid.nk, grid.nl]) : n
  m = m === nothing ? refinement * maximum([grid.nk, grid.nl]) : m

  # Calculate the shifted k and l
  lshift = range(-grid.nl/2, stop=grid.nl/2-1, length=grid.nl) * 2π/grid.Ly

  if size(fh)[1] == grid.nkr # conjugate symmetric form
    m = Int(m/2)                         # => half resolution in θ
    θ = range(-π/2, stop=π/2, length=m)  # θ-grid from k=0 to max(kr)
    fhshift = fftshift(fh, 2)            # shifted fh
    kshift = range(0, stop=grid.nkr-1, length=grid.nkr) * 2π/grid.Lx
  else # ordinary form
    θ = range(0, stop=2π, length=m)      # θ grid
    fhshift = fftshift(fh, [1, 2])       # shifted fh
    kshift = range(-grid.nk/2, stop=grid.nk/2-1, length=grid.nk) * 2π/grid.Lx
  end

  # Interpolator for fh
  itp = scale(interpolate(fhshift, BSpline(Linear())), kshift, lshift)

  # Get radial wavenumber vector
  ρmax = minimum([(grid.nk/2-1) * 2π/grid.Lx, (grid.nl/2-1) * 2π/grid.Ly])
  ρ = range(0, stop=ρmax, length=n)

  # Interpolate fh onto fine grid in (ρ, θ).
  fhρθ = zeros(eltype(fhshift), (n, m))

  for i₂=1:m, i₁=2:n # ignore zeroth mode; i₁≥2
    fhρθ[i₁, i₂] = itp(ρ[i₁] * cos(θ[i₂]), ρ[i₁] * sin(θ[i₂]))
  end

  # fhρ = ρ ∫ fh(ρ, θ) dθ  =>  Fh = ∫ fhρ dρ = ∫∫ fh dk dl
  dθ = θ[2] - θ[1]
  if size(fh)[1] == grid.nkr
    fhρ = 2ρ .* sum(fhρθ, dims=2) * dθ # multiply by 2 for conjugate symmetry
  else
    fhρ =  ρ .* sum(fhρθ, dims=2) * dθ
  end

  CUDA.@allowscalar fhρ[1] = fh[1, 1] # zeroth mode

  return ρ, fhρ
end

"""
    on_grid(func, grid)

Return an array, of the type compatible with the `device` of that the `grid` lives on,
that contains the values of function `func` evaluated on the `grid`.
"""
function on_grid(func, grid::OneDGrid{T}) where T
  f = zeros(grid.device, T, (grid.nx, ))

  @. f = func(grid.x)

  return f
end

function on_grid(func, grid::TwoDGrid{T}) where T
  f = zeros(grid.device, T, (grid.nx, grid.ny))

  x = reshape(grid.x, (grid.nx, 1,))
  y = reshape(grid.y, (1, grid.ny))
  
  @. f = func(x, y)

  return f
end

function on_grid(func, grid::ThreeDGrid{T}) where T
  f = zeros(grid.device, T, (grid.nx, grid.ny, grid.nz))

  x = reshape(grid.x, (grid.nx, 1, 1))
  y = reshape(grid.y, (1, grid.ny, 1))
  z = reshape(grid.z, (1, 1, grid.nz))
  
  @. f = func(x, y, z)

  return f
end

"""
    device_array(device::Device)
    device_array(device::Device, T, dim)

Return the proper array type according to the `device`, i.e., `Array` for CPU and
`CuArray` for GPU.
"""
device_array(device::CPU) = Array
device_array(device::CPU, T, dim) = Array{T, dim}
