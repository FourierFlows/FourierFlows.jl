# Discard `effort` argument for CuArrays
plan_flows_fft(a::Array, args...; kwargs...) = plan_fft(a, args...; kwargs...)
plan_flows_rfft(a::Array, args...; kwargs...) = plan_rfft(a, args...; kwargs...)
plan_flows_fft(a::CuArray, args...; flags=nothing, kwargs...) = plan_fft(a, args...; kwargs...)
plan_flows_rfft(a::CuArray, args...; flags=nothing, kwargs...) = plan_rfft(a, args...; kwargs...)

"""
    struct OneDGrid{T<:AbstractFloat, A, Tx, Tfft, Trfft, Talias} <: AbstractGrid{T, A, Talias, D}

A one-dimensional `grid`.

$(TYPEDFIELDS)
"""
struct OneDGrid{T<:AbstractFloat, A, R, Tfft, Trfft, Talias, D} <: AbstractGrid{T, A, Talias, D}
    "device which the grid lives on"
            device :: D
    "number of points in ``x``"
                nx :: Int
    "number of wavenumbers in ``x``"
                nk :: Int
    "number of positive wavenumbers in ``x`` (real Fourier transforms)"
               nkr :: Int
    "grid spacing in ``x``"
                dx :: T
    "domain extent in ``x``"
                Lx :: T
    "range with ``x``-grid-points"
                 x :: R
    "array with ``x``-wavenumbers"
                 k :: A
    "array with positive ``x``-wavenumbers (real Fourier transforms)"
                kr :: A
    "array with inverse squared ``k``-wavenumbers, ``1 / k²``"
            invksq :: A
    "array with inverse squared ``kᵣ``-wavenumbers, ``1 / kᵣ²``"
           invkrsq :: A
    "the FFT plan for complex-valued fields"
           fftplan :: Tfft
    "the FFT plan for real-valued fields"
          rfftplan :: Trfft
    "the fraction of wavenumbers that are aliased (e.g., 1/3 for quadradic nonlinearities)"
  aliased_fraction :: T
    "range of the indices of aliased ``x``-wavenumbers"
            kalias :: Talias
    "range of the indices of aliased positive ``x``-wavenumbers (real Fourier transforms)"
           kralias :: Talias
end

"""
    OneDGrid(dev::Device = CPU();
             nx, Lx,
             x0 = -Lx/2, nthreads = Sys.CPU_THREADS, effort = FFTW.MEASURE, 
             T = Float64, aliased_fraction = 1/3)


Construct a `OneDGrid` on `dev`ice; by default on `CPU()`. Grid size is `Lx`, resolution is `nx`, 
and leftmost position is `x0`. FFT plans are generated for `nthreads` CPUs using FFTW flag `effort`.
The float type is `T`. The `aliased_fraction` keyword determines the highest wavenubers that are
being zero-ed out by `dealias!` function; 1/3 is the nominal value for quadratic nonlinearities. 
"""
function OneDGrid(dev::Device = CPU();
                  nx, Lx,
                  x0 = -Lx/2, nthreads = Sys.CPU_THREADS, effort = FFTW.MEASURE, 
                  T = Float64, aliased_fraction = 1/3)

  mod(nx, 2) != 0 && error("nx must be even")

  dx = Lx/nx

   nk = nx
  nkr = Int(nx/2 + 1)

  # Physical grid
  x = range(T(x0), step=T(dx), length=nx)

  # Wavenubmer grid
   k = device_array(dev){T}( fftfreq(nx, 2π/Lx*nx))
  kr = device_array(dev){T}(rfftfreq(nx, 2π/Lx*nx))

   invksq = @. 1 / k^2
  invkrsq = @. 1 / kr^2
  CUDA.@allowscalar  invksq[1] = 0
  CUDA.@allowscalar invkrsq[1] = 0

  FFTW.set_num_threads(nthreads)
   fftplan = plan_flows_fft(device_array(dev){Complex{T}, 1}(undef, nx), flags=effort)
  rfftplan = plan_flows_rfft(device_array(dev){T, 1}(undef, nx), flags=effort)

  kalias, kralias = getaliasedwavenumbers(nk, nkr, aliased_fraction)

  R = typeof(x)
  A = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)
  Talias = typeof(kalias)
  D = typeof(dev)
   
  return OneDGrid{T, A, R, Tfft, Trfft, Talias, D}(dev, nx, nk, nkr, dx, Lx, x, k, kr, 
                                                   invksq, invkrsq, fftplan, rfftplan,
                                                   aliased_fraction, kalias, kralias)
end


"""
    struct TwoDGrid{T<:AbstractFloat, A, Tx, Tfft, Trfft, Talias, D} <: AbstractGrid{T, Tk, Talias, D}

A two-dimensional `grid`.

$(TYPEDFIELDS)
"""
struct TwoDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft, Talias, D} <: AbstractGrid{T, Tk, Talias, D}
    "device which the grid lives on"
           device :: D
    "number of points in ``x``"
               nx :: Int
    "number of points in ``y``"
               ny :: Int
    "number of wavenumbers in ``x``"
               nk :: Int
    "number of wavenumbers in ``y``"
               nl :: Int
    "number of positive wavenumers in ``x`` (real Fourier transforms)"
               nkr :: Int
    "grid spacing in ``x``"
                dx :: T
    "grid spacing in ``y``"
                dy :: T
    "domain extent in ``x``"
                Lx :: T
    "domain extent in ``y``"
                Ly :: T
    "range with ``x``-grid-points"
                 x :: Tx
    "range with ``y``-grid-points"
                 y :: Tx
    "array with ``x``-wavenumbers"
                 k :: Tk
    "array with ``y``-wavenumbers"
                 l :: Tk
    "array with positive ``x``-wavenumbers (real Fourier transforms)"
                kr :: Tk
    "array with squared total wavenumbers, ``k² + l²``"
               Ksq :: Tk
    "array with inverse squared total wavenumbers, ``1 / (k² + l²)``"
            invKsq :: Tk
    "array with squared total wavenumbers for real Fourier transforms, ``kᵣ² + l²``"
              Krsq :: Tk
    "array with inverse squared total wavenumbers for real Fourier transforms, ``1 / (kᵣ² + l²)``"
           invKrsq :: Tk
    "the FFT plan for complex-valued fields"
           fftplan :: Tfft
    "the FFT plan for real-valued fields"
          rfftplan :: Trfft
    "the fraction of wavenumbers that are aliased (e.g., 1/3 for quadradic nonlinearities)"
  aliased_fraction :: T
    "range of the indices of aliased ``x``-wavenumbers"
            kalias :: Talias
    "range of the indices of aliased positive ``x``-wavenumbers (real Fourier transforms)"
           kralias :: Talias
    "range of the indices of aliased ``y``-wavenumbers"
            lalias :: Talias
end

"""
    TwoDGrid(dev::Device=CPU(); nx, Lx, ny=nx, Ly=Lx,
             x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE,
             T=Float64, aliased_fraction=1/3)

Construct a `TwoDGrid` on `dev`ice; by default on `CPU()`. Grid size is `Lx`, `Ly`, resolution
is `nx`, `ny`, and leftmost positions are `x0`, `y0`. FFT plans are generated for `nthreads` CPUs using
FFTW flag `effort`. The float type is `T`. The `aliased_fraction` keyword determines the highest
wavenubers that are being zero-ed out by `dealias!` function; 1/3 is the nominal value for quadratic
nonlinearities. 
"""
function TwoDGrid(dev::Device=CPU(); nx, Lx, ny=nx, Ly=Lx,
                  x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE,
                  T=Float64, aliased_fraction=1/3)

  (mod(nx, 2) != 0 || mod(ny, 2) != 0) && error("nx and ny must be even")

  dx = Lx/nx
  dy = Ly/ny

   nk = nx
   nl = ny
  nkr = Int(nx/2 + 1)

  # Physical grid
  x = range(T(x0), step=T(dx), length=nx)
  y = range(T(y0), step=T(dy), length=ny)

  # Wavenubmer grid
   k = device_array(dev){T}(reshape( fftfreq(nx, 2π/Lx*nx), (nk, 1)))
   l = device_array(dev){T}(reshape( fftfreq(ny, 2π/Ly*ny), (1, nl)))
  kr = device_array(dev){T}(reshape(rfftfreq(nx, 2π/Lx*nx), (nkr, 1)))

     Ksq = @. k^2 + l^2
  invKsq = @. 1 / Ksq
  CUDA.@allowscalar invKsq[1, 1] = 0

     Krsq = @. kr^2 + l^2
  invKrsq = @. 1 / Krsq
  CUDA.@allowscalar invKrsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan = plan_flows_fft(device_array(dev){Complex{T}, 2}(undef, nx, ny), flags=effort)
  rfftplan = plan_flows_rfft(device_array(dev){T, 2}(undef, nx, ny), flags=effort)

  kalias, kralias = getaliasedwavenumbers(nk, nkr, aliased_fraction)
  lalias, _       = getaliasedwavenumbers(nl, nl,  aliased_fraction)
  
  R = typeof(x)
  A = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)
  Talias = typeof(kalias)
  D = typeof(dev)

  return TwoDGrid{T, A, R, Tfft, Trfft, Talias, D}(dev, nx, ny, nk, nl, nkr, dx, dy, Lx, Ly, x, y, k, l, kr, 
                                                   Ksq, invKsq, Krsq, invKrsq, fftplan, rfftplan,
                                                   aliased_fraction, kalias, kralias, lalias)
end

"""
    struct ThreeDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft, Talias} <: AbstractGrid{T, Tk, Talias}

A three-dimensional `grid`.

$(TYPEDFIELDS)
"""
struct ThreeDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft, Talias, D} <: AbstractGrid{T, Tk, Talias, D}
    "device which the grid lives on"
           device :: D
    "number of points in ``x``"
               nx :: Int
    "number of points in ``y``"
               ny :: Int
    "number of points in ``z``"
               nz :: Int
    "number of wavenumbers in ``x``"
               nk :: Int
    "number of wavenumbers in ``y``"
               nl :: Int
    "number of wavenumbers in ``z``"
               nm :: Int
    "number of positive wavenumers in ``x`` (real Fourier transforms)"
               nkr :: Int
    "grid spacing in ``x``"
                dx :: T
    "grid spacing in ``y``"
                dy :: T
    "grid spacing in ``z``"
                dz :: T
    "domain extent in ``x``"
                Lx :: T
    "domain extent in ``y``"
                Ly :: T
    "domain extent in ``z``"
                Lz :: T
    "range with ``x``-grid-points"
                 x :: Tx
    "range with ``y``-grid-points"
                 y :: Tx
    "range with ``z``-grid-points"
                 z :: Tx
    "array with ``x``-wavenumbers"
                 k :: Tk
    "array with ``y``-wavenumbers"
                 l :: Tk
    "array with ``z``-wavenumbers"
                 m :: Tk
    "array with positive ``x``-wavenumbers (real Fourier transforms)"
                kr :: Tk
    "array with squared total wavenumbers, ``k² + l² + m²``"
               Ksq :: Tk
    "array with inverse squared total wavenumbers, ``1 / (k² + l² + m²)``"
            invKsq :: Tk
    "array with squared total wavenumbers for real Fourier transforms, ``kᵣ² + l² + m²``"
              Krsq :: Tk
    "array with inverse squared total wavenumbers for real Fourier transforms, ``1 / (kᵣ² + l² + m²)``"
           invKrsq :: Tk
    "the FFT plan for complex-valued fields"
           fftplan :: Tfft
    "the FFT plan for real-valued fields"
          rfftplan :: Trfft
    "the fraction of wavenumbers that are aliased (e.g., 1/3 for quadradic nonlinearities)"
  aliased_fraction :: T
    "range of the indices of aliased ``x``-wavenumbers"
            kalias :: Talias
    "range of the indices of aliased positive ``x``-wavenumbers (real Fourier transforms)"
           kralias :: Talias
    "range of the indices of aliased ``y``-wavenumbers"
            lalias :: Talias
    "range of the indices of aliased ``m``-wavenumbers"
            malias :: Talias
end

"""
    ThreeDGrid(dev::Device=CPU(); nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx,
               x0=-Lx/2, y0=-Ly/2, z0=-Lz/2,
               nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, T=Float64,
               aliased_fraction=1/3)

Construct a `ThreeDGrid` on `dev`ice; by default on `CPU()`. Grid size is `Lx`, `Ly`, `Lz`, resolution
is `nx`, `ny`, `nz` and leftmost positions are `x0`, `y0`, `z0`. FFT plans are generated for `nthreads`
CPUs using FFTW flag `effort`. The float type is `T`. The `aliased_fraction` keyword determines the
highest wavenubers that are being zero-ed out by `dealias!` function; 1/3 is the nominal value for
quadratic nonlinearities. 
"""
function ThreeDGrid(dev::Device=CPU(); nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx,
                    x0=-Lx/2, y0=-Ly/2, z0=-Lz/2,
                    nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, T=Float64,
                    aliased_fraction=1/3)

  (mod(nx, 2) != 0 || mod(ny, 2) != 0 || mod(nz, 2) != 0) && error("nx, ny, and nz must be even")

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  nk = nx
  nl = ny
  nm = nz
  nkr = Int(nx/2 + 1)

  # Physical grid
  x = range(T(x0), step=T(dx), length=nx)
  y = range(T(y0), step=T(dy), length=ny)
  z = range(T(z0), step=T(dz), length=nz)

  # Wavenubmer grid
   k = device_array(dev){T}(reshape( fftfreq(nx, 2π/Lx*nx), (nk, 1, 1)))
   l = device_array(dev){T}(reshape( fftfreq(ny, 2π/Ly*ny), (1, nl, 1)))
   m = device_array(dev){T}(reshape( fftfreq(nz, 2π/Lz*nz), (1, 1, nm)))
  kr = device_array(dev){T}(reshape(rfftfreq(nx, 2π/Lx*nx), (nkr, 1, 1)))

     Ksq = @. k^2 + l^2 + m^2
  invKsq = @. 1 / Ksq
  CUDA.@allowscalar invKsq[1, 1, 1] = 0

     Krsq = @. kr^2 + l^2 + m^2
  invKrsq = @. 1 / Krsq
  CUDA.@allowscalar invKrsq[1, 1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan = plan_flows_fft(device_array(dev){Complex{T}, 3}(undef, nx, ny, nz), flags=effort)
  rfftplan = plan_flows_rfft(device_array(dev){T, 3}(undef, nx, ny, nz), flags=effort)

  kalias, kralias = getaliasedwavenumbers(nk, nkr,         aliased_fraction)
  lalias, _       = getaliasedwavenumbers(nl, Int(nl/2+1), aliased_fraction)
  malias, _       = getaliasedwavenumbers(nm, Int(nm/2+1), aliased_fraction)
  
  R = typeof(x)
  A = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)
  Talias = typeof(kalias)
  D = typeof(dev)

  return ThreeDGrid{T, A, R, Tfft, Trfft, Talias, D}(dev, nx, ny, nz, nk, nl, nm, nkr,
                                                    dx, dy, dz, Lx, Ly, Lz, x, y, z, k, l, m, kr,
                                                    Ksq, invKsq, Krsq, invKrsq, fftplan, rfftplan, 
                                                    aliased_fraction, kalias, kralias, lalias, malias)
end

Base.eltype(grid::OneDGrid) = eltype(grid.x)
Base.eltype(grid::TwoDGrid) = eltype(grid.x)
Base.eltype(grid::ThreeDGrid) = eltype(grid.x)

"""
    gridpoints(grid::OneDDGrid)
    gridpoints(grid::TwoDGrid)
    gridpoints(grid::ThreeDGrid)

Return the collocation points of the `grid` in 1D (`X`),  2D (`X, Y`) or 3D arrays (`X, Y, Z`).
"""
function gridpoints(grid::OneDGrid{T, A}) where {T, A}
  X = [ grid.x[i] for i=1:grid.nx ]

  return A(X)
end

function gridpoints(grid::TwoDGrid{T, A}) where {T, A}
  X = [ grid.x[i₁] for i₁=1:grid.nx, i₂=1:grid.ny ]
  Y = [ grid.y[i₂] for i₁=1:grid.nx, i₂=1:grid.ny ] 
  
  return A(X), A(Y)
end

function gridpoints(grid::ThreeDGrid{T, A}) where {T, A}
  X = [ grid.x[i₁] for i₁=1:grid.nx, i₂=1:grid.ny, i₃=1:grid.nz ]
  Y = [ grid.y[i₂] for i₁=1:grid.nx, i₂=1:grid.ny, i₃=1:grid.nz ]
  Z = [ grid.z[i₃] for i₁=1:grid.nx, i₂=1:grid.ny, i₃=1:grid.nz ]
  
  return A(X), A(Y), A(Z)
end

"""
    getaliasedwavenumbers(nk, nkr, aliased_fraction)

Return the top `aliased_fraction` highest wavenumbers, both for and real FFTs, `kalias` and 
`kralias` respectively. For example, `aliased_fraction=1/3` should return the indices of the 
top-most 1/6-th (in absolute value) for both positive and negative wavenumbers (i.e., 1/3 total) 
that should be set to zero after performing an FFT. 
"""
function getaliasedwavenumbers(nk, nkr, aliased_fraction)
  L = (1 - aliased_fraction)/2 # e.g., (1 - 1/3) / 2 + 1 = 1/3.
  R = (1 + aliased_fraction)/2 # e.g., (1 + 1/3) / 2 - 1 = 2/3.
  
  iL = floor(Int, L * nk) + 1
  iR =  ceil(Int, R * nk)

  aliased_fraction < 1 || error("`aliased_fraction` must be less than 1")
 
   kalias = (aliased_fraction > 0) ? (iL:iR) : nothing
  kralias = (aliased_fraction > 0) ? (iL:nkr) : nothing

  return kalias, kralias
end

"""
    dealias!(fh, grid)

Dealias array `fh` on the `grid` based on the `grids`'s `aliased_fraction`.
"""
function dealias!(fh, grid)
   _dealias!(fh, grid)
   
   return nothing
end

dealias!(::Any, ::AbstractGrid{T, A, Nothing}) where {T, A} = nothing

function _dealias!(fh, grid::OneDGrid)
  kalias = size(fh, 1) == grid.nkr ? grid.kralias : grid.kalias
  _dealias!(fh, grid, kalias)
  
  return nothing
end

function _dealias!(fh, grid::OneDGrid, kalias)
  @views @. fh[kalias, :] = 0
  
  return nothing
end

function _dealias!(fh, grid::TwoDGrid)
  kalias = size(fh, 1) == grid.nkr ? grid.kralias : grid.kalias
  _dealias!(fh, grid, kalias)
  
  return nothing
end

function _dealias!(fh, grid::TwoDGrid, kalias)
  @views @. fh[kalias, :, :] = 0
  @views @. fh[:, grid.lalias, :] = 0
  
  return nothing
end

function _dealias!(fh, grid::ThreeDGrid)
  kalias = size(fh, 1) == grid.nkr ? grid.kralias : grid.kalias
  _dealias!(fh, grid, kalias)
  
  return nothing
end

function _dealias!(fh, grid::ThreeDGrid, kalias)
  @views @. fh[kalias, :, :, :] = 0
  @views @. fh[:, grid.lalias, :, :] = 0
  @views @. fh[:, :, grid.malias, :] = 0
  
  return nothing
end

"""
    makefilter(K; order=4, innerK=0.65, outerK=1)

Return a filter acting on the non-dimensional wavenumber `K` that decays exponentially
for `K > innerK`, thus removing high-wavenumber content from a spectrum it is multiplied
with. The decay rate is determined by order and `outerK` determines the outer wavenumber
at which the filter is smaller than `Float64` machine precision.
"""
function makefilter(K::Array; order=4, innerK=0.65, outerK=1)
  TK = typeof(K)
  K = Array(K)
  decay = 15*log(10) / (outerK - innerK)^order # decay rate for filtering function
  filter = @. exp(-decay*(K - innerK)^order)
  filter[K .< innerK] .= 1
  
  return TK(filter)
end

function makefilter(g::OneDGrid; realvars=true, kwargs...)
  K = realvars ? g.kr*g.dx/π : @.(abs(g.k*g.dx/π))
  
  return makefilter(K; kwargs...)
end

function makefilter(g::TwoDGrid; realvars=true, kwargs...)
  K = realvars ?
      @.(sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2)) : @.(sqrt((g.k*g.dx/π)^2 + (g.l*g.dy/π)^2))
  
  return makefilter(K; kwargs...)
end

function makefilter(g::ThreeDGrid; realvars=true, kwargs...)
  K = realvars ?
      @.(sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2)) : @.(sqrt((g.k*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2))
  
  return makefilter(K; kwargs...)
end

makefilter(g, T, sz; kwargs...) = ones(T, sz) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
makefilter(eq; kwargs...) = makefilter(eq.grid, fltype(eq.T), eq.dims; kwargs...)

makefilter(K::CuArray; kwargs...) = CuArray(makefilter(Array(K); kwargs...))

makefilter(g::AbstractGrid{Tg, <:CuArray}, T, sz; kwargs...) where Tg =
    CuArray(ones(T, sz)) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)


show(io::IO, g::OneDGrid{T}) where T =
     print(io, "OneDimensionalGrid\n",
               "  ├─────────── Device: ", typeof(g.device), "\n",
               "  ├──────── FloatType: $T", "\n",
               "  ├────────── size Lx: ", g.Lx, "\n",
               "  ├──── resolution nx: ", g.nx, "\n",
               "  ├── grid spacing dx: ", g.dx, "\n",
               "  ├─────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]", "\n",
               "  └─ aliased fraction: ", g.aliased_fraction)

show(io::IO, g::TwoDGrid{T}) where T =
     print(io, "TwoDimensionalGrid\n",
               "  ├───────────────── Device: ", typeof(g.device), "\n",
               "  ├────────────── FloatType: $T", "\n",
               "  ├────────── size (Lx, Ly): ", (g.Lx, g.Ly), "\n",
               "  ├──── resolution (nx, ny): ", (g.nx, g.ny), "\n",
               "  ├── grid spacing (dx, dy): ", (g.dx, g.dy), "\n",
               "  ├───────────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]", "\n",
               "  |                          y ∈ [$(g.y[1]), $(g.y[end])]", "\n",
               "  └─ aliased fraction: ", g.aliased_fraction)

show(io::IO, g::ThreeDGrid{T}) where T =
     print(io, "ThreeDimensionalGrid\n",
               "  ├───────────────────── Device: ", typeof(g.device), "\n",
               "  ├────────────────── FloatType: $T", "\n",
               "  ├────────── size (Lx, Ly, Lz): ", (g.Lx, g.Ly, g.Lz), "\n",
               "  ├──── resolution (nx, ny, nz): ", (g.nx, g.ny, g.nz), "\n",
               "  ├── grid spacing (dx, dy, dz): ", (g.dx, g.dy, g.dz), "\n",
               "  ├────────────────────  domain: x ∈ [$(g.x[1]), $(g.x[end])]", "\n",
               "  |                              y ∈ [$(g.y[1]), $(g.y[end])]", "\n",
               "  |                              z ∈ [$(g.z[1]), $(g.z[end])]", "\n",
               "  └─ aliased fraction: ", g.aliased_fraction)
