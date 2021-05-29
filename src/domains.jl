plan_flows_fft(a::Array, args...; kwargs...) = plan_fft(a, args...; kwargs...)
plan_flows_rfft(a::Array, args...; kwargs...) = plan_rfft(a, args...; kwargs...)

"""
A placeholder grid object for `0D` problems (in other words, systems of ODEs).
"""
struct ZeroDGrid{T, A} <: AbstractGrid{T, A} end

"""
A one-dimensional `grid` object.
"""
struct OneDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft} <: AbstractGrid{T, Tk}
        nx :: Int
        nk :: Int
       nkr :: Int

        dx :: T
        Lx :: T

         x :: Tx
         
         k :: Tk
        kr :: Tk
    invksq :: Tk
   invkrsq :: Tk

   fftplan :: Tfft
  rfftplan :: Trfft

  # Range objects that access the aliased part of the wavenumber range
    kalias :: UnitRange{Int}
   kralias :: UnitRange{Int}
end

"""
    OneDGrid(nx, Lx; x0=-Lx/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, 
                      T=Float64, dealias=1/3, ArrayType=Array)

Constructs a OneDGrid object with size `Lx`, resolution `nx`, and leftmost position `x0`. 
FFT plans are generated for `nthreads` CPUs using FFTW flag `effort`. The float type is `T` 
and the array types is `ArrayType`. The `dealias` keyword determines the highest  wavenubers 
that are being zero-ed out by `dealias()` function; 1/3 is the nominal value for quadratic 
nonlinearities. 
"""
function OneDGrid(nx, Lx; x0=-Lx/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, 
                  T=Float64, dealias=1/3, ArrayType=Array)

  dx = Lx/nx
    
   nk = nx
  nkr = Int(nx/2 + 1)

  # Physical grid
  x = range(T(x0), step=T(dx), length=nx)

  # Wavenubmer grid
   k = ArrayType{T}( fftfreq(nx, 2π/Lx*nx))
  kr = ArrayType{T}(rfftfreq(nx, 2π/Lx*nx))

   invksq = @. 1 / k^2
  invkrsq = @. 1 / kr^2
  CUDA.@allowscalar  invksq[1] = 0
  CUDA.@allowscalar invkrsq[1] = 0

  FFTW.set_num_threads(nthreads)
   fftplan = plan_flows_fft(ArrayType{Complex{T}, 1}(undef, nx), flags=effort)
  rfftplan = plan_flows_rfft(ArrayType{T, 1}(undef, nx), flags=effort)

  kalias, kralias = getaliasedwavenumbers(nk, nkr, dealias)

  Tx = typeof(x)
  Tk = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)

  return OneDGrid{T, Tk, Tx, Tfft, Trfft}(nx, nk, nkr, dx, Lx, x, k, kr, 
                                      invksq, invkrsq, fftplan, rfftplan, kalias, kralias)
end

"A two-dimensional `grid` object."
struct TwoDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft} <: AbstractGrid{T, Tk}
        nx :: Int
        ny :: Int
        nk :: Int
        nl :: Int
       nkr :: Int

        dx :: T
        dy :: T
        Lx :: T
        Ly :: T

         x :: Tx
         y :: Tx
         
         k :: Tk
         l :: Tk
        kr :: Tk
       Ksq :: Tk
    invKsq :: Tk
      Krsq :: Tk
   invKrsq :: Tk

   fftplan :: Tfft
  rfftplan :: Trfft

  # Range objects that access the aliased part of the wavenumber range
    kalias :: UnitRange{Int}
   kralias :: UnitRange{Int}
    lalias :: UnitRange{Int}
end

"""
    TwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, 
                      effort=FFTW.MEASURE, T=Float64, dealias=1/3, ArrayType=Array)

Constructs a TwoDGrid object with size `Lx`, `Ly`, resolution `nx`, `ny`, and leftmost
positions `x0`, `y0`. FFT plans are generated for `nthreads` CPUs using
FFTW flag `effort`. The float type is `T` and the array types is `ArrayType`. 
The `dealias` keyword determines the highest wavenubers that are being zero-ed
out by `dealias()` function; 1/3 is the nominal value for quadratic nonlinearities. 
"""
function TwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, 
                  effort=FFTW.MEASURE, T=Float64, dealias=1/3, ArrayType=Array)

  dx = Lx/nx
  dy = Ly/ny

   nk = nx
   nl = ny
  nkr = Int(nx/2 + 1)

  # Physical grid
  x = range(T(x0), step=T(dx), length=nx)
  y = range(T(y0), step=T(dy), length=ny)

  # Wavenubmer grid
   k = ArrayType{T}(reshape( fftfreq(nx, 2π/Lx*nx), (nk, 1)))
   l = ArrayType{T}(reshape( fftfreq(ny, 2π/Ly*ny), (1, nl)))
  kr = ArrayType{T}(reshape(rfftfreq(nx, 2π/Lx*nx), (nkr, 1)))

     Ksq = @. k^2 + l^2
  invKsq = @. 1 / Ksq
  CUDA.@allowscalar invKsq[1, 1] = 0

     Krsq = @. kr^2 + l^2
  invKrsq = @. 1 / Krsq
  CUDA.@allowscalar invKrsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan = plan_flows_fft(ArrayType{Complex{T}, 2}(undef, nx, ny), flags=effort)
  rfftplan = plan_flows_rfft(ArrayType{T, 2}(undef, nx, ny), flags=effort)

  # Index endpoints for aliasfrac i, j wavenumbers
  kalias, kralias = getaliasedwavenumbers(nk, nkr, dealias)
  lalias, _ = getaliasedwavenumbers(nl, nl, dealias)
  
  Tx = typeof(x)
  Tk = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)

  return TwoDGrid{T, Tk, Tx, Tfft, Trfft}(nx, ny, nk, nl, nkr, dx, dy, Lx, Ly, x, y, k, l, kr, Ksq, invKsq, Krsq, invKrsq,
                                      fftplan, rfftplan, kalias, kralias, lalias)
end

"A three-dimensional `grid` object."
struct ThreeDGrid{T<:AbstractFloat, Tk, Tx, Tfft, Trfft} <: AbstractGrid{T, Tk}
        nx :: Int
        ny :: Int
        nz :: Int
        nk :: Int
        nl :: Int
        nm :: Int
       nkr :: Int

        dx :: T
        dy :: T
        dz :: T
        Lx :: T
        Ly :: T
        Lz :: T

         x :: Tx
         y :: Tx
         z :: Tx
         
         k :: Tk
         l :: Tk
         m :: Tk
        kr :: Tk
       Ksq :: Tk
    invKsq :: Tk
      Krsq :: Tk
   invKrsq :: Tk

   fftplan :: Tfft
  rfftplan :: Trfft

  # Range objects that access the aliased part of the wavenumber range
    kalias :: UnitRange{Int}
   kralias :: UnitRange{Int}
    lalias :: UnitRange{Int}
    malias :: UnitRange{Int}
end

"""
    ThreeDGrid(nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx; x0=-Lx/2, y0=-Ly/2, z0=-Lz/2, 
    nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, T=Float64, dealias=1/3, ArrayType=Array)

 Constructs a TwoDGrid object with size `Lx`, `Ly`, `Lz`, resolution `nx`, `ny`,
 `nz` and leftmost positions `x0`, `y0`, `z0`. FFT plans are generated for `nthreads` 
 CPUs using FFTW flag `effort`. The float type is `T` and the array types is `ArrayType`. 
 The `dealias` keyword determines the highest wavenubers that are being zero-ed
 out by `dealias()` function; 1/3 is the nominal value for quadratic nonlinearities. 
"""
function ThreeDGrid(nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx; x0=-Lx/2, y0=-Ly/2, z0=-Lz/2,
                  nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, T=Float64, dealias=1/3, ArrayType=Array)

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
   k = ArrayType{T}(reshape( fftfreq(nx, 2π/Lx*nx), (nk, 1, 1)))
   l = ArrayType{T}(reshape( fftfreq(ny, 2π/Ly*ny), (1, nl, 1)))
   m = ArrayType{T}(reshape( fftfreq(nz, 2π/Lz*nz), (1, 1, nm)))
  kr = ArrayType{T}(reshape(rfftfreq(nx, 2π/Lx*nx), (nkr, 1, 1)))

     Ksq = @. k^2 + l^2 + m^2
  invKsq = @. 1 / Ksq
  CUDA.@allowscalar invKsq[1, 1, 1] = 0

     Krsq = @. kr^2 + l^2 + m^2
  invKrsq = @. 1 / Krsq
  CUDA.@allowscalar invKrsq[1, 1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan = plan_flows_fft(ArrayType{Complex{T}, 3}(undef, nx, ny, nz), flags=effort)
  rfftplan = plan_flows_rfft(ArrayType{T, 3}(undef, nx, ny, nz), flags=effort)

  # Index endpoints for aliasfrac i, j wavenumbers
  kalias, kralias = getaliasedwavenumbers(nk, nkr, dealias)
  lalias, malias = getaliasedwavenumbers(nl, nm, dealias)
  
  Tx = typeof(x)
  Tk = typeof(k)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)

  return ThreeDGrid{T, Tk, Tx, Tfft, Trfft}(nx, ny, nz, nk, nl, nm, nkr, dx, dy, dz, Lx, Ly, Lz,
                                        x, y, z, k, l, m, kr, Ksq, invKsq, Krsq, invKrsq, fftplan, rfftplan, 
                                        kalias, kralias, lalias, malias)
end

Base.eltype(grid::OneDGrid) = eltype(grid.x)
Base.eltype(grid::TwoDGrid) = eltype(grid.x)
Base.eltype(grid::ThreeDGrid) = eltype(grid.x)

OneDGrid(dev::CPU, args...; kwargs...) = OneDGrid(args...; ArrayType=Array, kwargs...)
TwoDGrid(dev::CPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=Array, kwargs...)
ThreeDGrid(dev::CPU, args...; kwargs...) = ThreeDGrid(args...; ArrayType=Array, kwargs...)

"""
    gridpoints(grid::TwoDGrid)
    gridpoints(grid::ThreeDGrid)

Returns the collocation points of the `grid` in 2D or 3D arrays `X, Y` (and `Z`).
"""
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
    getaliasedwavenumbers(nk, nkr, aliasfraction)

Returns the top `aliasfraction` highest wavenumbers.
"""
function getaliasedwavenumbers(nk, nkr, aliasfraction)
  # Index endpoints for aliased i, j wavenumbers
  # 1/3 aliasfraction => upper 1/6 of +/- wavenumbers (1/3 total) are set to 0 after performing fft.
  # 1/2 aliasfraction => upper 1/4 of +/- wavenumbers (1/2 total) are set to 0 after performing fft.
  L = (1 - aliasfraction)/2 # (1 - 1/3) / 2 + 1 = 1/3.
  R = (1 + aliasfraction)/2 # (1 + 1/3) / 2 - 1 = 2/3.
  iL = floor(Int, L * nk) + 1
  iR =  ceil(Int, R * nk)

  aliasfraction < 1 || error("`aliasfraction` must be less than 1") # aliasfraction=1 is not sensible.
   kalias = (aliasfraction > 0) ? (iL:iR) : Int(nx/2 + 1)
  kralias = (aliasfraction > 0) ? (iL:nkr) : nkr

  return kalias, kralias
end

"""
    dealias!(a, grid, kalias)

Dealias array `a` on the `grid` with aliased x-wavenumbers `kalias`.
"""
function dealias!(a, grid::OneDGrid)
  kalias = size(a, 1) == grid.nkr ? grid.kralias : grid.kalias
  dealias!(a, grid, kalias)
  return nothing
end

function dealias!(a, grid::OneDGrid, kalias)
  @views @. a[kalias, :] = 0
  
  return nothing
end

function dealias!(a, grid::TwoDGrid)
  kalias = size(a, 1) == grid.nkr ? grid.kralias : grid.kalias
  dealias!(a, grid, kalias)
  
  return nothing
end

function dealias!(a, grid::TwoDGrid, kalias)
  @views @. a[kalias, :, :] = 0
  @views @. a[:, grid.lalias, :] = 0
  
  return nothing
end

function dealias!(a, grid::ThreeDGrid)
  kalias = size(a, 1) == grid.nkr ? grid.kralias : grid.kalias
  dealias!(a, grid, kalias)
  
  return nothing
end

function dealias!(a, grid::ThreeDGrid, kalias)
  @views @. a[kalias, :, :, :] = 0
  @views @. a[:, grid.lalias, :, :] = 0
  @views @. a[:, :, grid.malias, :] = 0
  
  return nothing
end

"""
    makefilter(K; order=4, innerK=0.65, outerK=1)

Returns a filter acting on the non-dimensional wavenumber K that decays exponentially
for K>innerK, thus removing high-wavenumber content from a spectrum it is multiplied
with. The decay rate is determined by order and outerK determines the outer wavenumber
at which the filter is smaller than Float64 machine precision.
"""
function makefilter(K::Array; order=4, innerK=0.65, outerK=1)
  TK = typeof(K)
  K = Array(K)
  decay = 15*log(10) / (outerK-innerK)^order # decay rate for filtering function
  filt = @. exp( -decay*(K-innerK)^order )
  filt[real.(K) .< innerK] .= 1
  
  return TK(filt)
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

"""
    griddevice(grid)

Returns the device on which the `grid` lives on.
"""
griddevice(grid::AbstractGrid{T, A}) where {T, A} = A<:Array ? "CPU" : "GPU"

show(io::IO, g::OneDGrid{T}) where T =
     print(io, "OneDimensionalGrid\n",
               "  ├─────────── Device: ", griddevice(g), '\n',
               "  ├──────── FloatType: $T", '\n', 
               "  ├────────── size Lx: ", g.Lx, '\n',
               "  ├──── resolution nx: ", g.nx, '\n',
               "  ├── grid spacing dx: ", g.dx, '\n',
               "  └─────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]")

show(io::IO, g::TwoDGrid{T}) where T =
     print(io, "TwoDimensionalGrid\n",
               "  ├───────────────── Device: ", griddevice(g), '\n',
               "  ├────────────── FloatType: $T", '\n', 
               "  ├────────── size (Lx, Ly): ", (g.Lx, g.Ly), '\n',
               "  ├──── resolution (nx, ny): ", (g.nx, g.ny), '\n',
               "  ├── grid spacing (dx, dy): ", (g.dx, g.dy), '\n',
               "  └───────────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]", '\n',
               "                             y ∈ [$(g.y[1]), $(g.y[end])]")

show(io::IO, g::ThreeDGrid{T}) where T =
     print(io, "ThreeDimensionalGrid\n",
               "  ├───────────────────── Device: ", griddevice(g), '\n',
               "  ├────────────────── FloatType: $T", '\n', 
               "  ├────────── size (Lx, Ly, Lz): ", (g.Lx, g.Ly, g.Ly), '\n',
               "  ├──── resolution (nx, ny, nz): ", (g.nx, g.ny, g.nz), '\n',
               "  ├── grid spacing (dx, dy, dz): ", (g.dx, g.dy, g.dz), '\n',
               "  └────────────────────  domain: x ∈ [$(g.x[1]), $(g.x[end])]", '\n',
               "                                 y ∈ [$(g.y[1]), $(g.y[end])]", '\n',
               "                                 z ∈ [$(g.z[1]), $(g.z[end])]")
