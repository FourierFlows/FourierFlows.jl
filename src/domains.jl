plan_flows_fft(a::Array, effort) = plan_fft(a; flags=effort)
plan_flows_rfft(a::Array, effort) = plan_rfft(a; flags=effort)

"""
    ZeroDGrid()

Constructs a placeholder grid object for "0D" problems (in other words, systems of ODEs).
"""
struct ZeroDGrid{T, Ta} <: AbstractGrid{T, Ta} end

function getaliasedwavenumbers(nk, nkr, aliasfraction)
  # Index endpoints for aliased i, j wavenumbers
  # 1/3 aliasfraction => upper 1/6 of +/- wavenumbers (1/3 total) are set to 0 after performing fft.
  # 1/2 aliasfraction => upper 1/4 of +/- wavenumbers (1/2 total) are set to 0 after performing fft.
  L = (1 - aliasfraction)/2 # (1 - 1/3) / 2 + 1 = 1/3.
  R = (1 + aliasfraction)/2 # (1 + 1/3) / 2 - 1 = 2/3.
  iL = floor(Int, L*nk) + 1
  iR =  ceil(Int, R*nk)

  aliasfraction < 1 || error("`aliasfraction` must be less than 1") # aliasfraction=1 is not sensible.
   kalias = (aliasfraction > 0) ? (iL:iR) : Int(nx/2+1)
  kralias = (aliasfraction > 0) ? (iL:nkr) : nkr

  kalias, kralias
end

"""
    OneDGrid(nx, Lx; x0=-Lx/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE)

Constructs a OneDGrid object with size `Lx`, resolution `nx`, and leftmost
position `x0`. FFT plans are generated for `nthreads` CPUs using
FFTW flag `effort`.
"""
struct OneDGrid{T<:AbstractFloat, Ta<:AbstractArray, Tfft, Trfft} <: AbstractGrid{T, Ta}
        nx :: Int
        nk :: Int
       nkr :: Int

        dx :: T
        Lx :: T

         x :: Ta
         k :: Ta
        kr :: Ta
    invksq :: Ta
   invkrsq :: Ta

   fftplan :: Tfft
  rfftplan :: Trfft

  # Range objects that access the aliased part of the wavenumber range
    kalias :: UnitRange{Int}
   kralias :: UnitRange{Int}
end

function OneDGrid(nx, Lx; x0=-Lx/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE, T=Float64, dealias=1/3,
                  ArrayType=Array)

  dx = Lx/nx
  x = ArrayType{T}(range(x0, step=dx, length=nx))

  nk = nx
  nkr = Int(nx/2+1)

  i₁ = 0:Int(nx/2)
  i₂ = Int(-nx/2+1):-1
   k = ArrayType{T}(2π/Lx*cat(i₁, i₂; dims=1))
  kr = ArrayType{T}(2π/Lx*cat(i₁; dims=1))

   invksq = @. 1/k^2
  invkrsq = @. 1/kr^2
   invksq[1] = 0
  invkrsq[1] = 0

  FFTW.set_num_threads(nthreads)
   fftplan = plan_flows_fft(ArrayType{Complex{T}, 1}(undef, nx), effort)
  rfftplan = plan_flows_rfft(ArrayType{T, 1}(undef, nx), effort)

  kalias, kralias = getaliasedwavenumbers(nk, nkr, dealias)

  Ta = typeof(x)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)

  OneDGrid{T, Ta, Tfft, Trfft}(nx, nk, nkr, dx, Lx, x, k, kr, 
                               invksq, invkrsq, fftplan, rfftplan, kalias, kralias)
end

"""
    TwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, effort=FFTW.MEASURE)

Constructs a TwoDGrid object.
"""
struct TwoDGrid{T<:AbstractFloat, Ta<:AbstractArray, Tfft, Trfft} <: AbstractGrid{T, Ta}
  nx::Int
  ny::Int
  nk::Int
  nl::Int
  nkr::Int

  dx::T
  dy::T
  Lx::T
  Ly::T

  x::Ta
  y::Ta
  k::Ta
  l::Ta
  kr::Ta
  Ksq::Ta
  invKsq::Ta
  Krsq::Ta
  invKrsq::Ta

  fftplan::Tfft
  rfftplan::Trfft

  # Range objects that access the aliased part of the wavenumber range
  kalias::UnitRange{Int}
  kralias::UnitRange{Int}
  lalias::UnitRange{Int}
end

function TwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_THREADS, 
                  effort=FFTW.MEASURE, T=Float64, dealias=1/3, ArrayType=Array)
                  
  dx = Lx/nx
  dy = Ly/ny

  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Physical grid
  x = ArrayType{T}(reshape(range(x0, step=dx, length=nx), (nx, 1)))
  y = ArrayType{T}(reshape(range(y0, step=dy, length=ny), (1, ny)))

  # Wavenubmer grid
  i₁ = 0:Int(nx/2)
  i₂ = Int(-nx/2+1):-1
  j₁ = 0:Int(ny/2)
  j₂ = Int(-ny/2+1):-1

   k = ArrayType{T}(reshape(2π/Lx*cat(i₁, i₂, dims=1), (nk, 1)))
   l = ArrayType{T}(reshape(2π/Ly*cat(j₁, j₂, dims=1), (1, nl)))
  kr = ArrayType{T}(reshape(2π/Lx*cat(i₁, dims=1), (nkr, 1)))

     Ksq = @. k^2 + l^2
  invKsq = @. 1/Ksq
  invKsq[1, 1] = 0

     Krsq = @. kr^2 + l^2
  invKrsq = @. 1/Krsq
  invKrsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan = plan_flows_fft(ArrayType{Complex{T}, 2}(undef, nx, ny), effort)
  rfftplan = plan_flows_rfft(ArrayType{T, 2}(undef, nx, ny), effort)

  # Index endpoints for aliasfrac i, j wavenumbers
  kalias, kralias = getaliasedwavenumbers(nk, nkr, dealias)
  lalias, _ = getaliasedwavenumbers(nl, nl, dealias)
  
  Ta = typeof(x)
  Tfft = typeof(fftplan)
  Trfft = typeof(rfftplan)

  TwoDGrid{T, Ta, Tfft, Trfft}(nx, ny, nk, nl, nkr, dx, dy, Lx, Ly, x, y, k, l, kr, Ksq, invKsq, Krsq, invKrsq,
           fftplan, rfftplan, kalias, kralias, lalias)
end

OneDGrid(dev::CPU, args...; kwargs...) = OneDGrid(args...; ArrayType=Array, kwargs...)
TwoDGrid(dev::CPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=Array, kwargs...)

"""
    gridpoints(g)

Returns the collocation points of the grid `g` in 2D arrays `X, Y`.
"""
function gridpoints(g::AbstractGrid{T, A}) where {T, A}
  X = [ g.x[i] for i=1:g.nx, j=1:g.ny]
  Y = [ g.y[j] for i=1:g.nx, j=1:g.ny]
  A(X), A(Y)
end

"""
    dealias!(a, g, kalias)

Dealias `a` on the grid `g` with aliased x-wavenumbers `kalias`.
"""
function dealias!(a, g::OneDGrid)
  kalias = size(a, 1) == g.nkr ? g.kralias : g.kalias
  dealias!(a, g, kalias)
  nothing
end

function dealias!(a, g::OneDGrid, kalias)
  @views @. a[kalias, :] = 0
  nothing
end

function dealias!(a, g::TwoDGrid)
  kalias = size(a, 1) == g.nkr ? g.kralias : g.kalias
  dealias!(a, g, kalias)
  nothing
end

function dealias!(a, g::TwoDGrid, kalias)
  @views @. a[kalias, g.lalias, :] = 0
  nothing
end

"""
    makefilter(K; order=4, innerK=0.65, outerK=1)

Returns a filter acting on the non-dimensional wavenumber K that decays exponentially
for K>innerK, thus removing high-wavenumber content from a spectrum it is multiplied with.
The decay rate is determined by order and outerK determines the outer wavenumber at which
the filter is smaller than Float64 machine precision.
"""
function makefilter(K::Array; order=4, innerK=0.65, outerK=1)
  TK = typeof(K)
  K = Array(K)
  decay = 15*log(10) / (outerK-innerK)^order # decay rate for filtering function
  filt = @. exp( -decay*(K-innerK)^order )
  filt[real.(K) .< innerK] .= 1
  TK(filt)
end

function makefilter(g::TwoDGrid; realvars=true, kwargs...)
  K = realvars ?
      @.(sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2)) : @.(sqrt((g.k*g.dx/π)^2 + (g.l*g.dy/π)^2))
  makefilter(K; kwargs...)
end

function makefilter(g::OneDGrid; realvars=true, kwargs...)
  K = realvars ? g.kr*g.dx/π : @.(abs(g.k*g.dx/π))
  makefilter(K; kwargs...)
end

makefilter(g, T, sz; kwargs...) = ones(T, sz) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
makefilter(eq) = makefilter(eq.grid, fltype(eq.T), eq.dims)
