export OneDGrid, TwoDGrid, dealias!

abstract type AbstractTwoDGrid <: AbstractGrid end
abstract type AbstractOneDGrid <: AbstractGrid end

"""
    ZeroDGrid(nvars)

Constructs a placeholder grid object for "0D" problems (in other words, systems of
ODEs) with nvars variables.
"""
struct ZeroDGrid <: AbstractGrid
  nvars::Int
end

"""
    OneDGrid(nx, Lx)

    OneDGrid(nx, Lx; x0=-Lx/2, nthreads=Sys.CPU_CORES, effort=FFTW.MEASURE)

Constrcut a OneDGrid object. The one-dimensional domain has size Lx, 
resolution nx, and leftmost grid point x0. FFT plans are generated 
which use nthreads threads with the specified planning effort. 
"""
struct OneDGrid{T} <: AbstractOneDGrid
  nx::Int
  nk::Int
  nkr::Int
  Lx::T
  dx::T
  x::Array{T,1}
  k::Array{T,1}
  kr::Array{T,1}
  invksq::Array{T,1}
  invkrsq::Array{T,1}

  fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{T},-1,false,1}
  ifftplan::Base.DFT.ScaledPlan{Complex{T},Base.DFT.FFTW.cFFTWPlan{Complex{T},1,false,1},T}
  rfftplan::Base.DFT.FFTW.rFFTWPlan{T,-1,false,1}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},Base.DFT.FFTW.rFFTWPlan{Complex{T},1,false, 1},T}

  # Range objects that access the aliased part of the wavenumber range
  ialias::UnitRange{Int}
  iralias::UnitRange{Int}
end

function OneDGrid(nx, Lx; x0=-0.5*Lx, nthreads=Sys.CPU_CORES, effort=FFTW.MEASURE)
  T = typeof(Lx)
  dx = Lx/nx
  x = Array{T}(linspace(x0, x0+Lx-dx, nx))

  nk = nx
  nkr = Int(nx/2+1)

  i1 = 0:Int(nx/2)
  i2 = Int(-nx/2+1):-1
  k = Array{T}(2π/Lx*cat(1, i1, i2))
  kr = Array{T}(2π/Lx*cat(1, i1))

  invksq = @. 1/k^2
  invksq[1] = 0
  invkrsq = @. 1/kr^2
  invkrsq[1] = 0

  FFTW.set_num_threads(nthreads)
  fftplan   = plan_fft(Array{Complex{T},1}(nx); flags=effort)
  ifftplan  = plan_ifft(Array{Complex{T},1}(nk); flags=effort)
  rfftplan  = plan_rfft(Array{T,1}(nx); flags=effort)
  irfftplan = plan_irfft(Array{Complex{T},1}(nkr), nx; flags=effort)

  # Index endpoints for aliased i, j wavenumbers
  iaL, iaR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  ialias  = iaL:iaR
  iralias = iaL:nkr

  OneDGrid(nx, nk, nkr, Lx, dx, x, k, kr, invksq, invkrsq,
           fftplan, ifftplan, rfftplan, irfftplan, ialias, iralias)
end

"""
    TwoDGrid(nx, Lx)
    TwoDGrid(nx, Lx, ny, Ly; x0=-Lx/2, y0=-Ly/2, nthreads=Sys.CPU_CORES, effort=FFTW.MEASURE)

Constrcut a TwoDGrid object. The two-dimensional domain has size (Lx, Ly), 
resolution (nx, ny) and bottom left corner at (x0, y0). FFT plans are generated 
which use nthreads threads with the specified planning effort. 
"""
struct TwoDGrid{T} <: AbstractTwoDGrid
  nx::Int
  ny::Int
  nk::Int
  nl::Int
  nkr::Int
  Lx::T
  Ly::T
  dx::T
  dy::T
  x::Array{T,2}
  y::Array{T,2}
  X::Array{T,2}
  Y::Array{T,2}
  k::Array{T,2}
  l::Array{T,2}
  kr::Array{T,2}
  K::Array{T,2}
  L::Array{T,2}
  Kr::Array{T,2}
  Lr::Array{T,2}
  KKsq::Array{T,2}      # K^2 + L^2
  invKKsq::Array{T,2}   # 1/KKsq, invKKsq[1, 1]=0
  KKrsq::Array{T,2}     # Kr^2 + Lr^2
  invKKrsq::Array{T,2}  # 1/KKrsq, invKKrsq[1, 1]=0

  fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{T},-1,false,2}
  ifftplan::Base.DFT.ScaledPlan{Complex{T},Base.DFT.FFTW.cFFTWPlan{Complex{T},1,false,2},T}
  rfftplan::Base.DFT.FFTW.rFFTWPlan{T,-1,false,2}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},Base.DFT.FFTW.rFFTWPlan{Complex{T},1,false,2},T}

  # Range objects that access the aliased part of the wavenumber range
  ialias::UnitRange{Int}
  iralias::UnitRange{Int}
  jalias::UnitRange{Int}
end

function TwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-0.5*Lx, y0=-0.5*Ly, nthreads=Sys.CPU_CORES, effort=FFTW.MEASURE)
  T = typeof(Lx)
  dx = Lx/nx
  dy = Ly/ny
  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Physical grid
  x = Array{T}(reshape(linspace(x0, x0+Lx-dx, nx), (nx, 1)))
  y = Array{T}(reshape(linspace(y0, y0+Ly-dy, ny), (1, ny)))
  X = [ x[i] for i = 1:nx, j = 1:ny]
  Y = [ y[j] for i = 1:nx, j = 1:ny]

  # Wavenubmer grid
  i1 = 0:Int(nx/2)
  i2 = Int(-nx/2+1):-1
  j1 = 0:Int(ny/2)
  j2 = Int(-ny/2+1):-1

  k  = reshape(2π/Lx*cat(1, i1, i2), (nk, 1))
  kr = reshape(2π/Lx*cat(1, i1), (nkr, 1))
  l  = reshape(2π/Ly*cat(1, j1, j2), (1, nl))

  K = [ k[i] for i = 1:nk, j = 1:nl]
  L = [ l[j] for i = 1:nk, j = 1:nl]
  Kr = [ kr[i] for i = 1:nkr, j = 1:nl]
  Lr = [ l[j]  for i = 1:nkr, j = 1:nl]

  KKsq  = @. K^2 + L^2
  invKKsq = 1./KKsq
  invKKsq[1, 1] = 0

  KKrsq = @. Kr^2 + Lr^2
  invKKrsq = 1./KKrsq
  invKKrsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan   = plan_fft(Array{Complex{T},2}(nx, ny); flags=effort)
  ifftplan  = plan_ifft(Array{Complex{T},2}(nk, nl); flags=effort)
  rfftplan  = plan_rfft(Array{T,2}(nx, ny); flags=effort)
  irfftplan = plan_irfft(Array{Complex{T},2}(nkr, nl), nx; flags=effort)

  # Index endpoints for aliased i, j wavenumbers
  iaL, iaR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  jaL, jaR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  ialias  = iaL:iaR
  iralias = iaL:nkr
  jalias  = iaL:iaR

  TwoDGrid(nx, ny, nk, nl, nkr, Lx, Ly, dx, dy, x, y, X, Y,
           k, l, kr, K, L, Kr, Lr, KKsq, invKKsq, KKrsq, invKKrsq,
           fftplan, ifftplan, rfftplan, irfftplan, ialias, iralias, jalias)
end

function dealias!(a, g::OneDGrid)
  if size(a)[1] == g.nkr; @views @. a[g.iralias, :] = 0
  else;                   @views @. a[g.ialias, :] = 0
  end
  nothing
end

function dealias!(a, g::AbstractTwoDGrid)
  if size(a)[1] == g.nkr; @views @. a[g.iralias, g.jalias, :] = 0
  else;                   @views @. a[g.ialias, g.jalias, :] = 0
  end
  nothing
end


"""
    makefilter(K; order=4, innerK=0.65, outerK=1)

Returns a filter acting on the non-dimensional wavenumber K that decays exponentially
for K>innerK, thus removing high-wavenumber content from a spectrum it is multiplied with. 
The decay rate is determined by order and outerK determines the outer wavenumber at which 
the filter is smaller than machine precision.
"""
function makefilter(K; order=4, innerK=0.65, outerK=1)
  decay = 15*log(10) / (outerK-innerK)^order # decay rate for filtering function
  filt = @. exp( -decay*(KK-innerK)^order )
  filt[real.(KK) .< innerK] .= 1
  filt
end

function makefilter(g::AbstractTwoDGrid; realvars=true, kwargs...)
  K = realvars ? 
      @.(sqrt((g.Kr*g.dx/π)^2 + (g.Lr*g.dy/π)^2)) : @.(sqrt((g.K*g.dx/π)^2 + (g.L*g.dy/π)^2))
  makefilter(K; kwargs...)
end

function makefilter(g::AbstractOneDGrid; realvars=true, kwargs...)
  K = realvars ? g.kr*g.dx/π : @.(abs(g.k*g.dx/π))
  makefilter(K; kwargs...)
end

makefilter(g, T, sz; kwargs...) = ones(T, sz).*makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
