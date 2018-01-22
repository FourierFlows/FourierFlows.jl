export TwoDGrid, dealias!


"""
Doc.
"""
struct ZeroDGrid <: AbstractGrid
  nvars::Int
end


"""
Doc.
"""
struct TwoDGrid <: AbstractGrid
  nx::Int
  ny::Int
  nk::Int
  nl::Int
  nkr::Int

  Lx::Float64
  Ly::Float64
  dx::Float64
  dy::Float64

  x::Array{Float64,2}
  y::Array{Float64,2}
  X::Array{Float64,2}
  Y::Array{Float64,2}

  k::Array{Float64,2}
  l::Array{Float64,2}
  kr::Array{Float64,2}

  K::Array{Float64,2}
  L::Array{Float64,2}
  Kr::Array{Float64,2}
  Lr::Array{Float64,2}

  # k^2 + l^2 and 1/(k^2+l^2) with zero wavenumber omitted
  KKsq::Array{Float64,2}
  invKKsq::Array{Float64,2}
  KKrsq::Array{Float64,2}
  invKKrsq::Array{Float64,2}

  # FFT plans
  fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{Float64}, -1, false, 2}

  ifftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.cFFTWPlan{Complex{Float64}, 1, false, 2}, Float64}

  rfftplan::Base.DFT.FFTW.rFFTWPlan{Float64, -1, false, 2}

  irfftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.rFFTWPlan{Complex{Float64}, 1, false, 2}, Float64}

  # Range objects that access the non-aliased part of the wavenumber range
  ialias::UnitRange{Int64}
  iralias::UnitRange{Int64}
  jalias::UnitRange{Int64}
end


"""
Construct a rectangular grid.
"""
function TwoDGrid(nx::Int, Lx::Float64, ny::Int=nx, Ly::Float64=Lx;
  x0=-0.5*Lx, y0=-0.5*Ly, nthreads=Sys.CPU_CORES, effort=FFTW.MEASURE)

  # Size attributes
  dx = Lx/nx
  dy = Ly/ny

  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Physical grid
  x = reshape(linspace(x0, x0+Lx-dx, nx), (nx, 1))
  y = reshape(linspace(y0, y0+Ly-dy, ny), (1, ny))

  X = [ x[i] for i = 1:nx, j = 1:ny]
  Y = [ y[j] for i = 1:nx, j = 1:ny]

  # Wavenubmer grid
  i1 = 0:Int(nx/2)
  i2 = Int(-nx/2+1):-1

  j1 = 0:Int(ny/2)
  j2 = Int(-ny/2+1):-1

  k  = reshape(2π/Lx * cat(1, i1, i2), (nk, 1))
  kr = reshape(2π/Lx * cat(1, i1), (nkr, 1))
  l  = reshape(2π/Ly * cat(1, j1, j2), (1, nl))

  K = [ k[i] for i = 1:nk, j = 1:nl]
  L = [ l[j] for i = 1:nk, j = 1:nl]

  Kr = [ kr[i] for i = 1:nkr, j = 1:nl]
  Lr = [ l[j]  for i = 1:nkr, j = 1:nl]

  KKsq  = K.^2 + L.^2
  invKKsq = 1./KKsq
  invKKsq[1, 1] = 0

  KKrsq = Kr.^2 + Lr.^2
  invKKrsq = 1./KKrsq
  invKKrsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  fftplan   = plan_fft(Array{Complex{Float64},2}(nx, ny); flags=effort)
  ifftplan  = plan_ifft(Array{Complex{Float64},2}(nk, nl); flags=effort)
  rfftplan  = plan_rfft(Array{Float64,2}(nx, ny); flags=effort)
  irfftplan = plan_irfft(Array{Complex{Float64},2}(nkr, nl), nx; flags=effort)

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

# Grid constructor for tupled arguments
function TwoDGrid(nxy::Tuple{Int, Int}, Lxy::Tuple{Float64, Float64};
  nthreads=Sys.CPU_CORES)
  nx, ny = nxy
  Lx, Ly = Lxy
  TwoDGrid(nx, Lx, ny, Ly; nthreads=nthreads)
end

function dealias!(a::Array{Complex{Float64},2}, g)
  if size(a)[1] == g.nkr
    a[g.iralias, g.jalias] = 0
  else
    a[g.ialias, g.jalias] = 0
  end
  nothing
end

function dealias!(a::Array{Complex{Float64},3}, g)
  if size(a)[1] == g.nkr
    @views @. a[g.iralias, g.jalias, :] = 0
  else
    @views @. a[g.ialias, g.jalias, :] = 0
  end
  nothing
end

"""
Returns an filter with an exponentially-decaying profile that, when multiplied
removes high-wavenumber content from a spectrum.
"""
function makefilter(g::TwoDGrid; order=4.0, innerK=0.65, outerK=1.0,
  realvars=true)

  # Get decay rate for filter
  decay = 15.0*log(10.0) / (outerK-innerK)^order

  # Non-dimensional square wavenumbers
  if realvars
    KK = sqrt.( (g.Kr*g.dx/π).^2 + (g.Lr*g.dy/π).^2 )
  else
    KK = sqrt.( (g.K*g.dx/π).^2  + (g.L*g.dy/π).^2  )
  end

  filt = exp.( -decay*(KK-innerK).^order )
  filt[ real.(KK) .< innerK ] = 1

  return filt
end

function makefilter(g, T, sz; kwargs...)
  if sz[1] == g.nkr; realvars = true
  else;              realvars = false
  end
  ones(T, sz) .* makefilter(g; realvars=realvars, kwargs...)
end
