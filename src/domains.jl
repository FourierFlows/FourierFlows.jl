__precompile__()

export TwoDGrid, dealias!, cubicdealias!
# Grid types and constructors

type TwoDGrid <: AbstractGrid
  nx::Int
  ny::Int

  Lx::Float64
  Ly::Float64

  nk::Int
  nl::Int
  nkr::Int

  dx::Float64
  dy::Float64

  x::Array{Float64, 1}
  y::Array{Float64, 1}

  # Range objects that access the non-aliased part of the wavenumber range
  ialias::UnitRange{Int64}
  iralias::UnitRange{Int64}
  jalias::UnitRange{Int64}

  ialias3::UnitRange{Int64}
  iralias3::UnitRange{Int64}
  jalias3::UnitRange{Int64}

  k::Array{Complex{Float64}, 1}
  l::Array{Complex{Float64}, 1}
  kr::Array{Complex{Float64}, 1}

  ksq::Array{Complex{Float64}, 1}
  lsq::Array{Complex{Float64}, 1}
  krsq::Array{Complex{Float64}, 1}

  ik::Array{Complex{Float64}, 1}
  il::Array{Complex{Float64}, 1}
  ikr::Array{Complex{Float64}, 1}

  X::Array{Float64, 2}
  Y::Array{Float64, 2}

  K::Array{Complex{Float64}, 2}
  L::Array{Complex{Float64}, 2}
  Kr::Array{Complex{Float64}, 2}
  Lr::Array{Complex{Float64}, 2}

  # Convenience arrays
  K2::Array{Complex{Float64}, 2}
  L2::Array{Complex{Float64}, 2}
  Kr2::Array{Complex{Float64}, 2}
  Lr2::Array{Complex{Float64}, 2}

  # k^2 + l^2 and 1/(k^2+l^2) with zero wavenumber omitted
  KKsq::Array{Complex{Float64}, 2}
  invKKsq::Array{Complex{Float64}, 2}
  KL::Array{Complex{Float64}, 2}
  KLr::Array{Complex{Float64}, 2}

  KKrsq::Array{Complex{Float64}, 2}
  invKKrsq::Array{Complex{Float64}, 2}

  # FFT plans!
  fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{Float64}, -1, false, 2}

  ifftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.cFFTWPlan{Complex{Float64}, 1, false, 2}, Float64}

  rfftplan::Base.DFT.FFTW.rFFTWPlan{Float64, -1, false, 2}

  irfftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.rFFTWPlan{Complex{Float64}, 1, false, 2}, Float64}
end


# Initializer for rectangular grids
function TwoDGrid(nx::Int, Lx::Float64, ny::Int=nx, Ly::Float64=Lx;
  x0=-0.5*Lx, y0=-0.5*Ly, nthreads=Sys.CPU_CORES)

  # Size attributes
  dx = Lx/nx
  dy = Ly/ny

  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Physical grid allocatio
  x = Array{Float64}(nx)
  y = Array{Float64}(ny)

  X = Array{Float64,2}(nx, ny)
  Y = Array{Float64,2}(nx, ny)

  # Wavenumber grid
  k  = Array{Complex{Float64}}(nk)
  l  = Array{Complex{Float64}}(nl)
  kr = Array{Complex{Float64}}(nkr)

  ksq  = Array{Complex{Float64}}(nk)
  lsq  = Array{Complex{Float64}}(nl)
  krsq = Array{Complex{Float64}}(nkr)

  K  = Array{Complex{Float64}}(nk, nl)
  L  = Array{Complex{Float64}}(nk, nl)
  Lr = Array{Complex{Float64}}(nkr, nl)
  Kr = Array{Complex{Float64}}(nkr, nl)

  K2 = Array{Complex{Float64}}(nk, nl)
  L2 = Array{Complex{Float64}}(nk, nl)
  Kr2 = Array{Complex{Float64}}(nkr, nl)
  Lr2 = Array{Complex{Float64}}(nkr, nl)

  KKsq    = Array{Complex{Float64}}(nk, nl)
  invKKsq = Array{Complex{Float64}}(nk, nl)
  KL      = Array{Complex{Float64}}(nk, nl)
  KLr     = Array{Complex{Float64}}(nkr, nl)

  KKrsq    = Array{Complex{Float64}}(nkr, nl)
  invKKrsq = Array{Complex{Float64}}(nkr, nl)

  ik  = Array{Complex{Float64}}(nk, nl)
  il  = Array{Complex{Float64}}(nk, nl)
  ikr = Array{Complex{Float64}}(nkr, nl)

  # 1D pre-constructors
  x = linspace(x0, x0+Lx-dx, nx)
  y = linspace(y0, y0+Ly-dy, ny)

  i1 = 0:1:Int(nx/2)
  i2 = Int(-nx/2+1):1:-1

  j1 = 0:1:Int(ny/2)
  j2 = Int(-ny/2+1):1:-1

  k  = 2π/Lx * cat(1, i1, i2)
  kr = 2π/Lx * cat(1, i1)
  l  = 2π/Ly * cat(1, j1, j2)

  ksq  = k.^2.0
  krsq = kr.^2.0
  lsq  = l.^2.0

  ik  = im*k
  ikr = im*kr
  il  = im*l

  X = [ x[i] for i = 1:nx, j = 1:ny]
  Y = [ y[j] for i = 1:nx, j = 1:ny]

  K = [ k[i] for i = 1:nk, j = 1:nl]
  L = [ l[j] for i = 1:nk, j = 1:nl]

  Kr = [ kr[i] for i = 1:nkr, j = 1:nl]
  Lr = [ l[j]  for i = 1:nkr, j = 1:nl]

  K2 = K.^2.0
  L2 = L.^2.0
  Kr2 = Kr.^2.0
  Lr2 = Lr.^2.0

  KKsq  = K.^2.0 + L.^2.0
  KL   = K.*L
  KLr  = Kr.*Lr
  invKKsq = 1.0./KKsq
  invKKsq[1, 1] = 0.0

  KKrsq = Kr.^2.0 + Lr.^2.0
  invKKrsq = 1.0./KKrsq
  invKKrsq[1, 1] = 0.0

  # Index endpoints for aliased i, j wavenumbers
  iaL, iaR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  jaL, jaR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  ialias  = iaL:iaR
  iralias = iaL:nkr
  jalias  = iaL:iaR

  ia3L, ia3R = Int(floor(nk/4))+1, 2*Int(ceil(nk/4))-1
  ja3L, ja3R = Int(floor(nl/4))+1, 2*Int(ceil(nl/4))-1

  ialias3 = ia3L:ia3R
  iralias3 = ia3L:nkr
  jalias3 = ja3L:ja3R


  # FFT plans; use grid vars.
  FFTW.set_num_threads(nthreads)
  effort = FFTW.MEASURE

  fftplan   = plan_fft(Array{Float64,2}(nx, ny); flags=effort)
  ifftplan  = plan_ifft(Array{Complex{Float64},2}(nk, nl); flags=effort)

  rfftplan  = plan_rfft(Array{Float64,2}(nx, ny); flags=effort)
  irfftplan = plan_irfft(Array{Complex{Float64},2}(nkr, nl), nx; flags=effort)

  return TwoDGrid(nx, ny, Lx, Ly, nk, nl, nkr, dx, dy, x, y,
          ialias, iralias, jalias, ialias3, iralias3, jalias3,
          k, l, kr, ksq, lsq, krsq, ik, il, ikr, X, Y, K, L, Kr, Lr,
          K2, L2, Kr2, Lr2, KKsq, invKKsq, KL, KLr, KKrsq, invKKrsq,
          fftplan, ifftplan, rfftplan, irfftplan)
end

# Grid constructor for tupled arguments
function TwoDGrid(nxy::Tuple{Int, Int}, Lxy::Tuple{Float64, Float64};
  nthreads=Sys.CPU_CORES)
  nx, ny = nxy
  Lx, Ly = Lxy
  TwoDGrid(nx, Lx, ny, Ly)
end

function dealias!(a::Array{Complex{Float64}, 2}, g)
  if size(a)[1] == g.nkr
    a[g.iralias, g.jalias] = 0im
  else
    a[g.ialias, g.jalias] = 0im
  end
  nothing
end

function dealias!(a::Array{Complex{Float64}, 3}, g)
  if size(a)[1] == g.nkr
    @views @. a[g.iralias, g.jalias, :] = 0im
  else
    @views @. a[g.ialias, g.jalias, :] = 0im
  end
  nothing
end


function cubicdealias!(a::Array{Complex{Float64}, 2}, g)
  if size(a)[1] == g.nkr
    a[g.iralias3, g.jalias3] = 0im
  else
    a[g.ialias3, g.jalias3] = 0im
  end
  nothing
end

function cubicdealias!(a::Array{Complex{Float64}, 3}, g)
  if size(a)[1] == g.nkr
    @views @. a[g.iralias3, g.jalias3, :] = 0im
  else
    @views @. a[g.ialias3, g.jalias3, :] = 0im
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
