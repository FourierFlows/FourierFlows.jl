__precompile__()

module Domain

export Grid
export dealias!

# -----------------------------------------------------------------------------
# Grid ------------------------------------------------------------------------
# -----------------------------------------------------------------------------
struct Grid
  nx::Int
  ny::Int

  Lx::Float64
  Ly::Float64

  nk::Int
  nl::Int
  nkr::Int

  # Range objects that access the non-aliased part of the wavenumber range
  krange::Array{Int64, 1}
  lrange::Array{Int64, 1}
  krrange::Array{Int64, 1}

  # Dealiasing "derange" objects are arrays, not UnitRanges.
  kderange::Array{Int64, 1}
  lderange::Array{Int64, 1}
  krderange::Array{Int64, 1}

  defilt::Array{Float64, 2}
  derfilt::Array{Float64, 2}

  dx::Float64
  dy::Float64

  x::Array{Float64, 1}
  y::Array{Float64, 1}

  k::Array{Complex128, 1}
  l::Array{Complex128, 1}
  kr::Array{Complex128, 1}

  ksq::Array{Complex128, 1}
  lsq::Array{Complex128, 1}
  krsq::Array{Complex128, 1}

  ik::Array{Complex128, 1}
  il::Array{Complex128, 1}
  ikr::Array{Complex128, 1}

  X::Array{Float64, 2}
  Y::Array{Float64, 2}

  K::Array{Complex128, 2}
  L::Array{Complex128, 2}
  Kr::Array{Complex128, 2}
  Lr::Array{Complex128, 2}

  # Convenience arrays
  K2::Array{Complex128, 2}
  L2::Array{Complex128, 2}

  # k^2 + l^2 and 1/(k^2+l^2) with zero wavenumber omitted
  KKsq::Array{Complex128, 2}
  invKKsq::Array{Complex128, 2}
  KL::Array{Complex128, 2}

  KKrsq::Array{Complex128, 2}
  invKKrsq::Array{Complex128, 2}

  # FFT plans!
  fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{Float64},-1,false,2}

  ifftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.cFFTWPlan{Complex{Float64},1,false,2},Float64}

  rfftplan::Base.DFT.FFTW.rFFTWPlan{Float64,-1,false,2}

  irfftplan::Base.DFT.ScaledPlan{Complex{Float64},
    Base.DFT.FFTW.rFFTWPlan{Complex{Float64},1,false,2},Float64}

end

# Initializer for square grids
function Grid(nx::Int, Lx::Float64)

  # Size attributes
  ny = nx
  Ly = Lx
  dx = Lx/nx
  dy = Ly/ny

  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Index ranges for physical, non-dealiased wavenumbers
  kcL, kcR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  lcL, lcR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  krange  = cat(1, 1:kcL, kcR:nk)
  lrange  = cat(1, 1:lcL, lcR:nl)
  krrange = cat(1, 1:kcL)

  kderange  = (kcL+1):(kcR-1)
  lderange  = (lcL+1):(lcR-1)
  krderange = (kcL+1):nkr

  # Dealiasing filters
  defilt  = Array{Complex128}(nk, nl)
  derfilt = Array{Complex128}(nkr, nl)

  # Physical grid allocatio
  x = Array{Float64}(nx)
  y = Array{Float64}(ny)

  X = Array{Float64,2}(nx, ny)
  Y = Array{Float64,2}(nx, ny)

  # Wavenumber grid
  k  = Array{Complex128}(nk)
  l  = Array{Complex128}(nl)
  kr = Array{Complex128}(nkr)

  ksq  = Array{Complex128}(nk)
  lsq  = Array{Complex128}(nl)
  krsq = Array{Complex128}(nkr)

  K  = Array{Complex128}(nk, nl)
  L  = Array{Complex128}(nk, nl)
  Lr = Array{Complex128}(nkr, nl)
  Kr = Array{Complex128}(nkr, nl)

  K2 = Array{Complex128}(nk, nl)
  L2 = Array{Complex128}(nk, nl)

  KKsq    = Array{Complex128}(nk, nl)
  invKKsq = Array{Complex128}(nk, nl)
  KL      = Array{Complex128}(nk, nl)
  E       = Array{Complex128}(nk, nl)

  KKrsq    = Array{Complex128}(nkr, nl)
  invKKrsq = Array{Complex128}(nkr, nl)

  ik  = Array{Complex128}(nk, nl)
  il  = Array{Complex128}(nk, nl)
  ikr = Array{Complex128}(nkr, nl)

  # 1D pre-constructors
  x = linspace(-Lx/2.0, Lx/2.0-dx, nx)
  y = linspace(-Ly/2.0, Ly/2.0-dy, ny)

  i1 = 0:1:Int(nx/2)
  i2 = Int(-nx/2+1):1:-1

  j1 = 0:1:Int(ny/2)
  j2 = Int(-ny/2+1):1:-1

  k  = 2.0*pi/Lx * cat(1, i1, i2)
  l  = 2.0*pi/Ly * cat(1, j1, j2)
  kr = 2.0*pi/Lx * cat(1, i1)

  ksq  = k.^2.0
  lsq  = l.^2.0
  krsq = kr.^2.0

  ik  = im*k
  il  = im*l
  ikr = im*kr

  # Build 2D physical arrays
  for j = 1:ny, i = 1:nx
    X[i, j] = x[i]
    Y[i, j] = y[j]
  end

  # Build 2D complex spectral arrays
  for j = 1:nl, i = 1:nk
    K[i, j] = k[i]
    L[i, j] = l[j]

    K2[i, j] = k[i]^2.0
    L2[i, j] = l[j]^2.0

    KKsq[i, j] = k[i]^2.0 + l[j]^2.0
    KL[i, j] = k[i]*l[j]

    if i == 1 && j == 1
      invKKsq[i, j] = 0.0
    else
      invKKsq[i, j] = 1.0/KKsq[i, j]
    end

    if i in krange && j in lrange
        defilt[i, j] = 1.0
    else
        defilt[i, j] = 0.0
    end

  end

  # Build 2D real spectral arrays
  for j = 1:nl, i = 1:nkr
    Kr[i, j] = kr[i]
    Lr[i, j] = l[j]

    KKrsq[i, j] = kr[i]^2.0 + l[j]^2.0

    if i == 1 && j == 1
      invKKrsq[i, j] = 0.0
    else
      invKKrsq[i, j] = 1.0/KKrsq[i, j]
    end

    if i in krrange && j in lrange
      derfilt[i, j] = 1.0
    else
      derfilt[i, j] = 0.0
    end


  end

  # FFT plans; use grid vars.
  effort = FFTW.MEASURE

  fftplan   =   plan_fft(K; flags=effort)
  ifftplan  =  plan_ifft(K; flags=effort)

  rfftplan  =  plan_rfft(X; flags=effort)
  irfftplan = plan_irfft(Kr, nx; flags=effort)


  return Grid(nx, ny, Lx, Ly, nk, nl, nkr, krange, lrange, krrange,
          kderange, lderange, krderange, defilt, derfilt, dx, dy, x, y,
          k, l, kr, ksq, lsq, krsq, ik, il, ikr, X, Y, K, L, Kr, Lr,
          K2, L2, KKsq, invKKsq,KL, KKrsq, invKKrsq,
          fftplan, ifftplan, rfftplan, irfftplan)
end









# Initializer for rectangular grids
function Grid(nx::Int, ny::Int, Lx::Float64, Ly::Float64)

  # Size attributes
  dx = Lx/nx
  dy = Ly/ny

  nk = nx
  nl = ny
  nkr = Int(nx/2+1)

  # Index ranges for physical, non-dealiased wavenumbers
  kcL, kcR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  lcL, lcR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  krange  = cat(1, 1:kcL, kcR:nk)
  lrange  = cat(1, 1:lcL, lcR:nl)
  krrange = cat(1, 1:kcL)

  kderange  = (kcL+1):(kcR-1)
  lderange  = (lcL+1):(lcR-1)
  krderange = (kcL+1):nkr

  # Dealiasing filters
  defilt  = Array{Complex128}(nk, nl)
  derfilt = Array{Complex128}(nkr, nl)

  # Physical grid allocatio
  x = Array{Float64}(nx)
  y = Array{Float64}(ny)

  X = Array{Float64,2}(nx, ny)
  Y = Array{Float64,2}(nx, ny)

  # Wavenumber grid
  k  = Array{Complex128}(nk)
  l  = Array{Complex128}(nl)
  kr = Array{Complex128}(nkr)

  ksq  = Array{Complex128}(nk)
  lsq  = Array{Complex128}(nl)
  krsq = Array{Complex128}(nkr)

  K  = Array{Complex128}(nk, nl)
  L  = Array{Complex128}(nk, nl)
  Lr = Array{Complex128}(nkr, nl)
  Kr = Array{Complex128}(nkr, nl)

  K2 = Array{Complex128}(nk, nl)
  L2 = Array{Complex128}(nk, nl)

  KKsq    = Array{Complex128}(nk, nl)
  invKKsq = Array{Complex128}(nk, nl)
  KL      = Array{Complex128}(nk, nl)
  E       = Array{Complex128}(nk, nl)

  KKrsq    = Array{Complex128}(nkr, nl)
  invKKrsq = Array{Complex128}(nkr, nl)

  ik  = Array{Complex128}(nk, nl)
  il  = Array{Complex128}(nk, nl)
  ikr = Array{Complex128}(nkr, nl)

  # 1D pre-constructors
  x = linspace(-Lx/2.0, Lx/2.0-dx, nx)
  y = linspace(-Ly/2.0, Ly/2.0-dy, ny)

  i1 = 0:1:Int(nx/2)
  i2 = Int(-nx/2+1):1:-1

  j1 = 0:1:Int(ny/2)
  j2 = Int(-ny/2+1):1:-1

  k  = 2.0*pi/Lx * cat(1, i1, i2)
  l  = 2.0*pi/Ly * cat(1, j1, j2)
  kr = 2.0*pi/Lx * cat(1, i1)

  ksq  = k.^2.0
  lsq  = l.^2.0
  krsq = kr.^2.0

  ik  = im*k
  il  = im*l
  ikr = im*kr

  # Build 2D physical arrays
  for j = 1:ny, i = 1:nx
    X[i, j] = x[i]
    Y[i, j] = y[j]
  end

  # Build 2D complex spectral arrays
  for j = 1:nl, i = 1:nk
    K[i, j] = k[i]
    L[i, j] = l[j]

    K2[i, j] = k[i]^2.0
    L2[i, j] = l[j]^2.0

    KKsq[i, j] = k[i]^2.0 + l[j]^2.0
    KL[i, j] = k[i]*l[j]

    if i == 1 && j == 1
      invKKsq[i, j] = 0.0
    else
      invKKsq[i, j] = 1.0/KKsq[i, j]
    end

    if i in krange && j in lrange
        defilt[i, j] = 1.0
    else
        defilt[i, j] = 0.0
    end

  end

  # Build 2D real spectral arrays
  for j = 1:nl, i = 1:nkr
    Kr[i, j] = kr[i]
    Lr[i, j] = l[j]

    KKrsq[i, j] = kr[i]^2.0 + l[j]^2.0

    if i == 1 && j == 1
      invKKrsq[i, j] = 0.0
    else
      invKKrsq[i, j] = 1.0/KKrsq[i, j]
    end

    if i in krrange && j in lrange
      derfilt[i, j] = 1.0
    else
      derfilt[i, j] = 0.0
    end


  end

  # FFT plans; use grid vars.
  effort = FFTW.MEASURE

  fftplan   =   plan_fft(K; flags=effort)
  ifftplan  =  plan_ifft(K; flags=effort)

  rfftplan  =  plan_rfft(X; flags=effort)
  irfftplan = plan_irfft(Kr, nx; flags=effort)


  return Grid(nx, ny, Lx, Ly, nk, nl, nkr, krange, lrange, krrange,
          kderange, lderange, krderange, defilt, derfilt, dx, dy, x, y,
          k, l, kr, ksq, lsq, krsq, ik, il, ikr, X, Y, K, L, Kr, Lr,
          K2, L2, KKsq, invKKsq,KL, KKrsq, invKKrsq,
          fftplan, ifftplan, rfftplan, irfftplan)
end





function dealias!(a::Array{Complex128, 2}, g::Grid)

  if size(a)[1] == g.nk       # Complex

    for j in g.lderange
      for i in g.kderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end

  else                        # Real

    for j in g.lderange
      for i in g.krderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end

  end
end


end
