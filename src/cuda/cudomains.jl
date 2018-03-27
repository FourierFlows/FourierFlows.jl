"""
    CuTwoDGrid(nx, Lx)
    CuTwoDGrid(nx, Lx, ny, Ly;  x0=-Lx/2, y0=-Ly/2)

Constrcut a CuTwoDGrid object. The two-dimensional domain has size (Lx, Ly), 
resolution (nx, ny) and bottom left corner at (x0, y0).
"""
struct CuTwoDGrid{T} <: AbstractTwoDGrid
  nx::Int
  ny::Int
  nk::Int
  nl::Int
  nkr::Int
  Lx::T
  Ly::T
  dx::T
  dy::T
  x::CuArray{T,2}
  y::CuArray{T,2}
  X::CuArray{T,2}
  Y::CuArray{T,2}
  k::CuArray{T,2}
  l::CuArray{T,2}
  kr::CuArray{T,2}
  K::CuArray{T,2}
  L::CuArray{T,2}
  Kr::CuArray{T,2}
  Lr::CuArray{T,2}
  KKsq::CuArray{T,2}      # K^2 + L^2
  invKKsq::CuArray{T,2}   # 1/KKsq, invKKsq[1, 1]=0
  KKrsq::CuArray{T,2}     # Kr^2 + Lr^2
  invKKrsq::CuArray{T,2}  # 1/KKrsq, invKKrsq[1, 1]=0

  fftplan::CuArrays.FFT.cCuFFTPlan{Complex{T},-1,false,2}
  ifftplan::Base.DFT.ScaledPlan{Complex{T},CuArrays.FFT.cCuFFTPlan{Complex{T},1,false,2},T}
  rfftplan::CuArrays.FFT.rCuFFTPlan{T,-1,false,2}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},CuCuArrays.FFT.rCuFFTPlan{Complex{T},1,false,2},T}

  # Range objects that access the aliased part of the wavenumber range
  ialias::UnitRange{Int}
  iralias::UnitRange{Int}
  jalias::UnitRange{Int}
end

function CuTwoDGrid(nx, Lx, ny=nx, Ly=Lx; x0=-0.5*Lx, y0=-0.5*Ly)
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

  # Convert Arrays to CuArrays
  @cuconvertarrays x y X Y k kr l K L Kr Lr KKsq invKKsq KKrsq invKKrsq

  # FFT plans
    fftplan = plan_fft(CuArray{Complex{T},2}(nx, ny); flags=effort)
   ifftplan = plan_ifft(CuArray{Complex{T},2}(nk, nl); flags=effort)
   rfftplan = plan_rfft(CuArray{T,2}(nx, ny); flags=effort)
  irfftplan = plan_irfft(CuArray{Complex{T},2}(nkr, nl), nx; flags=effort)

  # Index endpoints for aliased i, j wavenumbers
  iaL, iaR = Int(floor(nk/3))+1, 2*Int(ceil(nk/3))-1
  jaL, jaR = Int(floor(nl/3))+1, 2*Int(ceil(nl/3))-1

  ialias  = iaL:iaR
  iralias = iaL:nkr
  jalias  = iaL:iaR

  CuTwoDGrid(nx, ny, nk, nl, nkr, Lx, Ly, dx, dy, x, y, X, Y,
             k, l, kr, K, L, Kr, Lr, KKsq, invKKsq, KKrsq, invKKrsq,
             fftplan, ifftplan, rfftplan, irfftplan, ialias, iralias, jalias)
end

CuTwoDGrid(T::DataType, nx, Lx, args...; kwargs...) = CuTwoDGrid(nx, T(Lx), args...; kwargs...)
CuTwoDGrid(g::TwoDGrid) = CuTwoDGrid(g.nx, g.Lx, g.ny, g.Ly; x0=g.x[1], y0=g.y[1])
