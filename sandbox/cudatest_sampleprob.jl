#= twodturb.jl

Solves the 2D vorticity equation using a Fourier pseudospectral method
in a doubly-periodic box. Time-stepping is Forward Euler. Does not 
dealias.
=#

using CuArrays, BenchmarkTools

macro cuconvertarrays(vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = CuArray($(esc(var)));) for var in vars])
  expr 
end

struct Vars{T}
  k::Array{T,2}
  l::Array{T,2}
  Ksq::Array{T,2}
  invKsq::Array{T,2}
  q::Array{T,2}
  u::Array{T,2}
  v::Array{T,2}
  qh::Array{Complex{T},2}
  qsh::Array{Complex{T},2}
  rhs::Array{Complex{T},2}
  uh::Array{Complex{T},2}
  vh::Array{Complex{T},2}
  rfftplan::Base.DFT.FFTW.rFFTWPlan{T,-1,false,2}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},Base.DFT.FFTW.rFFTWPlan{Complex{T},1,false,2},T}
end

struct CuVars{T}
  k::CuArray{T,2}
  l::CuArray{T,2}
  Ksq::CuArray{T,2}
  invKsq::CuArray{T,2}
  q::CuArray{T,2}
  u::CuArray{T,2}
  v::CuArray{T,2}
  qh::CuArray{Complex{T},2}
  qsh::CuArray{Complex{T},2}
  rhs::CuArray{Complex{T},2}
  uh::CuArray{Complex{T},2}
  vh::CuArray{Complex{T},2}
  rfftplan::CuArrays.FFT.rCuFFTPlan{T,-1,false,2}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},CuArrays.FFT.rCuFFTPlan{Complex{T},1,false,2},T}
end

function Vars(T, nx, Lx; usegpu=false)
  # Construct the grid
  dx = Lx/nx
  nk, nl = Int(nx/2+1), nx
  k = Array{T,2}(reshape(2π/Lx*(0:Int(nx/2)), (nk, 1)))
  l = Array{T,2}(reshape(2π/Lx*cat(1, 0:nl/2, -nl/2+1:-1), (1, nl)))

  Ksq = @. k^2 + l^2
  invKsq = @. 1/Ksq
  invKsq[1, 1] = 0.0  # Will eliminate 0th mode during inversion

  # Random initial condition
  q = rand(nx, nx)
  qh = rfft(q)

  # Preallocate
  u = zeros(T, nx, nx)
  v = zeros(T, nx, nx)

  qsh  = zeros(Complex{T}, nk, nl)
  rhs  = zeros(Complex{T}, nk, nl)
  uh   = zeros(Complex{T}, nk, nl)
  vh   = zeros(Complex{T}, nk, nl)

  if usegpu
    @cuconvertarrays(k, l, Ksq, invKsq, q, u, v, qh, qsh, rhs, uh, vh)
    rfftplan  = plan_rfft(deepcopy(u))
    irfftplan = plan_irfft(deepcopy(uh), nx)
    return CuVars(k, l, Ksq, invKsq, q, u, v, qh, qsh, rhs, uh, vh, rfftplan, irfftplan)
  else
    FFTW.set_num_threads(Sys.CPU_CORES)
    rfftplan = plan_rfft(deepcopy(u); flags=FFTW.MEASURE)
    irfftplan = plan_irfft(deepcopy(uh), nx; flags=FFTW.MEASURE)
    return Vars(k, l, Ksq, invKsq, q, u, v, qh, qsh, rhs, uh, vh, rfftplan, irfftplan)
  end
end

function speedtest(vs, nu, dt, nsteps; dmsg=1000)
  k         = vs.k
  l         = vs.l
  q         = vs.q
  u         = vs.u
  v         = vs.v
  Ksq       = vs.Ksq
  invKsq    = vs.invKsq
  qh        = vs.qh
  uh        = vs.uh
  vh        = vs.vh
  qsh       = vs.qsh
  rhs       = vs.rhs
  rfftplan  = vs.rfftplan
  irfftplan = vs.irfftplan
  
  # Step forward
  t = 0.0
  for step = 1:nsteps
    if step % dmsg == 0 && step > 1
      @printf("step: %04d, t: %6.1f\n", step, t)
    end

    # Calculate right hand side of vorticity equation.
    #@. uh =  im * l / (k^2 + l^2) * qh
    #@. vh = -im * k / (k^2 + l^2) * qh
    #uh[1, 1] = 0
    #vh[1, 1] = 0
    @. uh =  im * l * invKsq * qh
    @. vh = -im * k * invKsq * qh
    qsh .= qh  # Necessary because irfft destroys its input.

    A_mul_B!(q, irfftplan, qsh)
    A_mul_B!(u, irfftplan, uh)
    A_mul_B!(v, irfftplan, vh)

    @. u *= q
    @. v *= q

    A_mul_B!(uh, rfftplan, u)
    A_mul_B!(vh, rfftplan, v)

    @. rhs = -im*k*uh - im*l*vh - nu*Ksq*qh

    # Step forward
    @. qh += dt*rhs
    t += dt
  end
  A_mul_B!(q, irfftplan, qh)
  q
end

# Parameters
nx = 128      # Resolution 
Lx = 2π       # Physical box size
nu = 8e-5     # Laplacian viscosity
dt = 1e-2     # Time step
nsteps = 100  # Number of time steps
T = Float64   # Computational precision

vcpu = Vars(T, nx, Lx; usegpu=false)
vgpu = Vars(T, nx, Lx; usegpu=true)
qc = speedtest(vcpu, nu, dt, 1)
qg = speedtest(vgpu, nu, dt, 1)

for nx in [64, 128, 256, 512, 1024, 2048, 4096]
  vcpu = Vars(T, nx, Lx; usegpu=false)
  vgpu = Vars(T, nx, Lx; usegpu=true)

  @printf "n=%d, CPU" nx
  @btime qc = speedtest(vcpu, nu, dt, nsteps)

  @printf "n=%d, GPU" nx
  @btime qg = speedtest(vgpu, nu, dt, nsteps)
end
