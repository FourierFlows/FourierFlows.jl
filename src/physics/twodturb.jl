module TwoDTurb

using FourierFlows
export InitialValueProblem, Params, Vars, Equation, set_q!, updatevars!

# Problem ---------------------------------------------------------------------
"""
Construct an initial value problem.
"""
function InitialValueProblem(;
   nx = 256,
   Lx = 2π,
   ny = nothing,
   Ly = nothing,
    ν = nothing,
   nu = nothing,
   nν = 2,
  nnu = nothing,
   dt = 0.01,
  withfilter = false
  )

  # Defaults
  if  nu != nothing;  ν = nu;  end
  if nnu != nothing; nν = nnu; end
  if  Ly == nothing; Ly = Lx;  end
  if  ny == nothing; ny = nx;  end

  if ν == nothing
    if withfilter; ν = 0.0
    else;          ν = 1e-1/(dt*(0.65π*nx/Lx)^nν)
    end
  end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.Params(ν, nν)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)

  if withfilter; ts = FilteredETDRK4TimeStepper(dt, eq.LC, g)
  else;          ts = ETDRK4TimeStepper(dt, eq.LC)
  end

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

function InitialValueProblem(n, L, ν, nν, dt, withfilter)
  InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, withfilter=withfilter)
end

function InitialValueProblem(nx, Lx, ny, Ly, ν, nν, dt, withfilter)
  InitialValueProblem(nx=nx, Lx=Lx, ny=ny, Ly=Ly, ν=ν, nν=nν, dt=dt,
                        withfilter=withfilter)
end

# Params
struct Params <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
end

struct ForcedParams <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
  μ::Float64        # Bottom drag
  calcF!::Function  # Function that calculates the forcing F
end

# Equations
function Equation(p::Params, g::TwoDGrid)
  LC = -p.ν * g.KKrsq.^(0.5*p.nν)
  FourierFlows.Equation{2}(LC, calcN!)
end

function Equation(p::ForcedParams, g::TwoDGrid)
  LC = -p.ν*g.KKrsq.^p.nν - μ
  FourierFlows.Equation{2}(LC, calcN!)
end




# Vars
physvars = [:q, :U, :V, :Uq, :Vq, :psi]
transvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih]

expr = FourierFlows.getexpr_varstype(:Vars, physvars, transvars)
eval(expr)

function Vars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih
  Vars(0.0, sol, q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih)
end


forcedtransvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih, :F]
expr = FourierFlows.getexpr_varstype(:ForcedVars, physvars, forcedtransvars)
eval(expr)

function ForcedVars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih F
  Vars(0.0, sol, q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih, F)
end




# Solvers
function calcN_advection!(N::Array{Complex{Float64}, 2}, 
                sol::Array{Complex{Float64}, 2},
                t::Float64, v::AbstractVars, p::AbstractParams, g::TwoDGrid)
  v.qh .= sol
  A_mul_B!(v.q, g.irfftplan, v.qh) # destroys qh when using fftw

  @. v.Uh =  im * g.l  * g.invKKrsq * sol
  @. v.Vh = -im * g.kr * g.invKKrsq * sol

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.Uq = v.U * v.q
  @. v.Vq = v.V * v.q

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  @. N = -im*g.kr*v.Uqh - im*g.l*v.Vqh
  nothing
end

function calcN!(N::Array{Complex{Float64}, 2}, 
                sol::Array{Complex{Float64}, 2}, t::Float64, 
                v::Vars, p::Params, g::TwoDGrid)
  calcN_advection!(N, sol, t, v, p, g)
  nothing
end

function calcN!(N::Array{Complex{Float64}, 2}, 
                sol::Array{Complex{Float64}, 2}, t::Float64, 
                v::ForcedVars, p::ForcedParams, g::TwoDGrid)

  calcN_advection!(N, sol, t, v, p, g)
  p.calcF!(v.F, sol, t, v, p, g)

  @. N += v.F
  nothing
end


# Helper functions
"""
Update solution variables using qh in v.sol.
"""
function updatevars!(v, g)
  v.qh .= v.sol
  @. v.psih = -g.invKKrsq * v.qh
  @. v.Uh = -im*g.l  * v.psih
  @. v.Vh =  im*g.kr * v.psih

  qh1 = deepcopy(v.qh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  A_mul_B!(v.q, g.irfftplan, qh1)
  A_mul_B!(v.U, g.irfftplan, Uh1)
  A_mul_B!(v.V, g.irfftplan, Vh1)

  nothing
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.grid)
end


"""
Set the vorticity field.
"""
function set_q!(v, g, q)
  A_mul_B!(v.sol, g.rfftplan, q)
  updatevars!(v, g)
end

function set_q!(prob::AbstractProblem, q)
  set_q!(prob.vars, prob.grid, q)
end


"""
Calculate the domain integrated kinetic energy.
"""
function energy(v, g)
  0.5*(FourierFlows.parsevalsum2(g.Kr.*g.invKKrsq.*v.sol, g)
        + FourierFlows.parsevalsum2(g.Lr.*g.invKKrsq.*v.sol, g))
end

function energy(prob)
  energy(prob.vars, prob.grid)
end


"""
Returns the domain-integrated enstrophy.
"""
function enstrophy(v, g)
  0.5*FourierFlows.parsevalsum2(v.sol, g)
end

function enstrophy(prob)
  enstrophy(prob.vars, prob.grid)
end




""" Make a field of mature turbulence on a square grid.

  Args:
    nx: grid resolution
    Lx: grid extent
    qf: final maximum vorticity
    q0: initial maximum vorticity
    nν: order of hyperviscosity
    maxsteps: maximum νmber of steps to take
    dt: time step
    ν: hyperviscosity
    k0: initial waveνmber
    E0: initial energy
    tf: final time
    plots: whether or not to plot field evolution

  Returns
    q: The vorticity field
"""
function makematureturb(nx::Int, Lx::Real; qf=0.1, q0=0.2, nν=4,
  maxsteps=10000, dt=nothing, ν=nothing, k0=nx/2,
  E0=nothing, tf=nothing, plots=false, loginterval=5)

  g  = TwoDGrid(nx, Lx)
  vs = TwoDTurb.Vars(g)

  if E0 != nothing # set initial energy rather than vorticity

    # Closely following the formulation in Rocha, Wagner, Young
    modk = sqrt(g.KKsq)

    psik = zeros(g.nk, g.nl)
    psik =  (modk .* (1 + (modk/k0).^4)).^(-0.5)
    psik[1, 1] = 0.0
    C = real(sqrt(E0/sum(g.KKsq.*abs2.(psik))))

    psi = zeros(g.nx, g.ny)
    for i = 1:128
      for j = 1:128
        psi .+= real.(C*psik[i, j]*cos.(
          g.k[i]*g.X + g.l[j]*g.Y + 2*pi*rand(1)[1]))
      end
    end

    psih = rfft(psi)
    qi = -irfft(g.KKrsq.*psih, g.nx)
    set_q!(vs, g, qi)
    E0 = FourierFlows.parsevalsum(g.KKrsq.*abs2.(psih), g)

  else
    qi = FourierFlows.peaked_isotropic_spectrum(nx, k0; maxval=q0)
    set_q!(vs, g, qi)
    E0 = energy(vs, g)
  end

  maxq = q0 = maximum(abs.(vs.q))

  # Defaults
  if dt == nothing; dt = 0.1*g.dx/maximum([vs.U; vs.V]);  end
  if ν == nothing; ν = 0.1/(dt*(0.65*nx/Lx)^nν);          end
  if tf != nothing; maxsteps = ceil(Int, tf/dt); qf=0.0   end

  # Number of substeps between vorticity-checking
  substeps = ceil(Int, loginterval/(maxq*dt))

  pr = TwoDTurb.Params(ν, nν)
  eq = TwoDTurb.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LC)

  if plots
    fig, axs = subplots()
    imshow(vs.q)
    pause(0.01)
  end

  @printf("\nMaking a mature turbulence field...\n")
  starttime = time()
  while maxq > qf && ts.step < maxsteps

    stepforward!(vs, ts, eq, pr, g; nsteps=substeps)
    TwoDTurb.updatevars!(vs, g)
    maxq = maximum(abs.(vs.q))

    if plots
      imshow(vs.q)
      pause(0.01)
    end

    log1 = @sprintf("τ: %.3f s, step: %d, t*q0: %.1e, max q: %.3e, ",
      time()-starttime, ts.step, vs.t*q0, maxq)

    log2 = @sprintf("ΔE: %.3f, CFL: %.3f",
      energy(vs, g)/E0, maximum([vs.U; vs.V])*ts.dt/g.dx)

    println(log1*log2)

  end

  @printf("... done.")

  return vs.q
end







end
# E N D   T W O D T U R B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
