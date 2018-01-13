module TwoDTurb

using FourierFlows
Grid = TwoDGrid

# Problem ---------------------------------------------------------------------
"""
Construct an initial-value 2D turbulence problem.
"""
function InitialValueProblem(;
     nx = 256,
     Lx = 2π,
     ny = nx,
     Ly = Lx,
      ν = 0.0,
     nν = 2,
     dt = 0.01,
stepper = "ETDRK4"
  )

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.Params(ν, nν)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  
  FourierFlows.Problem(g, vs, pr, eq, ts)
end

function InitialValueProblem(n, L, ν, nν, dt, withfilter)
  InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, withfilter=withfilter)
end


"""
Construct a forced 2D turbulence problem.
"""
function ForcedProblem(;
        nx = 256,
        Lx = 2π,
        ny = nx,
        Ly = Lx,
         ν = 0.0,
        nν = 2,
         μ = 0.0,
        dt = 0.01,
   stepper = "RK4",
     calcF = nothing
  )

  if calcF == nothing; _calcF(F, sol, t, s, v, p, g) = nothing
  else;                _calcF = calcF
  end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.ForcedParams(ν, nν, μ, _calcF)
  vs = TwoDTurb.ForcedVars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
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
  FourierFlows.Equation{2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g::TwoDGrid)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ
  FourierFlows.Equation{2}(LC, calcN_forced!)
end




# Vars
physvars = [:q, :U, :V, :Uq, :Vq, :psi]
transvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih]

expr = FourierFlows.getexpr_varstype(:Vars, physvars, transvars)
eval(expr)

function Vars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih
  Vars(q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih)
end


forcedtransvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih, :F]
expr = FourierFlows.getexpr_varstype(:ForcedVars, physvars, forcedtransvars)
eval(expr)

function ForcedVars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih F
  ForcedVars(q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih, F)
end




# Solvers
function calcN_advection!(
  N::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, s::State, v::AbstractVars, p::AbstractParams, g::TwoDGrid)

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

function calcN_forced!(N::Array{Complex{Float64}, 2}, 
                sol::Array{Complex{Float64}, 2}, t::Float64, 
                s::State, v::ForcedVars, p::ForcedParams, g::TwoDGrid)

  calcN_advection!(N, sol, t, s, v, p, g)
  p.calcF!(v.F, sol, t, s, v, p, g)

  @. N += v.F
  nothing
end


# Helper functions
"""
Update solution variables using qh in v.sol.
"""
function updatevars!(s, v, g)
  v.qh .= s.sol
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
  updatevars!(prob.state, prob.vars, prob.grid)
end


"""
Set the vorticity field.
"""
function set_q!(s, v, g, q)
  A_mul_B!(s.sol, g.rfftplan, q)
  updatevars!(s, v, g)
end

function set_q!(prob::AbstractProblem, q)
  set_q!(prob.state, prob.vars, prob.grid, q)
end


"""
Calculate the domain integrated kinetic energy.
"""
function energy(s, v, g)
  @. v.Uh =  im * g.l  * g.invKKrsq * s.sol
  @. v.Vh = -im * g.kr * g.invKKrsq * s.sol
  0.5*(FourierFlows.parsevalsum2(v.Uh, g)+FourierFlows.parsevalsum2(v.Vh, g))
end

function energy(prob)
  energy(prob.state, prob.vars, prob.grid)
end


"""
Returns the domain-integrated enstrophy.
"""
function enstrophy(s, g)
  0.5*FourierFlows.parsevalsum2(s.sol, g)
end

function enstrophy(prob)
  enstrophy(prob.state, prob.grid)
end


function injection(s, v, g)
  @. v.psih = -g.invKKrsq * s.sol
  @. v.Uq = -real(v.psih*conj(v.F))
  FourierFlows.parsevalsum(v.Uq)
end




#=
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
=#


end # module
