# This example solves the 1D rotating shallow water equations
# for the depth-integrated velocities u and v and for the
# surface elevation η in a fluid with constant depth h,
# gravitational acceleration g and Coriolis parameter f:
#
# ∂u/∂t + u∂u/∂x - fv = -gh∂η/∂x - D,   (1)
# ∂v/∂t + u∂v/∂x + fu = - D,            (2)
# ∂η/∂t + ∂u/∂x = - D.                  (3)
#
# where D indicates a hyperviscous operator of
# the form ν*Laplacian^nν, where ν is the viscosity
# coefficient and nν is the order of the operator.
# The code time-steps Equations (1-3) simultaneously.
#
# The initial conditions prescribe a Gaussian bump larger
# than the deformation radius in a fluid at rest, with
# some random noise smaller than the deformation radius
# added on top. The system develops geostrophically-balanced
# jets around the Gaussian bump, while the smaller-scale noise
# propagates around as inertia-gravity waves.

using FourierFlows, PyPlot
using FFTW: rfft, irfft, fft, ifft
using LinearAlgebra: mul!, ldiv!
using Printf
using Random


function calcN!(N, sol, t, s, v, p, g)
  mul!(v.uh, g.rfftplan, deepcopy(v.u))
  mul!(v.vh, g.rfftplan, deepcopy(v.v))
  mul!(v.ηh, g.rfftplan, deepcopy(v.η))
  ldiv!(v.ux, g.rfftplan, -im*g.kr.*v.uh)
  ldiv!(v.vx, g.rfftplan, -im*g.kr.*v.vh)
  mul!(v.uuxh, g.rfftplan, v.u.*v.ux)
  mul!(v.uvxh, g.rfftplan, v.u.*v.vx)
  rhsu = p.f*v.vh + p.g*p.h*im*g.kr.*v.ηh - v.uuxh
  rhsv = - p.f*v.uh - v.uvxh
  rhsη = im*g.kr.*v.uh
  N[:, 1] .= rhsu
  N[:, 2] .= rhsv
  N[:, 3] .= rhsη
  dealias!(N, g, g.kralias)
  nothing
end

struct Vars{Aphys, Atrans} <: AbstractVars
        u :: Aphys
        v :: Aphys
        η :: Aphys
       ux :: Aphys
       vx :: Aphys
       uh :: Atrans
       vh :: Atrans
       ηh :: Atrans
     uuxh :: Atrans
     uvxh :: Atrans
end

function Vars(::Dev, g::AbstractGrid{T}) where {Dev, T}
  @devzeros Dev T g.nx u v η ux vx
  @devzeros Dev Complex{T} (g.nkr) uh vh ηh uuxh uvxh
  mul!(uh, g.rfftplan, u)
  mul!(vh, g.rfftplan, v)
  mul!(ηh, g.rfftplan, η)
  Vars(u, v, η, ux, vx, uh, vh, ηh, uuxh, uvxh)
end

struct Params{T} <: AbstractParams
       ν :: T         # Hyperviscosity coefficient
      nν :: Int       # Order of the hyperviscous operator
      g  :: T         # Gravitational acceleration
      h  :: T         # Fluid depth
      f  :: T         # Coriolis parameter
end
Params(ν, nν, g, h, f) = Params(ν, nν, g, h, f)

function Equation(p, g::AbstractGrid{T}) where T
  LC = zeros(Float64, g.nkr, 3)
  D = - p.ν*g.kr.^(2*p.nν)
  LC[:, 1] .= D # u
  LC[:, 2] .= D # v
  LC[:, 3] .= D # η
  FourierFlows.Equation(LC, calcN!, g)
end

"""
    updatevars!(prob)
Update the vars in v on the grid g with the solution in sol.
"""
function updatevars!(prob)
  v, g, sol, p = prob.vars, prob.grid, prob.sol, prob.params
  ldiv!(v.u, g.rfftplan, deepcopy(sol[:, 1]))
  ldiv!(v.v, g.rfftplan, deepcopy(sol[:, 2]))
  ldiv!(v.η, g.rfftplan, deepcopy(sol[:, 3]))
  nothing
end

function set_uvη!(u0, v0, η0, prob)
    g, sol, p, v = prob.grid, prob.sol, prob.params, prob.vars
    mul!(v.uh, g.rfftplan, u0)
    mul!(v.vh, g.rfftplan, v0)
    mul!(v.ηh, g.rfftplan, η0)
    sol[:, 1] = v.uh
    sol[:, 2] = v.vh
    sol[:, 3] = v.ηh
    updatevars!(prob)
end

function updateplot!(axs, prob)
  sca(axs[1])
  cla()
  plot(prob.grid.x, prob.vars.η)
  xlabel(L"$x$ [m]")
  ylabel(L"$\eta$ [m]")
  sca(axs[2])
  cla()
  plot(prob.grid.x, prob.vars.u/prob.params.h) # Plot depth-averaged u.
  xlabel(L"$x$ [m]")
  ylabel(L"$\bar{u}$ [m/s]")
  sca(axs[3])
  cla()
  plot(prob.grid.x, prob.vars.v/prob.params.h) # Plot depth-averaged v.
  xlabel(L"$x$ [m]")
  ylabel(L"$\bar{v}$ [m/s]")
  xlim(prob.grid.x[1], prob.grid.x[end])
  subplots_adjust(hspace=0)
  nothing
end

#---
nx = 256     # Resolution
Lx = 20000π  # Domain length
ν = 1e-32    # Viscosity
nν = 1       # Viscosity order (nν = 1 means Laplacian)
dt = 0.01
g = 10.0     # Gravitational acceleration
h = 200.0    # Fluid depth
f = 1e-2     # Coriolis parameter
dev = CPU()

nt = 100000

grd = OneDGrid(nx, Lx)
vars = Vars(dev, grd)
params = Params(ν, nν, g, h, f)
eqn = Equation(params, grd)
stepper = "FilteredRK4"
iplt = 100

# Set up problem.
prob = FourierFlows.Problem(eqn, stepper, dt, grd, vars, params, dev)

# Set up the initial conditions, in this case a large Gaussian bump
# with small-scale perturbations superimposed.
u0 = zeros(grd.nx)
v0 = zeros(grd.nx)
noise = deepcopy(u0)
noise = Random.randn!(noise)/10
s = 1e4
η0 = @. 3*exp(-grd.x.^2/s^2) + noise
set_uvη!(u0, v0, η0, prob)

# Set up figure.
fig, axs = subplots(nrows=3, figsize=(8, 4), sharex=true)
Ld = @sprintf "%d" sqrt(g*h)/f
suptitle(L"Deformation radius $L_d = \sqrt{gh}/f = $"*string(Ld)*" m")

sca(axs[1])
plot(prob.grid.x, prob.vars.η)
xlabel(L"$x$ [m]")
ylabel(L"$\eta$ [m]")

sca(axs[2])
plot(prob.grid.x, prob.vars.u/prob.params.h) # Plot depth-averaged u.
xlabel(L"$x$ [m]")
ylabel(L"$\bar{u}$ [m/s]")

sca(axs[3])
plot(prob.grid.x, prob.vars.v/prob.params.h) # Plot depth-averaged v.
xlabel(L"$x$ [m]")
ylabel(L"$\bar{v}$ [m/s]")
xlim(prob.grid.x[1], prob.grid.x[end])

subplots_adjust(hspace=0)
show()

# Time-step the problem.
for i=1:nt
    println(i)
    if rem(i, iplt)==0
        updateplot!(axs, prob)
    end
    stepforward!(prob)
    updatevars!(prob)
end
