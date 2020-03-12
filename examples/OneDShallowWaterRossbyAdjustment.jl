using FourierFlows, PyPlot
using FFTW: rfft, irfft, fft, ifft
using LinearAlgebra: mul!, ldiv!
using Printf
using Random

nothingfunction(args...) = nothing

function calcN!(N, sol, t, s, v, p, g)
  mul!(v.uh, g.rfftplan, deepcopy(v.u))
  mul!(v.vh, g.rfftplan, deepcopy(v.v))
  mul!(v.etah, g.rfftplan, deepcopy(v.eta))
  ldiv!(v.ux, g.rfftplan, -im*g.kr.*v.uh)
  ldiv!(v.vx, g.rfftplan, -im*g.kr.*v.vh)
  mul!(v.uuxh, g.rfftplan, v.u.*v.ux)
  mul!(v.uvxh, g.rfftplan, v.u.*v.vx)
  rhsu = p.f*v.vh + p.g*im*g.kr.*v.etah - v.uuxh
  rhsv = - p.f*v.uh - v.uvxh
  rhseta = p.h*im*g.kr.*v.uh
  N[:,1] .= rhsu
  N[:,2] .= rhsv
  N[:,3] .= rhseta
  dealias!(N, g, g.kralias)
  nothing
end

struct Vars{Aphys, Atrans} <: AbstractVars
        u :: Aphys
        v :: Aphys
      eta :: Aphys
       ux :: Aphys
       vx :: Aphys
       uh :: Atrans
       vh :: Atrans
     etah :: Atrans
     uuxh :: Atrans
     uvxh :: Atrans
end

function Vars(::Dev, g::AbstractGrid{T}) where {Dev, T}
  @devzeros Dev T g.nx u v eta ux vx
  @devzeros Dev Complex{T} (g.nkr) uh vh etah uuxh uvxh
  mul!(uh, g.rfftplan, u)
  mul!(vh, g.rfftplan, v)
  mul!(etah, g.rfftplan, eta)
  Vars(u, v, eta, ux, vx, uh, vh, etah, uuxh, uvxh)
end

struct Params{T} <: AbstractParams
       ν :: T         # Hyperviscosity coefficient
      nν :: Int       # Order of the hyperviscous operator
      g  :: T         # Gravitational acceleration
      h  :: T         # Fluid depth
      f  :: T         # Coriolis parameter
end
Params(ν, nν, g, h, f) = Params(ν, nν, g, h, f)

function Equation(p, g)
  LC = zeros(Float64, g.nkr, 3)
  D = - p.ν*g.kr.^p.nν
  LC[:,1] .= D # u
  LC[:,2] .= D # v
  LC[:,3] .= D # eta
  FourierFlows.Equation(LC, calcN!, g)
end

"""
    updatevars!(prob)
Update the vars in v on the grid g with the solution in sol.
"""
function updatevars!(prob)
  v, g, sol, p = prob.vars, prob.grid, prob.sol, prob.params
  ldiv!(v.u, g.rfftplan, deepcopy(sol[:,1]))
  ldiv!(v.v, g.rfftplan, deepcopy(sol[:,2]))
  ldiv!(v.eta, g.rfftplan, deepcopy(sol[:,3]))
  nothing
end

function set_uveta!(u0, v0, eta0, prob)
    g, sol, p, v = prob.grid, prob.sol, prob.params, prob.vars
    mul!(v.uh, g.rfftplan, u0)
    mul!(v.vh, g.rfftplan, v0)
    mul!(v.etah, g.rfftplan, eta0)
    sol[:,1] = rfft(u0)
    sol[:,2] = rfft(v0)
    sol[:,3] = rfft(eta0)
    updatevars!(prob)
end

function updateplot!(axs, prob)
  sca(axs[1])
  cla()
  plot(prob.grid.x, prob.vars.eta)
  xlabel(L"$x$ [m]")
  ylabel(L"$\eta$ [m]")
  sca(axs[2])
  cla()
  plot(prob.grid.x, prob.vars.u)
  xlabel(L"$x$ [m]")
  ylabel(L"$u$ [m/s]")
  sca(axs[3])
  cla()
  plot(prob.grid.x, prob.vars.v)
  xlabel(L"$x$ [m]")
  ylabel(L"$v$ [m/s]")
  xlim(prob.grid.x[1], prob.grid.x[end])
  subplots_adjust(hspace=0)
  nothing
end

#---
nx = 256
Lx = 20000π
ν = 1e-32
nν = 1
dt = 0.01
g = 10.0
h = 200.0
f = 1e-2
dev = CPU()

nt = 500000

grd = OneDGrid(nx, Lx)
vars = Vars(dev, grd)
params = Params(ν, nν, g, h, f)
eqn = Equation(params, grd)
stepper = "FilteredRK4"
iplt = 100

# Set up problem.
prob = FourierFlows.Problem(eqn, stepper, dt, grd, vars, params, dev)

u0 = zeros(grd.nx)
v0 = zeros(grd.nx)
noise = deepcopy(u0)
noise = Random.randn!(noise)/10
s = 1e4
eta0 = @. 5*exp(-grd.x.^2/s^2) + noise
set_uveta!(u0, v0, eta0, prob)

# Set up figure.
fig, axs = subplots(nrows=3, figsize=(8, 4), sharex=true)
Ld = @sprintf "%d" sqrt(g*h)/f
suptitle(L"$L_d = \sqrt{gh}/f = $"*string(Ld)*" m")

sca(axs[1])
plot(prob.grid.x, prob.vars.eta)
xlabel(L"$x$ [m]")
ylabel(L"$\eta$ [m]")

sca(axs[2])
plot(prob.grid.x, prob.vars.u)
xlabel(L"$x$ [m]")
ylabel(L"$u$ [m/s]")

sca(axs[3])
plot(prob.grid.x, prob.vars.v)
xlabel(L"$x$ [m]")
ylabel(L"$v$ [m/s]")
xlim(prob.grid.x[1], prob.grid.x[end])

subplots_adjust(hspace=0)
show()

# Time-step.
for i=1:nt
    println(i)
    if rem(i, iplt)==0
        updateplot!(axs, prob)
    end
    stepforward!(prob)
    updatevars!(prob)
end
