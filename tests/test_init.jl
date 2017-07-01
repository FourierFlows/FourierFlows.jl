include("../src/solver.jl")

using PyPlot
using Solver

f0   = 1e-4     # Inertial frequency
sig  = 2e-4     # Wave frequency
kap  = 1e-5     # Modewise wavenumber
nuq  = 1e-2     # Vorticity hyperviscosity
nuqn = 4        # Vorticity hyperviscosity order
nua  = 1e-2     # Wave hyperviscosity
nuan = 4        # Wave hyperviscosity order

Lx = 2.0*pi
nx = 256
dt = 0.01 * 2.0*pi/f0

FFTW.set_num_threads(Sys.CPU_CORES)

g = Grid(nx, Lx)
p = Params(f0, sig, kap, nuq, nua, nuqn, nuan, g)
v = Vars(g, p)
qts = ETDRK4TimeStepper(dt, p.LCq)
ats = ETDRK4TimeStepper(dt, p.LCa)

println("Initialize grid, vars, params, time stepper:")
@time g = Grid(nx, Lx)
@time p = Params(f0, sig, kap, nuq, nua, nuqn, nuan, g)
@time v = Vars(g, p)
@time qts = ETDRK4TimeStepper(dt, p.LCq)
@time ats = ETDRK4TimeStepper(dt, p.LCa)
