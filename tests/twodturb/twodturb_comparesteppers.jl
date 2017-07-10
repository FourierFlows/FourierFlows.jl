include("../../src/physics/twodturb.jl")

# twodturb_compileandspeedtests.jl
#
# This script constructs a sample two-dimensional turbulent problem, and then
# tests the integrity and speed of functions used to compile and solve the
# problem.

using TwoDTurb, PyPlot

# Problem parameters
nx  = 128     # Numerical problem size
dt  = 5e-3    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order

nsteps = 20000

# Initial condition
q0 = rand(nx, nx)



# Prepare problems for four time-steppers, and run.
g = Grid(nx, Lx)
p = Params(nu, nun, g)

vFE     = Vars(p, g)
vAB3    = Vars(p, g)
vRK4    = Vars(p, g)
vETDRK4 = Vars(p, g)

tsFE     = ForwardEulerTimeStepper(dt, p.LC)
tsAB3    = AB3TimeStepper(dt, p.LC)
tsRK4    = RK4TimeStepper(dt, p.LC)
tsETDRK4 = ETDRK4TimeStepper(dt, p.LC)

set_q!(vFE,     g, q0)
set_q!(vAB3,    g, q0)
set_q!(vRK4,    g, q0)
set_q!(vETDRK4, g, q0)

stepforward!(nsteps, tsFE, vFE, p, g)
stepforward!(nsteps, tsAB3, vAB3, p, g)
stepforward!(nsteps, tsRK4, vRK4, p, g)
stepforward!(nsteps, tsETDRK4, vETDRK4, p, g)

updatevars!(vFE,     g)
updatevars!(vAB3,    g)
updatevars!(vRK4,    g)
updatevars!(vETDRK4, g)




# Plot
fig, axs = subplots(nrows=2, ncols=2, sharex=true, sharey=true)

axes(axs[1, 1])
pcolormesh(g.x, g.y, vFE.q)
axis("tight")
title("Forward Euler")

axes(axs[1, 2])
pcolormesh(g.x, g.y, vAB3.q)
axis("tight")
title("3rd order Adams-Bashforth")

axes(axs[2, 1])
pcolormesh(g.x, g.y, vRK4.q)
axis("tight")
title("4th order Runge-Kutta")

axes(axs[2, 2])
pcolormesh(g.x, g.y, vETDRK4.q)
axis("tight")
title("ETDRK4")

savefig("fourtimesteppers.png"; dpi=240, facecolor="w")
