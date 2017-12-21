include("../../src/FourierFlows.jl")

# twodturb_compileandspeedtests.jl
#
# This script constructs a sample two-dimensional turbulent problem, and then
# tests the integrity and speed of functions used to compile and solve the
# problem.

using FourierFlows,
      PyPlot

import TwoDTurb

# Problem parameters
nx  = 128     # Numerical problem size
dt  = 5e-3    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order

nsteps = 20000

# Initial condition
q0 = rand(nx, nx)


# Prepare problems for four time-steppers
g  = TwoDTurb.Grid(nx, Lx)
p  = TwoDTurb.Params(nu, nun)
eq = TwoDTurb.Equation(p, g)

vFE     = TwoDTurb.Vars(g)
vAB3    = TwoDTurb.Vars(g)
vRK4    = TwoDTurb.Vars(g)
vETDRK4 = TwoDTurb.Vars(g)

tsFE     = ForwardEulerTimeStepper(dt, eq.LC)
tsAB3    = AB3TimeStepper(dt, eq.LC)
tsRK4    = RK4TimeStepper(dt, eq.LC)
tsETDRK4 = ETDRK4TimeStepper(dt, eq.LC)

TwoDTurb.set_q!(vFE,     g, q0)
TwoDTurb.set_q!(vAB3,    g, q0)
TwoDTurb.set_q!(vRK4,    g, q0)
TwoDTurb.set_q!(vETDRK4, g, q0)


# Run four problems
@printf("Stepping forward with forward Euler... ")
@time stepforward!(vFE,     nsteps, tsFE,     eq, p, g)

@printf("Stepping forward with 3rd-order Adams-Bashforth... ")
@time stepforward!(vAB3,    nsteps, tsAB3,    eq, p, g)

@printf("Stepping forward with 4th-order Runge-Kutta... ")
@time stepforward!(vRK4,    nsteps, tsRK4,    eq, p, g)

@printf("Stepping forward with ETDRK4... ")
@time stepforward!(vETDRK4, nsteps, tsETDRK4, eq, p, g)


TwoDTurb.updatevars!(vFE,     g)
TwoDTurb.updatevars!(vAB3,    g)
TwoDTurb.updatevars!(vRK4,    g)
TwoDTurb.updatevars!(vETDRK4, g)




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
