include("../../src/fourierflows.jl")

# twodturb_compileandspeedtests.jl
#
# This script constructs a sample two-dimensional turbulent problem, and then
# tests the integrity and speed of functions used to compile and solve the
# problem.

using FourierFlows
using FourierFlows.TwoDTurb


# Problem parameters
nx  = 256     # Numerical problem size
dt  = 1e-4    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order

# Compile functions related to problem construction
g  = Grid(nx, Lx)
p  = Params(nu, nun)
v  = Vars(g)
eq = Equation(p, g)
tsFE     = ForwardEulerTimeStepper(dt, eq.LC)
tsAB3    = AB3TimeStepper(dt, eq.LC)
tsRK4    = RK4TimeStepper(dt, eq.LC)
tsETDRK4 = ETDRK4TimeStepper(dt, eq.LC)

# Evaluate grid, param, var generation
println()
@printf "%34s" "Init grid:"
@time g = Grid(nx, Lx)

@printf "%34s" "Init params:"
@time p = Params(nu, nun)

@printf "%34s" "Init vars:"
@time v = Vars(g)

@printf "%34s" "Init equation:"
@time eq = Equation(p, g)

@printf "%34s" "Init ForwEuler time stepper:"
@time ts = ForwardEulerTimeStepper(dt, eq.LC)

@printf "%34s" "Init AB3 time stepper:"
@time ts = AB3TimeStepper(dt, eq.LC)

@printf "%34s" "Init RK4 time stepper:"
@time ts = RK4TimeStepper(dt, eq.LC)

@printf "%34s" "Init ETDRK4 time stepper:"
@time ts = ETDRK4TimeStepper(dt, eq.LC)





# Evaluate var updating -------------------------------------------------------

function manyupdates!(nloops::Int)
  for i = 1:nloops; updatevars!(v, p, g); end
end

function manyupdates!(nloops::Int, withloops::Bool)
  for i = 1:nloops; updatevars!(v, p, g, withloops); end
end

# Compile step
updatevars!(v, p, g)
updatevars!(v, p, g, true)
manyupdates!(1)
manyupdates!(1, true)

nloops = 100
println()
@printf "%22s %4d %s" "Loop updatevars!" nloops "times:"
@time manyupdates!(nloops, true)

@printf "%22s %4d %s" "Fused updatevars!" nloops "times:"
@time manyupdates!(nloops)







# Evaluate stepping forward ---------------------------------------------------

stepforward!(v, 3, tsFE, eq, p, g)
stepforward!(v, 3, tsAB3, eq, p, g)
stepforward!(v, 3, tsRK4, eq, p, g)
stepforward!(v, 3, tsETDRK4, eq, p, g)

nsteps = 1000
println()
@printf "%22s %4d %s" "ForwEuler stepforward" nsteps "times:"
@time stepforward!(v, nsteps, tsFE, eq, p, g)

@printf "%22s %4d %s" "AB3 stepforward" nsteps "times:"
@time stepforward!(v, nsteps, tsAB3, eq, p, g)

@printf "%22s %4d %s" "RK4 stepforward" nsteps "times:"
@time stepforward!(v, nsteps, tsRK4, eq, p, g)

@printf "%22s %4d %s" "ETDRK4 stepforward" nsteps "times:"
@time stepforward!(v, nsteps, tsETDRK4, eq, p, g)
