include("../src/physics/twodturb.jl")

# twodturb_compileandspeedtests.jl
#
# This script constructs a sample two-dimensional turbulent problem, and then
# tests the integrity and speed of functions used to compile and solve the
# problem.

using TwoDTurb

# Problem parameters
nx  = 256     # Numerical problem size
dt  = 1e-4    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order

# Compile functions related to problem construction
g = Grid(nx, Lx)
p = Params(nu, nun, g)
v = Vars(p, g)
tsFE = ForwardEulerTimeStepper(dt, p.LC)
tsAB3 = AB3TimeStepper(dt, p.LC)
tsRK4 = RK4TimeStepper(dt, p.LC)
tsETDRK4 = ETDRK4TimeStepper(dt, p.LC)

# Evaluate grid, param, var generation
println()
@printf "%34s" "Init grid:"
@time g = Grid(nx, Lx)

@printf "%34s" "Init params:"
@time p = Params(nu, nun, g)

@printf "%34s" "Init vars:"
@time v = Vars(p, g)

@printf "%34s" "Init ForwEuler time stepper:"
@time ts = ForwardEulerTimeStepper(dt, p.LC)

@printf "%34s" "Init AB3 time stepper:"
@time ts = AB3TimeStepper(dt, p.LC)

@printf "%34s" "Init RK4 time stepper:"
@time ts = RK4TimeStepper(dt, p.LC)

@printf "%34s" "Init ETDRK4 time stepper:"
@time ts = ETDRK4TimeStepper(dt, p.LC)





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

stepforward!(3, tsFE, v, p, g)
stepforward!(3, tsAB3, v, p, g)
stepforward!(3, tsRK4, v, p, g)
stepforward!(3, tsETDRK4, v, p, g)

nsteps = 1000
println()
@printf "%22s %4d %s" "ForwEuler stepforward" nsteps "times:"
@time stepforward!(nsteps, tsFE, v, p, g)

@printf "%22s %4d %s" "AB3 stepforward" nsteps "times:"
@time stepforward!(nsteps, tsAB3, v, p, g)

@printf "%22s %4d %s" "RK4 stepforward" nsteps "times:"
@time stepforward!(nsteps, tsRK4, v, p, g)

@printf "%22s %4d %s" "ETDRK4 stepforward" nsteps "times:"
@time stepforward!(nsteps, tsETDRK4, v, p, g)
