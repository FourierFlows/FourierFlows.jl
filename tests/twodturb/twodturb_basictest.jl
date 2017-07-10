include("../../src/physics/twodturb.jl")

using TwoDTurb

# Problem parameters
nx  = 256     # Numerical problem size
dt  = 1e-3    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order


@printf "Generating grid, params, and vars..."
g  = Grid(nx, Lx)
p  = Params(nu, nun)
v  = Vars(g)
eq = Equation(p, g)
@printf " well, that seemed to go well.\n"


@printf "Initializing four time-steppers..."
tsFE     = ForwardEulerTimeStepper(dt, eq.LC)
tsAB3    = AB3TimeStepper(dt, eq.LC)
tsRK4    = RK4TimeStepper(dt, eq.LC)
tsETDRK4 = ETDRK4TimeStepper(dt, eq.LC)
@printf " well, that seemed to go well.\n"


@printf "Testing step-forward methods for each time-stepper..."
stepforward!(v, 3, tsFE,     eq, p, g)
stepforward!(v, 3, tsAB3,    eq, p, g)
stepforward!(v, 3, tsRK4,    eq, p, g)
stepforward!(v, 3, tsETDRK4, eq, p, g)
@printf " well, that seemed to go well.\n"


@printf "Testing var updating..."
updatevars!(v, p, g)
@printf " well, that seemed to go well.\n"
