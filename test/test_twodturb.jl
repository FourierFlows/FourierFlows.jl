include("../src/fourierflows.jl")

using FourierFlows
import FourierFlows.TwoDTurb

# Problem parameters
nx  = 256     # Numerical problem size
dt  = 1e-3    # Time step
Lx  = 2.0*pi  # Physical domain size
nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order


@printf "Generating grid, params, and vars..."
g  = TwoDTurb.Grid(nx, Lx)
p  = TwoDTurb.Params(nu, nun)
v  = TwoDTurb.Vars(g)
eq = TwoDTurb.Equation(p, g)
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
TwoDTurb.updatevars!(v, p, g)
@printf " well, that seemed to go well.\n"


# Initial condition with two ellipsoid vortices for comparison.
ampl = 1.131562576275490e-04
qh = ampl*rfft( 200.0*exp.(-((g.X-1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y-1).^2./.5^2)
  - 100.0* exp.(-((g.X+1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y+1).^2./.5^2) )

qh[1, 1] = 0

g = Grid(nx, Lx)
p = Params(f0, nuq, nuqn, g)
v = Vars(p, g)
qts = ForwardEulerTimeStepper(dt, p.LC)

Solver.updatevars!(v, p, g)

for n in 1:3
  println("step ", n)
  nsteps = 50

  Solver.stepforward!(nsteps, qts, v, p, g)
  Solver.updatevars!(v, p, g)
end
