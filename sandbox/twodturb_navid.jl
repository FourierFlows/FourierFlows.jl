include("../src/fourierflows.jl")

using FourierFlows,
      PyPlot
import FourierFlows.TwoDTurb

TwoDTurb.InitialValueProblem()

# g  = TwoDGrid(nx, Lx, ny, Ly)
# pr = TwoDTurb.Params(ν, nν)
# vs = TwoDTurb.Vars(g)
# eq = TwoDTurb.Equation(pr, g)
# ts = ETDRK4TimeStepper(dt, eq.LC)

# while maxq > qf && ts.step < maxsteps

for n in 1:300
  stepforward!(vs, ts, eq, pr, g; nsteps=100)
  TwoDTurb.updatevars!(vs, g)
  maxq = maximum(abs.(vs.q))

  imshow(vs.q)
  pause(0.01)
end


asdfasdf


nu  = 1e-6    # Vorticity hyperviscosity
nun = 4       # Vorticity hyperviscosity order
Lx  = 2.0*pi
nx  = 128
dt  = 1e-4;


println("Initialize grid, vars, params, time stepper:")
@time g = TwoDTurb.TwoDGrid(nx, Lx)
@time p = TwoDTurb.Params(nu, nun)
@time v  = TwoDTurb.Vars(g)
@time eq = TwoDTurb.Equation(p, g)

@time ts = ForwardEulerTimeStepper(dt, eq.LC);
# @time qts = ETDRK4TimeStepper(dt, p.LCq);

TwoDTurb.updatevars!(v, p, g)

# println(v.q[10,14])

# Lz = qts.LC.*v.qh;
# println(Lz[4,5])


figure(1);
pcolor(g.X,g.Y,v.q)
colorbar()
axis(:equal);xlim(-pi,pi);ylim(-pi,pi)
pause(2)
asdfasdf

for n in 1:3
  println("step ",n)
  nsteps = 50;
  @time Solver.stepforward!(nsteps, qts, v, p, g)
  Solver.updatevars!(v,p,g);
  # figure(1);clf()
  # pcolor(g.X,g.Y,v.q)
  # axis(:equal);xlim(-pi,pi);ylim(-pi,pi);
  # colorbar()
end
