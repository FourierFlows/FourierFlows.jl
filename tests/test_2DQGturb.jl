include("../src/physics/2DQGturb.jl")

# using PyPlot
using Solver

f0   = 1e-4     # Inertial frequency
nuq  = 1e-6    # Vorticity hyperviscosity
nuqn = 4        # Vorticity hyperviscosity order

Lx = 2.0*pi
nx = 128
dt = 1e-4;

FFTW.set_num_threads(Sys.CPU_CORES)


println("Initialize grid, vars, params, time stepper:")
@time g = Grid(nx, Lx);
@time p = Params(f0, nuq, nuqn, g);
@time v = Vars(p, g);
@time qts = ForwardEulerTimeStepper(dt, p.LC);
# @time qts = ETDRK4TimeStepper(dt, p.LCq);

Solver.updatevars!(v,p,g);

# println(v.q[10,14])

# Lz = qts.LC.*v.qh;
# println(Lz[4,5])


# figure(1);
# pcolor(g.X,g.Y,v.q)
# colorbar();
# axis(:equal);xlim(-pi,pi);ylim(-pi,pi);
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
