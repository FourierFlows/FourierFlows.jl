include("../../src/fourierflows.jl")

using  FourierFlows, PyPlot
import FourierFlows, FourierFlows.TracerAdvDiff
       

# Domain
L     = 2.0*pi
nx    = 128
u0    = 1.0
k = m = 4 * 2.0*pi/Lx
kap   = 1e-2
dt    = 0.5*(Lx/nx) / u0

nsteps = 10000
nsubs  = 100

h(x) = 1.0 + 0.4*sin(2*pi*x/L)


# psi = u0/m * sin(kx) * sin(mz). u = psi_z, w = -psi_x
u(x, z, t) = u0*sin(k*x)*cos(m*z)
w(x, z, t) = -u0*k/m*cos(k*x)*sin(m*z)

g  = FourierFlows.TwoDGrid(nx, Lx)
pr = TracerAdvDiff.ConstDiffParams(kap, u, w)
vs = TracerAdvDiff.Vars(g)
eq = TracerAdvDiff.Equation(pr, g)
ts = FourierFlows.ETDRK4TimeStepper(dt, eq.LC) 

# Initial condition
dc = 0.05*Lx
c0(x, z) = exp( -(x^2+z^2)/(2*dc^2) ) / (2*pi*dc^2)
TracerAdvDiff.set_c!(vs, pr, g, c0)


fig, axs = subplots()
fig[:show]()

while ts.step < nsteps

  stepforward!(vs, ts, eq, pr, g; nsteps=nsubs)
  TracerAdvDiff.updatevars!(vs, pr, g)

  println(ts.step)

  imshow(vs.c)
  pause(0.01)
  #fig[:canvas][:draw]()

end
