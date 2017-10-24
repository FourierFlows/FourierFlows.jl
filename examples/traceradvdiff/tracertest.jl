include("../../src/fourierflows.jl")

using FourierFlows
import FourierFlows, FourierFlows.TracerAdvDiff

include("./tracerutils.jl")
       

nx, kap, Lx, Ly = 128, 1e-3, 2*pi, 1.0    # Domain
nsteps, nsubs = 1000, 10                  # Time integration


# Barotropic flow
U, ep, k1 = 1.0, 0.6, 2*pi/Lx
h(x)  = Ly*(1.0 - ep*sin(k1*x))
hx(x) = -Ly*k1*ep*cos(k1*x) 

ub(x, y, t) = U/h(x)
wb(x, y, t) = y*ub(x, y, t)*hx(x)/h(x)

# Time stepping
CFL = 0.1
del, umax = Lx/nx, U/(Ly*(1-ep))

# Initial condition
x0, y0, dc = 0.0, -0.5*Ly, 0.001*Lx
c0(x, y) = exp( -((x-x0)^2+(y-y0)^2)/(2*dc^2) ) / (2*pi*dc^2)


# Initialize problem
g  = FourierFlows.TwoDGrid(nx, Lx, nx, Ly+ep; y0=-Ly-ep)
pr = TracerAdvDiff.ConstDiffParams(kap, ub, wb)
vs = TracerAdvDiff.Vars(g)
eq = TracerAdvDiff.Equation(pr, g)
ts = FourierFlows.ETDRK4TimeStepper(CFL*del/umax, eq.LC) 
TracerAdvDiff.set_c!(vs, pr, g, c0)


# Diagnostics
ctot(vs, pr, g) = Myn(vs.c, g, 0)
ycen(vs, pr, g) = Myn(vs.c, g, 1)/Myn(vs.c, g, 0)
xcen(vs, pr, g) = Mxn(vs.c, g, 1)/Mxn(vs.c, g, 0)
yvar(vs, pr, g) = Cy2(vs.c, g)

ntdiag = Int(nsteps/nsubs)+1
totalc   = Diagnostic(ctot, nsubs, vs, pr, ts, g)
ycenter  = Diagnostic(ycen, nsubs, vs, pr, ts, g)
xcenter  = Diagnostic(xcen, nsubs, vs, pr, ts, g)
variance = Diagnostic(yvar, nsubs, vs, pr, ts, g)
diags = [totalc, xcenter, ycenter, variance]


mask = Array{Bool}(g.nx, g.ny)
for j=1:g.ny, i=1:g.nx
  mask[i, j] = g.y[j] < -h(g.x[i]) ? true : false
end


function makeplot(vs, g, fig, axs)
  cplot = NullableArray(vs.c, mask)

  axes(axs[1])
  pcolormesh(g.X, g.Y, PyObject(cplot))

  axes(axs[2])
  plot(vs.t, variance.value, "c0.")

  pause(0.01)
  nothing
end

fig, axs = subplots(nrows=2)

# Integrate
while ts.step < nsteps

  stepforward!(vs, ts, eq, pr, g; nsteps=nsubs)
  TracerAdvDiff.updatevars!(vs, pr, g)
  increment!(diags)

  @printf("step: %04d, xbar: %.2e, ybar: %.2e, ysig: %.2e\n", 
    ts.step, xcenter.value, ycenter.value, variance.value)

  makeplot(vs, g, fig, axs)

end
