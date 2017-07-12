include("../../src/physics/barotropicqg.jl")

using BarotropicQG, PyPlot

# Problem parameters
nx       = 128
betastar = 1.0 
fdt      = 0.1 

g, p, v, eq = southern_gaussian_mountain(nx, betastar)
ts = ETDRK4TimeStepper(fdt*2.0*pi/p.f0, eq.LC)

stepforward!(v, 10, ts, eq, p, g)

fig = figure()
imshow(v.zet)
#pcolormesh(g.X, g.Y, v.zet)



