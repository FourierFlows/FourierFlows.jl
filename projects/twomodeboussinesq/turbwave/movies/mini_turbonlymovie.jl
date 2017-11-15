include("../../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoDTurb

@pyimport mpl_toolkits.axes_grid1 as pltgrid

plotpath = "./plots"
plotname = "miniturbonly"
filename = "twodturb_n512_Ro20_nnu4_nu1e+12_fwd.jld2"
filepath = "../data"
fullfilename = joinpath(filepath, filename)


# Parameters
f = 1e-4
Ro0 = 0.10

# Recreate problem
jldopen(fullfilename, "r") do file

  #ν, nν = file["params"]["ν"], file["params"]["nν"]
  ν, nν = 1e12, 4
  n, L = 512, 2π*1600e3
  dt = 1.0

   #n = file["grid"]["nx"]
   #L = file["grid"]["Lx"]
  #dt = file["timestepper"]["dt"]

  steps = parse.(Int, keys(file["timeseries/t"]))

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt)

  x, y = prob.grid.X, prob.grid.Y

  for (istep, step) in enumerate(steps)
    t = file["timeseries/t/$step"]
    Qh = file["timeseries/sol/$step"]

    Q = irfft(Qh, prob.grid.nx) 
    Ro = Q/f

    #TwoDTurb.set_q!(prob, Q) 
    #TwoDTurb.updatevars!(prob)
    
    close("all")
    fig, ax = subplots()
    qplot = pcolormesh(x, y, Ro, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

    axis("equal")
    ax[:set_adjustable]("box-forced")
    ax[:axis]("off")

    tight_layout(rect=(0.0, 0.0, 1.0, 0.95))

    savename = @sprintf("%s_%06d.png", joinpath(plotpath, plotname), istep)
    println(savename)
    savefig(savename, dpi=240)
  end

end
