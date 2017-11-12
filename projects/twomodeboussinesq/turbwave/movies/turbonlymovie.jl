include("../../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: mode1u, mode1w 

@pyimport mpl_toolkits.axes_grid1 as pltgrid

filename = "twodturb_n512_Ro20_nnu4_nu1e+12.jld2"
filepath = "../data"
fullfilename = joinpath(filepath, filename)

plotname = "turbonly"
if !isdir(joinpath(".", plotname)); mkdir(plotname); end

# Parameters
f = 1e-4
Ro0 = 0.2

# Recreate problem
jldopen(filename, "r") do file

  ν, nν = file["params"]["ν"], file["params"]["nν"]

   n = file["grid"]["nx"]
   L = file["grid"]["Lx"]
  dt = file["timestepper"]["dt"]

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
    ax[:tick_params](axis="both", which="both", length=0)

    divider = pltgrid.make_axes_locatable(ax)
    cax = divider[:append_axes]("right", size="5%", pad="8%")
    cb = colorbar(qplot, cax=cax, orientation="vertical")

    ax[:set_xlabel](L"x/R")
    ax[:set_ylabel](L"y/R")

    cb[:ax][:xaxis][:set_ticks_position]("right")
    cb[:ax][:xaxis][:set_label_position]("right")

    cb[:ax][:tick_params](axis="both", which="both", length=0)
    cb[:set_ticks]([-Ro0, 0.0, Ro0])

    labelpad = 16.0
    cb[:set_label]("Barotropic PV \$= Q/f\$", labelpad=labelpad)

    tight_layout()

    savename = joinpath(".", plotname, 
      @sprintf("%s_%06d.png", plotname, istep))
    println(savename)
    savefig(savename, dpi=240)
  end

end
