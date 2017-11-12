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

# Plot parameters
Ro0 = 0.2
ε0  = 0.1


# Recreate problem
jldopen(filename, "r") do file

  ν0, nν0 = file["params"]["ν0"], file["params"]["nν0"]
  ν1, nν1 = file["params"]["ν1"], file["params"]["nν1"]

   n = file["grid"]["nx"]
   L = file["grid"]["Lx"]
   f = file["params"]["f"]
   N = file["params"]["N"]
   m = file["params"]["m"]
  dt = file["timestepper"]["dt"]

  steps = parse.(Int, keys(file["timeseries/t"]))

  prob = TwoModeBoussinesq.PrognosticAPVInitialValueProblem(
    nx=n, Lx=L, ν0=ν0, nν0=nν0, ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt)

  # Non-retrievable parameters
  nkw = 16
   kw = 2π*nkw/L
    α = (N*kw)^2/(f*m)^2
    σ = f*sqrt(1+α)
   tσ = 2π/σ
   Lk = 2π/kw

  x, y = prob.grid.X/Lk, prob.grid.Y/Lk

  iscolorbar = false
  for (istep, step) in enumerate(steps)
    t = file["timeseries/t/$step"]
    Qh = file["timeseries/solr/$step"]
    solc = file["timeseries/solc/$step"]

    Q = irfft(Qh, prob.grid.nx) 
    u = ifft(solc[:, :, 1])
    v = ifft(solc[:, :, 2])
    p = ifft(solc[:, :, 3])

    TwoModeBoussinesq.set_Q!(prob, Q)
    TwoModeBoussinesq.set_uvp!(prob, u, v, p)
    
    Umax = maximum(sqrt.(prob.vars.U.^2.0+prob.vars.V.^2.0))

    #u = mode1u(prob)/Umax
    Ro = Q/f
     ε = mode1u(prob)*kw/σ


    close("all")
    fig, axs = subplots(ncols=2, nrows=1, sharex=true, sharey=true, 
      figsize=(8, 5))

    axes(axs[1]); axis("equal")
    plot_Ro = pcolormesh(x, y, Ro, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

    axes(axs[2]); axis("equal")
    plot_ε = pcolormesh(x, y, ε, cmap="RdBu_r", vmin=-ε0, vmax=ε0)


    plots = [plot_Ro, plot_ε]
    cbs = []
    for (i, ax) in enumerate(axs)
      ax[:set_adjustable]("box-forced")
      ax[:tick_params](axis="both", which="both", length=0)

      divider = pltgrid.make_axes_locatable(ax)
      cax = divider[:append_axes]("top", size="5%", pad="8%")
      cb = colorbar(plots[i], cax=cax, orientation="horizontal")

      cb[:ax][:xaxis][:set_ticks_position]("top")
      cb[:ax][:xaxis][:set_label_position]("top")
      cb[:ax][:tick_params](axis="x", which="both", length=0)

      push!(cbs, cb)
      iscolorbar = true
    end

    axs[1][:set_xlabel](L"x/R")
    axs[2][:set_xlabel](L"x/R")
    axs[1][:set_ylabel](L"y/R")

    ticks = [-12, -6, 0, 6, 12]
    axs[1][:xaxis][:set_ticks](ticks)
    axs[1][:yaxis][:set_ticks](ticks)
    axs[2][:xaxis][:set_ticks](ticks)

    cbs[1][:ax][:xaxis][:set_ticks_position]("top")
    cbs[1][:ax][:xaxis][:set_label_position]("top")

    cbs[1][:set_ticks]([-Ro0, 0.0, Ro0])
    cbs[2][:set_ticks]([-ε0, 0.0, ε0])

    labelpad = 16.0
    cbs[1][:set_label]("Barotropic PV \$= Q/f\$", labelpad=labelpad)
    cbs[2][:set_label]("Normalized \$\\hat u(z=0) = (u + u^*)k/\\sigma\$",
      labelpad=labelpad)

    msg = @sprintf("\$t = %02d\$ wave periods", t/tσ)
    figtext(0.51, 0.95, msg, horizontalalignment="center", fontsize=14)

    tight_layout(rect=(0.00, 0.00, 0.95, 0.90), h_pad=0.05)

    savename = @sprintf("%s_%06d.png", plotname, istep)
    println(savename)
    savefig(savename, dpi=240)
  end

end
