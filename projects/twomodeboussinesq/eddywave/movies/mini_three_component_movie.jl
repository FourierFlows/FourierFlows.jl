include("../../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1w, mode1u 

@pyimport mpl_toolkits.axes_grid1 as pltgrid

plotname = "./plots/miniAPV"
filename = "../data/example_nk64_nu2e+32_122f_n1024_ep20_Ro05_nkw64.jld2"

# Plot parameters
eddylim = 10
Ro0 = 0.04
u0  = 4.00

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

  prob = TwoModeBoussinesq.InitialValueProblem(
    nx=n, Lx=L, ν0=ν0, nν0=nν0, ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt)

  # Non-retrievable parameters
  nkw = 64
   kw = 2π*nkw/L
    α = (N*kw)^2/(f*m)^2
    σ = f*sqrt(1+α)
   tσ = 2π/σ
   Re = L/80

  x, y = prob.grid.X/Re, prob.grid.Y/Re

  for (istep, step) in enumerate(steps)
    t = file["timeseries/t/$step"]
    Zh = file["timeseries/solr/$step"]
    solc = file["timeseries/solc/$step"]

    Z = irfft(Zh, prob.grid.nx) 
    u = ifft(solc[:, :, 1])
    v = ifft(solc[:, :, 2])
    p = ifft(solc[:, :, 3])

    TwoModeBoussinesq.set_Z!(prob, Z)
    TwoModeBoussinesq.set_uvp!(prob, u, v, p)

    Umax = maximum(sqrt.(prob.vars.U.^2.0+prob.vars.V.^2.0))

    w = mode1w(prob)
    u = mode1u(prob)/Umax
    Q = mode0apv(prob)

    println(step)

    close("all")
    fig, axs = subplots(ncols=1, nrows=3, sharex=true, sharey=true, 
      figsize=(5, 12))


    axes(axs[1]); axis("equal")
    Qplot = pcolormesh(x, y, Q/f, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

    axes(axs[2]); axis("equal")
    uplot = pcolormesh(x, y, u, cmap="RdBu_r", vmin=-u0, vmax=u0)

    axes(axs[3]); axis("equal")
    Zplot = pcolormesh(x, y, Z/f, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)


    plots = [Zplot, Qplot, uplot]
    for (i, ax) in enumerate(axs)
      ax[:axis]("off")
      ax[:set_adjustable]("box-forced")
      ax[:set_xlim](-eddylim, eddylim)
      ax[:set_ylim](-eddylim, eddylim)
      ax[:tick_params](axis="both", which="both", length=0)
      ax[:xaxis][:set_ticks]([])
      ax[:yaxis][:set_ticks]([])
   end

    tight_layout(rect=(0.00, 0.00, 1.0, 1.00))

    savename = @sprintf("%s_%06d.png", plotname, istep)
    savefig(savename, dpi=240)
  end

end
