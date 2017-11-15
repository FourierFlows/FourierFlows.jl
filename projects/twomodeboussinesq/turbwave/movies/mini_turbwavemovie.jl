include("../../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: mode1u, mode1w 

@pyimport mpl_toolkits.axes_grid1 as pltgrid


plotpath = "./plots"
plotname = "turbwave_141f"
filename = "plane_example_141f_n512_ep10_Ro20_nk16_nnu4_nu1e+12.jld2"
filepath = "../data"
fullfilename = joinpath(filepath, filename)

# Plot parameters
Ro0 = 0.1 
ε0  = 2.0


# Recreate problem
jldopen(fullfilename, "r") do file

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

  x, y = kw*prob.grid.X, kw*prob.grid.Y

  lim = maximum(x)

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
     #ε = mode1u(prob)*kw/σ
     ε = mode1u(prob)/Umax


    close("all")
    fig, axs = subplots(ncols=1, nrows=2, sharex=true, sharey=true, 
      figsize=(5, 10))

    axes(axs[1]); axis("equal")
    plot_Ro = pcolormesh(x, y, Ro, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

    axes(axs[2]); axis("equal")
    plot_ε = pcolormesh(x, y, ε, cmap="RdBu_r", vmin=-ε0, vmax=ε0)


    plots = [plot_Ro, plot_ε]
    cbs = []
    for (i, ax) in enumerate(axs)
      ax[:set_adjustable]("box-forced")
      ax[:set_xlim](-lim, lim)
      ax[:set_ylim](-lim, lim)
      ax[:axis]("off")
    end

    tight_layout()

    savename = @sprintf("%s_%06d.png", joinpath("./plots", plotname), istep)
    println(savename)
    savefig(savename, dpi=240)
  end

end
