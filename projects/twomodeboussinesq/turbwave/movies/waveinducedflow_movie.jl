include("../../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1u, mode1w, 
  chi, UVchi
  

@pyimport mpl_toolkits.axes_grid1 as pltgrid

plotname = "./plots/wif"
filename = "../data/plane_example_141f_n512_ep10_Ro20_nk16_nnu4_nu1e+12.jld2" 

# -- Plot parameters --
  eddylim = 10
      Ro0 = 0.1
       u0 = 2.00
      sp0 = 0.8
ncontours = 10


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
   Re = L/80

  #x, y = prob.grid.X/Re, prob.grid.Y/Re
  
  x, y = kw*prob.grid.X, kw*prob.grid.Y
  eddylim = maximum(x)

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

    Umax = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2))
    
    u = mode1u(prob)/Umax
    PsiQ = irfft(-prob.grid.invKKrsq.*Qh, prob.grid.nx)

    Chi  = chi(prob)
    uc, vc = UVchi(prob)
    spc = sqrt.(uc.^2.0+vc.^2.0)/Umax


    close("all")
    fig, axs = subplots(ncols=3, nrows=1, sharex=true, sharey=true, 
      figsize=(12, 5))


    axes(axs[1]); axis("equal")
    plot1 = pcolormesh(x, y, Q/f, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)
    contour(x, y, PsiQ, 2ncontours, colors="k", linewidths=0.5, alpha=0.2)


    axes(axs[2]); axis("equal")
    plot2 = pcolormesh(x, y, spc, cmap="YlGnBu_r", vmin=0.0, vmax=sp0)
    contour(x, y, Chi, 2ncontours, colors="w", linewidths=0.3, alpha=0.2)


    axes(axs[3]); axis("equal")
    plot3 = pcolormesh(x, y, u, cmap="RdBu_r", vmin=-u0, vmax=u0)


    plots = [plot1, plot2, plot3]
    cbs = []
    for (i, ax) in enumerate(axs)
      ax[:set_adjustable]("box-forced")
      ax[:set_xlim](-eddylim, eddylim)
      ax[:set_ylim](-eddylim, eddylim)
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

    axs[1][:set_xlabel](L"kx")
    axs[2][:set_xlabel](L"kx")
    axs[3][:set_xlabel](L"kx")
    axs[1][:set_ylabel](L"ky")

    ticks = [-40, -20, 0, 20, 40]
    axs[1][:xaxis][:set_ticks](ticks)
    axs[1][:yaxis][:set_ticks](ticks)
    axs[2][:xaxis][:set_ticks](ticks)
    axs[3][:xaxis][:set_ticks](ticks)

    cbs[1][:set_ticks]([-Ro0, 0.0, Ro0])
    cbs[2][:set_ticks]([0.0, sp0])
    cbs[3][:set_ticks]([-u0, 0.0, u0])

    labelpad = 16.0
    cbs[1][:set_label]("\$Q/f\$ and \$\\Phi\$", labelpad=labelpad)
    cbs[2][:set_label]("\$| \\nabla \\chi |\$ and \$\\chi\$", 
      labelpad=labelpad)
    cbs[3][:set_label](L"(u + u^*)/\max(U)", labelpad=labelpad)

    msg = @sprintf("\$t = %02d\$ wave periods", t/tσ)
    figtext(0.51, 0.95, msg, horizontalalignment="center", fontsize=14)

    tight_layout(rect=(0.10, 0.05, 0.90, 0.85), h_pad=0.05)

    #pause(0.1)
    
    savename = @sprintf("%s_%06d.png", plotname, istep)
    println(savename)
    savefig(savename, dpi=240)
  end

end
