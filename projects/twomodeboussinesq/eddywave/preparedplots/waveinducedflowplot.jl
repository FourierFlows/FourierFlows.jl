include("./plottools")

using PyPlot, PyCall

import SLSPlotTools

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w, 
  wave_induced_speed, wave_induced_psi, wave_induced_uv, lagrangian_mean_uv,
  calc_chi, calc_chi_uv, totalenergy, mode0energy, mode1energy, CFL,
  lagrangian_mean_psi, apv_induced_psi

@pyimport mpl_toolkits.axes_grid1 as pltgrid

# Names
   simname = "nk32_nu1e+32_122f_n512_ep20_Ro05_nkw32"
plotprefix = "./plots/mitsls_"
  filename = "../data/$simname.jld2"

prob, steps = SLSPlotTools.load_boussinesq_ivp(filename)

# Plot parameters
eddylim = 10
    Ro0 = 0.04
    w0  = 0.08
    sp0 = 1.0
  istep = 40

# Non-retrievable parameters
nkw = 32
 kw = 2π*nkw/prob.grid.Lx
  α = (prob.params.N*kw)^2/(prob.params.f*prob.params.m)^2
  σ = prob.params.f*sqrt(1+α)
 tσ = 2π/σ
 Re = prob.grid.Lx/40


jldopen(filename, "r") do file

    x, y = prob.grid.X/Re, prob.grid.Y/Re
    step = steps[istep]
    t = file["timeseries/t/$step"]
    Zh = file["timeseries/solr/$step"]
    solc = file["timeseries/solc/$step"]

    Z = irfft(Zh, prob.grid.nx) 
    u = ifft(solc[:, :, 1])
    v = ifft(solc[:, :, 2])
    p = ifft(solc[:, :, 3])

    TwoModeBoussinesq.set_Z!(prob, Z)
    TwoModeBoussinesq.set_uvp!(prob, u, v, p)

    Umax = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2))
    
    w = mode1w(prob)
    Q = mode0apv(prob)
    sp = wave_induced_speed(σ, prob)/Umax
    Psiw = wave_induced_psi(σ, prob)
    PsiL = lagrangian_mean_psi(σ, prob)
    PsiQ = apv_induced_psi(prob)
    Chi  = calc_chi(prob)
    uL, vL = lagrangian_mean_uv(σ, prob)
    uc, vc = calc_chi_uv(prob)

    spc = sqrt.(uc.^2.0+vc.^2.0)/Umax

    close("all")
    fig, axs = subplots(ncols=3, nrows=1, sharex=true, sharey=true, 
      figsize=(12, 5))

    axes(axs[1]); axis("equal")
    plot1 = pcolormesh(x, y, Q/f, 
      cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

    #contour(x, y, Psiw, 20, colors="k", linewidths=0.2, alpha=0.5)
    #contour(x, y, PsiQ, 20, colors="m", linewidths=0.5, alpha=0.5)
    #contour(x, y, Chi, 8, colors="k", linewidths=0.5, alpha=0.2)

    contour(x, y, prob.vars.Psi, 15, colors="k", linewidths=0.5, alpha=0.4)

    axes(axs[2]); axis("equal")
    plot2 = pcolormesh(x, y, spc,
      cmap="YlGnBu_r", vmin=0.0, vmax=sp0)
    contour(x, y, Chi, 15, colors="w", linewidths=0.5, alpha=0.3)

    axes(axs[3]); axis("equal")
    plot3 = pcolormesh(x, y, w, 
      cmap="RdBu_r", vmin=-w0, vmax=w0)


    plots = [plot1, plot2, plot3]
    cbs = []
    for (i, ax) in enumerate(axs)
      ax[:set_adjustable]("box-forced")
      ax[:set_xlim](-eddylim, eddylim)
      ax[:set_ylim](-eddylim, eddylim)
      ax[:tick_params](axis="both", which="both", length=0)

      divider = pltgrid.make_axes_locatable(ax)
      cax = divider[:append_axes]("top", size="5%", pad="5%")
      cb = colorbar(plots[i], cax=cax, orientation="horizontal")

      cb[:ax][:xaxis][:set_ticks_position]("top")
      cb[:ax][:xaxis][:set_label_position]("top")
      cb[:ax][:tick_params](axis="x", which="both", length=0)

      push!(cbs, cb)
      iscolorbar = true
    end

    axs[1][:set_xlabel](L"x/R")
    axs[2][:set_xlabel](L"x/R")
    axs[3][:set_xlabel](L"x/R")
    axs[1][:set_ylabel](L"y/R")

    cbs[1][:set_ticks]([-Ro0, 0.0, Ro0])
    cbs[2][:set_ticks]([-sp0, 0.0, sp0])
    cbs[3][:set_ticks]([-w0, 0.0, w0])

    cbs[1][:set_label](L"Q/f", labelpad=12.0)
    cbs[2][:set_label](L"| \nabla \Psi^w |^2", labelpad=12.0)
    cbs[3][:set_label](L"\hat w(z=0) = w + w^*", labelpad=12.0)

    msg = @sprintf("\$t = %02d\$ wave periods", t/tσ)
    figtext(0.51, 0.95, msg, horizontalalignment="center", fontsize=14)

    tight_layout(rect=(0.00, 0.00, 0.95, 0.90), h_pad=0.05)
end
