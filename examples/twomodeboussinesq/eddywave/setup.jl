__precompile__()

include("../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w, 
  wave_induced_speed, wave_induced_psi, wave_induced_uv, lagrangian_mean_uv,
  calc_chi, calc_chi_uv, totalenergy, mode0energy, mode1energy, CFL,
  lagrangian_mean_psih

@pyimport mpl_toolkits.axes_grid1 as pltgrid




""" The EddyWave type organizes parameters in the eddy/wave simulation."""
struct EddyWave
  name
  L
  f
  N
  nkw
  α
  ε
  Ro 
  Reddy
  σ
  kw
  m
  twave
  nν0
  nν1
  dtfrac
  ν0frac
  ν1frac
  dt
  ν0
  ν1
  uw
  nsteps
  nsubs
end




function EddyWave(name, L, α, ε, Ro, Reddy;
  f=1e-4, N=5e-3, nkw=16, nν0=8, nν1=8, dtfrac=5e-2, ν0frac=1e-1, 
  ν1frac=1e-1, ν0=nothing, ν1=nothing, nperiods=400, nsubperiods=1)

  σ = f*sqrt(1+α)
  kw = 2π*nkw/L
  m = N*kw/(f*sqrt(α))
  twave = 2π/σ                      

  dt = dtfrac * twave
  if ν0 == nothing; ν0 = ν0frac/(dt*(0.65π*n/L)^nν0); end
  if ν1 == nothing; ν1 = ν1frac/(dt*(0.65π*n/L)^nν1); end   

  uw = minimum([ε*Reddy*σ, ε*σ/kw])

  nsteps = round(Int, nperiods*twave/dt)
  nsubs = round(Int, nsubperiods*twave/dt)

  msg = "\n"
  msg *= @sprintf("% 12s: %.3f\n",        "ε",        ε             )
  msg *= @sprintf("% 12s: %.3f\n",        "Ro",       Ro            )
  msg *= @sprintf("% 12s: %.3f \n",       "σ/f",      σ/f           )
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "f",        f             )
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "N",        N             )
  msg *= @sprintf("% 12s: %.1f m\n",      "m^-1",     1/m           )
  msg *= @sprintf("% 12s: %.2e km\n",     "N/fm",     1e-3*N/(f*m)  )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν0",       ν0, nν0       )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν1",       ν1, nν1       )

  println(msg)

  EddyWave(name, L, f, N, nkw, α, ε, Ro, Reddy, σ, kw, m, twave, 
    nν0, nν1, dtfrac, ν0frac, ν1frac, dt, ν0, ν1, uw, nsteps, nsubs)
end
 



""" Save EddyWave parameters. """
function saveeddywave(ew::EddyWave, filename::String)
  jldopen(filename, "a+") do file
      for name in fieldnames(ew)
        field = getfield(ew, name)
        file["eddywaveparams/$name"] = field
      end
  end
  nothing
end




function eddywavesetup(n, ew::EddyWave; perturbwavefield=false, 
  passiveapv=false)

  if passiveapv
    prob = TwoModeBoussinesq.PassiveAPVInitialValueProblem(
      nx=n, Lx=ew.L, ν0=ew.ν0, nν0=ew.nν0, ν1=ew.ν1, nν1=ew.nν1, 
      f=ew.f, N=ew.N, m=ew.m, dt=ew.dt, σ=ew.σ)
  else
    prob = TwoModeBoussinesq.InitialValueProblem(
      nx=n, Lx=ew.L, ν0=ew.ν0, nν0=ew.nν0, ν1=ew.ν1, nν1=ew.nν1, 
      f=ew.f, N=ew.N, m=ew.m, dt=ew.dt)
  end

  # Initial condition
  x, y = prob.grid.X, prob.grid.Y
  Z0 = ew.f*ew.Ro * exp.(-(x.^2+y.^2)/(2*ew.Reddy^2))
  TwoModeBoussinesq.set_Z!(prob, Z0)

  # ε = U/(Reddy*σ) or U*kw/σ
  TwoModeBoussinesq.set_planewave!(prob, ew.uw, ew.nkw)

  # Use a regular perturbation expansion solution to the mode-1
  # PV equation to reduce the initial mode-1 APV
  if perturbwavefield && ew.α > 0.0
    g = prob.grid

    p1 = eigenperturbation!(prob.vars.p, ew, prob.vars.Z, g)
    p = prob.vars.p + p1

    # Determine u, v from linear equations. Will fail if α=0.
    ph = fft(p)
    uh = -1.0/(ew.α*ew.f^2) * (-ew.σ*g.K .- ew.f*im*g.L).*ph 
    vh = -1.0/(ew.α*ew.f^2) * (-ew.σ*g.L .+ ew.f*im*g.K).*ph 

    u = ifft(uh)
    v = ifft(vh)

    TwoModeBoussinesq.set_uvp!(prob, u, v, p)
  end

  etot = Diagnostic(totalenergy, prob; nsteps=ew.nsteps)
  e0   = Diagnostic(mode0energy, prob; nsteps=ew.nsteps)
  e1   = Diagnostic(mode1energy, prob; nsteps=ew.nsteps) 
  diags = [etot, e0, e1]

# Output
  fileprefix = @sprintf("./data/%s_nu%.0e, %df_n%d_ep%02d_Ro%02d_nkw%02d",
    ew.name, ew.ν0, 100*ew.σ/ew.f, n, 100ε, 100Ro, ew.nkw)

  i, testprefix = 0, fileprefix
  while isfile(testprefix*".jld2"); i+=1; testprefix=fileprefix*"-$i"; end
  filename = testprefix * ".jld2"

  saveproblem(prob, filename)
  saveeddywave(ew, filename)

  getsolr(prob) = prob.vars.solr
  getsolc(prob) = prob.vars.solc

  output = Output(prob, filename, (:solr, getsolr), (:solc, getsolc))
    
  prob, diags, output
end




""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makefourplot(prob, ew; eddylim=nothing, 
  message=nothing, save=false, show=false, passiveapv=false)

  plotpath = "./plots"
  savename = @sprintf("%s_nu%.0e_n%d_%02df_ep%02d_Ro%02d_%06d.png", 
    joinpath(plotpath, ew.name), ew.ν0, prob.grid.nx, floor(Int, 100ew.σ/ew.f), 
    floor(Int, 100*ew.ε), floor(Int, 100*ew.Ro), prob.step)

  close("all")
  fig, axs = subplots(ncols=2, nrows=2, figsize=(8, 8)) 

  x, y = prob.grid.X/ew.Reddy, prob.grid.Y/ew.Reddy

  if eddylim == nothing
    eddylim = maximum(x)
  end

  # Some limits
  Z00 = ew.Ro*ew.f
  w00 = ew.uw*ew.kw/(2*ew.m)
  U00 = ew.Ro*ew.Reddy*ew.f

  # Quantities to plot
  Q      = mode0apv(prob)/ew.f
  w      = mode1w(prob)
  spw    = wave_induced_speed(ew.σ, prob)
  psiw   = wave_induced_psi(ew.σ, prob)
  uL, vL = lagrangian_mean_uv(ew.σ, prob)
  PsiLh  = lagrangian_mean_psih(ew.σ, prob)
  PsiL   = irfft(PsiLh, prob.grid.nx)

  # Plot
  axes(axs[1, 1])
  axis("equal")
  pcolormesh(x, y, Q, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)
    

  axes(axs[1, 2])
  axis("equal")
  pcolormesh(x, y, w, cmap="RdBu_r", vmin=-4w00, vmax=4w00)


  axes(axs[2, 1])
  axis("equal")
  if passiveapv
    pcolormesh(x, y, prob.vars.Q/ew.f, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)

    contour(x, y, PsiL, 20, colors="k", linewidths=0.2, alpha=0.5)

    nquiv = 32
    iquiv = floor(Int, prob.grid.nx/nquiv)
    quiverplot = quiver(
      x[1:iquiv:end, 1:iquiv:end], y[1:iquiv:end, 1:iquiv:end],
      uL[1:iquiv:end, 1:iquiv:end], vL[1:iquiv:end, 1:iquiv:end], 
      units="x", alpha=0.2, scale=2.0, scale_units="x")
  else
    pcolormesh(x, y, prob.vars.Z/ew.f, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)
  end


  axes(axs[2, 2])
  axis("equal")
  pcolormesh(x, y, spw, cmap="YlGnBu_r", vmin=0.0, vmax=U00)
  contour(x, y, psiw, 10, colors="w", linewidths=0.2, alpha=0.5)


  for ax in axs
    ax[:set_adjustable]("box-forced")
    ax[:set_xlim](-eddylim, eddylim)
    ax[:set_ylim](-eddylim, eddylim)
    ax[:tick_params](axis="both", which="both", length=0)
  end

  axs[2, 1][:set_xlabel](L"x/R")
  axs[2, 2][:set_xlabel](L"x/R")

  axs[1, 1][:set_ylabel](L"y/R")
  axs[2, 1][:set_ylabel](L"y/R")

  axs[1, 1][:set_xticks]([])
  axs[1, 2][:set_xticks]([])

  axs[1, 2][:set_yticks]([])
  axs[2, 2][:set_yticks]([])

  if message != nothing
    text(0.00, 1.03, message, transform=axs[1, 1][:transAxes], fontsize=14)
  end

  tight_layout(rect=(0.05, 0.05, 0.95, 0.95))

  if show; pause(0.1); end
  if save; savefig(savename, dpi=240); end

  nothing
end




""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makethreeplot!(axs, prob, ew, x, y, savename; eddylim=nothing, 
  message=nothing, save=false, show=false, makecolorbar=false, 
  remakecolorbar=nothing) 
  
  if eddylim == nothing
    eddylim = maximum(x)
  end


  # Some limits
  Z00 = ew.Ro*ew.f
  w00 = ew.uw*ew.kw/(2*ew.m)
  U00 = ew.Ro*ew.Reddy*ew.f

  # Quantities to plot
  Q      = mode0apv(prob)/ew.f
  w      = mode1w(prob)
  spw    = wave_induced_speed(ew.σ, prob)
  psiw   = wave_induced_psi(ew.σ, prob)


  # Plot
  axs[1][:cla]()
  axs[2][:cla]()
  axs[3][:cla]()

  if remakecolorbar != nothing
    for cb in cbs; cb[:ax][:cla](); end
    makecolorbar = true
  end


  axes(axs[1])
  axis("equal")
  Qplot = pcolormesh(x, y, Q, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)
    

  axes(axs[2])
  axis("equal")
  wplot = pcolormesh(x, y, w, cmap="RdBu_r", vmin=-6w00, vmax=6w00)


  axes(axs[3])
  axis("equal")
  spplot = pcolormesh(x, y, spw, cmap="YlGnBu_r", vmin=0.0, vmax=U00)
  contour(x, y, psiw, 10, colors="w", linewidths=0.2, alpha=0.5)


  plots = [Qplot, wplot, spplot]
  cbs = []
  for (i, ax) in enumerate(axs)

    ax[:set_adjustable]("box-forced")
    ax[:set_xlim](-eddylim, eddylim)
    ax[:set_ylim](-eddylim, eddylim)
    ax[:tick_params](axis="both", which="both", length=0)

    ax[:set_xlabel](L"x/R")

    if makecolorbar
      divider = pltgrid.make_axes_locatable(ax)
      cax = divider[:append_axes]("top", size="5%", pad="10%")
      push!(cbs, colorbar(plots[i], cax=cax, orientation="horizontal"))

      cbs[i][:ax][:xaxis][:set_ticks_position]("top")
      cbs[i][:ax][:xaxis][:set_label_position]("top")
      cbs[i][:ax][:tick_params](axis="x", which="both", length=0)

      locator_params(tight=true, nbins=5)
    end
  end


  axs[1][:set_ylabel](L"y/R")
  axs[2][:set_yticks]([])
  axs[3][:set_yticks]([])

  cbs[1][:set_ylabel](L"f^{-1}Q")
  cbs[2][:set_ylabel](L"w + w^*")
  cbs[3][:set_ylabel](L"| \nabla \Psi^w |^2")


  if message != nothing
    msg = text(0.00, 1.30, message, transform=axs[1][:transAxes], fontsize=14)
  else
    msg = nothing
  end

  tight_layout(rect=(0.00, 0.05, 1.00, 0.90))

  if show
    pause(0.1)
  end

  if save
    savefig(savename, dpi=240)
  end

  cbs, msg
end





""" Efficiently solve one iteration of the non-constant coeff,
APV-neutralizing A-eqn. """
function eigenperturbation!(p0, ew, Z, g)
  c = ew.α*ew.f*ew.m^2/ew.N^2

  # Invert in Fourier space
  p1h = c * fft(Z.*p0) ./ (-g.KKsq + ew.f*c)
  p1h[.!isfinite.(p1h)] = 0im  # Zero-out singular wavenumbers

  ifft(p1h)
end
