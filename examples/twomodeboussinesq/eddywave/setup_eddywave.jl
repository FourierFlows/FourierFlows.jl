__precompile__()

include("../../../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w, 
  wave_induced_speed, wave_induced_psi, wave_induced_uv, lagrangian_mean_uv,
  calc_chi, calc_chi_uv, totalenergy, mode0energy, mode1energy, CFL




struct EddyWave
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




function EddyWave(L, α, ε, Ro, Reddy;
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
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "f",        f             )
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "N",        N             )
  msg *= @sprintf("% 12s: %.2f \n",       "σ/f",      σ/f       )
  msg *= @sprintf("% 12s: %.2e km\n",     "m^-1",     1e-3/m        )
  msg *= @sprintf("% 12s: %.2e km\n",     "N/fm",     1e-3*N/(f*m)  )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν0",       ν0, nν0     )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν1",       ν1, nν1     )

  println(msg)

  EddyWave(L, f, N, nkw, α, ε, Ro, Reddy, σ, kw, m, twave, 
    nν0, nν1, dtfrac, ν0frac, ν1frac, dt, ν0, ν1, uw, nsteps, nsubs)
end
 



function eddywavesetup(n, ew::EddyWave; perturbwavefield=false)
  prob = TwoModeBoussinesq.InitialValueProblem(
    nx=n, Lx=ew.L, ν0=ew.ν0, nν0=ew.nν0, ν1=ew.ν1, nν1=ew.nν1, 
    f=ew.f, N=ew.N, m=ew.m, dt=ew.dt)

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
  
  prob, diags
end




""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makeplot!(axs, prob, ew, x, y, savename; eddylim=nothing, 
  message=nothing, save=false, show=false)

  if eddylim == nothing
    eddylim = maximum(x)
  end

  # Some limits
  Z00 = ew.Ro*ew.f
  w00 = ew.uw*ew.kw/(2*ew.m)
  U00 = ew.Ro*ew.Reddy*ew.f

  # Quantities to plot
  qc     = mode1apv(prob)/ew.f
  q      = real.(qc+conj.(qc))
  Q      = mode0apv(prob)/ew.f
  w      = mode1w(prob)
  spw    = wave_induced_speed(ew.σ, prob)
  uw, vw = wave_induced_uv(ew.σ, prob)
  uL, vL = lagrangian_mean_uv(ew.σ, prob)
  psiw   = wave_induced_psi(ew.σ, prob)
  nlresw = (real.(prob.vars.zeta - ew.m^2*ew.f/ew.N^2*prob.vars.p)
            /maximum(abs.(prob.vars.zeta)))

  chi = calc_chi(prob)
  uchi, vchi = calc_chi_uv(prob)
  spchi = sqrt.(uchi.^2+vchi.^2)

  spL = sqrt.(uL.^2+vL.^2)
  uL, vL = uL/maximum(spL), vL/maximum(spL)

  # Plot
  axs[1, 1][:cla]()
  axs[1, 2][:cla]()
  axs[2, 1][:cla]()
  axs[2, 2][:cla]()

  axes(axs[1, 1])
  axis("equal")
  pcolormesh(x, y, Q, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)
    

  axes(axs[1, 2])
  axis("equal")
  pcolormesh(x, y, w, cmap="RdBu_r", vmin=-4w00, vmax=4w00)


  axes(axs[2, 1])
  axis("equal")
  pcolormesh(x, y, prob.vars.Z, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)
  contour(x, y, psiw, 10, colors="k", linewidths=0.2, α=0.5)


  axes(axs[2, 2])
  axis("equal")
  pcolormesh(x, y, spw, cmap="YlGnBu_r", vmin=0.0, vmax=U00)
  contour(x, y, psiw, 10, colors="w", linewidths=0.2, α=0.5)




  for ax in axs
    ax[:set_xlim](-eddylim, eddylim)
    ax[:set_ylim](-eddylim, eddylim)
    ax[:tick_params](axis="both", which="both", length=0)
  end

  axs[2, 1][:set_xlabel](L"x/R")
  axs[2, 2][:set_xlabel](L"x/R")

  axs[1, 1][:set_ylabel](L"x/R")
  axs[2, 1][:set_ylabel](L"x/R")

  axs[1, 1][:set_xticks]([])
  axs[1, 2][:set_xticks]([])

  axs[1, 2][:set_yticks]([])
  axs[2, 2][:set_yticks]([])

  if message != nothing
    text(0.00, 1.03, message, transform=axs[1, 1][:transAxes], fontsize=14)
  end

  tight_layout(rect=(0.00, 0.00, 0.95, 0.95))

  if show
    pause(0.1)
  end

  if save
    savefig(savename, dpi=240)
  end

  nothing
end




""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makeplotwithquiver!(axs, prob, ew, x, y, savename; eddylim=nothing, 
  message=nothing, save=false, show=false)

  if eddylim == nothing
    eddylim = maximum(x)
  end

  # Some limits
  Z00 = ew.Ro*ew.f
  w00 = ew.uw*ew.kw/(2*ew.m)
  U00 = ew.Ro*ew.Reddy*ew.f

  # Quantities to plot
  qc     = mode1apv(prob)/ew.f
  q      = real.(qc+conj.(qc))
  Q      = mode0apv(prob)/ew.f
  w      = mode1w(prob)
  spw    = wave_induced_speed(ew.σ, prob)
  uw, vw = wave_induced_uv(ew.σ, prob)
  uL, vL = lagrangian_mean_uv(ew.σ, prob)
  psiw   = wave_induced_psi(ew.σ, prob)
  nlresw = (real.(prob.vars.zeta - ew.m^2*ew.f/ew.N^2*prob.vars.p)
            /maximum(abs.(prob.vars.zeta)))

  spL = sqrt.(uL.^2+vL.^2)
  uL, vL = uL/maximum(spL), vL/maximum(spL)

  # Plot
  axs[1, 1][:cla]()
  axs[1, 2][:cla]()
  axs[2, 1][:cla]()
  axs[2, 2][:cla]()

  axes(axs[1, 1])
  axis("equal")
  pcolormesh(x, y, Q, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)

  skip = floor(Int, prob.grid.n/32)
  quiverplot = quiver(
    xr[1:skip:end, 1:skip:end], yr[1:skip:end, 1:skip:end],
    uL[1:skip:end, 1:skip:end], vL[1:skip:end, 1:skip:end], 
    units="x", α=0.2, scale=2.0, scale_units="x")

    

  axes(axs[1, 2])
  axis("equal")
  pcolormesh(x, y, w, cmap="RdBu_r", vmin=-4w00, vmax=4w00)


  axes(axs[2, 1])
  axis("equal")
  pcolormesh(x, y, q, cmap="RdBu_r", vmin=-ew.Ro, vmax=ew.Ro)


  axes(axs[2, 2])
  axis("equal")
  pcolormesh(x, y, spw, cmap="YlGnBu_r", vmin=0.0, vmax=U00)

  contour(x, y, psiw, 10, colors="w", linewidths=0.2, α=0.5)




  for ax in axs
    ax[:set_xlim](-eddylim, eddylim)
    ax[:set_ylim](-eddylim, eddylim)
    ax[:tick_params](axis="both", which="both", length=0)
  end

  axs[2, 1][:set_xlabel](L"x/R")
  axs[2, 2][:set_xlabel](L"x/R")

  axs[1, 1][:set_ylabel](L"x/R")
  axs[2, 1][:set_ylabel](L"x/R")

  axs[1, 1][:set_xticks]([])
  axs[1, 2][:set_xticks]([])

  axs[1, 2][:set_yticks]([])
  axs[2, 2][:set_yticks]([])

  if message != nothing
    text(0.00, 1.03, message, transform=axs[1, 1][:transAxes], fontsize=14)
  end

  tight_layout(rect=(0.00, 0.00, 0.95, 0.95))

  if show
    pause(0.1)
  end

  if save
    savefig(savename, dpi=240)
  end

  nothing
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
