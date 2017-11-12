__precompile__()

include("../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoDTurb

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w, 
  chi, UVchi, totalenergy, mode0energy, mode1energy, CFL, mode1u, get_Z

@pyimport mpl_toolkits.axes_grid1 as pltgrid




struct TurbWave
  name
  n
  L
  f
  N
  m
  ε
  uw
  kw
  α
  σ
  Ro
  Lturb
  dt
  nsteps
  nsubs
  tσ
  initialwavefield
end






function turbwavesetup(name, n, L, α, ε, Ro;
  f=1e-4, N=5e-3, nkw=16, dkw=1, nν0=8, nν1=8, dtfrac=5e-2, ν0frac=1e-1, 
  ν1frac=1e-1, ν0=nothing, ν1=nothing, nperiods=400, nsubperiods=1, 
  k0turb=n/2, m=32N/(f*L), nowavecubics=true, wavefield="planar")


  # Initialize problem
  σ = f*sqrt(1+α)
  tσ = 2π/σ                      

  if α == 0.0 # NIW case
    nkw = 0 
    kw = (f*m/N)^2
  else
    kw = 2π*nkw/L
    m = N*kw/(f*sqrt(α))
  end

  dt = dtfrac * tσ
  if ν0 == nothing; ν0 = ν0frac/(dt*(0.65π*n/L)^nν0); end
  if ν1 == nothing; ν1 = ν1frac/(dt*(0.65π*n/L)^nν1); end   

  nsteps = round(Int, nperiods*tσ/dt)
  nsubs = round(Int, nsubperiods*tσ/dt)

  savename = @sprintf("./data/twodturb_n%d_Ro%02d_nnu%d_nu%.0e.jld2", 
    n, 100Ro, nν0, ν0)

  # Generate initial turbulence field
  if isfile(savename) 
    @load savename Z
  else
    qf = f*Ro
    Z = TwoDTurb.makematureturb(n, L; nnu=nν0, nu=ν0, k0=k0turb, 
      qf=qf, q0=1.2*qf, tf=40/qf)
    titlname = @sprintf("twodturb: \$n=%d\$, \$n_{\\nu}=%d\$, max(Z)=%.2e", 
      n, nν0, maximum(Z))
    plotname = @sprintf("./turbplots/twodturb_n%04d_Ro%02d_nnu%d_nu%.0e.png", 
      n, 100Ro, nν0, ν0)

    fig, axs = subplots()
    imshow(Z)
    title(titlname) 
    savefig(plotname, dpi=240)

    println("\nTurbulence field generated.")
    println("plot: ", plotname)
    println("data: ", savename, "\n")

    @save savename Z
  end

  prob = TwoModeBoussinesq.PrognosticAPVInitialValueProblem(nx=n, Lx=L, 
    ν0=ν0, nν0=nν0, ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt,
    nowavecubics=nowavecubics)

  TwoModeBoussinesq.set_Q!(prob, Z)

  # Initial wave field  
  # ε = U/(Reddy*σ) or U*kw/σ
  Lturb = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2)/f)

  Umax  = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2))
  Romax = maximum(abs.(prob.vars.Q))/f

  uw = Umax*ε/Romax

  tw = TurbWave(name, n, L, f, N, m, ε, uw, kw, α, σ, Ro, Lturb, 
    dt, nsteps, nsubs, tσ, wavefield)

  if wavefield == "planar"
    TwoModeBoussinesq.set_planewave!(prob, uw, nkw)
  elseif wavefield == "isotropic"
    kin  = (nkw - dkw/2)*2π/L
    kout = (nkw + dkw/2)*2π/L

    amplitude(k, l) = amplitude(sqrt(k^2+l^2))
    amplitude(K) = (abs(K) >= kin && abs(K) <= kout) ? 1.0 : 0.0

    TwoModeBoussinesq.set_isotropicwavefield!(prob, amplitude; maxspeed=uw)
  else
    throw("""Wavefield must be "planar" or "isotropic", not $wavefield.""")
  end

  # Diagnostics
  etot = Diagnostic(totalenergy, prob; nsteps=nsteps)
  e0   = Diagnostic(mode0energy, prob; nsteps=nsteps)
  e1   = Diagnostic(mode1energy, prob; nsteps=nsteps) 
  diags = [etot, e0, e1]

  # Prepare output
  fileprefix = @sprintf("./data/%s_%df_n%d_ep%02d_Ro%02d_nk%02d_nnu%d_nu%.0e",
    name, 100*sqrt(α+1), n, 100ε, 100Ro, nkw, nν0, ν0)

  i, testprefix = 0, fileprefix
  while isfile(testprefix*".jld2"); i+=1; testprefix=fileprefix*"-$i"; end
  filename = testprefix * ".jld2"

  saveproblem(prob, filename)

  getsolr(prob) = prob.vars.solr
  getsolc(prob) = prob.vars.solc

  output = Output(prob, filename, (:solr, getsolr), (:solc, getsolc))

  
  # Print message
  msg = "\n"
  msg *= @sprintf("% 12s: %.3f\n",        "ε",        ε             )
  msg *= @sprintf("% 12s: %.3f\n",        "Ro",       Ro            )
  msg *= @sprintf("% 12s: %.3f\n",        "max u",    uw            )
  msg *= @sprintf("% 12s: %.3f\n",        "max U",    Umax          )
  msg *= @sprintf("% 12s: %.3f\n",        "Lturb",    Lturb         )
  msg *= @sprintf("% 12s: %.3f \n",       "σ/f",      σ/f           )
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "f",        f             )
  msg *= @sprintf("% 12s: %.1e s^-1\n",   "N",        N             )
  msg *= @sprintf("% 12s: %.1f m\n",      "m^-1",     1/m           )
  msg *= @sprintf("% 12s: %.2e km\n",     "N/fm",     1e-3*N/(f*m)  )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν0",       ν0, nν0       )
  msg *= @sprintf("% 12s: %.2e (n=%d)\n", "ν1",       ν1, nν1       )

  println(msg)

  tw, prob, diags, output
end




""" Efficiently solve one iteration of the non-constant coeff,
APV-neutralizing A-eqn. """
function eigenperturbation!(p0, tw, Z, g)

  α, f, m, N = tw.α, tw.f, tw.m, tw.N

  # Invert in Fourier space
  c = α*f*m^2/N^2
  p1h = fft(c*Z.*p0) ./ (-g.KKsq + f*c)
  p1h[.!isfinite.(p1h)] = 0im  # Zero-out singular wavenumbers

  ifft(p1h)
end



""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makefourplot(prob, tw; eddylim=nothing, 
  message=nothing, save=false, show=false)

  TwoModeBoussinesq.updatevars!(prob)

  close("all")
  fig, axs = subplots(ncols=2, nrows=2, figsize=(8, 8)) 

  x, y = tw.kw*prob.grid.X, tw.kw*prob.grid.Y

  savename = @sprintf("%s_n%d_%02df_ep%02d_Ro%02d_nnu%d_nu%.0e_%06d.png", 
    joinpath("./plots", tw.name), tw.n, floor(Int, 100tw.σ/tw.f), 
    floor(Int, 100*tw.ε), floor(Int, 100*tw.Ro), prob.params.nν0, 
    prob.params.ν0, prob.step)

  if eddylim == nothing; eddylim = maximum(x); end

  # Color limits
  Z00 = tw.Ro*tw.f
  w00 = tw.uw*tw.kw/(2*tw.m)
  U00 = 2*tw.uw*tw.ε
  Ro0 = 0.8*tw.Ro
  if     tw.initialwavefield == "planar";    u00 = 2.0*tw.uw
  elseif tw.initialwavefield == "isotropic"; u00 = tw.uw
  end

  # Quantities to plot
  #Q      = mode0apv(prob)/tw.f
  Z      = get_Z(prob)/tw.f
  sp     = mode1speed(prob)
  Chi    = chi(prob)
  uc, vc = UVchi(prob)
  spc    = sqrt.(uc.^2+vc.^2)


  # Plot
  axes(axs[1, 1])
  axis("equal")
  pcolormesh(x, y, Z, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)
    

  axes(axs[1, 2])
  axis("equal")
  pcolormesh(x, y, sp, cmap="YlGnBu_r", vmin=0.0, vmax=2u00)


  axes(axs[2, 1])
  axis("equal")
  pcolormesh(x, y, prob.vars.Q/tw.f, cmap="RdBu_r", vmin=-Ro0, vmax=Ro0)

  axes(axs[2, 2])
  axis("equal")
  #pcolormesh(x, y, spc, cmap="YlGnBu_r", vmin=0.0, vmax=U00)
  #contour(x, y, Chi, 10, colors="w", linewidths=0.5, alpha=0.5)
  pcolormesh(x, y, uc, cmap="RdBu_r", vmin=-U00, vmax=U00)




  for ax in axs
    ax[:set_adjustable]("box-forced")
    ax[:set_xlim](-eddylim, eddylim)
    ax[:set_ylim](-eddylim, eddylim)
    ax[:tick_params](axis="both", which="both", length=0)
  end

  axs[2, 1][:set_xlabel](L"kx")
  axs[2, 2][:set_xlabel](L"kx")

  axs[1, 1][:set_ylabel](L"ky")
  axs[2, 1][:set_ylabel](L"ky")

  axs[1, 1][:set_xticks]([])
  axs[1, 2][:set_xticks]([])

  axs[1, 2][:set_yticks]([])
  axs[2, 2][:set_yticks]([])

  if message != nothing
    text(0.00, 1.03, message, transform=axs[1, 1][:transAxes], fontsize=14)
  end

  tight_layout(rect=(0.05, 0.05, 0.95, 0.95))

  if show
    pause(0.1)
  end

  if save
    savefig(savename, dpi=240)
  end

  nothing
end


