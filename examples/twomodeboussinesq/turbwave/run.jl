include("./setup.jl")

# -- Parameters --
  nkw = 16
    n = 384
    L = 2π*100e3*nkw
    α = 1.0             # Frequency parameter
    ε = 2e-1            # Wave amplitude
   Ro = 1e-1            # Eddy Rossby number
 name = "turbwave"


# Setup
tw, prob, diags, outs = turbwavesetup(name, n, L, α, ε, Ro; nkw=nkw,
  dtfrac=5e-3, nsubperiods=2, nν0=8, nν1=8, ν0=1e32, ν1=1e24,
  k0turb=floor(Int, 2n/3))

etot, e0, e1 = diags[1], diags[2], diags[3]


# Run
startwalltime = time()
while prob.step < tw.nsteps

  stepforward!(prob, diags; nsteps=tw.nsubs)
  TwoModeBoussinesq.updatevars!(prob)
  saveoutput(outs)

  log1 = @sprintf(
    "step: %04d, t: %d, max Ro: %.4f, ", 
    prob.step, prob.t/tw.twave, maximum(abs.(prob.vars.Z))/tw.f)

  log2 = @sprintf(
    "ΔE: %.3f, Δe: %.3f, Δ(E+e): %.6f, τ: %.2f min",
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1],
    (time()-startwalltime)/60)

  println(log1*log2)

  plotmsg1 = @sprintf(
    "\$t=% 3d\$ wave periods, \$\\Delta E=%.3f\$, ",
    round(Int, prob.t/tw.twave), e0.value/e0.data[1])
    
  plotmsg2 = @sprintf(
    "\$\\Delta e=%.3f\$, \$\\Delta (E+e)=%.6f\$",
    e1.value/e1.data[1], etot.value/etot.data[1])

  plotmsg = plotmsg1*plotmsg2

  makefourplot(prob, tw; save=true, message=plotmsg)

end
