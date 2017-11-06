include("./setup.jl")

# -- Parameters --
  nkw = 8
    n = 512
    L = 2π*100e3*nkw
    α = 0.2             # Frequency parameter
    ε = 1e-1            # Wave amplitude
   Ro = 1e-1            # Eddy Rossby number
 name = "turbwave"


# Setup
tw, prob, diags, outs = turbwavesetup(name, n, L, α, ε, Ro; nkw=nkw,
   dtfrac=1e-2, nsubperiods=4, nν0=6, nν1=6, ν0=1e20, ν1=1e6) 

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
    "\$t=% 3d\$ wave periods, \$E_0=%.3f\$, ",
    round(Int, prob.t/tw.twave), e0.value/e0.data[1])
    
  plotmsg2 = @sprintf(
    "\$E_1=%.3f\$, \$E_{\\mathrm{tot}}=%.6f\$",
    e1.value/e1.data[1], etot.value/etot.data[1])

  plotmsg = plotmsg1*plotmsg2

  makefourplot(prob, tw; save=true, message=plotmsg)

end
