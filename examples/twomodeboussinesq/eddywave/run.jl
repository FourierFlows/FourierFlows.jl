include("./setup.jl")



# -- Parameters --
  nkw = 64    
    n = 1024
    L = 2π*100e3*nkw
    α = 0.5             # Frequency parameter
    ε = 2e-1            # Wave amplitude
   Ro = 5e-2            # Eddy Rossby number
Reddy = L/80            # Eddy radius
   ν0 = 2e32
  nν0 = 8
 name = "example_nk64"

# Setup
ew = EddyWave(name, L, α, ε, Ro, Reddy; 
  nkw=nkw, dtfrac=1e-2, nsubperiods=1,
  nν0=nν0, ν0=ν0, nν1=8, ν1=1e16, nperiods=200) 

prob, diags, outs = eddywavesetup(n, ew, perturbwavefield=false) 
etot, e0, e1 = diags[1], diags[2], diags[3]

# Run
startwalltime = time()
while prob.step < ew.nsteps

  stepforward!(prob, diags; nsteps=ew.nsubs)
  TwoModeBoussinesq.updatevars!(prob)
  saveoutput(outs)

  walltime = (time()-startwalltime)/60

  log1 = @sprintf("step: %04d, t: %d, ", prob.step, prob.t/ew.twave)
  log2 = @sprintf("ΔE: %.3f, Δe: %.3f, Δ(E+e): %.6f, τ: %.2f min",
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1],
    walltime)

  println(log1*log2)

  plotmsg1 = @sprintf("\$t=% 3d\$ wave periods, \$\\Delta E=%.3f\$, ",
    round(Int, prob.t/ew.twave), e0.value/e0.data[1])
  plotmsg2 = @sprintf("\$\\Delta e=%.3f\$, \$\\Delta(E+e)=%.6f\$",
    e1.value/e1.data[1], etot.value/etot.data[1])

  makefourplot(prob, ew; message=plotmsg1*plotmsg2, save=true) 
end
