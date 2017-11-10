include("./setup.jl")

import FourierFlows.TwoModeBoussinesq: mode1speed

# -- Parameters --
  nkw = 16 
    n = 256
    L = 2π*1600e3
    α = 1.0             # Frequency parameter
    ε = 1e-1            # Wave amplitude
   Ro = 1e-1            # Eddy Rossby number
    m = 500/2π
 name = "isotropic"


# Setup
tw, prob, diags, outs = turbwavesetup(name, n, L, α, ε, Ro; nkw=nkw,
  dtfrac = 2e-2, nsubperiods = 1, k0turb = floor(Int, 2n/3), 
  nν0 = 6, nν1 = 6, ν0 = 1e20, ν1 = 1e20, wavecubics = false
)

# Spectral shape
k0, dk = 16*2π/L, 1.0*2π/L; 
kin, kout = k0-dk/2, k0+dk/2

amplitude(k, l) = amplitude(sqrt(k^2+l^2))
amplitude(K) = (abs(K) >= kin && abs(K) <= kout) ? 1.0 : 0.0

# Normalization by initial wavefield kinetic energy
σ = sqrt(prob.params.f^2 + prob.params.N^2*k0^2/prob.params.m^2)
KE = (ε*σ/k0)^2
maxspeed = ε*σ/k0

TwoModeBoussinesq.set_randomwavefield!(prob.vars, prob.params, prob.grid,
  amplitude; KE=KE, maxspeed=maxspeed)
update!(diags)

#fig, axs = subplots()
#imshow(real.(prob.vars.u))
#println(maximum(mode1speed(prob)*k0/σ))
#show()
 
etot, e0, e1 = diags[1], diags[2], diags[3]


# Run
startwalltime = time()
while prob.step < tw.nsteps

  stepforward!(prob, diags; nsteps=tw.nsubs)
  saveoutput(outs)

  log1 = @sprintf("step: %04d, t: %d, ", prob.step, prob.t/tw.tσ)
  log2 = @sprintf("ΔE: %.3f, Δe: %.3f, Δ(E+e): %.6f, τ: %.2f min",
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1],
    (time()-startwalltime)/60)

  println(log1*log2)

  plotmsg1 = @sprintf("\$t=% 3d\$ wave periods, \$\\Delta E=%.3f\$, ",
    round(Int, prob.t/tw.tσ), e0.value/e0.data[1])
  plotmsg2 = @sprintf("\$\\Delta e=%.3f\$, \$\\Delta (E+e)=%.6f\$",
    e1.value/e1.data[1], etot.value/etot.data[1])

  plotmsg = plotmsg1*plotmsg2

  if prob.step % 1000 == 0.0
    makefourplot(prob, tw; save=true, message=plotmsg)
  end

end
