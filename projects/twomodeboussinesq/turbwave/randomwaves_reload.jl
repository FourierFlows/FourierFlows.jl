include("./setup.jl")

import FourierFlows.TwoModeBoussinesq: mode1speed

# -- Parameters --
 name = "plane_example"
  nkw = 16
    ε = 1e-1            # Wave amplitude
   Ro = 2e-1            # Eddy Rossby number
    n = 512
    L = 2π*1600e3
    f = 1e-4
    N = 5e-3
   kw = nkw*2π/L
    α = 1.0 #(N*kw)^2/(f*m)^2
    m = 1/2000


# Setup
tw, prob, diags, _ = turbwavesetup(name, n, L, α, ε, Ro;
  nkw=nkw, f=f, N=N, m=m,
  dtfrac = 1e-2, nperiods = 400, nsubperiods = 1, k0turb = floor(Int, 2n/3), 
  nν0 = 4, nν1 = 6, ν0 = 1e12, ν1 = 1e20,
  nowavecubics = true, wavefield = "planar")

etot, e0, e1 = diags[1], diags[2], diags[3]

filename = "./data/plane_example_141f_n512_ep10_Ro20_nk16_nnu4_nu1e+12.jld2"
step = 9000

file = jldopen(filename, "r+")
   t = file["timeseries/t/$step"]
  Qh = file["timeseries/solr/$step"]
solc = file["timeseries/solc/$step"]
close(file)

u = ifft(solc[:, :, 1])
v = ifft(solc[:, :, 2])
p = ifft(solc[:, :, 3])
Q = irfft(Qh, n)
TwoModeBoussinesq.set_Q!(prob, Q)
TwoModeBoussinesq.set_uvp!(prob, u, v, p)
prob.vars.t = t
prob.ts.step = step

update!(diags)

getsolr(prob) = prob.vars.solr
getsolc(prob) = prob.vars.solc

outs = Output(prob, filename, (:solr, getsolr), (:solc, getsolc))


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

  if prob.step/tw.nsubs % 2 == 0
    makefourplot(prob, tw; save=true, message=plotmsg)
  end

end
