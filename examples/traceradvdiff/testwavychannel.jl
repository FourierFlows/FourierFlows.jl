include("./wavychannel.jl")


nx = 512
ep = 0.1
kap = 1e-4
dc = 0.02

try
  run(`rm ./'*'.jld2`)
end
prob, yvar = wavy_barotropic_tracer_run(ep, kap; nx=nx, dc=dc)

yvar_theory = 2*kap*yvar.time[end] + dc^2.0

println("enhancement: ", yvar.value/yvar_theory)
