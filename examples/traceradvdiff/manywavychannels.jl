include("./wavychannel.jl")

using JLD2

filename = "manyeps.jld2"

nx  = 1024
kap = 1e-4
dcx = 0.01
dcy = 0.02

epz = 0.0:0.1:0.9
nepz = length(epz)
enhancement = zeros(nepz)

for i = 1:nepz

  ep = epz[i]

  prob, yvar = wavy_barotropic_tracer_run(ep, kap; nx=nx, dcx=dcx, dcy=dcy)

  yvar_theory = 2*kap*yvar.time[end] + dcy^2.0

  groupname = @sprintf("ep%02d", Int(ep*10))
  jldopen(filename, "a+") do file
    file["$groupname/yvar"] = yvar
    file["$groupname/c"] = prob.vars.c
    file["$groupname/x"] = prob.grid.X
    file["$groupname/y"] = prob.grid.Y
  end

  enhancement[i] = yvar.value / yvar_theory
  println("enhancement: ", enhancement[i], "\n\n")

end


jldopen(filename, "a+") do file
  file["enhnacement"] = enhancement
end
