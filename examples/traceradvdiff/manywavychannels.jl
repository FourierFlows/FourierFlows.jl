include("./wavychannel.jl")

using JLD2

nx  = 512 
kap = 1e-4
eta = 1e-6
dcx = 0.01
dcy = 0.02


epz = 0.0:0.1:0.9
nepz = length(epz)

diffratio = zeros(nepz)

for i = 1:nepz

  ep = epz[i]

  prob, yvar = wavy_barotropic_tracer_run(ep, kap; nx=nx, dcx=dcx, dcy=dcy, 
    CFL=0.1, eta=eta)

  yvar_theory = 2*kap*yvar.time[end] + dcy^2.0

  diffratio[i] = yvar.value / yvar_theory
  println("diffusivity ratio: ", diffratio[i], "\n\n")

  time, variance = yvar.time, yvar.data
  c, X, Y = prob.vars.c, prob.grid.X, prob.grid.Y

  filename = @sprintf("wavychannel_nx%04d_ep%02d.jld2", nx, Int(ep*10))
  @save filename time variance c X Y

end

finalfilename = @sprintf("manyeps_nx%04d.jld2", nx)
@save finalfilename epz diffratio
