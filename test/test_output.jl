function test_withoutjld2()
  namewithjld2 = "blahblah.jld2"
  namewithoutjld2 = "blahblah"
  
  FourierFlows.withoutjld2(namewithjld2) == namewithoutjld2 && FourierFlows.withoutjld2(namewithoutjld2) == namewithoutjld2
end

function test_outputconstructor()
  prob = FourierFlows.Diffusion.Problem(nx=32, Lx=2π, kappa=1e-2, dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  get_sol(prob) = prob.sol
  get_c(prob) = prob.vars.c
  
  out1 = Output(prob, filename, (:sol, get_sol))
  out2 = Output(prob, filename, (:sol, get_sol), (:c, get_c))
  
  typeof(out1)<:Output && typeof(out2)<:Output
end