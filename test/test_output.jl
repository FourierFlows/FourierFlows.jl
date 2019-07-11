function test_withoutjld2()
  namewithjld2 = "blahblah.jld2"
  namewithoutjld2 = "blahblah"
  
  return FourierFlows.withoutjld2(namewithjld2) == namewithoutjld2 && FourierFlows.withoutjld2(namewithoutjld2) == namewithoutjld2
end

function test_uniquepath()
  path  = "foo.jld2"
  path1 = "foo_1.jld2"
  path2 = "foo_2.jld2"
  touch(path)
  test1 = FourierFlows.uniquepath(path) == path1
  touch(path1)
  test2 = FourierFlows.uniquepath(path) == path2
  return test1 && test2
end

function test_outputconstructor()
  prob = FourierFlows.Diffusion.Problem(nx=32, Lx=2π, kappa=1e-2, dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  get_sol(prob) = prob.sol
  get_c(prob) = prob.vars.c
  
  out1 = Output(prob, filename, (:sol, get_sol))
  out2 = Output(prob, filename, (:sol, get_sol), (:c, get_c))
  
  return  typeof(out1)<:Output && typeof(out2)<:Output
end

function test_getindex()
  prob = FourierFlows.Diffusion.Problem(nx=32, Lx=2π, kappa=1e-2, dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  
  ctest = zeros((prob.grid.nx, ))
  ctest[3] = π
  prob.vars.c .= ctest
  
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  
  return isapprox(ctest, getindex(out, :c), rtol=rtol_output)
end