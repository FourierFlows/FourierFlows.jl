function test_withoutjld2()
  namewithjld2 = "blahblah.jld2"
  namewithoutjld2 = "blahblah"
  
  return FourierFlows.withoutjld2(namewithjld2) == namewithoutjld2 && FourierFlows.withoutjld2(namewithoutjld2) == namewithoutjld2
end

function test_uniquepath()
  path  = "foo.jld2"
  path1 = "foo_1.jld2"
  path2 = "foo_2.jld2"
  for filename in (path, path1, path2)
    if isfile(filename); rm(filename); end
  end
  touch(path)
  test1 = FourierFlows.uniquepath(path) == path1
  touch(path1)
  test2 = FourierFlows.uniquepath(path) == path2
  return test1 && test2
end

function test_outputconstructor(dev::Device=CPU())
  prob = Problem(nx=32, Lx=2π, kappa=1e-2, dt=1e-7, stepper="ForwardEuler", dev=dev)
  filename = joinpath(".", "testoutput.jld2")
  get_sol(prob) = prob.sol
  get_c(prob) = prob.vars.c
  
  out1 = Output(prob, filename, (:sol, get_sol))
  out2 = Output(prob, filename, (:sol, get_sol), (:c, get_c))
  
  return  typeof(out1)<:Output && typeof(out2)<:Output
end

function test_getindex(dev::Device=CPU())
  prob = Problem(nx=32, Lx=2π, kappa=1e-2, dt=1e-7, stepper="ForwardEuler", dev=dev)
  filename = joinpath(".", "testoutput.jld2")
  
  ctest = devzeros(dev, Float64, (prob.grid.nx, ))
  ctest[3] = π
  prob.vars.c .= ctest
  
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  
  return isapprox(ctest, getindex(out, :c), rtol=rtol_output)
end

function test_saveproblem_saveoutput(dev::Device=CPU())
  nx = 32
  prob = Problem(nx=nx, Lx=2π, kappa=1e-2*ones(nx), dt=1e-7, stepper="ForwardEuler", dev=dev)
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); rm(filename); end
  
  ctest = devzeros(dev, Float64, (prob.grid.nx, ))
  ctest[3] = π
  prob.vars.c .= ctest
  
  get_c(prob) = collect(prob.vars.c)
  
  out = Output(prob, filename, (:c, get_c))
  
  saveproblem(out)
  
  saveoutput(out)
  
  file = jldopen(filename)
  
  return isfile(filename) && isapprox(file["snapshots"]["c"]["0"], collect(ctest), rtol=rtol_output) && isapprox(file["grid"]["Lx"], prob.grid.Lx, rtol=rtol_output) && isapprox(file["eqn"]["L"], prob.eqn.L, rtol=rtol_output)
end

function test_saveproblemTwoDGrid(dev::Device=CPU())
       nx = 32
       Lx = 2π
    kappa = 1e-2
       dt = 1e-7
  stepper = "ForwardEuler"

     grid = TwoDGrid(dev, nx, Lx)
   params = FourierFlows.Diffusion.Params(dev, kappa)
     vars = FourierFlows.Diffusion.Vars(dev, grid)
     
      # manually construct an Equation for a 2D grid
      L = zeros(dev, Float64, (grid.nkr, grid.nl))
      @. L = -kappa * grid.kr^2
      eqn = FourierFlows.Equation(L, FourierFlows.Diffusion.calcN!, grid)

     prob = FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, dev)
  
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); rm(filename); end
    
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  
  saveproblem(out)
    
  file = jldopen(filename)
  
  return isfile(filename) && isapprox(file["grid"]["Ly"], prob.grid.Ly, rtol=rtol_output) && isapprox(file["eqn"]["L"], collect(prob.eqn.L), rtol=rtol_output)
end

function test_savediagnostic(dev::Device=CPU())
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); rm(filename); end

  prob = Problem(nx=6, Lx=2π, dev=dev)
  getone(prob) = 1
  nsteps=100
  freq=1
  ndata=ceil(Int, (nsteps+1)/freq)
  d = Diagnostic(getone, prob, nsteps=nsteps, freq=freq, ndata=ndata)
  stepforward!(prob, d, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))  
  saveproblem(out)
  savediagnostic(d, "mydiagnostic", filename)
  
  file = jldopen(filename)
  
  return isfile(filename) && file["diags"]["mydiagnostic"]["steps"]==expectedsteps
end