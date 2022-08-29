function test_withoutjld2()
  namewithjld2 = "blahblah.jld2"
  namewithoutjld2 = "blahblah"
  
  return FourierFlows.withoutjld2(namewithjld2) == namewithoutjld2 &&
         FourierFlows.withoutjld2(namewithoutjld2) == namewithoutjld2
end

function test_uniquepath()
  path  = "foo.jld2"
  path₁ = "foo_1.jld2"
  path₂ = "foo_2.jld2"
  for filename in (path, path₁, path₂)
    if isfile(filename); rm(filename); end
  end
  touch(path)
  test₁ = FourierFlows.uniquepath(path) == path₁
  touch(path₁)
  test₂ = FourierFlows.uniquepath(path) == path₂
  return test₁ && test₂
end

function test_outputconstructor(dev::Device=CPU())
  prob = Problem(dev; nx=32, Lx=2π, κ=1e-2, dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  get_sol(prob) = prob.sol
  get_c(prob) = prob.vars.c
  
  out1 = Output(prob, filename, (:sol, get_sol))
  out2 = Output(prob, filename, (:sol, get_sol), (:c, get_c))
  
  return  typeof(out1)<:Output && typeof(out2)<:Output
end

function test_getindex(dev::Device=CPU())
  prob = Problem(dev; nx=32, Lx=2π, κ=1e-2, dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  
  ctest = devzeros(dev, Float64, (prob.grid.nx, ))
  CUDA.@allowscalar ctest[3] = π
  @. prob.vars.c = ctest
  
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  
  return isapprox(ctest, getindex(out, :c), rtol=rtol_output)
end

function test_saveproblem_saveoutput(dev::Device=CPU())
  nx = 32
  prob = Problem(dev; nx, Lx=2π, κ=1e-2*ones(nx), dt=1e-7, stepper="ForwardEuler")
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); rm(filename); end
  
  ctest = devzeros(dev, Float64, (prob.grid.nx, ))
  CUDA.@allowscalar ctest[3] = π
  prob.vars.c .= ctest
  
  get_c(prob) = collect(prob.vars.c)
  
  out = Output(prob, filename, (:c, get_c))
  
  saveproblem(out)
  
  saveoutput(out)
  
  file = jldopen(filename)
  
  return isfile(filename) &&
         isapprox(file["snapshots"]["c"]["0"], collect(ctest), rtol=rtol_output) &&
         isapprox(file["grid"]["Lx"], prob.grid.Lx, rtol=rtol_output) &&
         isapprox(file["eqn"]["L"], prob.eqn.L, rtol=rtol_output)
end

function test_saveproblemTwoDGrid(dev::Device=CPU())
       nx = 32
       Lx = 2π
        κ = 1e-2
       dt = 1e-7
  stepper = "ForwardEuler"

     grid = TwoDGrid(dev; nx, Lx)
   params = FourierFlows.Diffusion.Params(dev, κ)
     vars = FourierFlows.Diffusion.Vars(dev, grid)
     
  # manually construct an Equation for a 2D grid
     L = zeros(dev, Float64, (grid.nkr, grid.nl))
  @. L = - κ * grid.kr^2
  equation = FourierFlows.Equation(L, FourierFlows.Diffusion.calcN!, grid)

  prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
  
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); GC.gc(); rm(filename); end
    
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  
  saveproblem(out)
    
  file = jldopen(filename)
  
  return isfile(filename) &&
         isapprox(file["grid"]["Ly"], prob.grid.Ly, rtol=rtol_output) &&
         isapprox(file["eqn"]["L"], collect(prob.eqn.L), rtol=rtol_output)
end

function test_savediagnostic(dev::Device=CPU())
  filename = joinpath(".", "testoutput.jld2")
  if isfile(filename); GC.gc(); rm(filename); end
  
  prob = Problem(dev; nx=6, Lx=2π)
  getone(prob) = 1
  nsteps = 100
  freq = 2
  ndata = ceil(Int, (nsteps+1)/freq)
  d = Diagnostic(getone, prob; nsteps, freq, ndata)
  
  nsteps += 5 #to check what happens if we time-step a bit longer than the size of the diagnostic
  stepforward!(prob, d, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  
  get_c(prob) = prob.vars.c
  out = Output(prob, filename, (:c, get_c))
  saveproblem(out)
  savediagnostic(d, "mydiagnostic", filename)
  
  file = jldopen(filename)
  
  return isfile(filename) &&
         file["diagnostics"]["mydiagnostic"]["steps"] == expectedsteps &&
         file["diagnostics"]["mydiagnostic"]["data"] == d.data[1:d.i]
end

struct TestbedParams{T, Trfft} <: AbstractParams
  parameter :: T
       func :: Function
       plan :: Trfft
end

function test_savefields(dev::Device=CPU(); parameter=2.2)  
  func(x) = sin(x^2)
  
  grid = TwoDGrid(dev; nx=6, Lx=2π)
  plan = grid.rfftplan
  
  params = TestbedParams(parameter, func, plan)
  
  filename = joinpath(".", "testsavefields.jld2")
  if isfile(filename); GC.gc(); rm(filename); end
  
  file = jldopen(filename, "a+")
  FourierFlows.savefields(file, params)
  
  return file
end
