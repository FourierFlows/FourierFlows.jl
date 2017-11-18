__precompile__()

include("../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, NullableArrays
import FourierFlows
import FourierFlows.TracerAdvDiff

@pyimport warnings
warnings.filterwarnings("ignore")

@pyimport numpy.ma as ma
PyObject(a::NullableArray) = pycall(ma.array, Any, a.values, mask=a.isnull)




function wavybarotropicrun(ε;
         n = 512,
         L = 2π,
         H = 1.0,
         κ = 1e-4,
         η = nothing,
         U = 1.0,
       CFL = 0.1,
        δx = 0.02,
        δy = 0.02,
  savefreq = 100,
  plotfreq = 10,
      name = "test")


   δ = L/n
  um = U / (L*(1-ε))
  tf = H*L/U

  Δt₁ = CFL*δ/um
  totsteps₁ = tf/Δt₁
  filename = @sprintf("%s_n%04d_ep%02d.jld2", name, n, 100ε)

  if isfile(filename); rm(filename); end
  if η == nothing; η = κ; end

  # Total steps, intermediate steps.
  intsteps = ceil(Int, totsteps₁/savefreq) # Number of snapshots and diags
  totsteps = intsteps*savefreq
  Δt = tf/totsteps

  # Sinsuidal topography and barotropic velocity
  u, v = wavychannelbarotropicflow(ε, H, L, U)
  #u, v = leewaveflow(U, ε, 2π/L, 2π/H; evanescent=false) 
  
  # Initial condition
  x₀, y₀ = 0.0, -H/2
  ci(x, y) = exp( -(x-x₀)^2/(2*δx^2) - (y-y₀)^2/(2*δy^2) ) / (2π*δx*δy)

  # Initialize the problem
  g = FourierFlows.TwoDGrid(n, L, n, H+ε; y0=-H-ε, effort=FFTW.PATIENT)
  prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(grid=g, κ=κ, η=η, u=u, v=v) 
  TracerAdvDiff.set_c!(prob, ci)

  # Output
  getc(prob) = irfft(prob.vars.sol, prob.grid.nx)
  output = Output(prob, filename, (:c, getc))

  # Diagnostics
  variance(prob) = FourierFlows.cumulant_2y(prob.vars.c, prob.grid)
  variancediag = Diagnostic(variance, prob)
  diags = [variancediag]


  # Plotting
  fig, ax = subplots(figsize=(8, 4))
  mask = wavychannelmask(ε, H, L, prob.grid)
  
  # Integrate
  startwalltime = time()
  println("Running with ε = $ε, n = $n...")
  while prob.step < totsteps

    stepforward!(prob, diags; nsteps=intsteps)
    saveoutput(output)

    if prob.step % plotfreq == 0.0
      makeplot(ax, prob, output[:c], mask, filename, save=true)
    end

    @printf("frac step: %.2f, walltime = %.2f\n", 
      prob.step/totsteps, time()-startwalltime)

    
  end

  nothing
end




function makeplot(ax, prob, c, mask, filename; save=false, show=false)
  axes(ax)
  pcolormesh(prob.grid.X, prob.grid.Y, PyObject(NullableArray(c, mask)))

  if save
    savename = @sprintf("./plots/%s_%06d.png", filename[1:end-5], prob.step)
    savefig(savename, dpi=240)
  end

  if show
    pause(0.1)
  end

  nothing
end




"""
Returns u, v functions that define barotropic flow in a wavy channel with 
wave height ε, channel height H, and length L.
"""
function wavychannelbarotropicflow(ε, H, L, U) 
  k₁ = 2π/L      

   h(x) = H*(1 - ε*sin(k₁*x))       # Topography
  ∂h(x) = -H*k₁*ε*cos(k₁*x)        # Topographic x-gradient

  u(x, y) = U/h(x)                 # Barotropic x-velocity
  v(x, y) = y*u(x, y)*∂h(x)/h(x)   # Topography-following z-velocity

  u, v
end




function wavychannelmask(ε, H, L, grid)

  k₁ = 2π/L      
  h(x) = H*(1 - ε*sin(k₁*x))       # Topography

  mask = Array{Bool}(grid.nx, grid.ny)
  for j=1:grid.ny, i=1:grid.nx;
    mask[i, j] = grid.y[j] < -h(grid.x[i]) ? true : false 
  end

  mask
end




""" 
Returns u, v for a lee wave forming in steady flow with velocity U, 
wavenumbers k, m, and amplitude h, so that

  psi = Uz - U*h*cos(kx+mz); u = psi_y, v = -psi_x

"""
function leewaveflow(U, h₀, k, m; evanesent=false)
  if evanescent
    u(x, y) = U + m*U*h₀*cos(k*x)*exp(-m*z)
    v(x, y) = -k*U*h₀*sin(k*x)*exp(-m*z)
  else
    u(x, y) = U + m*U*h₀*sin(k*x + m*y)
    v(x, y) = k*U*h₀*sin(k*x + m*y)
  end

  u, v
end
