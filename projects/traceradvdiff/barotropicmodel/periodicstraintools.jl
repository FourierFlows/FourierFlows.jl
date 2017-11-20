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
        δx = 0.05,
        δy = 0.05,
    nsaves = 100,
    nplots = 100,
      name = "test")

    dx = L/n              # grid spacing
    dy = H/n              # grid spacing
    tf = H*L/U            # final time
  umax = U / (L*(1-ε))    # maximum u-velocity
  vmax = ε*U/H            # maximum v-velocity

  Δt₁ = minimum([CFL*dx/umax, CFL*dy/vmax]) # first guess at timestep
  totsteps₁ = tf/Δt₁
  filename = @sprintf("%s_n%04d_ep%02d.jld2", name, n, 100ε)

  if isfile(filename); rm(filename); end
  if η == nothing; η = κ; end

  # Total steps, intermediate steps.
  intsteps = ceil(Int, totsteps₁/nsaves)    # Number of snapshots and diags
  totsteps = intsteps*nsaves
  plotfreq = floor(Int, totsteps/nplots)
  Δt = tf/totsteps

  # Sinsuidal topography and barotropic velocity
  u, v = wavychannelbarotropicflow(ε, H, L, U)
  #u, v = leewaveflow(U, ε, 2π/L, 2π/H; evanescent=false) 
  
  # Initial condition
  x₀, y₀ = 0.0, -H/2
  ci(x, y) = exp( -(x-x₀)^2/(2*δx^2) - (y-y₀)^2/(2*δy^2) ) / (2π*δx*δy)

  # Initialize the problem
  g = FourierFlows.TwoDGrid(n, L, n, H+ε; y0=-H-ε, effort=FFTW.MEASURE)
  prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(grid=g, κ=κ, η=η, u=u, v=v) 
  TracerAdvDiff.set_c!(prob, ci)

  # Output
  getc(prob) = irfft(prob.vars.sol, prob.grid.nx)
  output = Output(prob, filename, (:c, getc))

  # Diagnostics
  """ Return a vector of the strain tensor calculated at the center of
  mass of the tracer blob.  Note that since u = U/h, we have

    ux = - U hₓ / h^2        vx = y U ( hₓₓ / h^2 - 2*hₓ^2/h^3 )
    uy = 0                   vy = u hₓ / h^2  

  where hₓ = ∂h/∂x and hₓₓ = ∂²h/∂x².
  """
  h, hₓ, hₓₓ = sinusoidaltopo(ε, H, L)
  function straintensor(prob)
    prob.vars.ch .= prob.vars.sol
    A_mul_B!(prob.vars.c, g.irfftplan, prob.vars.ch)

    c = prob.vars.c
    g = prob.grid

    x = FourierFlows.cumulant_1x(c, g)
    y = FourierFlows.cumulant_1y(c, g)

    # Each strain component is a scalar
    ux = - U*hₓ(x) / h(x)^2
    vx = y*U * ( hₓₓ(x) / h(x)^2 - 2*hₓ(x)^2/h(x)^3 )
    vy = U*hₓ(x) / h(x)^2
    uy = 0.0

    [ux, vx, uy, vy]
  end

  #variance(prob) = FourierFlows.cumulant_2y(prob.vars.c, prob.grid)

  function variance(prob)
    prob.vars.ch .= prob.vars.sol
    A_mul_B!(prob.vars.c, g.irfftplan, prob.vars.ch)

    # Domain average
    C = g.dx*g.dy*sum(prob.vars.c)

    # Ybar
    @. prob.vars.cv = g.Y*prob.vars.c
    c1y = g.dx*g.dy*sum(prob.vars.cv) / C 

    # 2nd cumulant
    @. prob.vars.cv = (g.Y - c1y)^2.0 * prob.vars.c
    g.dx*g.dy*sum(prob.vars.cv) / C
  end

  variance(prob) = 1

  variancediag = Diagnostic(variance, prob; nsteps=totsteps)
    straindiag = Diagnostic(straintensor, prob; nsteps=totsteps)

  #diags = [variancediag, straindiag]
  #diags = [straindiag]
  diags = [variancediag]


  # Plotting
  mask = wavychannelmask(ε, H, L, prob.grid)
  

  # Integrate
  msg = "Wavy barotropic run with \n"
  msg *= @sprintf("   ε: %0.4f\n", ε)
  msg *= @sprintf("   n: %d\n", n)
  msg *= @sprintf("   κ: %.2e\n", κ)
  msg *= @sprintf("   η: %.2e\n", η)
  msg *= @sprintf("  δx: %.2e\n", δx)
  msg *= @sprintf("  δy: %.2e\n", δy)

  println(msg)

  startwalltime = time()
  println("Running with ε = $ε, n = $n...")
  while prob.step < totsteps
    @time stepforward!(prob, diags; nsteps=intsteps)
    #@time stepforward!(prob, nsteps=ntsteps)
    #saveoutput(output)

    #if prob.step % plotfreq == 0.0
    #  makeplot(prob, output[:c], mask, filename, save=true)
    #end

    #@printf("prog: %.2f, wall: %.2f min, c2: %.2f\n",
    #  prob.step/totsteps, (time()-startwalltime)/60, 
    #  variancediag.value/variancediag.data[1])
  end

  #savediagnostic(variancediag, "variance", output.filename)
  #savediagnostic(straindiag, "strain", output.filename)

  nothing
end




function makeplot(prob, c, mask, filename; save=false, show=false)
  close("all")
  fig, ax = subplots(figsize=(8, 4))
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

   h(x) = H*(1 - ε*sin(k₁*x))      # Topography
  hₓ(x) = -H*k₁*ε*cos(k₁*x)        # Topographic x-gradient

  u(x, y) = U/h(x)                 # Barotropic x-velocity
  v(x, y) = y*u(x, y)*hₓ(x)/h(x)   # Topography-following z-velocity

  u, v
end




function sinusoidaltopo(ε, H, L)
  k₁ = 2π/L

  # Sinusoidal topo and its derivatives
     h(x) = H*(1 - ε*sin(k₁*x))      
    hₓ(x) = -H*k₁*ε*cos(k₁*x)       
   hₓₓ(x) = H*k₁^2*ε*sin(k₁*x)     

  h, hₓ, hₓₓ 
end




"""
Returns an array of booleans that masks sinusoidal topography.
"""
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
