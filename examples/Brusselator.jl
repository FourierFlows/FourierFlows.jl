# # Linear rotating shallow water dynamics
#
#md # This example can be run online via [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/OneDShallowWaterGeostrophicAdjustment.ipynb). 
#md # Also, it can be viewed as a Jupyter notebook via [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/OneDShallowWaterGeostrophicAdjustment.ipynb).
# 
#
# This example solves the linear 1D rotating shallow water equations
# for the $u(x, t)$, $v(x, t)$ and the surface surface elevation $\eta(x, t)$, 
# for a fluid with constant rest-depth $H$. That is, the total fluid's depth 
# is $H+\eta(x,t)$ with $|\eta| \ll H$.
# 
# The linearized equations for the evolution of $u$, $v$, $\eta$ are:
#
# ```math
# \frac{\partial u}{\partial t} - f v = - g \frac{\partial \eta}{\partial x} - \mathrm{D} u, \\
# \frac{\partial v}{\partial t} + f u = - \mathrm{D} v, \\
# \frac{\partial \eta}{\partial t} + H \frac{\partial u}{\partial x} = - \mathrm{D} \eta.
# ```
#
# Above, $g$ is the gravitational acceleration, $f$ is the  Coriolis parameter,
# and $\mathrm{D}$ indicates a hyperviscous linear operator of the
# form $(-1)^{n_\nu} \nu \nabla^{2 n_\nu}$, with $\nu$ the viscosity
# coefficient and $n_\nu$ the order of the operator.
#
# Rotation introduces the deformation length scale, $L_d = \sqrt{g H}/f$. 
# Disturbances with length scales much smaller than $L_d$ don't "feel" 
# the rotation and propagate as inertia-gravity waves. Disturbances with 
# length scales comparable or larger than $L_d$ should be approximately 
# in geostrophic balance, i.e., the Coriolis acceleration
# $f\widehat{\boldsymbol{z}} \times \boldsymbol{u}$ should be in approximate 
# balance with the pressure gradient $-g \boldsymbol{\nabla} \eta$.


using FourierFlows, Plots
using FFTW: rfft, irfft
using LinearAlgebra: mul!, ldiv!
using Printf
using Random

# ## Coding up the equations
# ### A demonstration of FourierFlows.jl framework
#
# What follows is a step-by-step tutorial showin how you can create your own
# solver for an equation of your liking.

# The basic building blocks for a `FourierFlows.Problem()` are:
# - `Grid` struct containining the physical and wavenumber grid for the problem,
# - `Params` struct containining all the parameters of the problem,
# - `Vars` struct containining arrays with the variables used in the problem,
# - `Equation` struct containining the coefficients of the linear operator $\mathcal{L}$ and the function that computes the nonlinear terms, usually named `calcN!()`.
# 
# The `Grid` structure is provided by FourierFlows.jl. One simply has to call one of
# the `OneDGrid()`,  `TwoDGrid()`, or `ThreeDGrid()` grid constructors, depending
# on the dimensionality of the problem. All other structs mentioned above are problem-specific
# and need to be constructed for every set of equations we want to solve.

# First lets construct the `Params` struct that contains all parameters of the problem 

struct Params{T} <: AbstractParams
   B :: T         # B parameter
   D :: T         # D parameter
   E :: T         # E parameter
end
nothing #hide

# Now the `Vars` struct that contains all variables used in this problem. For this
# problem `Vars` includes the represenations of the flow fields in physical space 
# `u` and `v` and their Fourier transforms `uh` and `vh`.

struct Vars{Aphys, Atrans} <: AbstractVars
   u :: Aphys
   v :: Aphys
  uh :: Atrans
  vh :: Atrans
end
nothing #hide

# A constructor populates empty arrays based on the dimension of the `grid`
# and then creates `Vars` struct.
"""
    Vars!(dev, grid)
Constructs Vars based on the dimensions of arrays of the `grid`.
"""
function Vars(::Dev, grid::AbstractGrid) where Dev
  T = eltype(grid)
  @devzeros Dev T grid.nx u v
  @devzeros Dev Complex{T} grid.nkr uh vh
  return Vars(u, v, uh, vh)
end
nothing #hide


# In Fourier space, the 1D linear shallow water dynamics are:
#
# ```math
# \frac{\partial \hat{u}}{\partial t} = \underbrace{ f \hat{v} - \mathrm{i} k g \hat{\eta} }_{\mathcal{N}_u} \; \underbrace{- \nu |\boldsymbol{k}|^2 }_{\mathcal{L}_u} \hat{u}, \\
# \frac{\partial \hat{v}}{\partial t} = \underbrace{ - f \hat{u} }_{\mathcal{N}_v} \; \underbrace{- \nu |\boldsymbol{k}|^2 }_{\mathcal{L}_v} \hat{v}, \\
# \frac{\partial \hat{\eta}}{\partial t} = \underbrace{ - \mathrm{i} k H \hat{u} }_{\mathcal{N}_{\eta}} \; \underbrace{- \nu |\boldsymbol{k}|^2 }_{\mathcal{L}_{\eta}} \hat{\eta}.
# ```
# Although, e.g., terms involving the Coriolis accelaration are, in principle, 
# linear we include them in the nonlinear term $\mathcal{N}$ because they render 
# the linear operator $\mathcal{L}$ non-diagonal.
#
# With these in mind, we construct function `calcN!` that computes the nonlinear terms.
#
"""
    calcN!(N, sol, t, clock, vars, params, grid)
The function that computes the nonlinear terms for our problem.
"""
function calcN!(N, sol, t, clock, vars, params, grid)
  @. vars.uh = sol[:, 1]
  @. vars.vh = sol[:, 2]
  ldiv!(vars.u, grid.rfftplan, deepcopy(vars.uh))
  ldiv!(vars.v, grid.rfftplan, vars.vh)
  @. vars.v *= vars.u^2                     # v->v*u²
  mul!(vars.vh, grid.rfftplan, vars.v)  # vh -> FFT(v*u²)
  rhsu = @. + vars.vh + params.E             # + v*u² + E
  rhsv = @. - vars.vh + params.B * vars.uh   # - v*u² + B*u
  N[:, 1] .= rhsu
  N[:, 2] .= rhsv
  return nothing
end
nothing #hide
 
# Next we construct the `Equation` struct:

"""
    Equation!(prob)
Construct the equation: the linear part, in this case the hyperviscous dissipation,
and the nonlinear part, which is computed by `caclN!` function.
"""
function Equation(dev, params, grid::AbstractGrid)
  T = eltype(grid)
  L = zeros(dev, T, (grid.nkr, 2))
  diffusion = @. - grid.kr^2
  @. L[:, 1] = params.D * diffusion  - (params.B + 1) # D ∂²u/∂x² - (B+1) u
  @. L[:, 2] = diffusion # ∂²v/∂x²
  return FourierFlows.Equation(L, calcN!, grid)
end
nothing #hide

# We now have all necessary building blocks to construct a `FourierFlows.Problem`. 
# It would be useful, however, to define some more "helper functions". For example,
# a function that updates all variables given the solution `sol` which comprises $\hat{u}$,
# $\hat{v}$ and $\hat{\eta}$:

"""
    updatevars!(prob)
Update the variables in `prob.vars` using the solution in `prob.sol`.
"""
function updatevars!(prob)
  vars, grid, sol = prob.vars, prob.grid, prob.sol
  
  @. vars.uh = sol[:, 1]
  @. vars.vh = sol[:, 2]
  
  ldiv!(vars.u, grid.rfftplan, deepcopy(sol[:, 1])) # use deepcopy() because irfft destroys its input
  ldiv!(vars.v, grid.rfftplan, deepcopy(sol[:, 2])) # use deepcopy() because irfft destroys its input
  return nothing
end
nothing #hide

# Another useful function is one that prescribes an initial condition to the state variable `sol`.

"""
    set_uv!(prob, u0, v0)
Sets the state variable `prob.sol` as the Fourier transforms of `u0` and `v0`
and update all variables in `prob.vars`.
"""
function set_uv!(prob, u0, v0)
  vars, grid, sol = prob.vars, prob.grid, prob.sol
  
  A = typeof(vars.u) # determine the type of vars.u
  
  mul!(vars.uh, grid.rfftplan, A(u0)) # A(u0) converts u0 to the same type as vars expects (useful if u0 is a CPU array while working on the GPU)
  mul!(vars.vh, grid.rfftplan, A(v0)) # A(v0) converts u0 to the same type as vars expects (useful if v0 is a CPU array while working on the GPU)

  @. sol[:, 1] = vars.uh
  @. sol[:, 2] = vars.vh
    
  updatevars!(prob)
  return nothing
end
nothing #hide

# ## Let's prescibe parameter values and solve the PDE
#
# We are now ready to write up a program that sets up parameter values, constructs 
# the problem `prob`, # time steps the solutions `prob.sol` and plots it.

# ## Choosing a device: CPU or GPU

dev = CPU()    # Device (CPU/GPU)
nothing # hide

# ## Numerical parameters and time-stepping parameters

     nx = 512            # grid resolution
stepper = "RK4"          # timestepper
     dt = 0.01           # timestep
 nsteps = 2000           # total number of time-steps
nothing # hide


# ## Physical parameters

E = 1.4

D_c = ( (sqrt(1+E^2)-1) / E )^2
B_H = (1 + E*sqrt(D))^2 
   
kc = sqrt( E / (sqrt(1+E^2) - 1) )

L = 137.37     # Domain length
 
ϵ = 0.1
μ = 25
ρ = 0.178

B = B_H + ϵ^2 * μ
D = D_c + ϵ^2 * ρ

nothing # hide


# ## Construct the `struct`s and you are ready to go!
# Create a `grid` and also `params`, `vars`, and the `equation` structs. Then 
# give them all as input to the `FourierFlows.Problem()` constructor to get a
# problem struct, `prob`, that contains all of the above.

grid = OneDGrid(dev, nx, L)
params = Params(B, D, E)
vars = Vars(dev, grid)
equation = Equation(dev, params, grid)

prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
nothing #hide

get_u(prob) = irfft(prob.sol[:, 1], prob.grid.nx)
get_v(prob) = irfft(prob.sol[:, 2], prob.grid.nx)

u = Diagnostic(get_u, prob; nsteps=nsteps)
v = Diagnostic(get_v, prob; nsteps=nsteps)
diags = [u]
# ## Setting initial conditions

# For initial condition we take the fluid at rest ($u=v=0$). The free surface elevation
# is perturbed from its rest position ($\eta=0$); the disturbance we impose a Gaussian 
# bump with half-width greater than the deformation radius and on top of that we 
# superimpose some random noise with scales smaller than the deformation radius. 
# We mask the small-scale perturbations so that it only applies in the central part 
# of the domain by applying
#
# The system develops geostrophically-balanced jets around the Gaussian bump, 
# while the smaller-scale noise propagates away as inertia-gravity waves. 

# First let's construct the Gaussian bump.
# Next the noisy perturbation. The `mask` is simply a product of $\tanh$ functions.

l = 30;
θ = @. (grid.x > -l/2) & (grid.x < l/2)


au, av =  (E^2 + kc^2) / B, 1
cu, cv = -E * (E + im) / B, 1

u0 = @.  E  + ϵ * real( au * exp(im * kc * grid.x) * θ + cu * (1-θ) )
v0 = @. B/E + ϵ * real( av * exp(im * kc * grid.x) * θ + cv * (1-θ) )

plot_u0 = plot(grid.x, u0,
               color = :black,
              legend = :false,
           linewidth = [3 2],
               alpha = 0.7,
               xlims = (-L/2, L/2),
               # ylims = (-0.3, 0.3),
              xlabel = "x",
              ylabel = "u(x, 0)")

plot_v0 = plot(grid.x, v0,
               color = :black,
              legend = :false,
           linewidth = [3 2],
               alpha = 0.7,
               xlims = (-L/2, L/2),
               # ylims = (-0.3, 0.3),
              xlabel = "x",
              ylabel = "v(x, 0)")


title = plot(title = "initial conditions", grid = false, showaxis = false, bottom_margin = -20Plots.px)

plot(title, plot_u0, plot_v0,
           layout = @layout([A{0.01h}; [B; C]]),
             size = (600, 400))

# Sum the Gaussian bump and the noise and then call `set_uvη!()` to set the initial condition to the problem `prob`.

set_uv!(prob, u0, v0)

# ## Visualizing the simulation

# We define a function that plots the surface elevation $\eta$ and the 
# depth-integrated velocities $u$ and $v$. 

function plot_output(prob)
  plot_u = plot(grid.x, vars.u,
                 color = :blue,
                legend = false,
             linewidth = 2,
                 alpha = 0.7,
                 xlims = (-L/2, L/2),
                 # ylims = (-0.3, 0.3),
                xlabel = "x",
                ylabel = "u")

  plot_v = plot(grid.x, vars.v,
                 color = :red,
                legend = false,
             linewidth = 2,
                 alpha = 0.7,
                 xlims = (-L/2, L/2),
                 # ylims = (-0.3, 0.3),
                xlabel = "x",
                ylabel = "v")

  
  title = plot(title = "Brusselator", grid = false, showaxis = false, bottom_margin = -20Plots.px)
  
  return plot(title, plot_u, plot_v, 
           layout = @layout([A{0.01h}; [B; C]]),
             size = (600, 400))
end
nothing # hide


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time. We update variables by calling 
# `updatevars!()` and we also update the plot. We enclose the `for` loop in 
# an `@animate` macro to produce an animation of the solution.





 
 
stepforward!(prob, diags, nsteps)
U = zeros(grid.nx, nsteps+1)
[ U[:, j] = u.data[j] for j=1:nsteps+1 ]

heatmap(grid.x, u.t, U', c = :grayC)

asdfasdf
p = plot_output(prob)

anim = @animate for j=0:nsteps
  updatevars!(prob)
    
  p[2][1][:y] = vars.u    # updates the plot for u
  p[2][:title] = "t = "*@sprintf("%.1f", prob.clock.t / 60 )*" min" # updates time in the title
  p[3][1][:y] = vars.v    # updates the plot for v

  stepforward!(prob)
end

mp4(anim, "brusselator.mp4", fps=18)
