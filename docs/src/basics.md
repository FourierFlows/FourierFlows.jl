# Code Basics


## Basic Notation

The code solves partial differential equations of the general form:

$\partial_t u = \mathcal{L}u + \mathcal{N}(u)\ .$

(Note: ODEs are special cases of the above. Thus the code also solves ODEs.)

We decompose the right hand side of the above in a linear part ($\mathcal{L}u$) and a nonlinear part ($\mathcal{N}(u)$). The time steppers treat the linear and nonlinear parts differently.

Boundary conditions in all spatial dimensions are periodic. That allows us to expand all variables using a Fourier decomposition. For example, a variable $\phi(x, t)$ that depends in one spatial dimension is expanded as:

```math
\phi(x, t) = \sum_{k} \widehat{\phi}(k, t)\,e^{\mathrm{i} k x}\ ,
```
where wavenumbers $k$ take the values $\tfrac{2\pi}{L_x}[0,\pm 1,\pm 2,\dots]$. The equation is time-stepped forward in Fourier space. That way $u$ becomes the array with all Fourier coefficients of the solution.

The coefficients for the linear operator $\mathcal{L}$ are stored in an array called `LC`. The term $\mathcal{N}(u)$ is computed for by calling the function `calcN!`.





## Abstract SuperTypes

The code is divided along conceptual lines into problem-agnostic and problem-specific components. Files that contain problem-agnostic parts of the code are stored in `/src`. Files in `/src` define the domain, 'AbstractTypes' that supertype problem-specific types, and time-stepper types and routines. Problem-specific modules are stored in `/src/physics`.

Below is a list of all Abstract Supertypes used by the code:

- `AbstractGrid`: Includes all variables that have to do with the grid, both in physical space as well as in wavenumber space. Currently implemented are: `ZeroGrid` for ODEs, `OneGrid` for PDEs with one spatial dimension, and `TwoGrid` for PDEs with two spatial dimensions. Grids are generic and work for any problem of that particular dimension.

- `AbstractParams`: Includes all parameters or functions related with the problem do not vary throughout the integration.

- `AbstractVars`: Includes all variables of the problem that change along the integration.

- `AbstractEquation`: Includes the array with the coefficients of the linear part of the equation, `LC` as well as function `calcN!` that computes the nonlinear part of the equation.

- `AbstractState`: Includes the solution `sol` at current time-step as well as the time-step `dt`, the time `t`, and `step` which counts the number of time-steps taken.

- `AbstractTimeStepper`: Includes all details for the time-stepper (e.g., `dt`, various coefficients, `sol` at intermediate time-step values). Time-steppers are generic and work for any problem.

- `AbstractProblem`: A super-supertype that includes all of the above. That is `problem` includes `grid`, `vars`, `params`, `eqn`, `ts`, `state`, and also `t` and `step`.

Grids and time-steppers are generic and work for any problem of that particular dimension. State and Problem just gathers things together. Thus, to write a solver for a new physical problem you only need to prescribe `params`, `vars`, the coefficients of the linear part, `LC`, and function `calcN!`.


## Source code organization

The code is divided along conceptual lines into problem-agnostic and problem-specific components. Files that contain problem-agnostic parts of the code are stored in `/src`. Files in `/src` define the domain, 'AbstractTypes' that supertype problem-specific types, and time-stepper types and routines. Problem-specific modules are stores in `/src/physics`.

Here's an overview of the code structure:

- `/src/`
    - `FourierFlows.jl`
        - Defines supertyping AbstractParams, AbstractGrid, etc.
        - Defines a `Problem` type to organize the grid, vars, params, equation, and timestepper into a single structure.
        - Includes all sources files and physics files.
   - `timesteppers.jl`: defines modules and `stepforward!` routines for
        various time-steppers. Current implemented time-steppers are:
        - Forward Euler
        - 3rd-order Adams-Bashforth (AB3)
        - 4th-order Runge-Kutta (RK4)
        - 4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)
        - 4th-order Dual Runge-Kutta (DualRK4)
        - 4th-order Dual Runge-Kutta Exponential Time Differencing (DualETDRK4)

        For each time-stepper exists also a "filtered" version that filters out high-wavenumber spectral components of the solution. The `Dual` time-steppers evolve a state variable that comprises both of real valued         and complex valued fields.

    - `physics/`
        - `twodturb.jl`: Defines a `TwoDTurb` module that provides a solver for the two-dimensional vorticity equation.
        - `barotropicqg.jl`: Defines a `BarotropicQG` module that provides several solvers for the barotropic QG model that permit beta, topography, beta + topography, and forcing.
        - `kuramotosivashinsky.jl`: Defines a `KuramotoSivashinsky` module that solves the Kuramoto-Sivashinsky.
        - `verticallyfourierboussinesq.jl`: Defines a `VerticallyFourierBoussinesq` module that solves the two-mode truncation of the Fourier series thin-layer approximation to the hydrostatic Boussinesq equations.
        - `verticallycosinerboussinesq.jl`: Defines a `VerticallyCosineBoussinesq` module that solves the two-mode truncation of the Sin/Cos series thin-layer approximation to the hydrostatic Boussinesq equations.
        - `traceradvdiff.jl`: Defines a `TracerAdvDiff` module that provides a solver for a two-dimensional and periodic tracer field in a given 2D flow (u, w), which can be an arbitrary function of x, z, and t.



## Basic steps for solving a problem: step through an example script

To illustrate the basic steps for solving a problem consider the 1D Kuramoto-Sivashinsky equation for $u(x, t)$:

$$\partial_t u + \partial_x^4 u + \partial_x^2 u + u\partial_x u = 0\ ,$$

which in Fourier base reads:

$$\partial_t \widehat{u} = \underbrace{(- k_x^4 + k_x^2) \widehat{u}}_{\mathcal{L}\widehat{u}}
+ \underbrace{\widehat{ -u\partial_x u }}_{\mathcal{N}(\widehat{u})}\ .$$


The steps to construct an `AbstractProblem` for the above are:
- Construct an `AbstractGrid`; for this problem we use the `OneGrid`.
- Construct an `AbstractParams`; for this problem `params` is be empty as there are no parameters in the equation. (Note that e.g., the domain size `Lx` and the number of gridpoints `nx` belong to the grid.)
- Construct an `AbstractVars`; for this problem `vars` includes $u$, $\partial_x u$, $u\partial_x u$ and their Fourier transforms $\widehat{u}$, $\widehat{\partial_x u}$, $\widehat{u\partial_xu}$.
- Construct the equations by prescribing coefficients for the linear part as an array `LC` and a function `calcN!` that computes $\mathcal{N}(\widehat{u})$.
- Construct the time-stepper which includes function `stepforward!` that time-steps the solution.
- Construct the `state` and gather everything as an `AbstractProblem`.


The example script found in  `examples/kuramotosivashinsky/trefethenexample.jl` demonstrates the above steps needed to construct an `AbstractProblem`. The `prob` is constructed by calling `prob = InitialValueProblem(nx=nx, Lx=Lx, dt=dt, stepper="ETDRK4")`. Looking into the  `InitialValueProblem` function we can see the above steps:
```julia
function InitialValueProblem(;
     nx = 256,
     Lx = 2π,
     dt = 0.01,
stepper = "RK4"
)

g  = OneDGrid(nx, Lx)
pr = Params()
vs = Vars(g)
eq = Equation(pr, g)
ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

FourierFlows.Problem(g, vs, pr, eq, ts)
end
```

The `OneDGrid` function is called for the grid. Within grid the wavenumber array is constructed:
```julia
i1 = 0:Int(nx/2)
i2 = Int(-nx/2+1):-1
k = Array{T}(2π/Lx*cat(1, i1, i2))
kr = Array{T}(2π/Lx*cat(1, i1))
```

For real-valued fields we use `rfft` and thus only positive wavenumbers are involved: array `kr`. E.g., for `nx=8` and `Lx=2π` the wavenumber grids are: `k = [0, 1, 2, 3, 4, -3, -2, -1]` and `kr = [0, 1, 2, 3, 4]`.

The construction of the grids only works for *even* number of grid points. Moreover, since the code relies on the $\mathrm{FFT}$ algorithm, we suggest you use a power of 2 as the number of grid points, since then $\mathrm{FFT}$ is most efficient. 


Function `Vars(g)` initialize variables `u`, `ux`, and `uux` as real valued arrays of length `nx` and variables `uh`, `uxh`, and `uuxh` as complex valued arrays of length `nkr = Int(nx/2+1)` (the same length as `kr`). As a general convention variable names with `h` denote the Fourier transforms of the corresponding variable (`h` stands for 'hat').

The array `LC` is constructed by `Equation` function
```julia
function Equation(p, g)
  LC = @. g.kr^2 - g.kr^4
  FourierFlows.Equation(LC, calcN!)
end
```
Also `eq` includes function `calcN!` which computes the nonlinear term $\mathcal{N}(\widehat{u})$:
```julia
function calcN!(N, sol, t, s, v, p, g)
  @. v.uh = sol
  @. v.uxh = im*g.kr*sol
  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.ux, g.irfftplan, v.uxh)
  @. v.uux = v.u*v.ux
  A_mul_B!(v.uuxh, g.rfftplan, v.uux)
  @. N = -v.uuxh
  dealias!(N, g)
  nothing
end
```

The time-stepper is constructed and stored as `ts`. Finally, all supertypes are gathered together as an `AbstractProblem`.



## Tutorials

```@contents
Pages = [
    "modules/kuramotosivashinsky.md",
    "modules/twodturb.md",
    "modules/barotropicqg.md",
    "modules/traceradvdiff.md"
        ]
Depth = 1
```
