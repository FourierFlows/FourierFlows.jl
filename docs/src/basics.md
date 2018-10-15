# Code Basics

## Notation

The code solves differential equations of the form

$\partial_t u = \mathcal{L}u + \mathcal{N}(u)\ ,$

using Fourier transforms on periodic domains. The right side term $\mathcal{L}u$ is a 'linear' part of the equation.
The term $\mathcal{N}(u)$ is, in general, a 'nonlinear' part. In FourierFlows, $\mathcal{L}u$ is specified by 
physical modules in Fourier space as an array with the same size as $\hat u$, the Fourier transform of $u$.
The nonlinear term $\mathcal{N}u$ is specified by a function.

Boundary conditions in all spatial dimensions are periodic. That allows us to expand all variables using a Fourier 
decomposition. For example, a variable $\phi(x, t)$ that depends in one spatial dimension is expanded as:


```math
\phi(x, t) = \sum_{k} \widehat{\phi}(k, t)\,e^{\mathrm{i} k x}\ ,
```
where wavenumbers $k$ take the values $\tfrac{2\pi}{L_x}[0,\pm 1,\pm 2,\dots]$. The equation is time-stepped 
forward in Fourier space. That way $u$ becomes the array with all Fourier coefficients of the solution.

The coefficients for the linear operator $\mathcal{L}$ are stored in an array called `LC`. The term $\mathcal{N}(u)$ 
is computed for by calling the function `calcN!`.

## Basic steps for solving a problem: step through an example script

To illustrate the basic steps for solving a problem consider the 1D Kuramoto-Sivashinsky equation for $u(x, t)$:

$$\partial_t u + \partial_x^4 u + \partial_x^2 u + u\partial_x u = 0 \,$$

which in Fourier space reads:

$$\partial_t \widehat{u} = \underbrace{(- k^4 + k^2) \widehat{u}}_{\mathcal{L}\widehat{u}}
+ \underbrace{\widehat{ -u\partial_x u }}_{\mathcal{N}(\widehat{u})}\ .$$

A `FourierFlows.Problem` is composed of the following types:
- `Grid` (`OneDGrid` in this case)
- `Params` (empty in this case)
- `Vars`, which holds $u$, $\partial_x u$, $u\partial_x u$ and their Fourier transforms $\widehat{u}$, $\widehat{\partial_x u}$, $\widehat{u\partial_xu}$.
- `Equation`, which holds the linear coefficients `LC` and a function `calcN!` that computes $\mathcal{N}(\widehat{u})$.
- `TimeStepper` for stepping the solution forward,
- `State`, which holds the solution `sol` and current time `t`.


The example script found in  `examples/kuramotosivashinsky/trefethenexample.jl` demonstrates the above steps 
needed to construct a KuramotoSivashinsky `Problem`. For this we call 
`prob = Problem(nx=nx, Lx=Lx, dt=dt, stepper="ETDRK4")`. Looking into the  `Problem` 
function we can see the above steps:
```julia
function Problem(; nx=256, Lx=2π, dt=0.01, stepper="RK4")
    g  = OneDGrid(nx, Lx)
    pr = Params()
    vs = Vars(g)
    eq = Equation(pr, g)
    ts = TimeStepper(stepper, dt, eq.LC, g)
    FourierFlows.Problem(g, vs, pr, eq, ts)
end
```

`OneDGrid(nx, Lx)` builds a one-dimensional grid with a Fourier wavenumber array `kr`:
```julia
i1 = 0:Int(nx/2)
i2 = Int(-nx/2+1):-1
kr = Array{T}(2π/Lx*cat(1, i1))
```

For real-valued fields we use `rfft` and thus only positive wavenumbers are involved: array `kr`. 
Foe example, with `nx=8` and `Lx=2π` the wavenumber grids are: `k = [0, 1, 2, 3, 4, -3, -2, -1]` and 
`kr = [0, 1, 2, 3, 4]`.

The construction of the grids only works for an *even* number of grid points. Moreover, the Fourier transforms are 
most efficient when the number of grid points is the product of powers of 2 and 3. For example: 
$2^7=128$, $2^6 3^1=192$, or $2^8=256$.

`Vars(g)` initializes variables `u`, `ux`, and `uux` as real valued arrays of length `nx` and variables 
`uh`, `uxh`, and `uuxh` as complex valued arrays of length `nkr = Int(nx/2+1)` (the same length as `kr`). 
We use the convention that the Fourier transform of a variable is appended with an `h`, which stands for 'hat'.
For example, the transform of `phi` is `phih`.

The array `LC` is constructed by the `Equation` constructor
```julia
function Equation(p, g)
  LC = @. g.kr^2 - g.kr^4
  FourierFlows.Equation(LC, calcN!)
end
```

One of the fields of `Equation` is the function `calcN!`, 
which computes the nonlinear term $\mathcal{N}(\widehat{u})$, storing the result in `N`:
```julia
function calcN!(N, sol, t, s, v, p, g)
  @. v.uh = sol
  @. v.uxh = im*g.kr*sol
  ldiv!(v.u, g.irfftplan, v.uh) # irfft
  ldiv!(v.ux, g.irfftplan, v.uxh)
  @. v.uux = v.u*v.ux
  mul!(v.uuxh, g.rfftplan, v.uux)
  @. N = -v.uuxh
  dealias!(N, g)
  nothing
end
```

The time-stepper is constructed and stored as `ts`. Finally, all supertypes are gathered together as a 
`FourierFlows.Problem`.
