# GPU

FourierFlows allows you to easily construct and run problems on GPUs.

Upon calling

```julia
using FourierFlows
```

FourierFlows will check whether any CUDA enabled device is present. If such a device is 
found then FourierFlows makes sure that CUDA related packages are loaded and also it will 
overload all methods to work with `GPU()` device as their argument (instead of the standard 
`CPU()` device).

It's easy to construct a grid that lives on the GPU. Calling:

```julia
dev = GPU()
n, L = 16, 2.0
grid = OneDGrid(dev, n, L)

OneDimensionalGrid
  ├─────────── Device: GPU
  ├──────── FloatType: Float64
  ├────────── size Lx: 2.0
  ├──── resolution nx: 16
  ├── grid spacing dx: 0.125
  └─────────── domain: x ∈ [-1.0, 0.875]
```

gives out a grid whose arrays are `CuArrays`. (If you simply call `grid = OneDGrid(n, L)` it
defaults to `grid = OneDGrid(CPU(), n, L)`.)

When we construct the `Params`, `Vars`, and `Equation` for our problem we need to make sure
that we create arrays on the appropriate device, i.e., `Arrays` for `CPU` or `CuArrays` for
the `GPU`. Function `ArrayType` is useful in constructing appropriately chosen arrays.

```@docs
ArrayType
```

The `FourierFlows.Problem` constructor then takes an optional positional argument `dev::Device` If not provided anything, the default values for `dev=CPU()`.

```julia
problem = Problem(equation, stepper, dt, grid, vars, params, GPU())
```

The `FourierFlows.Diffusion` module is written in a way such that switching from CPU to GPU 
is only a matter of calling `FourierFlows.Diffusion.Problem()` with `dev=GPU()`. All physics modules in [GeophysicalFlows.jl](https://github.com/FourierFlows/GeophysicalFlows.jl) can 
also seamlessly run on a GPU with `dev=GPU()` argument.