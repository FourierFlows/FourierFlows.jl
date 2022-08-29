# GPU

FourierFlows allows you to easily construct and run problems on GPUs.

Upon calling

```julia
using FourierFlows
```

FourierFlows.jl will check whether any CUDA enabled device is present. If such a device is 
found then FourierFlows.jl makes sure that CUDA related packages are loaded and also it will 
overload all methods to work with `GPU()` device as their argument (instead of the standard 
`CPU()` device).

It's easy to construct a grid that lives on the GPU. Calling:

```julia
dev = GPU()
n, L = 16, 2.0
grid = OneDGrid(dev; n, L)

OneDimensionalGrid
  ├─────────── Device: GPU
  ├──────── FloatType: Float64
  ├────────── size Lx: 2.0
  ├──── resolution nx: 16
  ├── grid spacing dx: 0.125
  ├─────────── domain: x ∈ [-1.0, 0.875]
  └─ aliased fraction: 0.3333333333333333
```

gives out a grid whose arrays are `CuArrays`. (Calling `OneDGrid(; n, L)` defaults to CPU, i.e., 
`OneDGrid(CPU(); n, L)`.)

When we construct the `Params`, `Vars`, and `Equation` for our problem we need to make sure
that we create arrays on the appropriate device, i.e., `Arrays` for `CPU` or `CuArrays` for
the `GPU`. Function `ArrayType` is useful in constructing appropriately chosen arrays.

```@docs
ArrayType
```

The `FourierFlows.Problem` constructor then takes an optional positional argument 
`dev::Device`. If not provided anything, the default values for `dev=CPU()`.

```julia
problem = Problem(equation, stepper, dt, grid, vars, params, GPU())
```

The `FourierFlows.Diffusion` module is written in a way such that switching from CPU to GPU 
is only a matter of calling `FourierFlows.Diffusion.Problem()` with `dev=GPU()`. All physics 
modules in [GeophysicalFlows.jl](https://github.com/FourierFlows/GeophysicalFlows.jl) can 
also seamlessly run on a GPU with `dev=GPU()` argument.


## Selecting GPU device

`FourierFlows.jl` can only utilize a single GPU. If your machine has more than one GPU available, 
then functionality within `CUDA.jl` package enables the user to choose the GPU device that 
`FourierFlows.jl` should use. The user is referred to the [`CUDA.jl` Documentation](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#Device-Management); in particular, [`CUDA.devices`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.devices) and [`CUDA.CuDevice`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.CuDevice).
