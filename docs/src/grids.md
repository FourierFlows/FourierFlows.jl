# Grids


1D, 2D, and 3D grids are supported. We demonstrate here the construction of a 
one-dimensional grid and how one can use it to perform Fourier transforms and 
compute spatial derivatives.

A one-dimensional grid with ``n_x = 64`` grid points and length ``L_x = 2 \pi`` is 
constructed by

```@meta
DocTestSetup = quote
    using FourierFlows
end
```

```jldoctest
using FourierFlows

nx, Lx = 64, 2π

grid = OneDGrid(; nx, Lx, x0 = -Lx/3)

# output

OneDimensionalGrid
  ├─────────── Device: CPU
  ├──────── FloatType: Float64
  ├────────── size Lx: 6.283185307179586
  ├──── resolution nx: 64
  ├── grid spacing dx: 0.09817477042468103
  ├─────────── domain: x ∈ [-2.0943951023931953, 4.090615434361711]
  └─ aliased fraction: 0.3333333333333333
```

The grid domain is, by default, constructed symmetrically around ``x = 0``, but this can be 
altered using the `x0` keyword argument of the `OneDGrid` constructor. The grid spacing 
is ``L_x / n_x``. Note that the last point of the domain is a grid-spacing before ``L_x / 2``. 
This is because periodicity implies that the values of any field at the end-points of the 
domain are equal and, therefore, grid-point values at both these end-points are redundant.

We can define an array `u` that contains the values of a function ``u(x)`` on this 
grid as

```@setup 1
using FourierFlows
using LinearAlgebra: mul!, ldiv!
using CairoMakie
set_theme!(Theme(linewidth = 3, fontsize = 20))
nx, Lx = 64, 2π
grid = OneDGrid(; nx, Lx, x0 = -Lx/3)
```

```@example 1
u = @. sin(2 * grid.x) + 1/2 * cos(5 * grid.x)

nothing # hide
```

Note that we chose a function that *is* periodic on our domain. We can visualize
`u` by

```@example 1
using CairoMakie

lines(grid.x, u, label="u", axis = (xlabel="x",))
axislegend()

current_figure() # hide
```

Function ``u(x)`` can be expanded in Fourier series:

```math
u(x) = \sum_{k} \hat{u}(k) \, e^{i k x} ,
```

where ``\hat{u}(k)`` is Fourier transform of ``u(x)`` and ``k`` the discrete set of 
wavenumbers that fit within our finite domain. We can compute ``\hat{u}`` via a 
Fast Fourier Transform (FFT). Since our `u` array is real-valued then we should 
use the `real-FFT` algorithm. The real-valued FFT transform only saves the Fourier 
coefficients for ``k \ge 0``; the coefficients for negative wavenumbers can be 
obtained via ``\hat{u}(-k) = \hat{u}(k)^{*}``.

The wavenumbers used in FFT are contained in `grid.k` and they are ordered as:
```math
\frac{2\pi}{L_x} \{ 0, 1, \dots, n_x/2-1, -n_x/2, -n_x/2+1, \dots, -1 \} ,
```
while the wavenumbers for real FFT are in `grid.kr`:

```math
\frac{2\pi}{L_x} \{ 0, 1, \dots, n_x/2-1 \} .
```


The `grid` also includes the `FFT` plans for both real-valued and complex valued transforms:

```@example 1
grid.fftplan
```

```@example 1
grid.rfftplan
```

We use the convention that variables names with `h` at the end stand for variable-hat, i.e., 
``\hat{u}``  is the Fourier transform of ``u`` and is stored in array `uh`. Since `u` is of 
size ``n_x``, the real-Fourier transform should be of size ``n_{kr} = n_x/2 + 1``.

```@example 1
uh = Complex.(zeros(grid.nkr))

nothing # hide
```

The `FFT` transform is done as an in-place matrix multiplication using `mul!`.

```@example 1
using LinearAlgebra: mul!

mul!(uh, grid.rfftplan, u)

nothing # hide
```

The `FFT` algorithm does not output exactly the Fourier coefficients ``\hat{u}(k)`` but
rather, due to different normalization, `FFT` outputs something proportional to ``\hat{u}(k)``. 
To obtain ``\hat{u}`` we need to divide the `FFT` output by the length of the original
array and by ``e^{i k x_0}``, where ``x_0`` is the first point of our domain array.

```@example 1
uhat = @. uh / (nx * exp(im * grid.kr * grid.x[1])) # due to normalization of FFT

fig = Figure()
ax = Axis(fig[1, 1];
          xlabel = "k",
          limits = ((-0.5, 10.5), (-0.55, 0.55)),
          xticks = 0:10,
          yticks = -0.5:0.25:0.5)

scatterlines!(ax, grid.kr, real.(uhat);
              markersize = 16, label = "real(û)")

scatterlines!(ax, grid.kr, imag.(uhat);
              markersize = 22, marker = :diamond, label = "imag(û)")

axislegend()

current_figure() # hide
```

We can compute the derivative of ``u(x)`` via Fourier transforms. To do that we use the
`FFTW` plans that are constructed with the grid. First we allocate some empty arrays
where the values of the derivative will be stored, both in physical and in Fourier space:

```@example 1
∂ₓu  = similar(u)
∂ₓuh = similar(uh)

nothing # hide
```

The grid contains the wavenumbers (both for real-value functions `grid.kr` and 
for complex-valued functions `grid.k`). We populate array `∂ₓuh` is with ``\mathrm{i} k \hat{u}``:

```@example 1
@. ∂ₓuh = im * grid.kr * uh

nothing # hide
```

Then the derivative in physical space, `∂ₓu`, is obtained with an inverse Fourier 
transform. The latter is obtained again using the `FFTW` plans but now via `ldiv!`:

```@example 1
using LinearAlgebra: ldiv!

ldiv!(∂ₓu, grid.rfftplan, ∂ₓuh)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x")

lines!(ax, grid.x, u, label = "u")
lines!(ax, grid.x, ∂ₓu, label = "∂u/∂x")

axislegend()

current_figure() # hide
```
