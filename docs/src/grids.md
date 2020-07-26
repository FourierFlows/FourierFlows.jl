# Grids


1D, 2D, and 3D grids are supported.

A one-dimensional grid with $n_x = 64$ grid points and length $L_x = 2.0$ is constructed using

```@meta
DocTestSetup = quote
    using FourierFlows
end
```

```jldoctest
julia> nx, Lx = 64, 2π
julia> grid = OneDGrid(nx, Lx)
OneDimensionalGrid
                Device: CPU
             FloatType: Float64
               size Lx: 6.283185307179586
         resolution nx: 64
       grid spacing dx: 0.09817477042468103
                domain: x ∈ [-3.141592653589793, 3.0434178831651124])```
```

The grid spacing is simply $L_x/n_x$. Note that $x\in [-L_x/2, L_x/2 - L_x/n_x]$ since periodicity
implies that the value of any field at $x=L_x/2$ is the same as $x=-L_x/2$.

We can define an array $f$ on this grid as

```@setup 1
using FourierFlows
using LinearAlgebra: mul!, ldiv!
using Plots # hide
nx, Lx = 64, 2π
grid = OneDGrid(nx, Lx)
```


```@example 1
f = @. sin(2 * grid.x) + 1/2 * cos(5 * grid.x)
nothing # hide
```

and plot it

```@example 1
plot(grid.x, f, label="f", xlabel="x")
savefig("assets/plot1.png"); nothing # hide
```

![](assets/plot1.png)

The Fourier transform of $f$ is $\hat{f}$. Since our `f` array is real-valued then we should use the `real-FFT` algorithm. The real-valued FFT transform only saves the Fourier coefficients for $k\ge 0$; the coefficients for negative wavenumbers can be obtained via $\hat{f}(-k) = \hat{f}(k)^\star$.

The `grid` includes the `FFT` plans for both real-valued and complex valued transforms:

```@example 1
grid.fftplan
```

```@example 1
grid.rfftplan
```

We use the convention that variables names with `h` at the end stand for variable-hat, i.e., $\hat{f}$  is the Fourier transform of $f$ and is stored in array `fh`.

```@example 1
fh = zeros(FourierFlows.cxeltype(grid), grid.nkr)

mul!(fh, grid.rfftplan, f)

plot(grid.kr, [real.(fh / nx), imag.(fh / nx)],
          label = ["real\\( rfft\\(f\\) \\)" "imag\\( rfft\\(f\\) \\)"],
         xlabel = "k",
          xlims = (-0.5, 10.5),
         marker = :auto)

savefig("assets/plot2.png"); nothing # hide
```

![](assets/plot2.png)

We can compute its derivative via Fourier transforms. To do that we can use the
`FFTW` plans that are constructed with the grid. First we allocate an empty array
where the values of the derivative will be stored,

```@example 1
∂ₓf = similar(f)
∂ₓfh = similar(fh)

nothing # hide
```

Array `∂ₓfh` is populated with $\mathrm{i} k \hat{f}$

```@example 1
@. ∂ₓfh = im * grid.kr * fh

nothing # hide
```

Then `∂ₓf` is computed with an inverse Fourier tranform. The latter is obtained 
again using the `FFTW` plans but now via `ldiv!`:

```@example 1
ldiv!(∂ₓf, grid.rfftplan, ∂ₓfh)

plot(grid.x, [f ∂ₓf], label=["f" "∂f/∂x"], xlabel="x")

savefig("assets/plot3.png"); nothing # hide
```

![](assets/plot3.png)
