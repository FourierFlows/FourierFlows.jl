# Aliasing


```@setup 1
using FourierFlows
using Plots
Plots.scalefontsizes(1.25)
Plots.default(lw=3)
```

In pseudospectral methods, often aliasing errors come into play. These errors originate from
the discrete version of the grid. A grid discretized with ``n_x`` points can only resolve a 
total of ``n_x`` wavenumbers in Fourier space. 

On a grid with a total of ``n_x`` wavenumbers, both harmonics ``e^{2\pi i k x / L_x}`` and 
``e^{2\pi i (k+n_x) x / L_x}``, with ``k`` an integer, are indistinguishable when evaluated
on the discrete grid-points of this grid. When we compute nonlinear terms in physical space, 
we may end up with terms that project on higher wavenumbers, beyond those that our grid can 
represent. In that case, those wavenumbers will be erroneously projected onto some lower 
wavenumber mode that fits our domain.

Take, for example, functions ``\cos(4x)`` and ``\cos(6x)`` and let's see how they are represented 
on a grid ``x \in [-π, π)`` with ``n_x = 10``.

```@example 1
nx, Lx = 10, 2π
grid = OneDGrid(nx, Lx)

f1(x) = cos(4x)
f2(x) = cos(6x)

p = plot(grid.x, f1.(grid.x), lw=0, marker=:circle, c=:red, ms=8, ylims=(-1.6, 1.6), label="cos(4x)")
plot!(p, f1, lw=3, alpha=0.2, c=:red, xlims=(-Lx/2, Lx/2), label="")
plot!(p, grid.x, f2.(grid.x), lw=0, marker=:star5, ms=8.5, color=:blue, alpha=0.8, label="cos(6x)")
plot!(p, f2, lw=3, alpha=0.2, c=:blue, xlims=(-Lx/2, Lx/2), label="")

plot(p, xlabel="x", xlims=(-3.3, 3.3))

savefig("assets/plot4.svg"); nothing # hide
```

![](assets/plot4.svg)

The take home message is that on this grid we cannot distinguish harmonics of wavenumbers 4 and 6
and attempting to represent harmonics with wavenumber 6 on this grid will lead to aliasing errors.
For example, say that we are solving an equation on this grid and at some point we compute the product 
``\cos(2x) \cos(4x)``. The result is ``\frac1{2} \cos(2x) + \frac1{2} \cos(6x)``, but on this 
grid ``\cos(6x)`` is indistinguishable from ``\cos(4x)`` and, therefore, we get an answer 
which is the sum of ``\frac1{2} \cos(2x) + \frac1{2} \cos(4x)``!

There are two ways to avoid aliasing errors we either *(i)* need to discard some of the wavenumber 
components in Fourier space before we tranform to physical space, or *(ii)* pad our Fourier 
represresentation with more wavenumbers that will have zero power. In FourierFlows.jl the former
is implemented

!!! info "De-aliasing scheme"
    FourierFlows.jl curently implement dealiasing by zeroing out the top-`aliased_fraction` 
    wavenumber components on a `grid`.

How many wavenumber components we need to discard depends on the order of the nonlinearity. For
quadradic nonlinearities, one would intuitively say that we need to discard the top-1/2 of the 
wavenumber components. However, Orszag (1972) pointed out that simply only discarding the 
top-1/3 of wavenumber components is enough. Actally, with Orszag's so-called 2/3-rule for dealiasing, 
still some aliasing errors occur, but only into wavenumbers that will be zero-ed out next time 
we dealias.

When constructing a `grid` we can specify the `aliased_fraction` parameter. By default, this is 
set to ``1/3``, appropriate for quadratic nonlinearities. Then `dealias!(fh, grid)` will zero-out 
the top-`aliased_fraction` wavenumber components of `fh`. 

If we construct a grid with `aliased_fraction=0`, e.g.,

```@example 1
grid_nodealias = OneDGrid(nx, Lx; aliased_fraction=0)
```

then `dealias!(fh, grid_nodealias)` will have _no effect_ whatsoever on `fh`.
