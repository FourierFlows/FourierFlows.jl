# TwoDTurb Module

```math
\newcommand{\J}{\mathsf{J}}
```

### Basic Equations

This module solves two-dimensional incompressible turbulence. The flow is given
through a streamfunction $\psi$ as $(u,\upsilon) = (-\partial_y\psi, \partial_x\psi)$.
The dynamical variable used here is the component of the vorticity of the flow
normal to the plane of motion, $q=\partial_x \upsilon- \partial_y u = \nabla^2\psi$.
The equation solved by the module is:

$$\partial_t q + \J(\psi, q) = \underbrace{-\left[\mu(-1)^{n_\mu} \nabla^{2n_\mu}
+\nu(-1)^{n_\nu} \nabla^{2n_\nu}\right] q}_{\textrm{dissipation}} + f\ .$$

where $\J(a, b) = (\partial_x a)(\partial_y b)-(\partial_y a)(\partial_x b)$. On
the right hand side, $f(x,y,t)$ is forcing, $\mu$ is hypoviscosity, and $\nu$ is
hyperviscosity. Plain old linear drag corresponds to $n_{\mu}=0$, while normal
viscosity corresponds to $n_{\nu}=1$.

### Implementation

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{q} = - \widehat{J(\psi, q)} -\left[\mu k^{2n_\mu}
+\nu k^{2n_\nu}\right] \widehat{q}  + \widehat{f}\ .$$

In doing so the Jacobian is computed in the conservative form: $\J(f,g) =
\partial_y [ (\partial_x f) g] -\partial_x[ (\partial_y f) g]$.

Thus:

$$\mathcal{L} = -\mu k^{-2n_\mu} - \nu k^{2n_\nu}\ ,$$
$$\mathcal{N}(\widehat{q}) = - \mathrm{i}k_x \mathrm{FFT}(u q)-
	\mathrm{i}k_y \mathrm{FFT}(\upsilon q)\ .$$


## Examples

- `examples/twodturb/McWilliams.jl`: A script that simulates decaying two-dimensional turbulence reproducing the results of the paper by

  > McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. *J. Fluid Mech.*, **146**, 21-43

- `examples/twodturb/IsotropicRingForcing.jl`: A script that simulates stochastically forced two-dimensional turbulence. The forcing is temporally delta-corraleted and its spatial structure is isotropic with power in a narrow annulus of total radius `kf` in wavenumber space.
