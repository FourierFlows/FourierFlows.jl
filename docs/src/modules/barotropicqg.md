# BarotropicQG Module

```math
\newcommand{\J}{\mathsf{J}}
```

### Basic Equations

This module solves the quasi-geostrophic barotropic vorticity equation on a beta-plane of variable fluid depth $H-h(x,y)$. The flow is obtained through a streamfunction $\psi$ as $(u, \upsilon) = (-\partial_y\psi, \partial_x\psi)$. All flow fields can be obtained from the quasi-geostrophic potential vorticity (QGPV). Here the QGPV is

$$\underbrace{f_0 + \beta y}_{\text{planetary PV}} + \underbrace{(\partial_x \upsilon
	- \partial_y u)}_{\text{relative vorticity}} +
	\underbrace{\frac{f_0 h}{H}}_{\text{topographic PV}}.$$

The dynamical variable is the component of the vorticity of the flow normal to the plane of motion, $\zeta\equiv \partial_x \upsilon- \partial_y u = \nabla^2\psi$. Also, we denote the topographic PV with $\eta\equiv f_0 h/H$. Thus, the equation solved by the module is:

$$\partial_t \zeta + \J(\psi, \underbrace{\zeta + \eta}_{\equiv q}) +
\beta\partial_x\psi = \underbrace{-\left[\mu + \nu(-1)^{n_\nu} \nabla^{2n_\nu}
\right] \zeta }_{\textrm{dissipation}} + f\ .$$

where $\J(a, b) = (\partial_x a)(\partial_y b)-(\partial_y a)(\partial_x b)$. On the right hand side, $f(x,y,t)$ is forcing, $\mu$ is linear drag, and $\nu$ is hyperviscosity. Plain old viscosity corresponds to $n_{\nu}=1$. The sum of relative vorticity and topographic PV is denoted with $q\equiv\zeta+\eta$.

### Implementation

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{\zeta} = - \widehat{\J(\psi, q)} +\beta\frac{\mathrm{i}k_x}{k^2}\widehat{\zeta} -\left(\mu
+\nu k^{2n_\nu}\right) \widehat{\zeta}  + \widehat{f}\ .$$

In doing so the Jacobian is computed in the conservative form: $\J(f,g) =
\partial_y [ (\partial_x f) g] -\partial_x[ (\partial_y f) g]$.

Thus:

$$\mathcal{L} = \beta\frac{\mathrm{i}k_x}{k^2} - \mu - \nu k^{2n_\nu}\ ,$$
$$\mathcal{N}(\widehat{\zeta}) = - \mathrm{i}k_x \mathrm{FFT}(u q)-
	\mathrm{i}k_y \mathrm{FFT}(\upsilon q)\ .$$


## Examples

- `examples/barotropicqg/decayingbetaturb.jl`: An script that simulates decaying quasi-geostrophic flow on a beta-plane demonstrating zonation.

- `examples/barotropicqg/forcedbetaturb.jl`: An script that simulates forced-dissipative quasi-geostrophic flow on a beta-plane demonstrating zonation. The forcing is temporally delta-corraleted and its spatial structure is isotropic with power in a narrow annulus of total radius `kf` in wavenumber space.

- `examples/barotropicqg/ACConelayer.jl`: A script that simulates barotropic quasi-geostrophic flow above topography reproducing the results of the paper by

  > Constantinou, N. C. (2018). A barotropic model of eddy saturation. *J. Phys. Oceanogr.*, **48 (2)**, 397-411.
