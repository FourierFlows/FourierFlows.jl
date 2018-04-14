# TwoDTurb Module

This module solves two-dimensional incompressible turbulence. The flow is given
through a streamfunction $\psi$ as $(u,\upsilon) = (-\partial_y\psi, \partial_x\psi)$.
The dynamical variable used here is the component of the vorticity of the flow
normal to the plane of motion, $q=\partial_x \upsilon- \partial_y u = \nabla^2\psi$.
The equation solved by the module is:

$$\partial_t q + J(\psi, q) = \underbrace{-\left[\mu(-1)^{n_\mu} \nabla^{2n_\mu}
+\nu(-1)^{n_\nu} \nabla^{2n_\nu}\right] q}_{\textrm{dissipation}} + f\ .$$

where $J(f,g) = (\partial_xf)(\partial_y g)-(\partial_x g)(\partial_y f)$. On
the right hand side, $f(x,y,t)$ is forcing, $\mu$ is hypoviscosity, and $\nu$ is
hyperviscosity. Plain old linear drag corresponds to $n_{\mu}=0$, while normal
viscosity corresponds to $n_{\nu}=1$.

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{q} = - \widehat{J(\psi, q)} -\left[\mu k^{2n_\mu}
+\nu k^{2n_\nu}\right] \widehat{q}  + \widehat{f}\ .$$

In doing so the Jacobian is computed in the conservative form: $J(f,g) =
\partial_y [ (\partial_x f) g] -\partial_x[ (\partial_y f) g]$.

Thus:

$$\mathcal{L} = -\mu k^{-2n_\mu} - \nu k^{2n_\nu}\ ,$$
$$\mathcal{N}(\widehat{q}) = - \mathrm{i}k_x \mathrm{FFT}(u q)-
	\mathrm{i}k_y \mathrm{FFT}(\upsilon q)\ .$$
