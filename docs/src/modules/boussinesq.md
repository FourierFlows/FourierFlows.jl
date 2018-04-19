# Thin-layer Boussinesq modules

```math
\newcommand{\bcdot}{\boldsymbol \cdot}
\newcommand{\bnabla}{\boldsymbol \nabla}
\newcommand{\pnabla}{\bnabla_{\! \perp}}

\newcommand{\com}{\, ,}
\newcommand{\per}{\, .}

\newcommand{\bu}{\boldsymbol u}
\newcommand{\bU}{\boldsymbol U}
\newcommand{\buu}{{\boldsymbol{u}}}
\newcommand{\bb}{{b}}
\newcommand{\pp}{{p}}
\newcommand{\ww}{{w}}
\newcommand{\uu}{{u}}
\newcommand{\v}{\upsilon}
\newcommand{\vv}{{\upsilon}}
\newcommand{\zzeta}{{\zeta}}
\newcommand{\oomega}{{\omega}}
\newcommand{\boomega}{\boldsymbol{\oomega}}

\newcommand{\bxh}{\widehat{\boldsymbol{x}}}
\newcommand{\byh}{\widehat{\boldsymbol{y}}}
\newcommand{\bzh}{\widehat{\boldsymbol{z}}}
\newcommand{\ii}{\mathrm{i}}
\newcommand{\ee}{\mathrm{e}}
\newcommand{\cc}{\mathrm{c.c.}}
\newcommand{\J}{\mathsf{J}}

\newcommand{\p}{\partial}
```

These modules solve various thin-layer approximations to the hydrostatic Boussinesq
equations. A thin-layer approximation is one that is appropriate for dynamics with
small aspect ratios, or small vertical scales and large horizontal scales.
Thin layer approximations include the shallow-water system, layered system, and
spectral approximations that apply a Fourier or Sin/Cos eigenfunction expansion in the
vertical coordinate to the Boussinesq equations, and truncate the expansion at just
two or three modes. Approximations of this last flavor are described here.

The three-dimensional rotating, stratified, hydrostatic Boussinesq equations are

```math
\p_t\buu + \left ( \buu \bcdot \bnabla \right ) \buu + f \bzh \times \buu + \bnabla \pp = D^{\buu} \com \\
\p_z \pp = \bb \com \\
\p_t\bb + \ww N^2 = D^\bb \com \\
\bnabla \bcdot \buu = 0 \com
```

where $\bu = (u, \v, w)$ is the three-dimensional velocity, $b$ is buoyancy, $p$ is pressure, $N^2$ is the
buoyancy frequency (constant), and $f$ is the rotation or Coriolis frequency. The operators $D^{\buu}$ and $D^{\bb}$ are arbitrary dissipation that we define only after projecting onto vertical Fourier or Sin/Cos modes.
Taking the curl of the horizontal momentum equation yields an evolution
equation for vertical vorticity, $\zzeta = \p_x \vv - \p_y \uu$:

```math
\p_t\zzeta + \buu \bcdot \bnabla \zzeta - \left (f \bzh + \boomega \right )
    \bcdot \bnabla \ww = D^{\zzeta} \per
```

## Vertically Fourier Boussinesq

The vertically-Fourier Boussinesq module solves the Boussinesq system obtained by expanding the hydrostatic
Boussinesq equations in a Fourier series. The horizontal velocity $\uu$, for example, is expanded with

```math
\uu(x, y, z, t) \mapsto U(x, y, t) + \ee^{\ii m z} u(x, y, t) + \ee^{-\ii m z} u^*(x, y, t) \com
```

The other variables $\vv$, $\bb$, $\pp$, $\zzeta$, and $\boomega$ are expanded identically. The barotropic
horizontal velocity is $V$ and the barotropic vertical vorticity is $Z = \p_x V - \p_y U$. The barotropic
vorticity obeys
```math
\p_t Z + \J \left ( \Psi, Z \right )
    + \bnabla \bcdot \left ( \bu \zeta^* \right ) + \ii m \pnabla \bcdot \left ( \bu w^* \right ) + \cc
    = D_0 Z \com
```

where $\cc$ denotes the complex conjugate and contraction with $\pnabla = -\p_y \bxh + \p_x \byh$
gives the vertical component of the curl.

The baroclinic components obey

```math
\p_t u - f \v + \p_x p = - \J \left ( \Psi, u \right ) - \bu \bcdot \bnabla U + D_1 u \com \\
\p_t \v + f u + \p_y p = - \J \left ( \Psi, \v \right ) - \bu \bcdot \bnabla V + D_1 \v \com \\
\p_t p - \tfrac{N^2}{m} w = - \J \left ( \Psi, p \right ) + D_1 p \per
```

The dissipation operators are defined

```math
D_0 = \nu_0 (-1)^{n_0} \nabla^{2n_0} + \mu_0 (-1)^{m_0} \nabla^{2m_0} \com \\
D_1 = \nu_1 (-1)^{n_1} \nabla^{2n_1} + \mu_1 (-1)^{m_1} \nabla^{2m_1}
```

where $U$ is the barotropic velocity and $u$ is the amplitude of the first baroclinic mode with periodic
vertical structure $\mathrm{e}^{\mathrm{i} m z}$.

### Implementation

Coming soon.


## Vertically Cosine Boussinesq

The vertically-Cosine Boussinesq module solves the Boussinesq system obtained by expanding the
hydrostatic Boussinesq equations in a Sin/Cos series. The horizontal velocity, for example, becomes

```math
\uu(x, y, z, t) \mapsto U(x, y, t) + \cos(mz) u(x, y, t) \per
```

The horizontal velocity $\vv$, pressure $\pp$, and vertical vorticity $\zzeta$ are also expanded in $\cos(mz)$,
where $Z = \p_x V - \p_y U$ denotes the barotropic component of the vertical vorticity. The vertical velocity $\ww$
and buoyancy $\bb$ are expanded with $\sin(mz)$.

### Basic governing equations

Projecting the vertical vorticity equation onto Sin/Cos modes an equation for the evolution of $Z$,

```math
\p_t Z + \J \left ( \Psi, Z \right )
    + \tfrac{1}{2} \bnabla \bcdot \left ( \bu \zeta \right ) + \tfrac{m}{2} \pnabla \bcdot \left ( \bu w \right )
    = D_0 Z \com
```

where $\J(a, b) = (\p_x a)(\p_y b) - (\p_y a)(\p_x b)$ is the Jacobian operator, contraction with $\pnabla = -\p_y \bxh + \p_x \byh$ gives the vertical component of the curl, and $\Psi$ is the barotropic streamfunction defined so that

```math
\bU = -\p_y\Psi \bxh + \p_x\Psi \byh \qquad \text{and} \qquad Z = \nabla^2 \Psi \per
```

The baroclinic components obey

```math
\p_t u - f \v + \p_x p = - \J \left ( \Psi, u \right ) - \bu \bcdot \bnabla U + D_1u \com \\
\p_t \v + f u + \p_y p = - \J \left ( \Psi, \v \right ) - \bu \bcdot \bnabla V + D_1\v \com \\
\p_t p - \tfrac{N^2}{m} w = - \J \left ( \Psi, p \right ) + D_1p \per
```

The dissipation operators are defined

```math
D_0 = \nu_0 (-1)^{n_0} \nabla^{2n_0} + \mu_0 (-1)^{m_0} \nabla^{2m_0} \com \\
D_1 = \nu_1 (-1)^{n_1} \nabla^{2n_1} + \mu_1 (-1)^{m_1} \nabla^{2m_1} \com
```

where $2n_0$ and $2m_0$ are the hyperviscous orders of the arbitrary barotropic dissipation operators
with coefficients $\nu_0$ and $\mu_0$, while $2n_1$ and $2m_1$ are the orders of the baroclinic
dissipation operators.

A passive tracer in the Vertically Cosine Boussinesq system is assumed to satisfy a no-flux condition
at the upper and lower boundaries, and thus expanded in cosine modes so that

```math
c(x, y, z, t) = C(x, y, t) + \cos(mz) c(x, y, t) \per
```

The barotropic and baroclinic passive tracer components then obey

```math
\p_t C + \J(\Psi, C) + \tfrac{1}{2} \bnabla \bcdot \left ( \bu c \right ) =
    \kappa (-1)^{n_{\kappa}} \nabla^{2{n_{\kappa}}} C \com \\
\p_t c + \J(\Psi, c) + \bu \bcdot \bnabla C = \kappa (-1)^{n_{\kappa}} \nabla^{2{n_{\kappa}}} c \com
```

where $\kappa$ and $n_{\kappa}$ are the tracer hyperdiffusivity and order of the hyperdiffusivity, respectively.
The choice $n_{\kappa} = 1$ corresponds to ordinary Fickian diffusivity.

### Implementation

Coming soon.
