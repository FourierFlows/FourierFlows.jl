```math
\newcommand{\sqr}{\mbox{sqr}}
\newcommand{\saw}{\mbox{saw}}
\newcommand{\ind}{\mbox{ind}}
\newcommand{\sgn}{\mbox{sgn}}
\newcommand{\erfc}{\mbox{erfc}}
\newcommand{\erf}{\mbox{erf}}

%% An average
\newcommand{\avg}[1]{\mathrm{avg}[ {#1} ]}
%% The right way to define new functions
\newcommand{\sech}{\mathop{\rm sech}\nolimits}
\newcommand{\cosech}{\mathop{\rm cosech}\nolimits}

%% A nice definition
\newcommand{\defn}{\stackrel{\mathrm{def}}{=}}
%%%%%%%%% %%%%

\newcommand{\ol}[1]{\overline{#1}}


%% Various boldsymbols
\newcommand{\bx}{\boldsymbol{x}}
\newcommand{\by}{\boldsymbol{y}}
\newcommand{\bq}{\boldsymbol{q}}
\newcommand{\bpsi}{\boldsymbol{\psi}}
\newcommand{\bu}{\boldsymbol{u}}
\newcommand{\bG}{\boldsymbol{\mathcal{G}}}
\newcommand{\G}{\mathcal{G}}
\newcommand{\ba}{\boldsymbol{a}}
\newcommand{\bb}{\boldsymbol{b}}
\newcommand{\bc}{\boldsymbol{c}}
\newcommand{\bv}{\boldsymbol{v}}
\newcommand{\bk}{\boldsymbol{k}}
\newcommand{\bX}{\boldsymbol{X}}
\newcommand{\br}{\boldsymbol{r}}
\newcommand{\J}{\mathsf{J}}
\newcommand{\D}{\mathsf{D}}
\renewcommand{\L}{\mathsf{L}}
\newcommand{\sL}{\mathsf{L}}
%\newcommand{\G}{\boldsymbol{\mathsf{G}}}
\newcommand{\bA}{\boldsymbol {A}}
\newcommand{\bU}{\boldsymbol {U}}
\newcommand{\bE}{\boldsymbol {E}}
\newcommand{\bJ}{\boldsymbol {J}}
\newcommand{\bXX}{\boldsymbol {\mathcal{X}}}
\newcommand{\bFF}{\ensuremath {\boldsymbol {F}}}
\newcommand{\bF}{\ensuremath {\boldsymbol {F}^{\sharp}}}
\newcommand{\bL}{\ensuremath {\boldsymbol {L}}}
\newcommand{\bI}{\ensuremath {\boldsymbol {I}}}
\newcommand{\bN}{\ensuremath {\boldsymbol {N}}}

\newcommand{\I}{\ensuremath {\mathsf{I}}}
\renewcommand{\L}{\ensuremath {\mathsf{L}}}
\renewcommand{\S}{\ensuremath {\mathsf{S}}}

\newcommand{\bSigma}{\ensuremath {\boldsymbol {\Sigma}}}
\newcommand{\kmax}{k_{\mathrm{max}}}
\newcommand\bnabla{\boldsymbol{\nabla}}
\newcommand\bcdot{\boldsymbol{\cdot}}

\def\ii{{\rm i}}
\def\dd{{\rm d}}
\def\ee{{\rm e}}
\def\DD{{\rm D}}
%%% Cals here %%%%%

%%%%%  Euler caligraphics %%%%%
\newcommand{\A}{\mathscr{A}}
% \newcommand{\B}{\mathscr{B}}
\newcommand{\B}{\mathcal{B}}
\newcommand{\E}{\mathscr{E}}
\newcommand{\F}{\mathscr{F}}
\newcommand{\K}{\mathscr{K}}
\newcommand{\N}{\mathscr{N}}
\newcommand{\U}{\mathscr{U}}
\newcommand{\LL}{\mathscr{L}}
\newcommand{\M}{\mathscr{M}}
\newcommand{\T}{\mathscr{T}}
\def\la{\langle}
\def\ra{\rangle}
\def\laa{\left \langle}
\def\raa{\right \rangle}
\def\Ek{\mathrm{Ek}}
\newcommand{\hzon}{h_{\mathrm{zon}}}
\newcommand{\lap}{\triangle}
\newcommand{\p}{\partial}
\newcommand{\half }{\tfrac{1}{2}}
\newcommand{\grad}{\boldsymbol {\nabla}}
\newcommand{\pde}{\textsc{pde}}
\newcommand{\ode}{\textsc{ode}}
\newcommand{\cc}{\textsc{cc}}
\newcommand{\dc}{\textsc{dc}}
\newcommand{\dbc}{\textsc{dbc}}
\newcommand{\byu}{\textsc{byu}}
\newcommand{\rhs}{\textsc{rhs}}
\newcommand{\lhs}{\textsc{lhs}}
\newcommand{\com}{\,,}
\newcommand{\per}{\,.}
\newcommand{\z}{\zeta}
\newcommand{\h}{\eta}
\renewcommand{\(}{\left(}
\renewcommand{\[}{\left[}
\renewcommand{\)}{\right)}
\renewcommand{\]}{\right]}
\newcommand{\<}{\left\langle}
\renewcommand{\>}{\right\rangle}
\renewcommand{\A}{\mathcal{A}}
\renewcommand{\N}{\mathcal{N}}
\newcommand{\C}{\mathcal{C}}
\newcommand{\transp}{\textrm{T}}
\newcommand{\zhat}{\widehat{\mathbf{z}}}

\newcommand{\bit}{\vphantom{\dot{W}}}
\newcommand{\sd}{b}
```

# Forcing

The code implements forcing in various modules (currently in `TwoDTurb` and `BarotropicQG`). Forcing can be either
deterministic or stochastic (random). For deterministic forcing the implementation is straightforward; for stochastic
forcing there are two main train of thoughts: Itô calculus and Stratonovich calculus.

Both stochastic calculi give the same results. But once we decide to use one of the two calculi we have to remain
consistent and use that calculus for everywhere. There is a lot of confusion and mostly the confusion stems from not
using the same stochastic calculus consistently throughout the computation but rather interchanging between the two.

`FourierFlows` uses Stratonovich calculus throughout the code. This choise was made because Stratonovich calculus
works the same with both stochastic and deterministic forcing, i.e. with Stratonovich calculus we have the same
chain rules for differentiation for stochastic functions as the chain rules we learn in normal-deterministic calculus).
Therefore, the code written as is does not really "care" of what forcing the user implements.

If you are interested in learning more regarding the two stochastic calculi and how they are numerically
implemented then read on; otherwise you can skip this section of the documentation and go to the Module Tutorials.

## Stochastic Differential Equations (SDEs)

A differential equation in the form:

```math
\begin{equation}
	\frac{\dd x}{\dd t} = f(x)\com\quad x(t_0)=0\com
\end{equation}
```

can also be written in an integral form:

```math
\begin{equation}
	x(t) = \int_{t_0}^{t} f(x(s))\,\dd s\per
\end{equation}
```

In a similar manner, a stochastic differential equation

```math
\begin{equation}
	\dd x = f(x)\,\dd t + g(x)\,\dd W_t\com\quad x(t_0)=0\com
\end{equation}
```

with $\dd W_t$ a white-noise process, can be written in an integral form as:

```math
\begin{equation}
	x(t) = \int_{t_0}^{t} f(x(s))\,\dd s + \int_{t_0}^{t} g(x(s))\,\dd W_s \per
\end{equation}
```

Of course now, the last integral is a stochastic integral and there is not a single straight-forward way of
computing it --- there are a lot of different ways we can approximate it as a Riemannian sum and each of them
leads to a different answer. The two most popular ways for computing such stochastic integrals are:

```math
{\color{Green}\text{Itô}: \int_{t_0}^{t} g(x(s))\,\dd W_s\approx\sum_{j} g\left(x(t_j)\right)(W_{j+1}-W_j)}\com\\
{\color{Magenta}\text{Stratonovich}: \int_{t_0}^{t} g(x(s))\,\dd W_s \approx \sum_{j} g\left(x\left(\half(t_j+t_{j+1})\right)\right)(W_{j+1}-W_j)}\per
```

Because the white noise process is not continuous the two definitions do not converge to the same result; the two
definitions give thoroughly different results. And to overcome that they come along with different chain rules,
i.e., chain rules that are not necessarily the same as those in plain old calculus.

An SDE can be written also in differential form. Because we cannot formally form $\dd W/\dd t$, since $W$ is
nowhere differentiable, we write an SDE in differential form as:

```math
{\color{Green}\text{Itô}: \dd x_t = f(x_t)\dd t + g(x_t)\dd W_t}\com\\
{\color{Magenta}\text{Stratonovich}: \dd x_t = f(x_t)\dd t + g(x_t)\circ\dd W_t\per}
```

The circle in $g(x_t)\circ\dd W_t$ is used to differentiate between Itô or Stratonovich calculus.

A variable change $y=G(x)$ is done as follows according to the two different calculi:

```math
{\color{Green}\text{Itô}: \dd y_t = \frac{\dd G}{\dd x}\dd x_t + \half g(x_t)^2 \frac{\dd^2 G}{\dd x^2}\dd t =\left[ \frac{\dd G}{\dd x}f(x_t) + \half g(x_t)^2 \frac{\dd^2 G}{\dd x^2}\right]\dd t + \frac{\dd G}{\dd x}g(x_t)\dd W_t}\com\\
{\color{Magenta}\text{Stratonovich}: \dd y_t  = \frac{\dd G}{\dd x}\dd x_t =\frac{\dd G}{\dd x} f(x_t) \dd t + \frac{\dd G}{\dd x}g(x_t)\circ\dd W_t}\per
```

The above are the so called *stochastic chain rules*. All derivatives of $G$ are evaluated at $x_t$.

It's easy to see that the extra drift-term in Itô's interpretation of the stochastic integral,
i.e., ${\color{Green}\half g^2\, \dd^2G/\dd x^2}$  is *exactly* equal to the ensemble mean of the
Stratonovich stochastic integral. That's because the Itô stochastic integral has, by construction,
zero ensemble mean since at every instant the noise is multiplied with $g$ evaluated before the action of the
noise; $g$ and $\dd W$ are uncorrelated and thus:

```math
\begin{equation}
{\color{Green}\laa g(x_t)\dd W_t \raa =0}\quad\text{while}\quad {\color{Magenta}\laa g(x_t)\circ\dd W_t \raa \ne 0}\per
\end{equation}
```

The above is demonstrated by evaluating the simple stochastic integral:

```math
{\color{Green}\text{Itô}: \laa \int_{t_0}^{t} W_s\,\dd W_s \raa \approx\sum_{j} \laa W_j(W_{j+1}-W_j)\raa}\\
{\color{Green}\hspace{7.3em} = \sum_{j} \laa W_j W_{j+1}\raa - \laa W_jW_j\raa \sim \sum_{j} t_j - t_j = 0} \com\\
{\color{Magenta}\text{Stratonovich}: \laa\int_{t_0}^{t} W_s\circ\dd W_s\raa \approx \sum_{j} \laa \frac1{2}(W_{j} + W_{j+1}) (W_{j+1}-W_j)\raa }\\
{\color{Magenta}\hspace{7.3em} = \frac1{2}\sum_{j} \laa W_{j+1} W_{j+1}\raa - \laa W_{j} W_{j}\raa  \sim \frac1{2}\sum_{j} t_{j+1} - t_j = \frac{t}{2}}\per
```

SDEs rarely can be solved in closed form; most often numerical solution of SDEs is brought to the rescue.
Itô calculus has the advantage that is very easily implemented numerically.
On the other hand, the chain rule in Stratonovich calculus coincides with that in normal calculus. This stems from the fact
that in the Stratonovich interpretation the white noise process is as a series of colored noise processes with the de-correlation time tending to zero. This made Stratonovich calculus more popular in the physics community.
A nice discussion on the differences and similarities between the two calculi is given by [van Kampen](https://doi.org/10.1007/BF01007642).

## A simple Stochastic Differential Equation (SDE): the Ornstein--Uhlenbeck process

One of the simpler SDEs is the Ornstein--Uhlenbeck process. A variation of which is:

```math
\begin{equation}
	x(t) = \int_{t_0}^{t} -\mu x(s)\,\dd s + \int_{t_0}^{t} \sqrt{\sigma}\,\dd W_s \per\label{eq:OU}
\end{equation}
```

Note that in differential form this is:

```math
\begin{equation}
	\dd x_t = -\mu x_t \,\dd t + \sqrt{\sigma}\,\dd W_s \per\label{eq:1}
\end{equation}
```

Luckily, here there is no need to distinguish between Itô and Stratonovich. This is because $g$ is independent of $x(t)$. But we stress that  this is often not the case; it is only a fortuitous coincident here.

How do we time-step this SDE numerically? Let us assume a discretization of time into time-steps
of $\tau$: $t_j=(j-1)\tau$. (What follows can be easily transfer to non-uniform time discretization.)
With that, we denote $x_j\defn x(t_j)$. Then the Euler--Mayorama time-step scheme for \eqref{eq:1} is

```math
\begin{equation}
	x_{j+1} = x_j + (-\mu x_j)\tau + \sqrt{\sigma}(W_{j+1}-W_j)\per
\end{equation}
```

Now let us ask the following question: How can we compute the work done by the noise?
In other words, if we are interested in the evolution of the "energy" $E\defn \half x^2$, how is the noise term
attributing in the growth of $E$? To answer that we first have to find the SDE that energy $E$ obeys.
But, in doing so, it is important to adopt a single interpretation for computing stochastic integrals as now a
transformation of variables is needed. That is, depending on whether we choose to interpret the stochastic integrals
according to Itô or to Stratonovich calculus, $E$ evolves as:

```math
\begin{equation}
{\color{Green}\text{Itô}:  \dd E_t  = \left( -2\mu E_t + \half \sigma \right)\dd t  + x_t \sqrt{\sigma}\dd W_t}\com\label{eq:E_ito}
\end{equation}
```

```math
\begin{equation}
{\color{Magenta}\text{Stratonovich}: \dd E_t  = -2\mu E_t  \dd t + x_t\circ \sqrt{\sigma}\dd W_t}\per\label{eq:E_str}
\end{equation}
```

How do we compute the work $P$ done by the noise? It is:

```math
{\color{Green}\text{Itô}: P_t = \half \sigma \dd t + \sqrt{\sigma} x_t \dd W_t \approx  \half \sigma + \sqrt{\sigma} x_j (W_{j+1}-W_j)\com}\\
{\color{Magenta}\text{Stratonovich}: P_t =  x_t \circ\sqrt{\sigma} \dd W_t \approx \sqrt{\sigma} x\left(\half(t_j+t_{j+1})\right)(W_{j+1}-W_j)}\per
```

Say we didn't know the rules for transforming Stratonovich to Itô and we were wondering what is the extra drift term
we have to include in the Itô formulations, i.e. the $\half\sigma$ term. We can compute the Itô's drift-term using
that it is exactly equal to $\la x_t\circ\sqrt{\sigma}\dd W_t\ra$; and for the latter we can use the "usual" calculus.
That is, rewrite \eqref{eq:OU} as:

```math
\begin{equation}
\dot{x} = -\mu x + \xi\com\label{eq:OUcont}
\end{equation}
```

where $\xi(t)$ is understood to be the "continuous" version of the white-noise process which is formally only
understood in terms of distributions. The forcing $\xi$ has the properties:

```math
\laa \xi(t)\raa = 0 \quad\text{and}\quad \laa \xi(t)\xi(t')\raa = \sigma \delta(t-t')\per
```

Thus we need to compute $\la P_t \ra = \la x(t) \xi(t) \ra$. But \eqref{eq:OUcont} has formally the solution:

```math
x(t) = \ee^{-\mu t} x(0) + \int_0^t \ee^{-\mu(t-s)}\xi(s)\,\dd s\per
```

and utilizing the above we get

```math
\la P_t \ra = \la x(t) \xi(t)  \ra
=  \ee^{-\mu t} \underbrace{\la x(0)\xi(t)\ra}_{=0} + \int_0^t \ee^{-\mu(t-s)}\la \xi(t)\xi(s)\ra\,\dd s
= \sigma \int_0^t \ee^{-\mu(t-s)} \delta(t-s)\,\dd s =  \frac{\sigma}{2} \per
```

Above we used that $\int_0^t\delta(t-s)\dd s = \half$, which is consistent with Stratonovich symmetric interpretation
of stochastic integrals.

### Numerical implementation

How do we time-step \eqref{eq:E_ito}? We use the Euler--Maruyama time-stepping scheme:

```math
	E_{j+1} = E_j + \left(-2\mu E_j + \frac{\sigma}{2}\right)\tau + \sqrt{\sigma}x_j(W_{j+1}-W_j)\per
```
However, we cannot use Euler--Maruyama for time-stepping \eqref{eq:E_str} since the Euler--Maruyama is "Itô"-thinking.
To time-step \eqref{eq:E_str} we have to approximate $g$ in the middle of the time-step. There are many ways to do
that, one of which is the, so called, Euler--Heun method:

```math
	\widetilde{E}_{j+1} = E_j + (-2\mu E_j)\tau + \sqrt{\sigma}x_j(W_{j+1}-W_j)\com\\
	E_{j+1} = E_j + \left(-2\mu \frac{E_j+\widetilde{E}_{j+1}}{2}\right)\tau + \sqrt{\sigma}\frac{x_j+x_{j+1}}{2}(W_{j+1}-W_j)\per
```

![energy_comparison](assets/energy_comparison.png)

Figure above shows a comparison of the energy evolution as done from:
- direct computation as $\half x_t^2$,
- time-integration of \eqref{eq:E_ito}, and
- time-integration of \eqref{eq:E_str}.

Figures below show the ensemble mean energy budgets (using 1000 ensemble members) as computed using Itô and
Stratonovich. For the energy budget to close we have to be consistent: if we time-step the energy equation based
on Stratonovich calculus then we must compute the work also according to Stratonovich.
(For more details see `examples/forcing/simpleSDEItoStratonovich.jl`.

![energy_budgets_Ito](assets/energy_budgets_Ito.png)
![energy_budgets_Stratonovich](assets/energy_budgets_Stratonovich.png)

## A simple Stochastic Partial Differential Equation (SPDE)

We want now to transfer all the knowledge we got from the previous sections to PDEs. In particular we'll focus on the simple SPDE:

\begin{equation}
\partial_t \nabla^2\psi(\bx, t) =  -\mu \nabla^2\psi(\bx, t) + \xi(\bx,t)\com\label{eq:PDEcont}
\end{equation}

which is also equivalently written as:

```math
\dd \nabla^2\psi_{t}(\bx) = -\mu \nabla^2 \psi_{t} (\bx) \dd t + \sqrt{\sigma} \dd W_{t} (\bx) \per
```

The form \eqref{eq:PDEcont} is the continuous version understood in the Stratonovich interpretation
(similar to \eqref{eq:OUcont}). Thus, forcing $\xi$ obeys now:

```math
\la\xi(\bx,t)\ra = 0 \quad\text{and}\quad\la\xi(\bx,t)\xi(\bx',t') \ra= Q(\bx-\bx')\delta(t-t')\com
```

that is the forcing is white in time but spatially correlated; its spatial correlation is prescribed by the
function $Q$ which is, necessarily, homogeneous in all its arguments
(see discussion by [Constantinou](http://arxiv.org/abs/1503.07644); Appendix A).

The above describes the vorticity evolution of a two-dimensional fluid $\nabla^2\psi$ which is stochastically
forced while dissipated by linear drag $\mu$. The energy of the fluid is:

```math
E = \half\overline{|\grad\psi|^2}^{x,y} = -\half\overline{\psi\nabla^2\psi}^{x,y}\com
```

where the overbar denotes average over $x$ and $y$. To obtain the energy equation we multiply
\eqref{eq:PDEcont} with $-\psi$ and average over the whole domain. Thus, the work done by the forcing is given by the term:

```math
P = -\,\overline{\psi\,\xi}^{x,y}\com
```

but the above is a stochastic integral and it is meaningless without a rule for computing the stochastic integral.

Numerically, the work done by the forcing can be obtained Stratonovich-wise as:

```math
\begin{align}
P_j = -\,\overline{\frac{\psi(\bx,t_j)+\psi(\bx,t_{j+1})}{2}  \xi(\bx,t_{j+1}) }^{x,y}\com
\end{align}
```

or Itô-wise

```math
\begin{align}
P_j = -\,\overline{ \psi(\bx,t_j) \xi(\bx,t_{j+1}) }^{x,y} + \text{drift}\com
\end{align}
```

But how much is the Itô drift term in this case? As in the previous section, the drift is *precisely* the ensemble
mean of the Stratonovich work, i.e.:

```math
\textrm{Ito drift}= -\overline{ \la\underbrace{\psi(\bx,t)\circ  \xi(\bx,t)}_{\textrm{Stratonovich}} \ra }^{x,y}\com
```

But again the above can be computed relatively easy if we use the "formal" solution of \eqref{eq:PDEcont}:

```math
\psi(\bx,t) = \ee^{-\mu t}\psi(\bx,0) + \int_0^t \ee^{-\mu(t-s)}\nabla^{-2}\xi(\bx,s)\,\dd s\com
```

which implies

```math
\text{drift} = -\overline{\ee^{-\mu t}\underbrace{\laa\psi(\bx,0) \xi(\bx,t)\raa}_{=0} }^{x,y} - \int_0^t \ee^{-\mu(t-s)}\overline{\nabla^{-2}\laa \xi(\bx,s)\xi(\bx,t)\raa}^{x,y}\,\dd s \\
= -\int_0^t \ee^{-\mu(t-s)}\overline{\underbrace{\left[\nabla^{-2} Q(\bx)\right]\big|_{\bx=0}}_{\text{independent of }x,y}\,\delta(t-s)}^{x,y}\,\dd s = -\frac1{2} \nabla^{-2} Q(\bx)\big|_{\bx=0} \\
= -\frac1{2} \left[\nabla^{-2} \int \frac{\dd^2\bk}{(2\pi)^2} \widehat{Q}(\bk)\,\ee^{\ii\bk\bcdot\bx} \right] _{\bx=0}
= \int \frac{\dd^2\bk}{(2\pi)^2} \frac{\widehat{Q}(\bk)}{2k^2}\per
```

Thus, the drift, or in this case the mean energy input rate by the stochastic forcing, is precisely determined
from the spatial correlation of the forcing. Let us denote:

```math
\varepsilon \defn \int \frac{\dd^2\bk}{(2\pi)^2} \frac{\widehat{Q}(\bk)}{2k^2}\per\label{eq:def_epsilon}
```

Therefore, work for a single forcing realization is computed numerically as:

```math
\begin{align}
{\color{Green}\text{Itô}}&: {\color{Green} P_j  =  -\overline{ \psi(\bx,t_j) \xi(\bx,t_{j+1}) }^{x,y}  + \varepsilon}\com\\
{\color{Magenta}\text{Stratonovich}} &: {\color{Magenta}P_j = -\overline{\frac{\psi(\bx,t_j)+\psi(\bx,t_{j+1})}{2}  \xi(\bx,t_{j+1}) }^{x,y}}\per \label{eq:PtStrat}
\end{align}
```

Remember, previously the work done by the stochastic forcing was:
```math
\dd P_t = {\color{Green}\frac{\sigma}{2}\dd t + \sqrt{\sigma}x_t\dd W_t} = {\color{Magenta}\sqrt{\sigma} x_t\circ\dd W_t}\com
```
and by sampling over various forcing realizations:
```math
\langle \dd P_t\rangle = \frac{\sigma}{2}\dd t = \langle\sqrt{\sigma} x_t\circ\dd W_t\rangle
```

The code uses Stratonovich. For example, the work done by the forcing in the `TwoDTurb` module is computed based
on \eqref{eq:PtStrat} with the function

```julia
@inline function work(s, v::ForcedVars, g)
  @. v.Uh = g.invKKrsq * (v.prevsol + s.sol)/2.0 * conj(v.Fh)
  1/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end
```

## A less-simple SPDE

It turns out that nothing changes if we include the nonlinear terms in the vorticity equation:

```math
\begin{equation}
\partial_t \nabla^2\psi(\bx, t) + \J(\psi,\nabla^2\psi) =  -\mu \nabla^2\psi(\bx, t) + \xi(\bx,t)\per\per\label{eq:PDEcont2}
\end{equation}
```

The nonlinearity does not alter the Itô drift; thus the ensemble mean energy input by the stochastic forcing,
remains the same. We can easily verify this from the "formal" solution of \eqref{eq:PDEcont2}:

```math
\psi(\bx,t) = \ee^{-\mu t}\psi(\bx,0) + \int_0^t \ee^{-\mu(t-s)}\nabla^{-2}\xi(\bx,s)\,\dd s - \int_0^t \nabla^{-2}\J\left(\psi(\bx,s),\nabla^2\psi(\bx,s)\right)\,\dd s\com
```

When multiplied with $\xi(\bx,t)$ the last term vanishes since its only non-zero contribution comes from the point
$s=t$ which is of measure zero (in the integrated sense).

Figure below shows the energy budgets for a numerical solution of \eqref{eq:PDEcont2}  starting from rest
($\psi(\bx,0)=0$) in a doubly periodic square domain of size $L$ (`examples/twodturb/IsotropicRingForcing.jl`).
The forcing was prescribed to have power in a narrow ring in wavenumber space:

```math
\widehat{Q}(\bk)\propto \ee^{-(|\bk|-k_f)^2/(2\delta_f^2)}\com
```

with $k_f L/(2\pi) = 12$ and $\delta_f L/(2\pi) = 2$. The mean energy input rate was set to $\varepsilon = 0.1$.

![energy_budgets_SPDE_Stratonovich](assets/energy_budgets_SPDE_Stratonovich.png)
