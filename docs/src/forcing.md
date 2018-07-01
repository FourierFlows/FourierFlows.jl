# Forcing

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
\newcommand{\zhat}{\hat{\mathbf{z}}}

\newcommand{\bit}{\vphantom{\dot{W}}}
\newcommand{\sd}{b}
```



## Stochastic Differential Equations (SDEs)

A differential equation in the form:
\begin{equation}
	\frac{\dd x}{\dd t} = f(x)\com
\end{equation}
can also be written in an integral form:
\begin{equation}
	x(t) = \int_{t_0}^{t} f(x(s))\,\dd s\per
\end{equation}
In a similar manner, a stochastic differential equation
\begin{equation}
	\dd x = f(x)\,\dd t + g(x)\,\dd W_t\com
\end{equation}
with $\dd W_t$ a white-noise process, can be written in an integral form as:
\begin{equation}
	x(t) = \int_{t_0}^{t} f(x(s))\,\dd s + \int_{t_0}^{t} g(x(s))\,\dd W_s \per
\end{equation}
Of course now, the last integral is a stochastic integral and there is not a single straight-forward way of computing it --- there are a lot of different ways we can approximate it as a Riemannian sum and each of them leads to a different answer. The two most popular ways for computing such stochastic integrals are:

```math
{\color{green}\text{Itô}: \int_{t_0}^{t} g(x(s))\,\dd W_s\approx\sum_{j} g\left(x(t_j)\right)(W_{j+1}-W_j)}\com\\
{\color{Blue}\text{Stranotovich}: \int_{t_0}^{t} g(x(s))\,\dd W_s \approx \sum_{j} g\left(x\left(\half(t_j+t_{j+1})\right)\right)(W_{j+1}-W_j)}\per
```

Although that the two definitions would converge to the very same result if the $x(t)$ was continuous, in SDE's that is not the case; the two different definitions give thoroughly different results. And to overcome that they come along with different chain rules, i.e. chain rules that are not necessarily the same as those in plain old calculus.

An SDE can be written also in differential form. Because we cannot formally form $\dd W/\dd t$, since $W$ is nowhere differentiable, we write an SDE in differential form as:

```math
{\color{green}\text{Itô}: \dd x_t = f(x_t)\dd t + g(x_t)\dd W_t}\com\\
{\color{Blue}\text{Stranotovich}: \dd x_t = f(x_t)\dd t + g(x_t)\circ\dd W_t\per}
```

The circle in $g(x_t)\circ\dd W_t$ is used to differentiate between Itô or Stratonovich calculus.

A variable change $y=G(x)$ is done as follows according to the two different calculi:

```math
{\color{green}\text{Itô}: \dd y_t = \frac{\dd G}{\dd x}\dd x_t + \half g(x_t)^2 \frac{\dd^2 G}{\dd x^2}\dd t =\left[ \frac{\dd G}{\dd x}f(x_t) + \half g(x_t)^2 \frac{\dd^2 G}{\dd x^2}\right]\dd t + \frac{\dd G}{\dd x}g(x_t)\dd W_t}\com\\
{\color{Blue}\text{Stranotovich}: \dd y_t  = \frac{\dd G}{\dd x}\dd x_t =\frac{\dd G}{\dd x} f(x_t) \dd t + \frac{\dd G}{\dd x}g(x_t)\dd W_t}\per
```

The above are the so called \textit{stochastic chain rules}. All derivatives of $G$ are evaluated at $x_t$.

It's easy to see that the extra drift-term in Itô's interpretation of the stochastic integral, i.e., ${\color{green}\half g^2\, \dd^2G/\dd x^2}$  is _exactly_ equal to the ensemble mean of the Stratonovich stochastic integral. This is that case because, by construction, the Itô stochastic integral has zero ensemble mean since at every instant the noise is multiplied with $g$ evaluated before the action of the noise occurs; $g$ and $\dd W$ are uncorrelated and thus:
\begin{equation}
{\color{green}\laa g(x_t)\dd W_t \raa =0}\quad\text{while}\quad {\color{Blue}\laa g(x_t)\circ\dd W_t \raa \ne 0}\per
\end{equation}
The above is demonstrated by evaluating the simple stochastic integral:

```math
{\color{green}\text{Itô}: \laa \int_{t_0}^{t} W_s\,\dd W_s \raa \approx\sum_{j} \laa W_j(W_{j+1}-W_j)\raa}\\
{\color{green}\hspace{7.3em} = \sum_{j} \laa W_j W_{j+1}\raa - \laa W_jW_j\raa \sim \sum_{j} t_j - t_j = 0} \com\\
{\color{Blue}\text{Stranotovich}: \laa\int_{t_0}^{t} W_s\circ\dd W_s\raa \approx \sum_{j} \laa \frac1{2}(W_{j} + W_{j+1}) (W_{j+1}-W_j)\raa }\\
{\color{Blue}\hspace{7.3em} = \frac1{2}\sum_{j} \laa W_{j+1} W_{j+1}\raa - \laa W_{j} W_{j}\raa  \sim \frac1{2}\sum_{j} t_{j+1} - t_j = \frac{t}{2}}\per
```

SDEs rarely can be solved in closed form; most often numerical solution of SDEs is brought to the rescue. Itô calculus has the advantage that is very easily implemented numerically. On the other hand, Stratonovich calculus coincides with that from normal calculus and this stems from the fact that it vies the white noise process as a series of colored noise processes with the de-correlation time tending to zero. This last fact is what made Stratonovich calculus more popular in the physics community. A nice discussion on the differences and similarities between the two calculi is done by [van Kampen](https://doi.org/10.1007/BF01007642).

## A simple SDE: the Ornstein--Uhlenbeck process

One of the simpler SDEs is the Ornstein--Uhlenbeck process. A variation of which is:
\begin{equation}
	x(t) = \int_{t_0}^{t} -\mu x(s)\,\dd s + \int_{t_0}^{t} \sqrt{\sigma}\,\dd W_s \per\label{eq:OU}
\end{equation}
Note that in differential form this is:
\begin{equation}
	\dd x_t = -\mu x_t \,\dd t + \sqrt{\sigma}\,\dd W_s \per\label{eq:1}
\end{equation}
Luckily, here there is no need to distinguish between Itô and Stratonovich. This is because $g$ is independent of $x(t)$. But we stress that this is only a fortuitous However, this is often not the case.

How do we time-step this SDE numerically? Let us assume a discretization of time into time-steps of $\tau$: $t_j=(j-1)\tau$. (What follows can be easily transfer to non-uniform time discretization.) With that, we denote $x_j\defn x(t_j)$. Then the Euler--Mayorama time-step scheme for \eqref{eq:1} is
\begin{equation}
	x_{j+1} = x_j + (-\mu x_j)\tau + \sqrt{\sigma}(W_{j+1}-W_j)\per
\end{equation}

Now let us ask the following question: How can we compute the work done by the noise? In other words, if we are interested in the evolution of the "energy" $E\defn \half x^2$, how is the noise term attributing in the growth of $E$? To answer that we first have to find the SDE that energy $E$ obeys. But, in doing so, it is important to adopt a single interpretation for computing stochastic integrals as now a transformation of variables is needed. That is, depending on whether we choose to interpret the stochastic integrals according to Itô or to Stratonovich calculus, $E$ evolves as:



```math
{\color{green}\text{Itô}:  \dd E_t  = \left( -2\mu E_t + \frac{1}{2} \sigma \right)\dd t  + x_t \sqrt{\sigma}\dd W_t}\com\\
{\color{Blue}\text{Stratonovich}: \dd E_t  = -2\mu E_t  \dd t + x_t\circ \sqrt{\sigma}\dd W_t}\per
```

How do we compute the work $P$ done by the noise? Thus the work done by the stochastic forcing is:

```math
{\color{green}\text{Itô}: P_t = \half \sigma \dd t + \sqrt{\sigma} x_t \dd W_t \approx  \half \sigma + \sqrt{\sigma} x_j (W_{j+1}-W_j)\com}\\
{\color{Blue}\text{Stratonovich}: P_t =  x_t \circ\sqrt{\sigma} \dd W_t \approx \sqrt{\sigma} x\left(\half(t_j+t_{j+1})\right)(W_{j+1}-W_j)}\per
```

Say we didn't know the rules for transforming Stratonovich to Itô and we were wondering what is the extra drift term we have to include in the Itô formulations, i.e. the $\half\sigma$ term. We can compute the Itô's drift-term using that it is exactly equal to $\la x_t\circ\sqrt{\sigma}\dd W_t\ra$; and for the latter we can use the "usual" calculus. That is, rewrite \eqref{eq:OU} as:

\begin{equation}
\dot{x} = -\mu x + \xi\com\label{eq:OUcont}
\end{equation}

where $\xi(t)$ is understood to be the "continuous" version of the white-noise process which is formally only understood in terms of distributions. The forcing $\xi$ has the properties:

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

Above we used that $\int_0^t\delta(t-s)\dd s = \half$, which is consistent with Stratonovich symmetric interpretation of stochastic integrals.
