# [Time-stepping](@id timestepping)

FourierFlows.jl includes several time-stepping algorithms.

Most of the time-stepping algorithms are fully explicit schemes: [`ForwardEulerTimeStepper`](@ref), [`AB3TimeStepper`](@ref), [`RK4TimeStepper`](@ref), and [`LSRK54TimeStepper`](@ref)
but also implemented is the [`ETDRK4TimeStepper`](@ref).

## High-wavenumber filtering

Most of the time steppers also come with their `Filtered` equivalents: [`FilteredForwardEulerTimeStepper`](@ref), [`FilteredAB3TimeStepper`](@ref), [`FilteredRK4TimeStepper`](@ref), and [`FilteredLSRK54TimeStepper`](@ref), and [`FilteredETDRK4TimeStepper`](@ref).

The filtered time steppers include a high-wavenumber filter that is applied to the solution state vector after each time step. The motivation behind filtering is to 

This filter is

```math
\mathrm{filter}(\boldsymbol{k}) = 
     \begin{cases}
       1 & \quad |\boldsymbol{k}| ≤ k_{\textrm{cutoff}} \, ,\\ 
       \exp{ \left [- \alpha (|\boldsymbol{k}| - k_{\textrm{cutoff}})^p \right]} & \quad |\boldsymbol{k}| > k_{\textrm{cutoff}} \, .
     \end{cases}
```

For fluid equations with quadratic non-linearities it makes sense to choose a cutoff wavenumber
at 2/3 of the highest wavenumber resolved in our domain, ``k_{\textrm{cutoff}} = \tfrac{2}{3} k_{\textrm{max}}`` (see discussion in [Aliasing section](@ref aliasing)).

Given the order ``p``, we calculate the coefficient ``\alpha`` so that the the filter value
that corresponds to the highest allowed wavenumber in our domain is a small value, ``\delta``,
taken to be close to machine precision.

That is:

```math
\alpha = \frac{- \log\delta}{(k_{\textrm{max}} - k_{\textrm{cutoff}})^p} \ .
```

The above filter originates from Canuto et al. (1988). In geophysical turbulence applications
it was used by LaCasce (1996) and later by Arbic & Flierl (2003).