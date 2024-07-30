# Installation instructions

To install the latest version of FourierFlows.jl use the Julia's built-in package manager
(accessed by pressing `]` in the Julia REPL command prompt):

```julia
julia> ]
(v1.10) pkg> add FourierFlows
```

We recommend installing FourierFlows.jl with the built-in Julia package manager, because
this installs a stable, tagged release. Later, you can update FourierFlows.jl to the
latest tagged release again via the package manager by

```julia
(v1.10) pkg> update FourierFlows
```

Note that some releases might induce breaking changes to certain modules. If after anything
happens or your code stops working, please open an [issue](https://github.com/FourierFlows/FourierFlows.jl/issues) 
or start a [discussion](https://github.com/FourierFlows/FourierFlows.jl/discussions). We're
more than happy to help with getting your simulations up and running.

!!! warn "Use Julia 1.10 or newer"
    The latest FourierFlows.jl requires at least Julia v1.10 (the current long-term-release).
    Installing FourierFlows with an older version of Julia will install an older version
    of FourierFlows.jl (the latest version compatible with your version of Julia).
