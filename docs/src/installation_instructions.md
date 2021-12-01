# Installation instructions

You can install the latest version of FourierFlows.jl via the built-in package manager 
(accessed by pressing `]` in the Julia REPL command prompt) to add the package and also to 
instantiate/build all the required dependencies

```julia
julia>]
(v1.6) pkg> add FourierFlows
(v1.6) pkg> instantiate
```

We recommend installing FourierFlows.jl with the built-in Julia package manager, because 
this installs a stable, tagged release. Later on, you can update FourierFlows.jl to the 
latest tagged release again via the package manager by typing

```julia
(v1.6) pkg> update FourierFlows
```

Note that some releases might induce breaking changes to certain modules. If after anything 
happens or your code stops working, please open an [issue](https://github.com/FourierFlows/FourierFlows.jl/issues) 
or start a [discussion](https://github.com/FourierFlows/FourierFlows.jl/discussions). We're 
more than happy to help with getting your simulations up and running.

!!! warn "Use Julia 1.6 or newer"
    The latest FourierFlows.jl requires at least Julia v1.6 (the current long-term-release).
    Installing FourierFlows with an older version of Julia will install an older version 
    of FourierFlows.jl (the latest version compatible with your version of Julia).

    Last version that works with Julia v1.5 is FourierFlows.jl v0.7.2.
