include("./rungaussianwaves.jl")

nx    = 256      # Resolution
ep    = 1e-2     # Wave nonlinearity
Ro    = 2e-1     # Eddy Rossby number

for nkw in [1, 4, 16]
  plotname = @sprintf("wave%02dn", nkw)
  runsimulation(-Ro, ep, nkw, nx, plotname)
end
