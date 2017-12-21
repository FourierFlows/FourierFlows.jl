fourierflowspath = ".."
include(joinpath(fourierflowspath, "src", "FourierFlows.jl"))

using FourierFlows, PyPlot
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

ns = 100
nx = 256

for nx in [64, 128, 256, 512, 1024]
  prob = TwoDTurb.InitialValueProblem(withfilter=true, nx=nx, Lx=2π, dt=2e-1)
  TwoDTurb.set_q!(prob, rand(nx, nx))
  stepforward!(prob; nsteps=1)

  walltime = @elapsed stepforward!(prob; nsteps=ns)
  @printf "nx: %04d, time per step: %.6f\n" nx walltime/ns
end
