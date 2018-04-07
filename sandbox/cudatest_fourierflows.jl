using CuArrays, BenchmarkTools
using FourierFlows
import FourierFlows.TwoDTurb

 L = 2π   # Domain
nu = 1e-6
dt = 1.0
stepper = "FilteredRK4"

function setup(n, withgpu)
  q0 = rand(n, n)
  if withgpu
    prob = TwoDTurb.CuProblem(nx=n, Lx=L, nu=nu, dt=dt, stepper=stepper)
    TwoDTurb.set_q!(prob, CuArray(q0))
  else
    prob = TwoDTurb.Problem(nx=n, Lx=L, nu=nu, dt=dt, stepper=stepper)
    TwoDTurb.set_q!(prob, q0)
  end
  prob
end

nt = 100
for n = [128, 256, 512, 1024, 2048]
  gpuprob = setup(n, true)
  cpuprob = setup(n, false)

  stepforward!(gpuprob, 1)
  stepforward!(cpuprob, 1)

  tgpu = @belapsed stepforward!($gpuprob, nt)
  tcpu = @belapsed stepforward!($cpuprob, nt)

  @printf "\n%d steps on the GPU at %d^2: %.6f s\n" nt n tgpu
  @printf "%d steps on the CPU at %d^2: %.6f s\n" nt n tcpu
  @printf "speedup: %.2f\n\n" tcpu/tgpu
end
