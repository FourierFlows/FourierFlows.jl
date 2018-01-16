using PyPlot, FourierFlows, Interpolations
import FourierFlows.TwoDTurb

 n = 256
 L = 2π
 ν = 1e-6  # Laplacian viscosity
nν = 1
dt = 1e0   # Time step
nt = 1000  # Number of time steps

prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt,
  stepper="FilteredRK4")
TwoDTurb.set_q!(prob, rand(n, n))

x, y = prob.grid.X, prob.grid.Y
dk, dl = (prob.grid.kr[2]-prob.grid.kr[1]), (prob.grid.l[2]-prob.grid.l[1]) 
KK = sqrt.(prob.grid.KKrsq)

# pyqg method
function radialinterp_pyqg(qh)
  kmax = maximum([maximum(prob.grid.kr), maximum(prob.grid.l)])

  dkr = sqrt(dk^2 + dl^2)
  nkr = round(Int, kmax/dkr) + 1
  kr = range(0.5*dkr, dkr, nkr)
  qhr = zeros(Complex{Float64}, nkr)

  for i = 1:nkr
    ir = (kr[i]-0.5*dkr) .<= KK .<= (kr[i]+0.5*dkr)
    dth = π / (sum(ir)-1)
    qhr[i] = sum(qh[ir]) * kr[i] * dth
  end

  kr, qhr
end

# Step forward
close("all")
fig, axs = subplots(ncols=2, figsize=(8, 4))
tic()
for i = 1:10
  stepforward!(prob, nt)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.grid.dx/prob.ts.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  krp, qhrp = radialinterp_pyqg(prob.vars.qh)
  kri, qhri = FourierFlows.radialspectrum(prob.vars.qh, prob.grid, 
    refinement=1)

  krc = 0.25*(kri[1:end-3]+kri[2:end-2]+kri[3:end-1]+kri[4:end])
  qhrc = 0.25*(qhri[1:end-3] + qhri[2:end-2] + qhri[3:end-1] + qhri[4:end])

  sca(axs[1]); cla()
  pcolormesh(x, y, prob.vars.q)

  sca(axs[2]); cla()
  plot(krp, abs.(qhrp), "-")
  plot(kri, abs.(qhri), "--")
  #plot(krc, abs.(qhrc), "-")
  axs[2][:set_yscale]("log")
  xlim(0, 50)
  ylim(ymin=1, ymax=1000)


  show()
  pause(1)
  readline(STDIN)

end

show()
