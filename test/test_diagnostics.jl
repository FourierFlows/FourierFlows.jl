import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, enstrophy

nx = 64
beta = 2.0
Lx = 2π
mu = 0.02
nu = 0.0

prob = BarotropicQG.InitialValueProblem(nx=nx, Lx=Lx, beta=beta, mu=mu, nu=nu, stepper="FilteredETDRK4", dt=0.01)
s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

# the Rossby wave initial condition
ampl = 1e-2
kwave = 3.0*2π/g.Lx
lwave = 2.0*2π/g.Ly
  ω = -p.beta*kwave/(kwave^2.0 + lwave^2.0)
 ζ0 = ampl*cos.(kwave*g.X).*cos.(lwave*g.Y)
ζ0h = rfft(ζ0)

BarotropicQG.set_zeta!(prob, ζ0)

nsteps = 200
extrasteps = 20

freqE = 2
E = Diagnostic(energy, prob; nsteps=nsteps, freq=freqE)

freqZ = 1
Z = Diagnostic(enstrophy, prob; nsteps=nsteps, freq=freqZ)

diags = [E, Z]


nstepstot = nsteps + extrasteps
while prob.step < nstepstot
  stepforward!(prob, diags, 1)
end

(
 isapprox(length(E.data), Int(round(nstepstot/freqE+1))) &&
 isapprox(length(Z.data), Int(round(nstepstot/freqZ+1))) &&
 isapprox(norm(E.data), norm(E.data[1]*exp.(-2*mu*E.time)), rtol=1e-13) &&
 isapprox(norm(Z.data), norm(Z.data[1]*exp.(-2*mu*Z.time)), rtol=1e-13)
)
