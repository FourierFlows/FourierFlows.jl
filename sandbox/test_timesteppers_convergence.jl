using PyPlot
using FourierFlows


function test_baroQG_RossbyWave(stepper, dt, nsteps, g, p, v, eq)

    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
    prob = FourierFlows.Problem(g, v, p, eq, ts)

    s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

    # the Rossby wave initial condition
     ampl = 1e-2
    kwave = 3.0*2π/g.Lx
    lwave = 2.0*2π/g.Ly
        ω = -p.beta*kwave/(kwave^2.0 + lwave^2.0)
       ζ0 = ampl*cos.(kwave*g.X).*cos.(lwave*g.Y)
      ζ0h = rfft(ζ0)

    BarotropicQG.set_zeta!(prob, ζ0)

    stepforward!(prob, nsteps)
    BarotropicQG.updatevars!(prob)

    ζ_theory = ampl*cos.(kwave*(g.X-ω/kwave*s.t)).*cos.(lwave*g.Y)
    # println(norm(ζ_theory - v.zeta)/norm(ζ_theory))
    # isapprox(ζ_theory, v.zeta, rtol=g.nx*g.ny*nsteps*1e-12)
    residual = norm(ζ_theory - v.zeta)/norm(ζ_theory)
    residual
end



# -----------------------------------------------------------------------------
# Running the tests

nx  = 64
ν  = 0.0
νn = 2

f0 = 1.0
 β = 2.0
Lx = 2π
 μ = 0.0

η(x,y) = zeros(nx, nx)
FU(t)  = 0

g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)

dt0 = 0.5
N = round.(logspace(1,4,21))

residual1 = zeros(1, length(N))
residual2 = zeros(1, length(N))
residual3 = zeros(1, length(N))
for i=1:length(N)
  println(i)
  nsteps = N[i]
  dt = dt0/N[i]
  residual1[i]=test_baroQG_RossbyWave("FilteredForwardEuler", dt, nsteps, g, p, v, eq)
  residual2[i]=test_baroQG_RossbyWave("FilteredAB3", dt, nsteps, g, p, v, eq)
  residual3[i]=test_baroQG_RossbyWave("FilteredRK4", dt, nsteps, g, p, v, eq)
end

figure()
loglog( N, residual1', "bo-", linewidth=2, label="Euler")
loglog( N, residual2', "ro-", linewidth=2, label="AB3")
loglog( N, residual3', "go-", linewidth=2, label="RK4")
loglog( N, 0.1*N.^(-1)', "b:", linewidth=2, label=L"$N^{-1}$")
loglog( N, 0.1*N.^(-2)', "r:", linewidth=2, label=L"$N^{-2}$")
loglog( N, 0.0001*N.^(-4)', "g:", linewidth=2, label=L"$N^{-4}$")
legend(fontsize=10)
xlabel("number of points")
ylabel("relative error")
