using PyPlot, Random

import Statistics.mean
import Random: randn

Random.seed!(1234)

stochastic = true

μ = 0.2
σ = 1/5
dt = 0.01
nsteps= 2001
T = 0:dt:(nsteps-1)*dt

# Theoretical results
if stochastic
  nens = 1000
  ΔW = sqrt(σ)*randn(nsteps, nens)/sqrt(dt)
  if μ == 0
       E_theory = 0.5*σ*T
    dEdt_theory = 0.5*σ*ones(size(T))
  else
       E_theory = @. σ/4μ * (1 - exp(-2μ*T))
    dEdt_theory = @. σ/2  * exp(-2μ*T)
  end
else
  nens = 1
  ΔW = sqrt(σ)*ones(nsteps, nens)
  if μ == 0
       E_theory = @. 0.5*σ*T^2
    dEdt_theory = σ*T
  else
       E_theory = @. σ/(2*μ^2) * (1 - exp(-μ*T))^2
    dEdt_theory = @. σ/μ * exp(-μ*T) * (1 - exp(-μ*T))
  end

end

# Numerical calculation
X = zeros(size(ΔW))
E_ito = zeros(size(ΔW))
E_str = zeros(size(ΔW))
E_numerical = zeros(size(ΔW))


for it = 1:nsteps-1 # time step the equation

    @views @.   X[it+1, :] = X[it, :] + (-μ*X[it, :] + ΔW[it, :])*dt

    @views @.   E_ito[it+1, :] = E_ito[it, :] + (-2*μ*E_ito[it, :]
                                    .+ σ/2)*dt + X[it, :]*ΔW[it, :]*dt

    Ebar = E_str[it, :] + (-2*μ*E_str[it, :])*dt + X[it, :].*ΔW[it, :]*dt
    @views @.   E_str[it+1, :] = E_str[it, :] + (-2*μ*(0.5*(E_str[it, :] +
                        Ebar)))*dt + (0.5*(X[it, :]+X[it+1, :]))*ΔW[it, :]*dt
end

# Energy
@views @. E_numerical = 0.5*X^2

# compute dE/dt numerically
dEdt_ito = mean((E_ito[2:nsteps, :] - E_ito[1:nsteps-1, :])/dt, dims=2)
dEdt_str = mean((E_str[2:nsteps, :] - E_str[1:nsteps-1, :])/dt, dims=2)

# compute the work and dissipation
work_ito = mean(ΔW[1:nsteps-1, :].* X[1:nsteps-1, :], dims=2) .+ σ/2
work_str = mean(ΔW[1:nsteps-1, :].*(X[1:nsteps-1, :]+X[2:nsteps, :])/2, dims=2)
diss_ito = 2*μ*(mean(E_ito[1:nsteps-1, :], dims=2))
diss_str = 2*μ*(mean(E_str[1:nsteps-1, :], dims=2))




# Make plots: compare E(t) evolution Ito, Stratonovich, direct 0.5*x^2

fig, axs = subplots(nrows=1, figsize=(8, 4))

sca(axs); cla()
plot(μ*T, E_numerical[:, 1],
                linewidth=3, label="\$\\frac{1}{2}x_t^2\$ ")
plot(μ*T, E_ito[:, 1],
                linestyle="--", linewidth=2, label="\$E_t\$ (Ito) ")
plot(μ*T, E_str[:, 1],
                linestyle="-.", linewidth=1, label="\$E_t\$ (Stratonovich) ")
xlabel(L"$\mu t$")
ylabel(L"$E$")
title(L"comparison of $E(t)$ for single realization")
legend()

savefig("energy_comparison.eps", dpi=240)
savefig("energy_comparison.png", dpi=240)


# Make plots: energy budgets for a realization of the Ito integration
fig, axs = subplots(nrows=2, figsize=(8, 6), sharex=true)

sca(axs[1]); cla()
axs[1][:plot](μ*T, E_theory,
                linewidth=3, label="theoretical \$\\langle E\\rangle\$")
axs[1][:plot](μ*T, mean(E_ito, dims=2), linestyle="--",
                linewidth=2, label="\$\\langle E\\rangle\$ from $nens ensemble member(s)")
xlabel(L"$\mu t$")
ylabel(L"$E$")
if stochastic; title(L"Ito: $\mathrm{d}E_t = (-2\mu E_t + \frac{\sigma}{2})\mathrm{d}t + \sqrt{\sigma}X_t\mathrm{d}W$")
else;          title(L"Ito: $\mathrm{d}X/\mathrm{d}t = -\mu X + \sqrt{\sigma}$")
end
legend()

sca(axs[2]); cla()
plot(μ*T[1:nsteps-1], dEdt_ito[1:nsteps-1, 1],
                linestyle="--", linewidth=2,
                label="numerical \$d\\langle E\\rangle/dt\$")
plot(μ*T[1:nsteps-1], work_ito[1:nsteps-1, 1] - diss_ito[1:nsteps-1, 1],
                linestyle="-.", linewidth=1,
                label="\$\\langle\$work \$-\$ dissipation\$\\rangle\$")
plot(μ*T[1:nsteps-1], dEdt_theory[1:nsteps-1],
                linewidth=3, label="theoretical \$d\\langle E\\rangle/dt\$")
xlabel(L"$\mu t$")
ylabel(L"$dE/dt$")

legend()

savefig("energy_budgets_Ito.eps", dpi=240)
savefig("energy_budgets_Ito.png", dpi=240)



# Make plots: energy budgets for a realization of the Stratonovich integration
fig, axs = subplots(nrows=2, figsize=(8, 6), sharex=true)

sca(axs[1]); cla()
plot(μ*T, E_theory,
                linewidth=3, label="theoretical \$\\langle E\\rangle\$")
plot(μ*T, mean(E_str, dims=2),
                linestyle="--", linewidth=2,
                    label="\$\\langle E\\rangle\$ from $nens ensemble member(s)")
xlabel(L"$\mu t$")
ylabel(L"$E$")
if stochastic; title(L"Stratonovich: $\mathrm{d}E_t = -2\mu E_t\mathrm{d}t + \sqrt{\sigma}x_t\circ\mathrm{d}W$")
else;          title(L"Stratonovich: $\mathrm{d}E/\mathrm{d}t = -2\mu E + \sqrt{\sigma} \dot{x}$")
end
legend()

sca(axs[2]); cla()
plot(μ*T[1:nsteps-1], dEdt_str[1:nsteps-1],
                linestyle="--", linewidth=2,
                label="numerical \$d\\langle E\\rangle/dt\$")
plot(μ*T[1:nsteps-1], work_str[1:nsteps-1] - diss_str[1:nsteps-1],
                linestyle="-.", linewidth=1,
                label="\$\\langle\$work \$-\$ dissipation\$\\rangle\$")
plot(μ*T[1:nsteps-1], dEdt_theory[1:nsteps-1],
                linewidth=3, label="theoretical \$d\\langle E\\rangle/dt\$")
xlabel(L"$\mu t$")
ylabel(L"$dE/dt$")

legend()

savefig("energy_budgets_Stratonovich.eps", dpi=240)
savefig("energy_budgets_Stratonovich.png", dpi=240)
