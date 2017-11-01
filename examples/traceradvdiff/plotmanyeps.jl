using PyPlot, JLD2

@load "./data/manyeps_nx1024.jld2" epz diffratio

epztheory = 0.0:0.01:0.9
diffratiotheory = 1 + epztheory.^2

close("all")

fig, axs = subplots()
plot(epz, diffratio, "o", label="Simulation")
plot(epztheory, diffratiotheory, linestyle="--", label=L"1 + \epsilon^2")

leg = legend(markerscale=0.9)

xlabel(L"\epsilon")
ylabel(L"r = \kappa_e / \bar \kappa")
title("Enhancement of tracer dispersion in barotropic squeezing flow")

savefig("firsttry.png", dpi=240)

