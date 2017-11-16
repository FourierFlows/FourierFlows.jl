using PyPlot, JLD2

@load "./data/manyeps_nx1024.jld2" epz diffratio

epztheory = 0.0:0.01:0.9
diffratiotheory = 1 + 0.5*epztheory.^2
diffratiocomp = zeros(epztheory)

dx = 2*pi/1000
x = dx:dx:2*pi
for i = 1:length(epztheory)
  ep = epztheory[i]
  diffratiocomp[i] = dx*sum(1./(1+ep*sin.(x)))/(2*pi)
end

close("all")

fig, axs = subplots()
plot(epz, diffratio, "o", label="Simulation")
plot(epztheory, diffratiotheory, linestyle="--", 
  label=L"1 + \frac{1}{2}\epsilon^2")

plot(epztheory, diffratiocomp, linestyle="-", 
  label=L"\frac{1}{2 \pi} \int_0^{2\pi} \frac{1}{1-\epsilon \sin(x)} \, \mathrm{d} x")

leg = legend(markerscale=0.9)

xlabel(L"\epsilon")
ylabel(L"r = \kappa_e / \bar \kappa")
title("Enhancement of tracer dispersion in barotropic squeezing flow")

savefig("firsttry.png", dpi=240)

