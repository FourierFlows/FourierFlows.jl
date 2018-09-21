using PyPlot

# sample problem
  n = [64, 128, 256, 512, 1024, 2048, 4096];
cpu = [63.376e-3, 103.985e-3, 292.813e-3, 872.886e-3, 3.957, 22.645, 79.344 ];
gpu = [18.923e-3, 18.869e-3, 24.897e-3, 71.766e-3, 238.528e-3, 925.671e-3,4.984   ];

# FourierFlows TwoDTurb problem
  n = [128, 256, 512, 1024, 2048, 4096]
cpu = [0.446, 1.130, 3.344, 19.654, 98.487, 376.201]
gpu = [0.104, 0.152, 0.414, 1.405, 5.459, 23.724]


fig, axs = subplots(ncols=1, nrows=2, figsize=(6, 7))

subplot(2,1,1)
loglog(n, cpu, "*:", label="CPU", linewidth=0.8,  markersize=8)
loglog(n, gpu, "o:", label="GPU", linewidth=0.8,  markersize=6)
loglog(n, 0.0001*n.^2, "-k", label="slope 2", linewidth=1)
ylabel("walltime for 100 time-steps [s]", fontsize=12)
title("TwoDTurb module benchmark", fontsize=16)
# xlim(70, 10000)
legend()

subplot(2,1,2)
semilogx(n, cpu./gpu, "d:k", linewidth=0.8, markersize=6)
ylabel("GPU over CPU speedup", fontsize=12)
xlabel("number of grid-points per spatial dimension", fontsize=12)
# xlim(20, 10000)
ylim(0, 20)
legend()

savename = string("CPUversusGPU")
savefig(savename, dpi=300)
