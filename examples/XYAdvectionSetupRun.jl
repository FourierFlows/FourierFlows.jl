using FourierFlows, JLD2
using FourierFlows:XYAdvection as xy
using LinearAlgebra: mul!, ldiv!
using Plots

dev = CPU()
Nx = 80        # grid resolution
Lx = 80        # box size
dt = 0.02      # timestep (s)
nsteps = 100   # total number of time-steps
printFreq = 1  # data print frequency

stepper = "ForwardEuler"  # timestepper_Scheme


J = 1.0  # diffusion coefficient
A = 0.0  # advection coefficient
Γ = 1.0  # GL coefficient 
α = 1.0 # average order parameter


grid   = xy.set_2Dgrid(Nx,Lx; dev,aliased_fraction = 0)
params = xy.set_Params(;J=J, A=A, Γ=Γ, α=α)
vars   = xy.set_Vars(grid)
prob   = xy.set_Problem(grid,params,vars,dt,stepper)


px = rand(Nx,Nx)
py = rand(Nx,Nx)
xy.set_initialConditions!(prob,px,py)
heatmap(grid.x,grid.y,atan.(py,px)',clims=(-π,π),c=:hsv,axis=true,grid=false,colorbar=true,ratio=1,size=(400,400))


out = xy.run_and_save(prob,nsteps,printFreq)
path=out.path

# plot order parameter time series
orderP = Float64[]
for titr in 0:printFreq:nsteps
   px,py = xy.get_data(titr,path,Nx)
   mag = ((sum(px)^2 + sum(py)^2)^0.5)/(Nx*Nx)
   push!(orderP,mag)
end
plot(orderP,ylims=(0,1),ylabel="m", xlabel="t", size=(400,400))

titr = nsteps
px,py = xy.get_data(titr,path,Nx)
heatmap(grid.x,grid.y,atan.(py,px)',clims=(-π,π),c=:hsv,axis=true,grid=false,colorbar=true,ratio=1,size=(400,400))