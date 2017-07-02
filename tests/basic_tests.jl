include("../src/domain.jl")

using Domain
using PyPlot


nx = 8;
Lx = 2.0*pi;

g = Grid(nx, Lx);

if maximum(abs.(g.X))<1e-10
  info("X grid is ALL zeros!")
else
  println("X grid seems OK.")
end


if maximum(abs.(g.Y))<1e-10
  info("Y grid is ALL zeros!")
else
  println("Y grid seems OK.")
end
