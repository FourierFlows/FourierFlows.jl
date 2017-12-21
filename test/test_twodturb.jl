include("../src/FourierFlows.jl")

using Base.Test
using FourierFlows
import FourierFlows.TwoDTurb


function makebasicturbproblem(n, L, ν, nν)
  g  = TwoDTurb.TwoDGrid(n, L)
  p  = TwoDTurb.Params(ν, nν)
  v  = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(p, g)

  g, p, v, eq
end

function teststepforward(g, p, v, eq; dt=1e-16, nsteps=10,
  stepper="ForwardEuler")

  if stepper == "ForwardEuler"; ts = ForwardEulerTimeStepper(dt, eq.LC)
  elseif stepper == "AB3";      ts = AB3TimeStepper(dt, eq.LC)
  elseif stepper == "RK4";      ts = RK4TimeStepper(dt, eq.LC)
  elseif stepper == "ETDRK4";   ts = ETDRK4TimeStepper(dt, eq.LC)
  end

  prob = Problem(g, v, p, eq, ts)
  TwoDTurb.set_q!(prob, rand(g.nx, g.ny))

  absq₀ = sum(abs.(prob.vars.sol))
  stepforward!(prob; nsteps=nsteps)
  absq₁ = sum(abs.(prob.vars.sol))

  isapprox(absq₀, absq₁, atol=nsteps*g.nx*g.ny*1e-15)
end

function teststepforward(n::Int, L, ν, nν::Int; stepper="ForwardEuler")
  g, p, v, eq = makebasicturbproblem(n, L, ν, nν)
  teststepforward(g, p, v, eq; stepper=stepper)
end

@testset "TwoDTurb and Timestepper Tests" begin
  @test teststepforward(128, 2π, 1e-2, 2; stepper="ForwardEuler")
  @test teststepforward(128, 2π, 1e-2, 2; stepper="AB3")
  @test teststepforward(128, 2π, 1e-2, 2; stepper="RK4")
  @test teststepforward(128, 2π, 1e-2, 2; stepper="ETDRK4")
end
