include("../src/fourierflows.jl")

import FourierFlows
import FourierFlows.TwoModeBoussinesq

function test_ivp()
  prob = TwoModeBoussinesq.InitialValueProblem()
end

function test_stepforward()
  prob = TwoModeBoussinesq.InitialValueProblem()
  FourierFlows.stepforward!(prob, nsteps=3)
end

# Run test
@printf("Test ivp...")
test_ivp()
@printf(" ok.\n")

@printf("Test stepforward...")
test_stepforward()
@printf(" ok.\n")
