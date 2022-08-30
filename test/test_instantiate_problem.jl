function instantiate_problem(dev, stepper)
    problem = Problem(dev; nx=4, stepper)
    return typeof(problem) <: FourierFlows.Problem
end

function instantiate_problem_with_filter_kwargs(dev, stepper)
    dummy_problem = Problem(dev; nx=16, stepper=stepper)

    real_problem = FourierFlows.Problem(dummy_problem.eqn,
                                        stepper,
                                        1.0,
                                        dummy_problem.grid,
                                        dummy_problem.vars,
                                        dummy_problem.params,
                                        innerK = 0.0,
                                        outerK = 1/16)

    stepper = real_problem.timestepper
    
    return CUDA.@allowscalar stepper.filter[3] < 1e-16
end
