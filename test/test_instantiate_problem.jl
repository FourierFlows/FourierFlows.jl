function instantiate_problem(dev, stepper)
    problem = Problem(nx=4, dev=dev, stepper=stepper)
    return typeof(problem) <: FourierFlows.Problem
end

function instantiate_problem_with_filter_kwargs(dev, stepper)
    dummy_problem = Problem(nx=16, dev=dev, stepper=stepper)

    real_problem = FourierFlows.Problem(dummy_problem.eqn,
                                        stepper,
                                        1.0,
                                        dummy_problem.grid,
                                        dummy_problem.vars,
                                        dummy_problem.params,
                                        dev,
                                        innerK=0.0,
                                        outerK=1/16)

    stepper = real_problem.timestepper
    
    return CUDA.@allowscalar stepper.filter[3] < 1e-16
end
