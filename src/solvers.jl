# solvers.jl

# This file defines functions that operate on types common to all
# FourierFlows problems: 
#     
#       stepforward! takes a time-stepping method
#         as input and advances the variable Vars.sol nsteps times,
#         where 'nsteps' is an integer.






# ----------------------------------------------------------------------------- 
# Forward Euler
# ----------------------------------------------------------------------------- 
function stepforward!(nsteps::Int, ts::ForwardEulerTimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps
    calcNL!(ts.NL, v.sol, v.t, v, p, g)
    #ts.NL .= calcNL(ts.NL, v.sol, v.t, v, p, g)
    v.sol .+= ts.dt .* (ts.NL .+ ts.LC.*v.sol)

    ts.step += 1
    v.t += ts.dt
  end
end


# ----------------------------------------------------------------------------- 
# ETDRK4
# ----------------------------------------------------------------------------- 
function stepforward!(nsteps::Int, ts::ETDRK4TimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps

    calcNL!(ts.NL1, v.sol, v.t, v, p, g)

    # Substep 1
    ts.ti = v.t + 0.5*ts.dt
    ts.sol1 .= ts.expLCdt2.*v.sol .+ ts.zeta.*ts.NL1
    calcNL!(ts.NL2, ts.sol1, ts.ti, v, p, g)

    # Substep 2
    ts.sol2 .= ts.expLCdt2.*v.sol .+ ts.zeta.*ts.NL2
    calcNL!(ts.NL3, ts.sol2, ts.ti, v, p, g)

    # Substep 3
    ts.ti = v.t + ts.dt
    ts.sol2 .= ts.expLCdt2.*ts.sol1 .+ ts.zeta.*(2.0.*ts.NL3 .- ts.NL1)
    calcNL!(ts.NL4, ts.sol2, ts.ti, v, p, g)

    # Substep 4 and update
    v.sol .= (ts.expLCdt.*v.sol .+      ts.alph .* ts.NL1 
                                .+ 2.0.*ts.beta .* (ts.NL2 .+ ts.NL3)
                                .+      ts.gamm .* ts.NL4 ) 

    ts.step += 1
    v.t += ts.dt

  end
end



# ----------------------------------------------------------------------------- 
# RK4
# ----------------------------------------------------------------------------- 
function stepforward!(nsteps::Int, ts::RK4TimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps
    
    calcNL!(ts.RHS1, v.sol, v.t, v, p, g)
    ts.RHS1 .+= ts.LC.*v.sol

    # Substep 1
    ts.ti = v.t + 0.5*ts.dt
    ts.sol1 .= v.sol .+ (0.5*ts.dt).*ts.RHS1
    calcNL!(ts.RHS2, ts.sol1, v.t, v, p, g)
    ts.RHS2 .+= ts.LC.*ts.sol1

    # Substep 2
    ts.sol1 .= v.sol .+ (0.5*ts.dt).*ts.RHS2
    calcNL!(ts.RHS3, ts.sol1, v.t, v, p, g)
    ts.RHS3 .+= ts.LC.*ts.sol1

    # Substep 3
    ts.ti = v.t + ts.dt
    ts.sol1 .= v.sol .+ ts.dt.*ts.RHS3
    calcNL!(ts.RHS4, ts.sol1, v.t, v, p, g)
    ts.RHS4 .+= ts.LC.*ts.sol1

    # Substep 4 and final step
    v.sol .+= ts.dt.*( 
         (1.0/6.0).*ts.RHS1 .+ (1.0/3.0).*ts.RHS2
      .+ (1.0/3.0).*ts.RHS3 .+ (1.0/6.0).*ts.RHS4 )

  end
end


# ----------------------------------------------------------------------------- 
# AB3
# ----------------------------------------------------------------------------- 
function stepforward!(nsteps::Int, ts::AB3TimeStepper,
  v::Vars, p::Params, g::Grid; initstepper=false)

  # Determine initsteps value
  initstepper ? initsteps=0 : initsteps=ts.step

  while initsteps < 2
    # Take forward Euler steps to initialize AB3
    calcNL!(ts.RHS, v.sol, v.t, v, p, g)
    ts.RHS .+= ts.LC.*v.sol

    # Update 
    ts.RHSm2 .= ts.RHSm1
    ts.RHSm1 .= ts.RHS

    v.sol     .+= ts.dt.*ts.RHS
    v.t        += ts.dt
    ts.step    += 1
    initsteps  += 1

    # Account for steps taken
    nsteps -= 1
  end

  # Loop
  for step = 1:nsteps
    calcNL!(ts.RHS, v.sol, v.t, v, p, g)
    ts.RHS .+= ts.LC.*v.sol
    
    v.sol .+= ts.dt .* ( 
      (23.0/12.0).*ts.RHS .- (16.0/12.0).*ts.RHSm1 .+ (5.0/12.0).*ts.RHSm2 )

    ts.RHSm2 .= ts.RHSm1
    ts.RHSm1 .= ts.RHS
    ts.step += 1
    v.t += ts.dt
  end
end




