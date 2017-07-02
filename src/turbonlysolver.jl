# -----------------------------------------------------------------------------
# Turb only solver -------------------------------------------------------------
# -----------------------------------------------------------------------------
module TurbOnlySolver

using Domain, Framework, TimeSteppers

export updatevars!, stepforward!, set_q!






# -----------------------------------------------------------------------------
# Solver ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
function calc_nl!(NLqh::Array{Complex128, 2}, qh::Array{Complex128, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  for j in 1:g.nl
    @simd for i in 1:g.nkr
      @inbounds v.Uh[i, j] =  im*g.Lr[i, j]*qh[i, j]*g.invKKrsq[i, j]
      @inbounds v.Vh[i, j] = -im*g.Kr[i, j]*qh[i, j]*g.invKKrsq[i, j]
    end
  end

  # Transforms
  A_mul_B!( v.q, g.irfftplan, qh   )
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  for j in 1:g.ny
    @simd for i in 1:g.nx
      @inbounds v.Uq[i, j] = v.U[i, j]*v.q[i, j]
      @inbounds v.Vq[i, j] = v.V[i, j]*v.q[i, j]
    end
  end

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  for j in 1:g.nl
    @simd for i in 1:g.nkr
      @inbounds NLqh[i, j] = -im*g.Kr[i, j]*v.Uqh[i, j] - im*g.Lr[i, j]*v.Vqh[i, j]
    end
  end

  dealias!(NLqh, g)

end



function stepforward!(nsteps::Int, qts::ForwardEulerTimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps

    calc_nl!(qts.NL, v.qh, v.t, v, p, g)

    # Update
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds v.qh[i, j] += qts.dt*(qts.NL[i, j] + qts.LC[i, j]*v.qh[i, j])
      end
    end

    v.t += qts.dt

  end

end




function stepforward!(nsteps::Int, qts::ETDRK4TimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps

    calc_nl!(qts.NL1, v.qh, v.t, v, p, g)

    # Substep 1
    qts.ti = v.t + 0.5*qts.dt
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol1[i, j] = (qts.expLCdt2[i, j]*v.qh[i, j]
          + qts.zeta[i, j]*qts.NL1[i, j])
      end
    end

    calc_nl!(qts.NL2, qts.sol1, qts.ti, v, p, g)

    # Substep 2
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol2[i, j] = (qts.expLCdt2[i, j]*v.qh[i, j]
          + qts.zeta[i, j]*qts.NL2[i, j])
      end
    end

    calc_nl!(qts.NL3, qts.sol2, qts.ti, v, p, g)

    # Substep 3
    qts.ti = v.t + qts.dt
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol2[i, j] = (qts.expLCdt2[i, j]*qts.sol1[i, j]
          + qts.zeta[i, j]*(2.0*qts.NL3[i, j] - qts.NL1[i, j]))
      end
    end

    # Final eval of nonlinear term
    calc_nl!(qts.NL4, qts.sol2, qts.ti, v, p, g)

    # Update
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds v.qh[i, j] = (qts.expLCdt[i, j]*v.qh[i, j]
          +     qts.alph[i, j] * qts.NL1[i, j]
          + 2.0*qts.beta[i, j] * (qts.NL2[i, j] + qts.NL3[i, j])
          +     qts.gamm[i, j] * qts.NL4[i, j] )
      end
    end

    v.t = v.t + qts.dt

  end

end


# -----------------------------------------------------------------------------
# Helper functions ------------------------------------------------------------
# -----------------------------------------------------------------------------

function updatevars!(v::Vars, p::Params, g::Grid)

  v.q = irfft(v.qh, g.nx)
  v.psih = -v.qh .* g.invKKrsq
  v.U = -irfft(im*g.Lr.*v.psih, g.nx)
  v.V =  irfft(im*g.Kr.*v.psih, g.nx)

end




end
