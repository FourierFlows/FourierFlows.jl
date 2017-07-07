include("physics/2DQGturb.jl")
include("timesteppers.jl")

module Solver
using Domain, Framework, TimeSteppers

export Grid, Vars, Params, ETDRK4TimeStepper, ForwardEulerTimeStepper
export updatevars!, calc_nl_2DQG!, stepforward!



# -----------------------------------------------------------------------------
# Solver ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
function calc_nl!(
  NLqh::Array{Complex128, 2},
  # NLAh::Array{Complex128, 2},
  qh::Array{Complex128, 2},
  # Ah::Array{Complex128, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  for j in 1:g.nl
    @simd for i in 1:g.nkr
      @inbounds v.Uh[i, j] =  im*g.Lr[i, j]*qh[i, j]*g.invKKrsq[i, j]
      @inbounds v.Vh[i, j] = -im*g.Kr[i, j]*qh[i, j]*g.invKKrsq[i, j]
    end
  end

    #@simd for i in 1:g.nk
    #  @inbounds v.Axh[i, j]  = im*g.K[i, j] * Ah[i, j]
    #  @inbounds v.Ayh[i, j]  = im*g.L[i, j] * Ah[i, j]
    #  @inbounds v.Axxh[i, j] = -g.K2[i, j]  * Ah[i, j]
    #  @inbounds v.Ayyh[i, j] = -g.L2[i, j]  * Ah[i, j]
    #  @inbounds v.Axyh[i, j] = -g.KL[i, j]  * Ah[i, j]
    #  @inbounds v.EAh[i, j]  = p.E[i, j]    * Ah[i, j]
    #end
  #end

  #v.Axh  = im*g.K .* Ah
  #v.Ayh  = im*g.L .* Ah
  #v.Axxh = - g.K2 .* Ah
  #v.Ayyh = - g.L2 .* Ah
  #v.Axyh = - g.KL .* Ah

  # v.EAh  =    p.E .* Ah

  # Transforms
  A_mul_B!( v.q, g.irfftplan, qh   )
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  #A_mul_B!( v.Ax,  g.ifftplan, v.Axh  )
  #A_mul_B!( v.Ay,  g.ifftplan, v.Ayh  )
  #A_mul_B!( v.Axx, g.ifftplan, v.Axxh )
  #A_mul_B!( v.Ayy, g.ifftplan, v.Ayyh )
  #A_mul_B!( v.Axy, g.ifftplan, v.Axyh )

  # A_mul_B!( v.EA,  g.ifftplan, v.EAh  )

  for j in 1:g.ny
    @simd for i in 1:g.nx
      @inbounds v.Uq[i, j] = v.U[i, j]*v.q[i, j]
      @inbounds v.Vq[i, j] = v.V[i, j]*v.q[i, j]
    end
  end

      #@inbounds v.UEA[i, j] = v.U[i, j]*v.EA[i, j]
      #@inbounds v.VEA[i, j] = v.V[i, j]*v.EA[i, j]

      #@inbounds v.qrefxA[i, j] = v.q[i, j]*(
      #  im*p.sig*v.Ax[i, j] - p.f0*v.Ay[i, j])

      #@inbounds v.qrefyA[i, j] = v.q[i, j]*(
      #  im*p.sig*v.Ay[i, j] + p.f0*v.Ax[i, j])

      #@inbounds v.UVjacxA[i, j] = (
      #    v.V[i, j]*(im*p.sig*v.Axy[i, j] - p.f0*v.Ayy[i, j])
      #  - v.U[i, j]*(im*p.sig*v.Ayy[i, j] + p.f0*v.Axy[i, j]) )

      #@inbounds v.UVjacyA[i, j] = (
      #    v.U[i, j]*(im*p.sig*v.Axy[i, j] + p.f0*v.Axx[i, j])
      #  - v.V[i, j]*(im*p.sig*v.Axx[i, j] - p.f0*v.Axy[i, j]) )
  #  end
  #end

  #v.qrefxA = v.q .* ( im*p.sig*v.Ax - p.f0*v.Ay )
  #v.qrefyA = v.q .* ( im*p.sig*v.Ay + p.f0*v.Ax )

  #v.UVjacxA = (
  #    v.V .* ( im*p.sig*v.Axy - p.f0*v.Ayy )
  #  - v.U .* ( im*p.sig*v.Ayy + p.f0*v.Axy ) )

  #v.UVjacyA = (
  #    v.U .* ( im*p.sig*v.Axy + p.f0*v.Axx )
  #  - v.V .* ( im*p.sig*v.Ayy - p.f0*v.Axy ) )

  #A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  #A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  # A_mul_B!( v.UEAh,     g.fftplan, v.UEA     )
  # A_mul_B!( v.VEAh,     g.fftplan, v.VEA     )

  #A_mul_B!( v.qrefxAh,  g.fftplan, v.qrefxA  )
  #A_mul_B!( v.qrefyAh,  g.fftplan, v.qrefyA  )
  #A_mul_B!( v.UVjacxAh, g.fftplan, v.UVjacxA )
  #A_mul_B!( v.UVjacyAh, g.fftplan, v.UVjacyA )

  for j in 1:g.nl
    @simd for i in 1:g.nkr
      @inbounds NLqh[i, j] = -im*g.Kr[i, j]*v.Uqh[i, j] - im*g.Lr[i, j]*v.Vqh[i, j]
    end
  end

    #@simd for i in 1:g.nk
    #  @inbounds NLAh[i, j] = (
    #      im*p.invEk[i, j]   *     v.UEAh[i, j]
    #    + im*p.invEl[i, j]   *     v.VEAh[i, j]
    #    + im*p.invEf0k[i, j] *  v.qrefxAh[i, j]
    #    + im*p.invEf0l[i, j] *  v.qrefyAh[i, j]
    #    +step p.invE2sigk[i, j]  * v.UVjacxAh[i, j]
    #    + p.invE2sigl[i, j]  * v.UVjacyAh[i, j] )
    #end
  #end

  # NLAh = im*p.invEk.*v.UEAh + im*p.invEl.*v.VEAh

  dealias!(NLqh, g)
  # dealias!(NLAh, g)

end



function calc_nl_2DQG!(
  NLqh::Array{Complex{Float64},2},
  qh::Array{Complex128, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!

  # A_mul_B!( v.q, g.irfftplan, v.qh )
  v.q = irfft(qh, g.nx)

  v.Uh = +im*real.(g.Lr).*qh.*real.(g.invKKrsq)
  v.Vh = -im*real.(g.Kr).*qh.*real.(g.invKKrsq)

  # A_mul_B!( v.U, g.irfftplan, v.Uh )
  # A_mul_B!( v.V, g.irfftplan, v.Vh )
  v.U = irfft(v.Uh, g.nx)
  v.V = irfft(v.Vh, g.nx)

  # v.Uq = v.U .* v.q
  # v.Vq = v.V .* v.q

  # A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  # A_mul_B!( v.Vqh, g.rfftplan, v.Vq )
  v.Uqh = rfft(v.U .* v.q)
  v.Vqh = rfft(v.V .* v.q)

  NLqh = -im.*g.Kr.*v.Uqh -im.*g.Lr.*v.Vqh
end


function stepforward!(nsteps::Int,
  qts::ForwardEulerTimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps
    println(v.qh[10,20])
    calc_nl_2DQG!(qts.NL, v.qh, v.t, v, p, g)
    println("why is qts.NL not updated?")
    println(qts.NL[10,20])

    v.qh = v.qh + qts.dt * (qts.NL + qts.LC.*v.qh)
  end
end

################################
### COMMENT
#
#     for j in 1:g.nl
#       @simd for i in 1:g.nkr
#         @inbounds v.Uh[i, j] =  im*g.Lr[i, j]*qh[i, j]*g.invKKrsq[i, j]
#         @inbounds v.Vh[i, j] = -im*g.Kr[i, j]*qh[i, j]*g.invKKrsq[i, j]
#       end
#     end
#
#     # v.EAh  =    p.E .* Ah
#
#     # Transforms
#     A_mul_B!( v.q, g.irfftplan, qh   )
#     A_mul_B!( v.U, g.irfftplan, v.Uh )
#     A_mul_B!( v.V, g.irfftplan, v.Vh )
#
#     # A_mul_B!( v.EA,  g.ifftplan, v.EAh  )
#
#     for j in 1:g.ny
#       @simd for i in 1:g.nx
#         @inbounds v.Uq[i, j] = v.U[i, j]*v.q[i, j]
#         @inbounds v.Vq[i, j] = v.V[i, j]*v.q[i, j]
#       end
#     end
#
#     # v.UEA = v.U .* v.EA
#     # v.VEA = v.V .* v.EA
#
#     A_mul_B!( v.Uqh,  g.rfftplan, v.Uq  )
#     A_mul_B!( v.Vqh,  g.rfftplan, v.Vq  )
#     # A_mul_B!( v.UEAh, g.fftplan,  v.UEA )
#     # A_mul_B!( v.VEAh, g.fftplan,  v.VEA )
#
#     for j in 1:g.nl
#       @simd for i in 1:g.nkr
#         @inbounds qts.NL[i, j] = -im*g.Kr[i, j]*v.Uqh[i, j] - im*g.Lr[i, j]*v.Vqh[i, j]
#       end
#     end
#
#     # ats.NL = im*p.invEk.*v.UEAh + im*p.invEl.*v.VEAh
#
#     #sum1 = sum(abs.(p.invEk))
#     #sum2 = sum(abs.(p.invEl))
#     #sum3 = sum(abs.(v.UEAh))
#     #sum4 = sum(abs.(v.VEAh))
#     # sum5 = ats.dt*sum(abs.(ats.NL))
#     # sum6 = sum(abs.(v.Ah))
#     #@printf("(%.2e, %.2e, %.2e, %.2e) ", sum1, sum2, sum3, sum4)
#     # @printf("(%.2e, %.2e) ", sum5, sum6)
#
#     for j = 1:g.nl
#       for i = 1:g.nkr
#         v.qh[i, j] = v.qh[i, j] + qts.dt*(qts.NL[i, j] + qts.LC[i, j]*v.qh[i, j])
#       end
#
#       # for i = 1:g.nk
#       #   v.Ah[i, j] = v.Ah[i, j] + ats.dt*ats.NL[i, j]
#       # end
#
#     end
#
#     #for j in 1:g.nl
#     #  @simd for i in 1:g.nkr
#     #    @inbounds v.qh[i, j] += qts.dt*(qts.NL[i, j] + qts.LC[i, j]*v.qh[i, j])
#     #  end
#     #end
#
#     #for j in 1:g.nl
#     #  @simd for i in 1:g.nk
#     #    @inbounds v.Ah[i, j] += ats.dt*(ats.NL[i, j] + ats.LC[i, j]*v.Ah[i, j])
#     #  end
#     #end
#
#   end
#
# end
################################

function stepforward!(nsteps::Int,
  qts::ETDRK4TimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps

    calc_nl!(qts.NL1, ats.NL1, v.qh, v.Ah, v.t, v, p, g)

    # Substep 1
    qts.ti = v.t + 0.5*qts.dt
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol1[i, j] = (qts.expLCdt2[i, j]*v.qh[i, j]
          + qts.zeta[i, j]*qts.NL1[i, j])
      end

      for i in 1:g.nk
        @inbounds ats.sol1[i, j] = (ats.expLCdt2[i, j]*v.Ah[i, j]
          + ats.zeta[i, j]*ats.NL1[i, j])
      end
    end

    calc_nl!(qts.NL2, ats.NL2, qts.sol1, ats.sol1, qts.ti, v, p, g)

    # Substep 2
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol2[i, j] = (qts.expLCdt2[i, j]*v.qh[i, j]
          + qts.zeta[i, j]*qts.NL2[i, j])
      end

      @simd for i in 1:g.nk
        @inbounds ats.sol2[i, j] = (ats.expLCdt2[i, j]*v.Ah[i, j]
          + ats.zeta[i, j]*ats.NL2[i, j])
      end
    end

    calc_nl!(qts.NL3, ats.NL3, qts.sol2, ats.sol2, qts.ti, v, p, g)

    # Substep 3
    qts.ti = v.t + qts.dt
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds qts.sol2[i, j] = (qts.expLCdt2[i, j]*qts.sol1[i, j]
          + qts.zeta[i, j]*(2.0*qts.NL3[i, j] - qts.NL1[i, j]))
      end

      @simd for i in 1:g.nkr
        @inbounds ats.sol2[i, j] = (ats.expLCdt2[i, j]*ats.sol1[i, j]
          + ats.zeta[i, j]*(2.0*ats.NL3[i, j] - ats.NL1[i, j]))
      end
    end

    # Final eval of nonlinear term
    calc_nl!(qts.NL4, ats.NL4, qts.sol2, ats.sol2, qts.ti, v, p, g)

    # Update
    for j in 1:g.nl
      @simd for i in 1:g.nkr
        @inbounds v.qh[i, j] = (qts.expLCdt[i, j]*v.qh[i, j]
          +     qts.alph[i, j] * qts.NL1[i, j]
          + 2.0*qts.beta[i, j] * (qts.NL2[i, j] + qts.NL3[i, j])
          +     qts.gamm[i, j] * qts.NL4[i, j] )
      end

      @simd for i in 1:g.nk
        @inbounds v.Ah[i, j] = (ats.expLCdt[i, j]*v.Ah[i, j]
          +     ats.alph[i, j] * ats.NL1[i, j]
          + 2.0*ats.beta[i, j] * (ats.NL2[i, j] + ats.NL3[i, j])
          +     ats.gamm[i, j] * ats.NL4[i, j] )
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

  # Turbulent parameters
  # v.psih = -v.qh .* g.invKKrsq
  # v.U = -irfft(im*g.Lr.*v.psih, g.nx)
  # v.V =  irfft(im*g.Kr.*v.psih, g.nx)

  # # Wave velocities
  # v.uh = 1.0/(p.alpha*p.f0) * (p.sig*g.K + im*p.f0*g.L) .* v.Ah
  # v.vh = 1.0/(p.alpha*p.f0) * (p.sig*g.L - im*p.f0*g.K) .* v.Ah
  #
  # v.u = real( ifft(v.uh) + conj(ifft(v.uh)) )
  # v.v = real( ifft(v.vh) + conj(ifft(v.vh)) )

end





end
