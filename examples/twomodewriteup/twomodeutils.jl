__precompile__()

using FourierFlows, 
      PyPlot

import FourierFlows.TwoModeBoussinesq,
       FourierFlows.NIWQG

# Abbrev's
Vars          = TwoModeBoussinesq.Vars
TwoModeParams = TwoModeBoussinesq.TwoModeParams

rms(a) = sqrt(mean(a.^2))

abstract type AbstractPlot end

""" Compute the transform of the Jacobian of two fields a, b on a grid g. """
function jach(a, b, g::TwoDGrid)
  # J(a, b) = dx(a b_y) - dy(a b_x)

  bh = fft(b)
  bx = ifft(im*g.K.*bh)
  by = ifft(im*g.L.*bh)

  return im*g.K.*fft(a.*bx) - im*g.L.*fft(a.*by)
end

""" Compute the Jacobian of two fields a, b on a grid g. """
function jac(a, b, g::TwoDGrid)
  ifft(jach(a, b, g))
end

""" Calculate the wave-induced streamfunction and velocity fields. """
function calc_uw(qw, g::TwoDGrid)

  qwh = rfft(qw)

  psiwh = g.invKKrsq.*qwh
  uwh   = -im*g.Lr.*psiwh
  vwh   =  im*g.Kr.*psiwh

  psiw = irfft(psiwh, g.nx)
  uw   = irfft(uwh, g.nx)
  vw   = irfft(vwh, g.nx)

  return uw, vw
end




""" Calculate the wave-induced streamfunction and velocity fields. """
function calc_uw(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid)
  calc_uw(calc_qw(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid), 
    g::TwoDGrid)
end




""" Calculate the wave contribution to PV, qw. """
function calc_qw(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid)

  usig, vsig = calc_usigvsig(sig, v, p, g)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*fft(abs2.(usig)))
  vsig2yy = ifft(-g.L.^2.0.*fft(abs2.(vsig)))
  usigvsigxy = ifft(-g.K.*g.L.*fft(usig.*conj.(vsig) + conj.(usig).*vsig))

  # Assemble contributions
  qw = -real.( 
       2.0*im/sig * (jac(conj.(usig), usig, g) + jac(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (jac(conj.(vsig), usig, g) + jac(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  return qw
end






""" Calculate the wave contribution to PV, qw. """
function calc_qw(usig::AbstractArray, vsig::AbstractArray, sig::Real, 
  p::TwoModeParams, g::TwoDGrid)

  usigh = fft(usig)
  vsigh = fft(vsig)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*fft(abs2.(usig)))
  vsig2yy = ifft(-g.L.^2.0.*fft(abs2.(vsig)))
  usigvsigxy = ifft(-g.K.*g.L.*fft(usig.*conj.(vsig) + conj.(usig).*vsig))

  # Assemble contributionsCalulate
  qw = -real.( 
       2.0*im/sig * (jac(conj.(usig), usig, g) + jac(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (jac(conj.(vsig), usig, g) + jac(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  return qw
end




""" Calculate usig and vsig, the complex, sigm-ified amplitudes 
of u_1 and v_1. """
function calc_usigvsig(sig, v::Vars, p::TwoModeParams, g::TwoDGrid)

  # Calculate u_t and v_t, and use them to find usig and vsig

  @views @. v.uh = v.solc[:, :, 1]
  @views @. v.vh = v.solc[:, :, 2]
  @views @. v.ph = v.solc[:, :, 3]

  # This copy is necessary because calling A_mul_B(v.Z, g.irfftplan, sol) 
  # a few lines below destroys solr
  v.Zh .= v.solr

  @. v.psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.Lr*v.psih
  @. v.Vh   =  im*g.Kr*v.psih
  @. v.Uxh  = im*g.Kr*v.Uh
  @. v.Vxh  = im*g.Kr*v.Vh
  @. v.Uyh  = im*g.Lr*v.Uh
  @. v.Vyh  = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  # Inverse transforms
  A_mul_B!(v.U,  g.irfftplan, v.Uh)
  A_mul_B!(v.V,  g.irfftplan, v.Vh)
  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  A_mul_B!(v.u,  g.ifftplan, v.uh)
  A_mul_B!(v.v,  g.ifftplan, v.vh)

  @. v.Uu =  v.U * v.u
  @. v.Vu =  v.V * v.u
  @. v.Uv =  v.U * v.v
  @. v.Vv =  v.V * v.v
  @. v.uUx = v.u * v.Ux
  @. v.uVx = v.u * v.Vx
  @. v.vUy = v.v * v.Uy
  @. v.vVy = v.v * v.Vy

  A_mul_B!(v.Uuh,  g.fftplan, v.Uu)
  A_mul_B!(v.Uvh,  g.fftplan, v.Uv)
  A_mul_B!(v.Vuh,  g.fftplan, v.Vu)
  A_mul_B!(v.Vvh,  g.fftplan, v.Vv)
  A_mul_B!(v.uUxh, g.fftplan, v.uUx)
  A_mul_B!(v.uVxh, g.fftplan, v.uVx)
  A_mul_B!(v.vUyh, g.fftplan, v.vUy)
  A_mul_B!(v.vVyh, g.fftplan, v.vVy)

  # First-mode nonlinear terms:
  # u
  uth = ( -p.nu1*g.KKsq.^(0.5*p.nnu1) .* v.uh
    + p.f*v.vh - im*g.K.*v.ph
    - im*g.K.*v.Uuh - im*g.L.*v.Vuh - v.uUxh - v.vUyh
  )

  # v
  vth = ( -p.nu1*g.KKsq.^(0.5*p.nnu1) .* v.vh
    - p.f*v.uh - im*g.L.*v.ph
    - im*g.K.*v.Uvh - im*g.L.*v.Vvh - v.uVxh - v.vVyh
  )

  ut = ifft(uth)
  vt = ifft(vth)

  # Calculate amplitudes
  usig = exp(im*sig*v.t) * (v.u + im/sig*ut)
  vsig = exp(im*sig*v.t) * (v.v + im/sig*vt)

  return usig, vsig
end


# Potential plot functions
rossbyq(vs, pr, g)      = TwoModeBoussinesq.calc_apv(vs, pr, g) / pr.f
rossbynum(vs, pr, g)    = vs.Z / pr.f
waveu(vs, pr, g)        = real.(vs.u + conj.(vs.u))
wavev(vs, pr, g)        = real.(vs.v + conj.(vs.v))
wavew(vs, pr, g)        = real.(vs.w + conj.(vs.w))
wavespeed(vs, pr, g)    = sqrt.(waveu(vs, pr, g).^2.0 + wavev(vs, pr, g).^2.0)
wavepressure(vs, pr, g) = real.(vs.p + conj.(vs.p)) 
wavebuoyancy(vs, pr, g) = real.(im*pr.m*vs.p - im*pr.m*conj.(vs.p))
meanspeed(vs, pr, g)    = sqrt.(vs.U.^2.0 + vs.V.^2.0)

function waveinducedflow(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return sqrt.(uw.^2.0 + vw.^2.0)
end

function waveinducedu(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return uw
end

function waveinducedv(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return vw
end

function apvinducedflow(vs, pr, g)
  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  psiqh = -g.invKKrsq.*rfft(q)
  uq = irfft(-im*g.Lr.*psiqh, g.nx)
  vq = irfft( im*g.Kr.*psiqh, g.nx)
  return sqrt.(uq.^2.0+vq.^2.0)
end
