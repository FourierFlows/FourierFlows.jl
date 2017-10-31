include("../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

# Abbrev's
TwoDGrid      = FourierFlows.TwoDGrid
Vars          = TwoModeBoussinesq.Vars
TwoModeParams = TwoModeBoussinesq.TwoModeParams




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

  psiwh = g.KKrsq.*qwh
  uwh   = -im*g.Lr.*psiwh
  vwh   =  im*g.Kr.*psiwh

  psiw = irfft(psiwh, g.nx)
  uw   = irfft(uwh, g.nx)
  vw   = irfft(vwh, g.nx)

  return psiw, uw, vw
end




""" Calculate the wave contribution to PV, qw. """
function calc_qw(usig, vsig, sig, p::TwoModeParams, g::TwoDGrid)

  usigh = fft(usig)
  vsigh = fft(vsig)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*abs2.(usigh))
  vsig2yy = ifft(-g.L.^2.0.*abs2.(vsigh))
  usigvsigxy = ifft(-g.K.*g.L.*(usigh.*conj.(vsigh) + conj.(usigh).*vsigh))

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

  # This copy is necessary because calling A_mul_B(v.q, g.irfftplan, sol) 
  # a few lines below destroys solr
  v.qh .= solr

  @. v.psih = -g.invKKrsq*v.qh
  @. v.Uh   = -im*g.Lr*v.psih
  @. v.Vh   =  im*g.Kr*v.psih
  @. v.Uxh  = im*g.Kr*v.Uh
  @. v.Vxh  = im*g.Kr*v.Vh
  @. v.Uyh  = im*g.Lr*v.Uh
  @. v.Vyh  = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  # Inverse transforms
  A_mul_B!(v.q,  g.irfftplan, solr)
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
