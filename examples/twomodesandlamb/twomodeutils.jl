__precompile__()

using FourierFlows, 
      PyPlot

import FourierFlows.TwoModeBoussinesq,
       FourierFlows.NIWQG

# Abbrev's
Vars          = TwoModeBoussinesq.Vars
TwoModeParams = TwoModeBoussinesq.TwoModeParams




rms(a) = sqrt(mean(a.^2))




type FourComponentPlot
  g::TwoDGrid
  v::Vars
  p::TwoModeParams

  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String

  comp2::Function
  title2::String
  limz2::Array{Float64, 1}
  colors2::String

  comp3::Function
  title3::String
  limz3::Array{Float64, 1}
  colors3::String

  comp4::Function
  title4::String
  limz4::Array{Float64, 1}
  colors4::String

  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String
end




type TwoComponentPlot
  g::TwoDGrid
  v::Vars
  p::TwoModeParams

  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String

  comp2::Function
  title2::String
  limz2::Array{Float64, 1}
  colors2::String

  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String
end




type ComparisonPlot
  g::TwoDGrid
  v_bo::Vars
  p_bo::TwoModeParams
  v_qg::NIWQG.Vars
  p_qg::NIWQG.Params
  comp1_bo::Function
  comp2_bo::Function
  comp1_qg::Function
  comp2_qg::Function
  Lnorm::Float64
  limz1::Array{Float64, 1}
  limz2::Array{Float64, 1}
  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  title1_bo::String
  title2_bo::String
  title1_qg::String
  title2_qg::String
  colors1::String
  colors2::String
  xlabel::String
  ylabel::String
end




function makeplot!(axs, pl::FourComponentPlot)

  c1 = pl.comp1(pl.v, pl.p, pl.g)
  c2 = pl.comp2(pl.v, pl.p, pl.g)
  c3 = pl.comp3(pl.v, pl.p, pl.g)
  c4 = pl.comp4(pl.v, pl.p, pl.g)


  axes(axs[1, 1])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c1, 
    cmap=pl.colors1,
    vmin=pl.limz1[1], vmax=pl.limz1[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title1)
  ylabel(pl.ylabel)


  axes(axs[1, 2])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c2, 
    cmap=pl.colors2,
    vmin=pl.limz2[1], vmax=pl.limz2[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title2)


  axes(axs[2, 1])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c3, 
    cmap=pl.colors3,
    vmin=pl.limz3[1], vmax=pl.limz3[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title3)
  xlabel(pl.xlabel)
  ylabel(pl.ylabel)


  axes(axs[2, 2])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c4, 
    cmap=pl.colors4,
    vmin=pl.limz4[1], vmax=pl.limz4[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title4)
  xlabel(pl.xlabel)

  pause(0.01)

end








function makeplot!(axs, pl::TwoComponentPlot)

  c1 = pl.comp1(pl.v, pl.p, pl.g)
  c2 = pl.comp2(pl.v, pl.p, pl.g)


  axes(axs[1])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c1, 
    cmap=pl.colors1,
    vmin=pl.limz1[1], vmax=pl.limz1[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title1)
  xlabel(pl.xlabel)
  ylabel(pl.ylabel)


  axes(axs[2])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c2, 
    cmap=pl.colors2,
    vmin=pl.limz2[1], vmax=pl.limz2[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title2)
  xlabel(pl.xlabel)

  pause(0.01)

end





function makeplot!(axs, pl::ComparisonPlot)

  c1_qg = pl.comp1_qg(pl.v_qg, pl.p_qg, pl.g)
  c2_qg = pl.comp2_qg(pl.v_qg, pl.p_qg, pl.g)

  c1_bo = pl.comp1_bo(pl.v_bo, pl.p_bo, pl.g)
  c2_bo = pl.comp2_bo(pl.v_bo, pl.p_bo, pl.g)

  axes(axs[1, 1])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c1_qg, 
    cmap=pl.colors1,
    vmin=pl.limz1[1], vmax=pl.limz1[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title1_qg)
  ylabel(pl.ylabel)


  axes(axs[1, 2])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c2_qg, 
    cmap=pl.colors2,
    vmin=pl.limz2[1], vmax=pl.limz2[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title2_qg)




  axes(axs[2, 1])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c1_bo, 
    cmap=pl.colors1,
    vmin=pl.limz1[1], vmax=pl.limz1[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title1_bo)
  xlabel(pl.xlabel)
  ylabel(pl.ylabel)


  axes(axs[2, 2])

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c2_bo, 
    cmap=pl.colors2,
    vmin=pl.limz2[1], vmax=pl.limz2[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  title(pl.title2_bo)
  xlabel(pl.xlabel)


  pause(0.01)

end






""" Returns the wave speed projected onto the zeroth mode. """
function wavespeed(vs::Vars)
  sqrt.(2*abs2.(vs.u) + 2*abs2.(vs.v))
end

function wavespeed(vs::NIWQG.Vars)
  abs.(vs.phi)
end





function niwqgplot(axs, vs, pr, g, Z0, Uw, R, tnd, E0) 

  TwoModeBoussinesq.updatevars!(vs, pr, g)

  domfrac = 5
  xl, xr = -g.Lx/domfrac, g.Lx/domfrac
  yl, yr = -g.Ly/domfrac, g.Ly/domfrac

  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  E = TwoModeBoussinesq.calc_energy(vs, pr, g)
  sp = sqrt.(2*abs2.(vs.u) + 2*abs2.(vs.v))

  # All lengths non-dimensionalized by R
  axes(axs[1])
  pcolormesh(g.X/R, g.Y/R, q*tnd, cmap="RdBu_r", vmin=-Z0, vmax=Z0)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  title(L"q"); xlabel(L"x/R"); ylabel(L"y/R")

  axes(axs[2])
  pcolormesh(g.X/R, g.Y/R, sp, cmap="YlGnBu_r", vmin=0, vmax=Uw)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  title(L"\sqrt{u^2+v^2}"); xlabel(L"x/R")

  @printf("max Ro: %.2e, max speed: %.3f, t: %.3f, total energy: %.3f\n",
     maximum(abs.(q))/pr.f, maximum(sp), vs.t/tnd, E/E0)

  pause(0.01)

  nothing
end





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
