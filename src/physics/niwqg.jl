__precompile__()


module NIWQG

using FourierFlows

import FourierFlows: jacobian

export Timestepper, 
       TwoDGrid,
       Params,
       Vars,
       Equation

export set_q!, updatevars!




# Problems --------------------------------------------------------------------- 
function NIWQGProblem(nx::Int, Lx::Real, dt::Real, kap::Real, 
  nkap::Int, nu::Real, nnu::Int, eta::Real, f::Real, Us::Real, Vs::Real)

  gr = FourierFlows.TwoDGrid(nx, Lx)
  pr = Params(kap, nkap, nu, nnu, eta, f, Us, Vs) 
  eq = Equation(pr, gr)
  vs = Vars(gr)
  ts = Timestepper(dt, eq)

  Problem(gr, vs, pr, eq, ts)
end

function InitialValueProblem(nx, Lx, dt, kap, nkap, nu, nnu, eta, f, Us, Vs)
  NIWQGProblem(nx, Lx, dt, kap, nkap, nu, nnu, eta, f, Us, Vs)
end

function InitialValueProblem(nx, Lx, dt, kap, nkap, nu, nnu, eta, f)
  NIWQGProblem(nx, Lx, dt, kap, nkap, nu, nnu, eta, f, 0.0, 0.0)
end



function NIWQGProblem(nx::Int, Lx::Real; kap=nothing,
  nkap=4, nu=nothing, nnu=4, eta=1.0, f=1.0, Us=0.0, Vs=0.0, dt=nothing)

  if dt  == nothing;  dt = 0.01*2*pi/f;                     end
  if kap == nothing; kap = 0.1/(dt*(0.65*pi*nx/Lx)^nkap);   end
  if nu  == nothing;  nu = 0.1/(dt*(0.65*pi*nx/Lx)^nnu);    end

  gr = FourierFlows.TwoDGrid(nx, Lx)
  pr = Params(kap, nkap, nu, nnu, eta, f, Us, Vs) 
  eq = Equation(pr, gr)
  vs = Vars(gr)
  ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

  Problem(gr, vs, pr, eq, ts)
end






# Params ---------------------------------------------------------------------- 
abstract type NIWQGParams <: AbstractParams end

type Params <: NIWQGParams
  kap::Float64  # Potential vorticity hyperdiffusivity
  nkap::Int     # Order of Potential vorticity hyperdiffusivity
  nu::Float64   # Wave hyperviscosity
  nnu::Int      # Order of wave hyperviscosity
  eta::Float64  # Wave dispersivity
  kw::Float64   # Wave vertical-mode horizontal wavenumber N*m/f
  f::Float64    # Planetary vorticity
  Us::Float64   # x-velocity of uniform mean flow
  Vs::Float64   # y-velocity of uniform mean flow
end


""" Construct a Params type with zero mean velocity and kw computed from
f and eta. """
function Params(kap, nkap, nu, nnu, eta, f)
  kw = sqrt(f/eta)
  Params(kap, nkap, nu, nnu, eta, kw, f, 0.0, 0.0)
end


""" Construct a Params type with kw computed from f and eta. """
function Params(kap::Real, nkap::Int, nu::Real, nnu::Int, eta::Real, f::Real, 
  Us::Real, Vs::Real)
  kw = sqrt(2*eta/f)
  Params(kap, nkap, nu, nnu, eta, kw, f, Us, Vs)
end






# Equations ------------------------------------------------------------------- 
type Equation <: AbstractEquation
  LCc::Array{Complex{Float64}, 2}  # Coeff of the complex eqn's linear part
  LCr::Array{Complex{Float64}, 2}  # Coeff of the real eqn's linear part
  calcNL!::Function                # Function to calculate eqn's nonlinear part
end




""" Construct the default, fully-coupled NIWQG equation type with
hyperdiffusivity and hyperviscosity. """
function Equation(p::NIWQGParams, g::TwoDGrid)
  # Linear term for phi:
  LCc = -p.nu  * g.KKsq.^(0.5*p.nnu) - 0.5*im*p.eta*g.KKsq 

  # Linear term for q:
  LCr = -p.kap * g.KKrsq.^(0.5*p.nkap)
  Equation(LCc, LCr, calcNL!)
end





# Vars ------------------------------------------------------------------------ 
type Vars <: AbstractVars

  t::Float64
  solc::Array{Complex128, 2} # solc = phih
  solr::Array{Complex128, 2} # solr = qh

  # Auxiliary vorticity vars
  q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  zeta::Array{Float64, 2}
  psi::Array{Float64, 2}
  Uq::Array{Float64, 2}
  Vq::Array{Float64, 2}
  modphi::Array{Float64, 2}
  jacphi::Array{Float64, 2}

  # Auxiliary wave vars
  phi::Array{Complex{Float64}, 2}
  Uphi::Array{Complex{Float64}, 2}
  Vphi::Array{Complex{Float64}, 2}
  phix::Array{Complex{Float64}, 2}
  phiy::Array{Complex{Float64}, 2}
  zetaphi::Array{Complex{Float64}, 2}
  
  # Vorticity transforms
  qh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  zetah::Array{Complex{Float64}, 2}
  psih::Array{Complex{Float64}, 2}
  Uqh::Array{Complex{Float64}, 2}
  Vqh::Array{Complex{Float64}, 2}
  modphih::Array{Complex{Float64}, 2}
  jacphih::Array{Complex{Float64}, 2}

  # Wave transforms
  phih::Array{Complex{Float64}, 2}
  Uphih::Array{Complex{Float64}, 2}
  Vphih::Array{Complex{Float64}, 2}
  phixh::Array{Complex{Float64}, 2}
  phiyh::Array{Complex{Float64}, 2}
  zetaphih::Array{Complex{Float64}, 2}

end




""" Construct a vars type for the default NIWQG problem. """
function Vars(g::TwoDGrid)

  # Initialize with t=0
  t = 0.0
  solc  = zeros(Complex{Float64}, g.nk, g.nl)
  solr  = zeros(Complex{Float64}, g.nkr, g.nl)

  # Auxiliary vorticity vars
  q      = zeros(Float64, g.nx, g.ny)
  U      = zeros(Float64, g.nx, g.ny)
  V      = zeros(Float64, g.nx, g.ny)
  zeta   = zeros(Float64, g.nx, g.ny)
  psi    = zeros(Float64, g.nx, g.ny)
  Uq     = zeros(Float64, g.nx, g.ny)
  Vq     = zeros(Float64, g.nx, g.ny)
  modphi = zeros(Float64, g.nx, g.ny)
  jacphi = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary wave vars
  phi     = zeros(Complex{Float64}, g.nx, g.ny)
  Uphi    = zeros(Complex{Float64}, g.nx, g.ny)
  Vphi    = zeros(Complex{Float64}, g.nx, g.ny)
  phix    = zeros(Complex{Float64}, g.nx, g.ny)
  phiy    = zeros(Complex{Float64}, g.nx, g.ny)
  zetaphi = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  qh       = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh       = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh       = zeros(Complex{Float64}, g.nkr, g.nl)
  zetah    = zeros(Complex{Float64}, g.nkr, g.nl)
  psih     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uqh      = zeros(Complex{Float64}, g.nkr, g.nl)
  Vqh      = zeros(Complex{Float64}, g.nkr, g.nl)
  modphih  = zeros(Complex{Float64}, g.nkr, g.nl)
  jacphih  = zeros(Complex{Float64}, g.nkr, g.nl)

  phih     = zeros(Complex{Float64}, g.nk, g.nl)
  Uphih    = zeros(Complex{Float64}, g.nk, g.nl)
  Vphih    = zeros(Complex{Float64}, g.nk, g.nl)
  phixh    = zeros(Complex{Float64}, g.nk, g.nl)
  phiyh    = zeros(Complex{Float64}, g.nk, g.nl)
  zetaphih = zeros(Complex{Float64}, g.nk, g.nl)
  
  return Vars(t, solc, solr, 
    q, U, V, zeta, psi, Uq, Vq, modphi, jacphi,
    phi, Uphi, Vphi, phix, phiy, zetaphi,
    qh, Uh, Vh, zetah, psih, Uqh, Vqh, modphih, jacphih,
    phih, Uphih, Vphih, phixh, phiyh, zetaphih,
    )
end




# Solvers ---------------------------------------------------------------------

""" Calculate the 'nonlinear' right hand side of the NIWQG system. The 
result for the complex and real parts of the system are stored in the first
two arguments, NLc, and NLr. """
function calcNL!(
  NLc::Array{Complex{Float64}, 2},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 2}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::Params, g::TwoDGrid)


  # Prelim calcs
  @. v.qh = solr # because irfft destroys input

  # Wave contribution to PV
  @. v.phixh = im*g.K*solc
  @. v.phiyh = im*g.L*solc

  A_mul_B!(v.phi,  g.ifftplan, solc)
  A_mul_B!(v.phix, g.ifftplan, v.phixh)
  A_mul_B!(v.phiy, g.ifftplan, v.phiyh)

  @. v.modphi = abs2(v.phi)
  @. v.jacphi = real(im*conj(v.phix)*v.phiy - im*conj(v.phiy)*v.phix)

  A_mul_B!(v.modphih, g.rfftplan, v.modphi)
  A_mul_B!(v.jacphih, g.rfftplan, v.jacphi)


  # Mean flow calcs
  @. v.zetah = v.qh + 0.25/p.f*g.KKrsq*v.modphih - 0.5/p.f*v.jacphih
  @. v.psih  = -g.invKKrsq*v.zetah
  @. v.Uh    = -im*g.Lr*v.psih
  @. v.Vh    =  im*g.Kr*v.psih

  # Add mean flow (with normalization for irfft)
  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny


  # Inverse transforms
  A_mul_B!(v.q,    g.irfftplan, solr)
  A_mul_B!(v.U,    g.irfftplan, v.Uh)
  A_mul_B!(v.V,    g.irfftplan, v.Vh)
  A_mul_B!(v.zeta, g.irfftplan, v.zetah)


  # Multiplies 
  @. v.Uq = v.U * v.q
  @. v.Vq = v.V * v.q

  @. v.Uphi    = v.U    * v.phi
  @. v.Vphi    = v.V    * v.phi
  @. v.zetaphi = v.zeta * v.phi


  # Forward transforms
  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  A_mul_B!(v.Uphih,    g.fftplan, v.Uphi)
  A_mul_B!(v.Vphih,    g.fftplan, v.Vphi)
  A_mul_B!(v.zetaphih, g.fftplan, v.zetaphi)


  # Nonlinear terms
  @. NLr = -im*g.Kr*v.Uqh  - im*g.Lr*v.Vqh
  @. NLc = -im*g.K*v.Uphih - im*g.L*v.Vphih - 0.5*im*v.zetaphih

  nothing
end




# Default timesteppers --------------------------------------------------------- 

""" Construct an ETDRK4TimeStepper for solving the NIWQG problem. """
function Timestepper(dt::Real, eq::Equation)
  FourierFlows.ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)
end




# Helper functions ------------------------------------------------------------ 

""" Update variables for current value of solc and solr. """
function updatevars!(v::Vars, p::NIWQGParams, g::TwoDGrid)

  @. v.phih = v.solc 
  @. v.qh = v.solr

  # Wave calcs and inverse transforms
  @. v.phixh = im*g.K*v.phih
  @. v.phiyh = im*g.L*v.phih

  A_mul_B!(v.phi,  g.ifftplan, v.phih)
  A_mul_B!(v.phix, g.ifftplan, v.phixh)
  A_mul_B!(v.phiy, g.ifftplan, v.phiyh)


  # Wave parts of PV. J(conj(phi), phi) = conj(phix)*phiy - conj(phiy)*phix
  @. v.modphi = abs2(v.phi)
  @. v.jacphi = real(im*conj(v.phix)*v.phiy - im*conj(v.phiy)*v.phix)

  A_mul_B!(v.modphih, g.rfftplan, v.modphi)
  A_mul_B!(v.jacphih, g.rfftplan, v.jacphi)


  # Vorticity calcs and inverse transforms
  @. v.zetah = v.qh + 0.25/p.f*g.KKrsq*v.modphih - 0.5/p.f*v.jacphih
  @. v.psih  = -g.invKKrsq*v.zetah
  @. v.Uh    = -im*g.Lr*v.psih
  @. v.Vh    =  im*g.Kr*v.psih


  # We don't use A_mul_B here because irfft destroys its input.
  v.q    = irfft(v.qh, g.nx)
  v.U    = irfft(v.Uh, g.nx) 
  v.V    = irfft(v.Vh, g.nx) 
  v.zeta = irfft(v.zetah, g.nx) 

  nothing
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.params, prob.grid)
end




""" Update vars v given a vorticity field q in physical space. """
function set_q!(v::Vars, p::NIWQGParams, g::TwoDGrid, q)
  A_mul_B!(v.solr, g.rfftplan, q)
  updatevars!(v, p, g)
  nothing
end

function set_q!(prob::AbstractProblem, q)
  set_q!(prob.vars, prob.params, prob.grid, q)
end




""" Update vars v given a complex wave amplitude phi in physical space. """
function set_phi!(vs::Vars, pr::NIWQGParams, g::TwoDGrid, phi)
  phi = convert(Array{Complex{Float64}, 2}, phi)
  A_mul_B!(vs.solc, g.fftplan, phi)
  updatevars!(vs, pr, g)
  nothing
end

function set_phi!(prob::AbstractProblem, phi)
  set_phi!(prob.vars, prob.params, prob.grid, phi)
end


function delphi(vs, pr, g)
  0.5*(vs.phix-im*vs.phiy)
end

  

""" Return the horizontal velocities u and v, vertical velocity w, and 
buoyancy b. """
function getprimitivefields(vs::Vars, pr::NIWQGParams, g::TwoDGrid; z=0.0, 
  N=10*pr.f)

  updatevars!(vs, pr, g)
  
  # k = N*m/f
  m = pr.f*pr.kw/N

  phase = m*z - pr.f*vs.t

  delphi = 0.5*(vs.phix - im*vs.phiy)

  wc = im/m * exp(im*phase) * delphi
  bc = m*pr.eta * exp(im*phase) * delphi

  u = real.(exp(im*phase) * vs.phi) 
  v = imag.(exp(im*phase) * vs.phi) 
  w = real.(wc + conj.(wc))
  b = real.(bc + conj.(bc))

  return u, v, w, b
end

function getprimitivefields(prob::AbstractProblem)
  getprimitivefields(prob.vars, prob.params, prob.grid)
end


""" Returns the net NIW kinetic energy. """
function niwke(v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  0.5*FourierFlows.parsevalsum2(v.solc, g)
end

function niwke(prob::AbstractProblem)
  niwke(prob.vars, prob.params, prob.grid)
end 


""" Returns the net NIW potential energy. """
function niwpe(v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  #0.25/p.kw^2.0*g.dx*g.dy*sum(abs2.(v.phix)+abs2.(v.phiy))
  0.25/p.kw^2.0*FourierFlows.parsevalsum2(g.KKsq.*v.solc, g)
end

function niwpe(prob::AbstractProblem)
  niwpe(prob.vars, prob.params, prob.grid)
end


""" Returns the net QG energy. """
function qgenergy(v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  (0.5*g.dx*g.dy*sum(v.U.^2.0+v.V.^2.0) 
    + 0.25/p.kw^2.0*g.dx*g.dy*sum(abs2.(v.phix)+abs2.(v.phiy)))
end

function qgenergy(prob::AbstractProblem)
  qgenergy(prob.vars, prob.params, prob.grid)
end


""" Returns the net kinetic energy dissipation. """
function niwke_dissipation(
  v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  -p.nu*FourierFlows.parsevalsum(g.KKsq.^(p.nnu/2).*abs2.(v.phih), g)
end

function niwke_dissipation(prob::AbstractProblem)
  niwke_dissipation(prob.vars, prob.params, prob.grid)
end


""" Returns the net potential energy dissipation. """
function niwpe_dissipation(
  v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  -0.5/p.kw^2.0*p.nu*(
    FourierFlows.parsevalsum2(im*g.K.*g.KKsq.^(p.nnu/4).*v.phih, g)
    + FourierFlows.parsevalsum2(im*g.L.*g.KKsq.^(p.nnu/4).*v.phih, g)
  )
end

function niwpe_dissipation(prob::AbstractProblem)
  niwpe_dissipation(prob.vars, prob.params, prob.grid)
end


""" Returns the net dissipation of balanced kinetic energy. """
function qg_dissipation(
  v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  delq = irfft(g.KKrsq.^(p.nkap/2).*v.qh, g.nx)
  p.kap*g.dx*g.dy*sum(v.psi.*delq)
end

function qg_dissipation(prob::AbstractProblem)
  qg_dissipation(prob.vars, prob.params, prob.grid)
end


""" Returns the wave-QG cross-term that contributes to the balanced energy 
budget. """
function waveqg_dissipation(
  v::AbstractVars, p::AbstractParams, g::AbstractGrid)

  delphi = (-1.0)^(p.nnu/2).*ifft(g.KKsq.^(p.nnu/2).*v.phih)

  chi = 0.25.*p.nu./p.f .* v.zeta .* real.(
    conj.(v.phi).*delphi .+ v.phi.*conj.(delphi))

  chi .+= 0.5.*p.nu./p.f .* v.psi .* real.(
              im.*jacobian(conj.(v.phi), delphi, g) 
           .- im.*jacobian(v.phi, conj.(delphi), g)
  )

  g.dx*g.dy*sum(chi)
end

function waveqg_dissipation(prob::AbstractProblem)
  waveqg_dissipation(prob.vars, prob.params, prob.grid)
end



# End module
end
