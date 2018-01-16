using Interpolations
import SpecialFunctions

export @createarrays

# Utility for generating time-steppers.

# Time-steppers lists
steppers = [
  "ForwardEuler",
  "FilteredForwardEuler",
  "AB3",
  "RK4",
  "FilteredRK4",
  "ETDRK4",
  "FilteredETDRK4",
]

filteredsteppers = [
  "FilteredForwardEuler",
  "FilteredRK4",
  "FilteredETDRK4",
]

"""
Returns a time-stepper type defined by the prefix 'stepper', timestep dt
solution sol (used to construct variables with identical type and size as
the solution vector), and grid g.
"""
function autoconstructtimestepper(stepper, dt, sol, 
                                  g::AbstractGrid=ZeroDGrid(1))
  fullsteppername = Symbol(stepper, :TimeStepper)
  if stepper ∈ filteredsteppers
    tsexpr = Expr(:call, fullsteppername, dt, sol, g)
  else
    tsexpr = Expr(:call, fullsteppername, dt, sol)
  end

  eval(tsexpr)
end

function autoconstructtimestepper(stepper, dt, solc, solr)
  fullsteppername = Symbol(stepper, :TimeStepper)
  tsexpr = Expr(:call, fullsteppername, dt, solc, solr)
  eval(tsexpr)
end



"""
    @createarrays T dims a b c 

Create arrays of all zeros with element type T, size dims, and global names
a, b, c (for example). An arbitrary number of arrays may be created.
"""
macro createarrays(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, 
    [:( $(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end


"""
This function returns an expression that defines a Composite Type
of the AbstractVars variety.
"""
function structvarsexpr(name, physfields, transfields; soldims=2, vardims=2, 
                          parent=:AbstractVars) 
  
  physexprs = [:( $fld::Array{Float64,$vardims} ) 
    for fld in physfields]
  transexprs = [:( $fld::Array{Complex{Float64},$vardims} ) 
    for fld in transfields]

  if parent != nothing
    expr = quote
      mutable struct $name <: $parent
        $(physexprs...)
        $(transexprs...)
      end
    end
  else
    expr = quote
      mutable struct $name
        $(physexprs...)
        $(transexprs...)
      end
    end
  end
    
  expr
end


"""
    fftwavenums(n; L=1)

Return the fftwavenumber vector with length n and domain size L.
"""
fftwavenums(n::Int; L=1.0) = 2.0*pi/L*cat(1, 0:n/2, -n/2+1:-1)


""" 
    rms(q)

Return the root-mean-square of an array. 
"""
rms(q) = sqrt(mean(q.^2))


"""
    peaked_isotropic_spectrum(nkl, kpeak; ord=4, rms=1, maxval=1)

Generate a real and random two-dimensional distribution phi(x, y) with
a Fourier spectrum peaked around a central non-dimensional wavenumber kpeak.
The spectrum is normalized either by setting the root-mean-square value of phi
with the keyword 'rms', or the maximum value of phi with the keyword 'maxval'.

    peaked_isotropic_spectrum(nx, npeak; ord=4, rms=1, maxval=1)
"""
function peaked_isotropic_spectrum(nkl::Tuple{Int, Int}, kpeak::Real;
  ord=4.0, rms=1.0, maxval=0.0)

  # Non-dimensional wavenumbers
  nk, nl = nkl
  k, l   = fftwavenums(nk), fftwavenums(nl)

  K = zeros(Float64, nk, nl)
  for j = 1:nl, i = 1:nk
    K[i, j] = sqrt(k[i]^2.0 + l[j]^2.0)
  end

  # Generate random spectrum and then normalize
  phih = exp.(2.0*im*pi*rand(nk, nl)) ./ (1.0 .+ K./kpeak).^ord

  # Normalize by maximum value if specified
  if maxval > 0
    phi = real.(ifft(phih))
    phi = maxval * phi / maximum(abs.(phi))
  else
    phih .*= rms ./ sqrt.(sum(abs.(phih).^2.0))
    phi = real.(ifft(phih))
  end

  return phi
end

function peaked_isotropic_spectrum(nx::Int, npeak::Real; ord=4.0, rms=1.0,
  maxval=0.0)
  peaked_isotropic_spectrum((nx, nx), npeak; ord=ord, rms=rms, maxval=maxval)
end


""" 
    lambdipole(Ue, R, g; center=(x0, y0))

Return a 2D vorticity field corresponding to the Lamb Dipole with
strength Ue, radius R, and centered around
(xc, yc)=center. The default value of 'center' is the middle of the grid.
"""
function lambdipole(Ue::Real, R::Real, g::TwoDGrid; center=(nothing, nothing))

  if center == (nothing, nothing)
    xc = mean(g.x)
    yc = mean(g.y)
  else
    xc = center[1]
    yc = center[2]
  end

  # Wavenumber corresponding to radius R and the first bessel func zero.
  k = 3.8317 / R
  q0 = -2*Ue*k/SpecialFunctions.besselj(0, k*R)

  r = sqrt.((g.X-xc).^2.0 + (g.Y-yc).^2.0)
  q = q0 * SpecialFunctions.besselj.(1, k*r) .* (g.Y-yc)./r

  q[r .== 0.0] = 0.0 # just in case.
  q[r .> R] = 0.0

  return q
end


""" 
    gaussianvortex(q0, R, G; center=(x0, y0))

Return a vorticity field with magnitude q0, radius R, and center at
center[1], center[2] on a TwoDGrid g corresponding to a 'Gaussian vortex' with
Gaussian streamfunction. 
"""
function gaussianvortex(q0::Real, R::Real, g::TwoDGrid;
  center=(nothing, nothing))

  if center == (nothing, nothing)
    xc = mean(g.x)
    yc = mean(g.y)
  else
    xc = center[1]
    yc = center[2]
  end

  ( q0/R^2.0 * ( (g.X-xc).^2.0 + (g.Y-yc).^2.0 - 2*R^2.0 )
        .* exp.( -((g.X-xc).^2.0 + (g.Y-yc).^2.0) / (2.0*R^2.0)) )
end


""" 
    rmsrand(g, rmsval)

Return an array of random numbers on a TwoDGrid normalized to have a
specifed rms value. 
"""
function rmsrand(g::TwoDGrid, rmsval::Real)
  q = rand(g.nx, g.ny)
  q .*= rmsval / rms(q)
  return q
end


""" 
    parsevalsum2(uh, g)

Calculate ∫u = Σ|uh|² on a 2D grid, where uh is the Fourier transform of u.
Accounts for DFT normalization, grid resolution, and whether or not uh
is the product of fft or rfft.
"""
function parsevalsum2(uh, g::TwoDGrid)
  norm = g.Lx*g.Ly/(g.nx^2*g.ny^2)    # weird normalization for dft

  if size(uh)[1] == g.nkr             # uh is conjugate symmetric
    U = sum(abs2, uh[1, :])           # k=0 modes
    U += 2*sum(abs2, uh[2:end, :])    # sum k>0 modes twice     
  else                                # count every mode once
    U = sum(abs2, uh)
  end

  norm*U
end

""" 
    parsevalsum(uh, g)

Calculate real(Σ uh) on a 2D grid.  Accounts for DFT normalization, 
grid resolution, and whether or not uh is in a conjugate-symmetric form to 
save memory.
""" 
function parsevalsum(uh, g::TwoDGrid)
  norm = g.Lx*g.Ly/(g.nx^2*g.ny^2) # weird normalization for dft

  if size(uh)[1] == g.nkr       # uh is conjugate symmetric
    U = sum(uh[1, :])           # k=0 modes
    U += 2*sum(uh[2:end, :])    # sum k>0 modes twice     
  else # count every mode once
    U = sum(uh)
  end

  norm*real(U)
end


"""
    jacobianh(a, b, g)

Returns the transform of the Jacobian of two fields a, b on the grid g.
"""
function jacobianh(a, b, g::TwoDGrid)
  if eltype(a) <: Real
    bh = rfft(b)
    bx = irfft(im*g.kr.*bh, g.nx)
    by = irfft(im*g.l.*bh, g.nx)
    return im*g.kr.*rfft(a.*by)-im*g.l.*rfft(a.*bx)
  else
    # J(a, b) = dx(a b_y) - dy(a b_x)
    bh = fft(b)
    bx = ifft(im*g.k.*bh)
    by = ifft(im*g.l.*bh)
    return im*g.k.*fft(a.*by).-im*g.l.*fft(a.*bx)
  end
end


"""
    jacobian(a, b, g)

Returns the Jacobian of a and b.
"""
function jacobian(a, b, g::TwoDGrid)
  if eltype(a) <: Real
   return irfft(jacobianh(a, b, g), g.nx)
  else
   return ifft(jacobianh(a, b, g))
  end
end


"""
    radialspectrum(ah, g; nr=nothing, nθ=nothing, refinement=4)

Returns aρ = ∫ ah(ρ,θ) ρ dρ dθ, the radial spectrum of ah known on the 
Cartesian wavenumber grid (k,l). 

aρ is found by intepolating ah onto a polar wavenumber grid (ρ,θ), and 
then integrating over θ to find aρ. The default resolution (n,m) for the 
polar wave number grid is n=refinement*maximum(nk,nl), 
m=refinement*maximum(nk,nl), where refinement=4 by default. If 
ah is in conjugate symmetric form only the upper half plane in θ is
represented on the polar grid.

"""
function radialspectrum(ah, g::TwoDGrid; n=nothing, m=nothing, refinement=4)

  if n == nothing; n = refinement*maximum([g.nk, g.nl]); end
  if m == nothing; m = refinement*maximum([g.nk, g.nl]); end

  if size(ah)[1] == g.nkr # conjugate symmetric form
    m = Int(m/2) # need half resolution for conjugate symmetric form
    θ = linspace(-π/2, π/2, m) # θ grid from k=0 to max(kr)
    ahsh = fftshift(ah, 2)     # shifted ah
    ksh = linspace(0, (g.nkr-1)*2π/g.Lx, g.nkr) 
  else                     
    θ = linspace(0, 2π, m) # θ grid
    ahsh = fftshift(ah, [1, 2]) # shifted ah
    ksh = linspace((-g.nk+2)*π/g.Lx, g.nk*π/g.Lx, g.nk)
  end

  lsh = linspace((-g.nl+2)*π/g.Ly, g.nl*π/g.Ly, g.nl)
  ρmax = maximum([maximum(g.k), maximum(g.l)])
  ρ = linspace(0, ρmax, n)

  itp = scale(interpolate(ahsh, BSpline(Linear()), OnGrid()), ksh, lsh)
  ahρθ = zeros(eltype(ahsh), (n, m))

  for i=2:n, j=1:m # ignore zeroth mode
    kk = ρ[i]*cos(θ[j])
    ll = ρ[i]*sin(θ[j])
    ahρθ[i, j] = itp[kk, ll]
  end

  # aρ = ∫ a(p, θ) r dθ
  if size(ah)[1] == g.nkr
    aρ = 2*sum(ρ.*ahρθ, 2)*(θ[2]-θ[1]) # multiply by 2 for conjugate symmetry
  else
    aρ = sum(ρ.*ahρθ, 2)*(θ[2]-θ[1])
  end
  aρ[1] = ah[1, 1] # zeroth mode

  ρ, aρ 
end

# Moments and cumulants
domainaverage(c, g) = g.dx*g.dy*sum(c)
moment_x(c, g, n) = g.dx*g.dy*sum(g.X.^n.*c)
moment_y(c, g, n) = g.dx*g.dy*sum(g.Y.^n.*c)

cumulant_1x(c, g) = g.dx*g.dy*sum(g.X.*c) / domainaverage(c, g)
cumulant_1y(c, g) = g.dx*g.dy*sum(g.Y.*c) / domainaverage(c, g)

cumulant_2x(c, g) = (g.dx*g.dy*sum((g.X-cumulant_1x(c, g)).^2.0.*c)
  / domainaverage(c, g))
cumulant_2y(c, g) = (g.dx*g.dy*sum((g.Y.-cumulant_1y(c, g)).^2.0.*c)
  / domainaverage(c, g))
