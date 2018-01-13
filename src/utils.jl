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
function getexpr_varstype(name, physfields, transfields; soldims=2, vardims=2, 
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






"Return the fftwavenumber vector with length n and domain size L."
fftwavenums(n::Int; L=1.0) = 2.0*pi/L*cat(1, 0:n/2, -n/2+1:-1)

""" Return the root-mean-square of an array. """
rms(q) = sqrt(mean(q.^2))


"""
Generate a real and random two-dimensional distribution phi(x, y) with
a Fourier spectrum peaked around a central non-dimensional wavenumber kpeak.
The spectrum is normalized either by setting the root-mean-square value of phi
with the keyword 'rms', or the maximum value of phi with the keyword 'maxval'.
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

"Alternative input form for peaked_isotropic_spectrum for square grids."
function peaked_isotropic_spectrum(nx::Int, npeak::Real; ord=4.0, rms=1.0,
  maxval=0.0)
  peaked_isotropic_spectrum((nx, nx), npeak; ord=ord, rms=rms, maxval=maxval)
end




""" Return a 2D vorticity field corresponding to the Lamb Dipole with
strength Ue, radius R, and wavenumber k, and centered around
(xc, yc)=center. The default value of 'center' is the middle of the grid."""
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





""" Return a vorticity field with magnitude q0, radius R, and center at
center[1], center[2] on a TwoDGrid g corresponding to a 'Gaussian vortex' with
Gaussian streamfunction. """
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



""" Return an array of random numbers on a TwoDGrid normalized to have a
specifed rms value. """
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
Returns the transform of the Jacobian of two fields a, b on the grid g.
"""
function jacobianh(a, b, g::TwoDGrid)
  # J(a, b) = dx(a b_y) - dy(a b_x)
  bh = fft(b)
  bx = ifft(im*g.K.*bh)
  by = ifft(im*g.L.*bh)
  im*g.K.*fft(a.*by)-im*g.L.*fft(a.*bx)
end



"""
Returns the Jacobian of a and b. Uses ifft (not irfft) so some remnant imaginary
part of order machine precision might remain in the end.
"""
function jacobian(a, b, g::TwoDGrid)
 ifft(jacobianh(a, b, g))
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


#=
M0(c, g)     = g.dx*g.dy*sum(c)
Mxn(c, g, n) = g.dx*g.dy*sum(g.X.^n.*c)
Myn(c, g, n) = g.dx*g.dy*sum(g.Y.^n.*c)

Cx1(c, g) = g.dx*g.dy*sum(g.X.*c) / M0(c, g)
Cy1(c, g) = g.dx*g.dy*sum(g.Y.*c) / M0(c, g)


Cx3(c, g) = g.dx*g.dy*sum((g.X-Cx1(c, g)).^3.0.*c) / M0(c, g)
Cy3(c, g) = g.dx*g.dy*sum((g.Y-Cy1(c, g)).^3.0.*c) / M0(c, g)

intx(II, g) = g.dx*sum(II)
inty(II, g) = g.dy*sum(II)

mxn(c, g, n) = g.dx*sum(g.Y.^n.*c, 1)
myn(c, g, n) = g.dy*sum(g.Y.^n.*c, 2)

cx1(c, g) = g.dx*sum(g.X.*c) ./ mxn(c, g, 0)
cy1(c, g) = g.dy*sum(g.Y.*c) ./ myn(c, g, 0)

delyc1(c, g) = g.Y - broadcast(*, cy1(c, g), ones(g.nx, g.ny))
cy2(c, g) = g.dy * sum(delyc1(c, g).*c, 2)./myn(c, g, 0)
=#
