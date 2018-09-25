"""
    cxeltype(a)

Returns Complex{eltype(a)} if eltype(a) <: Real; eltype(a) otherwise.
"""
cxeltype(a) = eltype(a) <: Real ? Complex{eltype(a)} : eltype(a)

const steppers = [
  "ForwardEuler",
  "FilteredForwardEuler",
  "AB3",
  "FilteredAB3",
  "RK4",
  "FilteredRK4",
  "ETDRK4",
  "FilteredETDRK4",
]

const filteredsteppers = [
  "FilteredForwardEuler",
  "FilteredAB3",
  "FilteredRK4",
  "FilteredETDRK4",
]

"""
    autoconstructtimestepper(stepper, dt, sol, g=ZeroDGrid())

Returns a time-stepper type defined by the prefix `stepper`, timestep `dt`
solution `sol` (used to construct variables 
with identical type and size as the solution vector), and grid `g`.
"""
function autoconstructtimestepper(stepper, dt, sol, g::AbstractGrid=ZeroDGrid())
  fullsteppername = Symbol(stepper, :TimeStepper)
  tsexpr = occursin("Filtered", stepper) ?
      Expr(:call, fullsteppername, dt, sol, g) : Expr(:call, fullsteppername, dt, sol)
  eval(tsexpr)
end

"""
    autoconstructtimestepper(stepper, dt, solc, solr, g=ZeroDGrid())

Returns a time-stepper type defined by the prefix `stepper`, timestep `dt`
complex solution `solc` and real solution `solr` (used to construct variables 
with identical type and size as the solution vector), and grid `g`.

"""
function autoconstructtimestepper(stepper, dt, solc, solr, g::AbstractGrid=ZeroDGrid())
  fullsteppername = Symbol(stepper, :TimeStepper)
  tsexpr = occursin("Filtered", stepper) ?
    Expr(:call, fullsteppername, dt, solc, solr, g) : Expr(:call, fullsteppername, dt, solc, solr)
  eval(tsexpr)
end

"""
    @createarrays T dims a b c...

Create arrays of all zeros with element type `T`, size `dims`, and global names
`a`, `b`, `c` (for example). An arbitrary number of arrays may be created.
"""
macro createarrays(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

"Returns an expression that defines a Composite Type of the AbstractVars variety."
function structvarsexpr(name, physfields, transfields; vardims=2, parent=:AbstractVars, T=Float64, arraytype=:Array)
  physexprs = [:($fld::$arraytype{T,$vardims}) for fld in physfields]
  transexprs = [:($fld::$arraytype{Complex{T},$vardims}) for fld in transfields]
  expr = :(struct $name{T} <: $parent; $(physexprs...); $(transexprs...); end)
end

"""
    structvarsexpr(name, fieldspecs; parent=nothing)

Returns an expression that defines a composite type whose fields are given by
the name::type pairs specifed by the tuples in fieldspecs. The convention is
name = fieldspecs[i][1] and type = fieldspecs[i][2] for the ith element of
fieldspecs.
"""
function structvarsexpr(name, fieldspecs; parent=nothing)
  # name = spec[1]; type = spec[2]
  # example: fieldspecs[1] = (:u, Array{Float64,2})
  fieldexprs = [ :( $(spec[1])::$(spec[2]) ) for spec in fieldspecs ]
  if parent == nothing
    expr = :(struct $name{T}; $(fieldexprs...); end)
  else
    expr = :(struct $name{T} <: $parent; $(fieldexprs...); end)
  end
  expr
end


"""
    getfieldspecs(fieldnames, fieldtype)

Returns an array of (fieldname[i], fieldtype) tuples that can be given to the
function getstructexpr. This function makes it convenient to construct
fieldspecs for lists of variables of the same type.
"""
getfieldspecs(fieldnames, fieldtype) = collect(zip(fieldnames, [ fieldtype for name in fieldnames ]))

"""
    fftwavenums(n; L=1)

Return the fftwavenumber vector with length n and domain size L.
"""
fftwavenums(n::Int; L=1) = 2π/L*cat(0:n/2, -n/2+1:-1, dims=1)

"""
    rms(q)

Return the root-mean-square of an array.
"""
rms(q) = sqrt(mean(q.^2))

"""
    peakedisotropicspectrum(g, kpeak, E0; mask=mask, allones=false)

Generate a real and random two-dimensional vorticity field q(x, y) with
a Fourier spectrum peaked around a central non-dimensional wavenumber kpeak and
normalized so that its total energy is E0.
"""
function peakedisotropicspectrum(g::TwoDGrid, kpeak::Real, E0::Real; mask=ones(size(g.Kr)), allones=false)
  if g.Lx !== g.Ly
      error("the domain is not square")
  else
    k0 = kpeak*2π/g.Lx
    modk = sqrt.(g.KKrsq)
    psik = zeros(g.nk, g.nl)
    psik =  (modk.^2 .* (1 .+ (modk/k0).^4)).^(-0.5)
    psik[1, 1] = 0.0
    psih = (randn(g.nkr, g.nl)+im*randn(g.nkr, g.nl)).*psik
    if allones; psih = psik; end
    psih = psih.*mask
    Ein = real(sum(g.KKrsq.*abs2.(psih)/(g.nx*g.ny)^2))
    psih = psih*sqrt(E0/Ein)
    q = -irfft(g.KKrsq.*psih, g.nx)
  end
end



"""
    lambdipole(Ue, R, g; center=(x0, y0))

Return a 2D vorticity field corresponding to the Lamb Dipole with
strength Ue, radius R, and centered around
(xc, yc)=center. The default value of 'center' is the middle of the grid.
"""
function lambdipole(Ue::Real, R::Real, g::TwoDGrid; center=(nothing, nothing))
  xc, yc = center == (nothing, nothing) ? (mean(g.x), mean(g.y)) : (center[1], center[2])

  k = 3.8317059702075123156 / R # dipole wavenumber for radius R in terms of first zero of besselj
  q0 = -2Ue*k/besselj(0, k*R)

  r = @. sqrt((g.x-xc)^2 + (g.y-yc)^2)
  q = @. q0*besselj(1, k*r)*(g.y-yc)/r

  @. q[r == 0.0] = 0.0 # just in case.
  @. q[r > R] = 0.0

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
    radialspectrum(ah, g; n=nothing, m=nothing, refinement=2)

Returns `aρ = ∫ ah(ρ,θ) ρ dρ dθ`, the radial spectrum of `ah` known on the
Cartesian wavenumber grid (k,l).

`aρ` is found by intepolating `ah` onto a polar wavenumber grid (ρ,θ), and
then integrating over `θ` to find `aρ`. The default resolution (n,m) for the
polar wave number grid is `n=refinement*maximum(nk, nl),
m=refinement*maximum(nk, nl)`, where `refinement=2` by default. If
`ah` is in conjugate symmetric form only the upper half plane in `θ` is
represented on the polar grid.
"""
function radialspectrum(ah, g::TwoDGrid; n=nothing, m=nothing, refinement=2)

  n = n == nothing ? refinement*maximum([g.nk, g.nl]) : n 
  m = m == nothing ? refinement*maximum([g.nk, g.nl]) : m 

  # Calcualte shifted k and l
  lshift = range(-g.nl/2+1, stop=g.nl/2, length=g.nl)*2π/g.Ly

  if size(ah)[1] == g.nkr # conjugate symmetric form
    m = Int(m/2)                         # => half resolution in θ
    θ = range(-π/2, stop=π/2, length=m)  # θ-grid from k=0 to max(kr)
    ahshift = fftshift(ah, 2)            # shifted ah
    kshift = range(0, stop=g.nkr-1, length=g.nkr)*2π/g.Lx
  else # ordinary form
    θ = range(0, stop=2π, length=m)      # θ grid
    ahshift = fftshift(ah, [1, 2])       # shifted ah
    kshift = range(-g.nk/2+1, stop=g.nk/2, length=g.nk)*2π/g.Lx
  end

  # Interpolator for ah
  itp = scale(interpolate(ahshift, BSpline(Linear())), kshift, lshift)

  # Get radial wavenumber vector
  ρmax = minimum([(g.nk/2-1)*2π/g.Lx, (g.nl/2-1)*2π/g.Ly])
  ρ = range(0, stop=ρmax, length=n)

  # Interpolate ah onto fine grid in (ρ,θ).
  ahρθ = zeros(eltype(ahshift), (n, m)) 

  for i=2:n, j=1:m # ignore zeroth mode
    kk = ρ[i]*cos(θ[j])
    ll = ρ[i]*sin(θ[j])
    ahρθ[i, j] = itp(kk, ll)
  end

  # ahρ = ρ ∫ ah(ρ,θ) dθ  =>  Ah = ∫ ahρ dρ = ∫∫ ah dk dl
  dθ = θ[2]-θ[1]
  if size(ah)[1] == g.nkr
    ahρ = 2ρ.*sum(ahρθ, dims=2)*dθ # multiply by 2 for conjugate symmetry
  else
    ahρ = ρ.*sum(ahρθ, dims=2)*dθ
  end

  ahρ[1] = ah[1, 1] # zeroth mode

  ρ, ahρ
end

# Moments and cumulants
"Compute the average of `c` on the grid `g`."
domainaverage(c, g) = g.dx*g.dy*sum(c)/(g.Lx*g.Ly)

"Compute the `n`th x-moment of `c` on the grid `g`."
xmoment(c, g, n=1) = sum(g.X.^n.*c)/sum(c)

"Compute the `n`th y-moment of `c` on the grid `g`."
ymoment(c, g, n=1) = sum(g.Y.^n.*c)/sum(c)
