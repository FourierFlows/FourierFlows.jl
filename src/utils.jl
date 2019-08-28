"""
    innereltype(x)

Recursively determine the 'innermost' type in by the collection `x` (which may be, for example,
a collection of a collection).
"""
function innereltype(x)
  T = eltype(x)
  T <: AbstractArray ? innereltype(T) : return T
end

"""
    cxtype(T)

Returns `T` when `T` is `Complex`, or `Complex{T}` when `T` is `Real`.
"""
cxtype(::Type{T}) where T<:Number = T
cxtype(::Type{T}) where T<:Real = Complex{T}

"""
    fltype(T)

Returns `T` when `T<:AbstractFloat` or `Tf` when `T<:Complex{Tf}`.
"""
fltype(::Type{T})          where T<:AbstractFloat = T
fltype(::Type{Complex{T}}) where T<:AbstractFloat = T
fltype(T::Tuple) = fltype(T[1])

cxeltype(x) = cxtype(innereltype(x))
fleltype(x) = fltype(innereltype(x))

"""
    superzeros(T, A)

Returns an array like `A`, but full of zeros. If `innereltype(A)` can be promoted to `T`, then
the innermost elements of the array will have type `T`.
"""
superzeros(T, A::AbstractArray) = T(0)*A
superzeros(A::AbstractArray) = superzeros(innereltype(A), A)
superzeros(T, dims::Tuple) = eltype(dims) <: Tuple ? [ superzeros(T, d) for d in dims ] : zeros(T, dims)
superzeros(dims::Tuple) = superzeros(Float64, dims) # default
superzeros(T::Tuple, dims::Tuple) = [ superzeros(T[i], dims[i]) for i=1:length(dims) ]

"""
    @superzeros T a b c d...
    @superzeros T dims b c d...

Generate arrays `b, c, d...` with the super-dimensions of `a` and innereltype `T`.
"""
macro superzeros(T, ad, vars...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(var)) = superzeros($(esc(T)), $(esc(ad))); ) for var in vars])
  expr
end
supersize(a) = Tuple([size(ai) for ai in a])
supersize(a::Array{T}) where T<:AbstractArray = Tuple([size(ai) for ai in a])
supersize(a::Array{T}) where T<:Number = size(a)

macro createarrays(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

"""
    @zeros T dims a b c...

Create arrays of all zeros with element type `T`, size `dims`, and global names
`a`, `b`, `c` (for example). An arbitrary number of arrays may be created.
"""
macro zeros(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

Base.zeros(::CPU, T, dims) = zeros(T, dims)

"""
    devzeros(dev, T, dims)

Returns an array like `A` of type `T`, but full of zeros.
"""
devzeros(dev, T, dims) = zeros(dev, T, dims)


"""
    @devzeros dev T dims a b c...

Create arrays of all zeros with element type `T`, size `dims`, and global names
`a`, `b`, `c` (for example) on device `dev`.
"""
macro devzeros(dev, T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(dev))(), $(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end


"""
    varsexpression(name, fieldspecs; parent=:AbstractVars, typeparams=nothing)

    varsexpression(name, fieldspecs...; parent=:AbstractVars, typeparams=nothing)

    varsexpression(name, physicalfields, fourierfields; physicaltype=:Tp, fouriertype=:Tf,
                   parent=:AbstractVars, typeparams=[physicaltype, fouriertype])

Returns an expression that defines an `AbstractVars` type.
"""
function varsexpression(name, fieldspecs; parent=:AbstractVars,
                        typeparams::Union{Nothing,Symbol,Array{Symbol,1}}=nothing)
  if typeparams == nothing
    signature = name
  else
    try
      signature = Expr(:curly, name, typeparams...)
    catch
      signature = Expr(:curly, name, typeparams) # only one typeparam given?
    end
  end

  fieldexprs = [ :( $(spec[1])::$(spec[2]); ) for spec in fieldspecs ]

  :(struct $signature <: $parent; $(fieldexprs...); end)
end

function varsexpression(name, physicalfields, fourierfields; parent=:AbstractVars,  physicaltype=:Tp, fouriertype=:Tf,
                        typeparams=[physicaltype, fouriertype])
  physicalfieldspecs = getfieldspecs(physicalfields, physicaltype)
   fourierfieldspecs = getfieldspecs(fourierfields, fouriertype)
  varsexpression(name, cat(physicalfieldspecs, fourierfieldspecs; dims=1); parent=parent, typeparams=typeparams)
end


"""
    getfieldspecs(fieldnames, fieldtype)

Returns an array of (fieldname[i], fieldtype) tuples.
"""
getfieldspecs(fieldnames::AbstractArray, fieldtype) = [ (name, fieldtype) for name in fieldnames ]

"""
    parsevalsum2(uh, g)

Returns `∫|u|² = Σ|uh|²` on the grid `g`, where `uh` is the Fourier transform of `u`.
"""
function parsevalsum2(uh, g::TwoDGrid)
  if size(uh, 1) == g.nkr # uh is in conjugate symmetric form
    U = @views sum(abs2, uh[1, :])           # k=0 modes
    U += @views 2*sum(abs2, uh[2:end, :])    # sum k>0 modes twice
  else # count every mode once
    U = sum(abs2, uh)
  end
  norm = g.Lx*g.Ly/(g.nx^2*g.ny^2)    # normalization for dft
  norm*U
end

function parsevalsum2(uh, g::OneDGrid)
  if size(uh, 1) == g.nkr # uh is conjugate symmetric
    U = sum(abs2, uh[1])                  # k=0 modes
    U += @views 2*sum(abs2, uh[2:end])    # sum k>0 modes twice
  else # count every mode once
    U = sum(abs2, uh)
  end
  norm = g.Lx/g.nx^2 # normalization for dft
  norm*U
end

"""
    parsevalsum(uh, g)

Returns `real(Σ uh)` on the grid `g`.
"""
function parsevalsum(uh, g::TwoDGrid)
  if size(uh, 1) == g.nkr       # uh is conjugate symmetric
    U = sum(uh[1, :])           # k=0 modes
    U += 2*sum(uh[2:end, :])    # sum k>0 modes twice
  else # count every mode once
    U = sum(uh)
  end

  norm = g.Lx*g.Ly/(g.nx^2*g.ny^2) # weird normalization for dft
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


#TODO: delete structvarsexpr() but first make sure that GeophysicalFlows.jl does not require it
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

ArrayType(::CPU, T, dim) = Array{T, dim}
ArrayType(::CPU) = Array
