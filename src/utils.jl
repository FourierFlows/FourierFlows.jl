__precompile__()

import SpecialFunctions

"Return the fftwavenumber vector with length n and domain size L."
fftwavenums(n::Int; L=1.0) = 2.0*pi/L*cat(1, 0:n/2, -n/2+1:-1)


"Fast loop-based dealiasing method for complex, spectral-space vars."
function dealias!(a::Array{Complex{Float64}, 2}, g::AbstractGrid)
  if size(a)[1] == g.nk       # Transform of a complex var
    for j in g.lderange
      @simd for i in g.kderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end
  else                        # Transform of a real var
    for j in g.lderange
      @simd for i in g.krderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end
  end
end


"Dealiasing method for 3-arrays with broadcasting for third dimension."
function dealias!(a::Array{Complex{Float64}, 3}, g::AbstractGrid)
  if size(a)[1] == g.nk       # Transform of a complex var
    for j in g.lderange
      @simd for i in g.kderange
        @inbounds @views @. a[i, j, :] = Complex{Float64}(0.0)
      end
    end
  else                        # Transform of a real var
    for j in g.lderange
      @simd for i in g.krderange
        @inbounds @views @. a[i, j, :] = Complex{Float64}(0.0)
      end
    end
  end
end


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
    phih .*= rms ./ sqrt.(sum(abs.(phih).^2.0)) 
    phi = real.(ifft(phih))
  else
    phi = real.(ifft(phih))
    phi = maxval * phi / maximum(abs(phi))
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
  q0 = 2*Ue*k/SpecialFunctions.besselj(0, k*R)
  println("q0", q0)

  r = sqrt.((g.X-xc).^2.0 + (g.Y-yc).^2.0)
  q = q0 * SpecialFunctions.besselj.(1, k*r) .* (g.Y-yc)./r

  q[r .> R] = 0.0

  return q
end
