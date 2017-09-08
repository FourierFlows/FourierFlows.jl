__precompile__()


fftwavenums(n::Int; L=1.0) = 2.0*pi/L*cat(1, 0:n/2, -n/2+1:-1)


# Fast loop-based dealiasing method for complex, spectral-space vars
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


# Dealiasing method for 3-arrays with broadcasting for third dimension.
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


function peaked_isotropic_spectrum(nkl::Tuple{Int, Int}, kpeak::Real; 
  ord=4.0, rms=1.0, maxval=0.0)
  # This function generates a real, 'random' two-dimensional 
  # distribution phi(x, y) with a Fourier spectrum that is peaked 
  # around the central non-dimensional wavenumber kpeak. The spectrum is
  # normalized either by setting the root-mean-square value of phi with the
  # keyword 'rms', or the maximum value of phi with the keyword 'maxval'.

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

function peaked_isotropic_spectrum(nx::Int, npeak::Real; ord=4.0, rms=1.0, 
  maxval=0.0)
  peaked_isotropic_spectrum((nx, nx), npeak; ord=ord, rms=rms, maxval=maxval)
end
