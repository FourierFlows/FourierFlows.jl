function create_testfuncs(g::OneDGrid{Tg,<:Array}) where Tg
    g.nx > 8 || error("nx must be > 8")
    m = 5
    k₀ = g.k[2] # Fundamental wavenumber

    # Analytic function
    φ = π/3
    f₁ = @. cos(m*k₀*g.x + φ)

    # Transform of f₁ three ways
     f₁h = fft(f₁)
    f₁hr = rfft(f₁)
    f₁hr_mul = zeros(Complex{eltype(g.x)}, g.nkr)
    mul!(f₁hr_mul, g.rfftplan, f₁)

    # Analytical values of the fft and rfft
     f₁h_analytical = zeros(Complex{Tg}, size(f₁h))
    f₁hr_analytical = zeros(Complex{Tg}, size(f₁hr))

    for i in 1:g.nk
      if abs(real(g.k[i])) == m*k₀
        f₁h_analytical[i] = -exp(sign(real(g.k[i]))*im*φ)*g.nx/2
      end
    end

    for i in 1:g.nkr
      if abs(real(g.k[i])) == m*k₀
        f₁hr_analytical[i] = -exp(sign(real(g.kr[i]))*im*φ)*g.nx/2
      end
    end

    f₁, f₁h, f₁hr, f₁hr_mul, f₁h_analytical, f₁hr_analytical
end


function create_testfuncs(g::TwoDGrid{Tg,<:Array}) where Tg
    g.nx > 8 || error("nx must be > 8")
    x, y = gridpoints(g)
    nx, ny = g.nx, g.ny
    m, n = 5, 2
    k₀, l₀ = g.k[2], g.l[2]

    # Analytic functions
    f₁ = @. cos(m*k₀*x) * cos(n*l₀*y)
    f₂ = @. sin(m*k₀*x + n*l₀*y)

     f₁h = fft(f₁)
     f₂h = fft(f₂)
    f₁hr = rfft(f₁)
    f₂hr = rfft(f₂)

    f₁hr_mul = zeros(Complex{Tg}, (g.nkr, g.nl))
    f₂hr_mul = zeros(Complex{Tg}, (g.nkr, g.nl))
    mul!(f₁hr_mul, g.rfftplan, f₁)
    mul!(f₂hr_mul, g.rfftplan, f₂)

    # Theoretical results
     f₁h_analytical = zeros(Complex{Tg}, size(f₁h))
     f₂h_analytical = zeros(Complex{Tg}, size(f₂h))
    f₁hr_analytical = zeros(Complex{Tg}, size(f₁hr))
    f₂hr_analytical = zeros(Complex{Tg}, size(f₂hr))

    for j in 1:g.nl, i in 1:g.nk
      if ( abs(g.k[i]) == m*k₀ && abs(g.l[j]) == n*l₀ )
        f₁h_analytical[i, j] = -nx*ny/4
      end
      if ( g.k[i] == m*k₀ && g.l[j] == n*l₀ )
        f₂h_analytical[i, j] = -nx*ny/2
      elseif ( g.k[i] == -m*k₀ && g.l[j] == -n*l₀ )
        f₂h_analytical[i, j] = nx*ny/2
      end
    end
    f₂h_analytical = -im*f₂h_analytical

    for j in 1:g.nl, i in 1:g.nkr
      if ( abs(g.kr[i])==m*k₀ && abs(g.l[j])==n*l₀ )
        f₁hr_analytical[i, j] = - nx*ny/4
      end
      if ( g.kr[i] == m*k₀ && g.l[j] == n*l₀ )
        f₂hr_analytical[i, j] = -nx*ny/2
      end
    end
    f₂hr_analytical = -im*f₂hr_analytical

    f₁, f₂, f₁h, f₂h, f₁hr, f₂hr, f₁hr_mul, f₂hr_mul, f₁h_analytical, f₁hr_analytical, f₂h_analytical, f₂hr_analytical
end


function create_testfuncs(g::ThreeDGrid{Tg,<:Array}) where Tg
    g.nx > 8 || error("nx must be > 8")
    x, y, z = gridpoints(g)
    nx, ny, nz = g.nx, g.ny, g.nz
    mx, my, mz = 5, 2, 3
    k₀, l₀, m₀ = g.k[2], g.l[2], g.m[2]

    # Analytic functions
    f₁ = @. cos(mx*k₀*x) * cos(my*l₀*y) * cos(mz*m₀*z)
    f₂ = @. sin(mx*k₀*x + my*l₀*y + mz*m₀*z)

     f₁h = fft(f₁)
     f₂h = fft(f₂)
    f₁hr = rfft(f₁)
    f₂hr = rfft(f₂)

    f₁hr_mul = zeros(Complex{Tg}, (g.nkr, g.nl, g.nm))
    f₂hr_mul = zeros(Complex{Tg}, (g.nkr, g.nl, g.nm))
    mul!(f₁hr_mul, g.rfftplan, f₁)
    mul!(f₂hr_mul, g.rfftplan, f₂)

    # Theoretical results
     f₁h_analytical = zeros(Complex{Tg}, size(f₁h))
     f₂h_analytical = zeros(Complex{Tg}, size(f₂h))
    f₁hr_analytical = zeros(Complex{Tg}, size(f₁hr))
    f₂hr_analytical = zeros(Complex{Tg}, size(f₂hr))

    for k in 1:g.nm, j in 1:g.nl, i in 1:g.nk
      if ( abs(g.k[i]) == mx*k₀ && abs(g.l[j]) == my*l₀ && abs(g.m[k]) == mz*m₀ )
        f₁h_analytical[i, j, k] = nx*ny*nz/8
      end
      if ( g.k[i] == mx*k₀ && g.l[j] == my*l₀ && g.m[k] == mz*m₀ )
        f₂h_analytical[i, j, k] = nx*ny*nz/2
      elseif ( g.k[i] == -mx*k₀ && g.l[j] == -my*l₀ && g.m[k] == -mz*m₀ )
        f₂h_analytical[i, j, k] = -nx*ny*nz/2
      end
    end
    f₂h_analytical = -im*f₂h_analytical;

    for k in 1:g.nm, j in 1:g.nl, i in 1:g.nkr
      if ( abs(g.kr[i]) == mx*k₀ && abs(g.l[j]) == my*l₀ && abs(g.m[k]) == mz*m₀ )
        f₁hr_analytical[i, j, k] = nx*ny*nz/8
      end
      if ( g.kr[i] == mx*k₀ && g.l[j] == my*l₀ && g.m[k] == mz*m₀ )
        f₂hr_analytical[i, j, k] = nx*ny*nz/2
      end
    end
    f₂hr_analytical = -im*f₂hr_analytical;

    f₁, f₂, f₁h, f₂h, f₁hr, f₂hr, f₁hr_mul, f₂hr_mul, f₁h_analytical, f₁hr_analytical, f₂h_analytical, f₂hr_analytical
end


if CUDA.has_cuda()
  function create_testfuncs(g::OneDGrid{Tg, <:CuArray}) where Tg
    cpugrid = OneDGrid(g.nx, g.Lx)
    out = create_testfuncs(cpugrid)
    return map(x->CuArray(x), out)
  end

  function create_testfuncs(g::TwoDGrid{Tg, <:CuArray}) where Tg
    cpugrid = TwoDGrid(g.nx, g.Lx, g.ny, g.Ly)
    out = create_testfuncs(cpugrid)
    return map(x->CuArray(x), out)
  end

  function create_testfuncs(g::ThreeDGrid{Tg, <:CuArray}) where Tg
    cpugrid = ThreeDGrid(g.nx, g.Lx, g.ny, g.Ly, g.nz, g.Lz)
    out = create_testfuncs(cpugrid)
    return map(x->CuArray(x), out)
  end
end