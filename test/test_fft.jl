function create_testfuncs(g)
    if nx <= 8; error("This test requires grids with nx > 8."); end

    m, n = 5, 2

    # Fundamental wavenumbers
    k0 = g.k[2]
    l0 = g.l[2]

    # Test functions
    f1 = @. cos(m*k0*g.X) * cos(n*l0*g.Y)
    f2 = @. sin(m*k0*g.X + n*l0*g.Y)

    # Their fft's and rfft's
    f1h = fft(f1)
    f2h = fft(f2)
    f1hr = rfft(f1)
    f2hr = rfft(f2)

    f1hr_mul = zeros(Complex{eltype(g.X)}, size(g.Kr))
    A_mul_B!(f1hr_mul, g.rfftplan, f1)

    f2hr_mul = zeros(Complex{eltype(g.X)}, size(g.Kr))
    A_mul_B!(f2hr_mul, g.rfftplan, f2)

    # Theoretical results for fft and rfft
    f1h_th = zeros(eltype(f1h), size(f1h))
    f2h_th = zeros(eltype(f2h), size(f2h))
    f1hr_th = zeros(eltype(f1hr), size(f1hr))
    f2hr_th = zeros(eltype(f2hr), size(f2hr))

    for j in 1:g.nl, i in 1:g.nk
      if ( abs(real(g.K[i, j])) == m*k0 && abs(real(g.L[i, j])) == n*l0 )
        f1h_th[i, j] = - g.nx*g.ny/4
      end
      if ( real(g.K[i, j]) == m*k0 && real(g.L[i, j]) == n*l0 )
        f2h_th[i, j] = -g.nx*g.ny/2
      elseif ( real(g.K[i, j]) == -m*k0 && real(g.L[i, j]) == -n*l0 )
        f2h_th[i, j] = g.nx*g.ny/2
      end
    end
    f2h_th = -im*f2h_th;

    for j in 1:g.nl, i in 1:g.nkr
      if ( abs(g.Kr[i, j])==m*k0 && abs(g.L[i, j])==n*l0 )
        f1hr_th[i, j] = - g.nx*g.ny/4
      end
      if ( real(g.Kr[i, j]) == m*k0 && real(g.L[i, j]) == n*l0 )
        f2hr_th[i, j] = -g.nx*g.ny/2
      elseif ( real(g.Kr[i, j]) == -m*k0 && real(g.L[i, j]) == -n*l0 )
        f2hr_th[i, j] = g.nx*g.ny/2
      end
    end
    f2hr_th = -im*f2hr_th;

    return f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th
end

# FFT test functions
function test_fft_cosmxcosny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1h-f1h_th)/norm(f1h_th) < 1e-12
end

function test_rfft_cosmxcosny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1hr-f1hr_th)/norm(f1hr_th) < 1e-12
end

function test_rfft_AmulB_cosmxcosny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1hr_mul-f1hr_th)/norm(f1hr_th) < 1e-12
end

function test_fft_sinmxny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2h-f2h_th)/norm(f2h_th) < 1e-12
end

function test_rfft_sinmxny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2hr-f2hr_th)/norm(f2hr_th) < 1e-12
end

function test_rfft_AmulB_sinmxny(g)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2hr_mul-f2hr_th)/norm(f2hr_th) < 1e-12
end

# Run tests on square two-dimensional grid
nx = 32               # Number of points
Lx = 2π               # Domain width
g = TwoDGrid(nx, Lx)

@test test_fft_cosmxcosny(g)
@test test_fft_cosmxcosny(g)
@test test_rfft_cosmxcosny(g)
@test test_rfft_AmulB_cosmxcosny(g)
@test test_fft_sinmxny(g)
@test test_rfft_sinmxny(g)
@test test_rfft_AmulB_sinmxny(g)
