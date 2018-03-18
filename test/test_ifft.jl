# Requires the create_testfuncs function defined in test_fft.jl


# -----------------------------------------------------------------------------
# IFFT's TEST FUNCTIONS

tolerance = 1e-12

function test_ifft_cosmx(g::OneDGrid)
    # test fft for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_irfft_cosmx(g::OneDGrid)
    # test irfft for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = irfft(f1hr_c, nx)
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_irfft_AmulB_cosmx(g::OneDGrid)
    # test irfft with A_mul_B for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = Array{Float64}(g.nx)
    A_mul_B!( f1b, g.irfftplan, f1hr_c )
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_ifft_cosmxcosny(g::TwoDGrid)
    # test fft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_irfft_cosmxcosny(g::TwoDGrid)
    # test irfft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = irfft(f1hr_c, nx)
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_irfft_AmulB_cosmxcosny(g::TwoDGrid)
    # test irfft with A_mul_B for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = Array{Float64}(g.nx, g.ny)
    A_mul_B!( f1b, g.irfftplan, f1hr_c )
    norm(f1-f1b)/norm(f1) < tolerance
end

function test_ifft_sinmxny(g::TwoDGrid)
    # test fft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b3 = irfft(f2hr_c, g.nx)
    norm(f2-f2b3)/norm(f2) < tolerance
end

function test_irfft_sinmxny(g::TwoDGrid)
    # test irfft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b4 = real(ifft(f2h))
    norm(f2-f2b4)/norm(f2) < tolerance
end

function test_irfft_AmulB_sinmxny(g::TwoDGrid)
    # test irfft with A_mul_B for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b = Array{Float64}(g.nx, g.ny)
    A_mul_B!( f2b, g.irfftplan, f2hr_c )
    norm(f2-f2b)/norm(f2) < tolerance
end



# -----------------------------------------------------------------------------
# Running the tests

# Test 1D grid
nx = 32             # number of points
Lx = 2π             # Domain width
g1 = OneDGrid(nx, Lx)

# Test 2D rectangular grid
nx, ny = 32, 64     # number of points
Lx, Ly = 2π, 3π     # Domain width
g2 = TwoDGrid(nx, Lx, ny, Ly)

@test test_ifft_cosmx(g1)
@test test_irfft_cosmx(g1)
@test test_irfft_AmulB_cosmx(g1)

@test test_ifft_cosmxcosny(g2)
@test test_irfft_cosmxcosny(g2)
@test test_irfft_AmulB_cosmxcosny(g2)
@test test_ifft_sinmxny(g2)
@test test_irfft_sinmxny(g2)
@test test_irfft_AmulB_sinmxny(g2)
