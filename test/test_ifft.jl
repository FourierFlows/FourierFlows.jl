# -----------------------------------------------------------------------------
# IFFT's TEST FUNCTIONS

function test_ifft_cosmxcosny(g)
    # test fft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    norm(f1-f1b)/norm(f1) < 1e-12
end

function test_irfft_cosmxcosny(g)
    # test irfft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = irfft(f1hr_c, nx)
    norm(f1-f1b)/norm(f1) < 1e-12
end

function test_irfft_AmulB_cosmxcosny(g)
    # test irfft with A_mul_B for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = Array{Float64}(g.nx, g.ny)
    A_mul_B!( f1b, g.irfftplan, f1hr_c )
    norm(f1-f1b)/norm(f1) < 1e-12
end

function test_ifft_sinmxny(g)
    # test fft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr);
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b3 = irfft(f2hr_c, g.nx);
    norm(f2-f2b3)/norm(f2) < 1e-12
end

function test_irfft_sinmxny(g)
    # test irfft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b4 = real(ifft(f2h));
    norm(f2-f2b4)/norm(f2) < 1e-12
end

function test_irfft_AmulB_sinmxny(g)
    # test irfft with A_mul_B for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b = Array{Float64}(g.nx, g.ny)
    A_mul_B!( f2b, g.irfftplan, f2hr_c )
    norm(f2-f2b)/norm(f2) < 1e-12
end



# -----------------------------------------------------------------------------
# Running the tests

# Test square grid
nx = 32               # number of points
Lx = 2.0π             # Domain width
g = TwoDGrid(nx, Lx)


@test test_ifft_cosmxcosny(g)
@test test_ifft_cosmxcosny(g)
@test test_irfft_cosmxcosny(g)
@test test_irfft_AmulB_cosmxcosny(g)
@test test_ifft_sinmxny(g)
@test test_irfft_sinmxny(g)
@test test_irfft_AmulB_sinmxny(g)
