# -----------------------------------------------------------------------------
# FFT's TEST FUNCTIONS

function test_fft_cosmxcosny(g)
    # test fft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1h-f1h_th)/norm(f1h_th) < 1e-12
end

function test_rfft_cosmxcosny(g)
    # test rfft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1hr-f1hr_th)/norm(f1hr_th) < 1e-12
end

function test_rfft_AmulB_cosmxcosny(g)
    # test rfft with A_mul_B for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f1hr_mul-f1hr_th)/norm(f1hr_th) < 1e-12
end

function test_fft_sinmxny(g)
    # test fft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2h-f2h_th)/norm(f2h_th) < 1e-12
end

function test_rfft_sinmxny(g)
    # test rfft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2hr-f2hr_th)/norm(f2hr_th) < 1e-12
end

function test_rfft_AmulB_sinmxny(g)
    # test rfft with A_mul_B for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    norm(f2hr_mul-f2hr_th)/norm(f2hr_th) < 1e-12
end

# -----------------------------------------------------------------------------
# Running the tests

# Test square grid
nx = 32               # number of points
Lx = 2.0π             # Domain width
g = TwoDGrid(nx, Lx)


@test test_fft_cosmxcosny(g)
@test test_fft_cosmxcosny(g)
@test test_rfft_cosmxcosny(g)
@test test_rfft_AmulB_cosmxcosny(g)
@test test_fft_sinmxny(g)
@test test_rfft_sinmxny(g)
@test test_rfft_AmulB_sinmxny(g)
