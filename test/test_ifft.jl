function test_ifft_cosmx(g::OneDGrid)
    # test fft for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_cosmx(g::OneDGrid)
    # test irfft for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = irfft(f1hr_c, nx)
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_mul_cosmx(g::OneDGrid)
    # test irfft with ldiv! for cos(mx+φ)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = Array{Float64}(undef, g.nx)
    ldiv!( f1b, g.rfftplan, f1hr_c )
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_ifft_cosmxcosny(g::TwoDGrid)
    # test fft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_cosmxcosny(g::TwoDGrid)
    # test irfft for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = irfft(f1hr_c, nx)
    isapprox(f1, f1b, rtol=rtol_fft)
end

        function test_irfft_mul_cosmxcosny(g::TwoDGrid)
    # test irfft with ldiv! for cos(mx)cos(ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1hr_c = deepcopy(f1hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f1b = Array{Float64}(undef, g.nx, g.ny)
    ldiv!( f1b, g.rfftplan, f1hr_c )
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_ifft_sinmxny(g::TwoDGrid)
    # test fft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b3 = irfft(f2hr_c, g.nx)
    isapprox(f2, f2b3, rtol=rtol_fft)
end

function test_irfft_sinmxny(g::TwoDGrid)
    # test irfft for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b4 = real(ifft(f2h))
    isapprox(f2, f2b4, rtol=rtol_fft)
end

function test_irfft_mul_sinmxny(g::TwoDGrid)
    # test irfft with ldiv! for sin(mx+ny)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2hr_c = deepcopy(f2hr)
    # deepcopy is needed because FFTW irfft messes up with input!
    f2b = Array{Float64}(undef, g.nx, g.ny)
    ldiv!(f2b, g.rfftplan, f2hr_c)
    isapprox(f2, f2b, rtol=rtol_fft)
end
