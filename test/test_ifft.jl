function test_ifft_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1b = irfft(deepcopy(f1hr), g.nx)
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_mul_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    f1b = Array{Float64}(undef, g.nx)
    ldiv!(f1b, g.rfftplan, deepcopy(f1hr))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_ifft_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = real(ifft(f1h))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = irfft(deepcopy(f1hr), g.nx)
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_irfft_mul_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f1b = Array{Float64}(undef, g.nx, g.ny)
    ldiv!(f1b, g.rfftplan, deepcopy(f1hr))
    isapprox(f1, f1b, rtol=rtol_fft)
end

function test_ifft_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b3 = irfft(deepcopy(f2hr), g.nx)
    isapprox(f2, f2b3, rtol=rtol_fft)
end

function test_irfft_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b4 = real(ifft(f2h))
    isapprox(f2, f2b4, rtol=rtol_fft)
end

function test_irfft_mul_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    f2b = Array{Float64}(undef, g.nx, g.ny)
    ldiv!(f2b, g.rfftplan, deepcopy(f2hr))
    isapprox(f2, f2b, rtol=rtol_fft)
end
