function test_fft_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    isapprox(f1h, f1h_th, rtol=rtol_fft)
end

function test_rfft_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    isapprox(f1hr, f1hr_th, rtol=rtol_fft)
end

function test_rfft_mul_cosmx(g::OneDGrid)
    f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th = create_testfuncs(g)
    isapprox(f1hr_mul, f1hr_th, rtol=rtol_fft)
end

function test_fft_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f1h, f1h_th, rtol=rtol_fft)
end

function test_rfft_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f1hr, f1hr_th, rtol=rtol_fft)
end


function test_fft_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mulPk, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f2h, f2h_th, rtol=rtol_fft)
end

function test_rfft_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f2hr, f2hr_th, rtol=rtol_fft)
end

function test_rfft_mul_cosmxcosny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f1hr_mul, f1hr_th, rtol=rtol_fft)
end

function test_rfft_mul_sinmxny(g::TwoDGrid)
    f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul, f1h_th, f1hr_th, f2h_th, f2hr_th = create_testfuncs(g)
    isapprox(f2hr_mul, f2hr_th, rtol=rtol_fft)
end
