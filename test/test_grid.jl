include("../src/fourierflows.jl")

using Base.Test
using FourierFlows


# -----------------------------------------------------------------------------
# X-Y GRID TEST FUNCTIONS

function testXgrid(g)
    # test if the X grid is actually created
    maximum(abs.(g.X)) > 1e-10
end

function testYgrid(g)
    # test if the Y grid is actually created
    maximum(abs.(g.Y)) > 1e-10
end

# -----------------------------------------------------------------------------
# K-L GRID TEST FUNCTIONS
function testKgrid(g)
    test = 0
    for i in 1:g.nk
      test += abs(maximum(g.K[i, :])-minimum(g.K[i, :]));
    end
    test<1e-12  #if true K grid is oriented correctly
end

function testKrgrid(g)
    test = 0
    for i in 1:g.nkr
      test += abs(maximum(g.Kr[i, :])-minimum(g.Kr[i, :]));
    end
    test<1e-12  #if true Kr grid is oriented correctly
end

function testLgrid(g)
    test = 0
    for i in 1:g.nl
      test += abs(maximum(g.L[:, i])-minimum(g.L[:, i]));
    end
    test<1e-12  #if true L grid is oriented correctly
end

# -----------------------------------------------------------------------------

function create_testfuncs(g)
    if nx/2<=4
        error("asdf")
    end

    m = 5;
    n = 2;

    # the fundumental wavenumbers for this particular grid
    k0 = g.k[2, 1]
    l0 = g.l[1, 2]

    # some test functions
    f1   = cos.(m*k0*g.X).*cos.(n*l0*g.Y);
    f2   = sin.(m*k0*g.X + n*l0*g.Y);
    # and their fft's and rfft's
    f1h  = fft(f1);
    f2h  = fft(f2);
    f1hr = rfft(f1);
    f2hr = rfft(f2);

    f1hr_mul = Array{Complex128}(g.nkr, g.nl)
    A_mul_B!( f1hr_mul, g.rfftplan, f1 )

    f2hr_mul = Array{Complex128}(g.nkr, g.nl)
    A_mul_B!( f2hr_mul, g.rfftplan, f2 )

    ############################################################################
    # the theoretical values of the fft's and rfft's
    f1h_th   = zeros(size(f1h));  f2h_th   = zeros(size(f2h));
    f1hr_th  = zeros(size(f1hr)); f2hr_th  = zeros(size(f2hr));

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
    ############################################################################

    return f1, f2, f1h, f2h, f1hr, f2hr, f1hr_mul, f2hr_mul,
            f1h_th, f1hr_th, f2h_th, f2hr_th
end

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


@testset "Grid Tests" begin
  @test testXgrid(g)
  @test testYgrid(g)
  @test testKgrid(g)
  @test testKrgrid(g)
  @test testLgrid(g)
end

@testset "FFT Tests" begin
  @test test_fft_cosmxcosny(g)
  @test test_fft_cosmxcosny(g)
  @test test_rfft_cosmxcosny(g)
  @test test_rfft_AmulB_cosmxcosny(g)
  @test test_fft_sinmxny(g)
  @test test_rfft_sinmxny(g)
  @test test_rfft_AmulB_sinmxny(g)
end

@testset "IFFT Tests" begin
  @test test_ifft_cosmxcosny(g)
  @test test_ifft_cosmxcosny(g)
  @test test_irfft_cosmxcosny(g)
  @test test_irfft_AmulB_cosmxcosny(g)
  @test test_ifft_sinmxny(g)
  @test test_irfft_sinmxny(g)
  @test test_irfft_AmulB_sinmxny(g)
end
