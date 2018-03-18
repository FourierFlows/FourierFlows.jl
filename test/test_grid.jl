# -----------------------------------------------------------------------------
# X-Y GRID TEST FUNCTIONS

"Tests that grid spacing is as specified."
function testdx(g)
  sum(abs.(g.x[2:end]-g.x[1:end-1])-(g.nx-1)*g.dx) < 1e-15*g.nx
end


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

function create_testfuncs1D(g::OneDGrid)
    if nx/2<=4
        error("please create a larger grid for the tests")
    end

    m = 5

    # the fundumental wavenumbers for this particular grid
    k0 = g.k[2]

    # a test function
    φ = π/3
    f1   = cos.(m*k0*g.x + φ)
    # and its fft's and rfft's
    f1h  = fft(f1)
    f1hr = rfft(f1)

    f1hr_mul = Array{Complex128}(g.nkr)
    A_mul_B!( f1hr_mul, g.rfftplan, f1 )

    ############################################################################
    # the theoretical values of the fft's and rfft's
    f1h_th   = zeros(Complex{Float64}, size(f1h))
    f1hr_th  = zeros(Complex{Float64}, size(f1hr))

    for i in 1:g.nk
      if abs(real(g.k[i])) == m*k0
        f1h_th[i] = -exp(sign(real(g.k[i]))*im*φ)*g.nx/2
      end
    end

    for i in 1:g.nkr
      if abs(real(g.k[i]))==m*k0
        f1hr_th[i] = -exp(sign(real(g.kr[i]))*im*φ)*g.nx/2
      end
    end
    ############################################################################

    return f1, f1h, f1hr, f1hr_mul, f1h_th, f1hr_th
end


function create_testfuncs2D(g::TwoDGrid)
    if nx/2<=4
        error("please create a larger grid for the tests")
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
    f1h_th   = zeros(Complex{Float64}, size(f1h))
    f2h_th   = zeros(Complex{Float64}, size(f2h))
    f1hr_th  = zeros(Complex{Float64}, size(f1hr))
    f2hr_th  = zeros(Complex{Float64}, size(f2hr))

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
# Running the tests

# Test 1d grid
nx = 32
Lx = 2π
g1 = OneDGrid(nx, Lx)

@test testdx(g1)

# Test 2d grid
nx = 32  # number of points
Lx = 2π  # Domain width
g2 = TwoDGrid(nx, Lx)

@test testXgrid(g2)
@test testYgrid(g2)
@test testKgrid(g2)
@test testKrgrid(g2)
@test testLgrid(g2)
