# import FourierFlows.TwoDTurb

function integralsquare(func, grid)
    sum(abs2.(func))*grid.dx*grid.dy
end

function testparsevalsum(func, grid; realvalued=true)
    # compute the integral in physical space
    integral = integralsquare(func, grid)

    # compute the integral in wavenumber space using Parseval's theorem
    if realvalued==true
        funch = rfft(func)
    elseif realvalued==false
        funch = fft(func)
    end
    parsevalsum = FourierFlows.parsevalsum(abs2.(funch), grid)

    # compare the two results
    isapprox(integral, parsevalsum; atol=1.0e-14)
end

function testparsevalsum2(func, grid; realvalued=true)
    # compute the integral in physical space
    integral = integralsquare(func, grid)

    # compute the integral in wavenumber space using Parseval's theorem
    if realvalued==true
        funch = rfft(func)
    elseif realvalued==false
        funch = fft(func)
    end
    parsevalsum2 = FourierFlows.parsevalsum2(funch, grid)

    # compare the two results
    isapprox(integral, parsevalsum2; atol=1.0e-14)
end


"""
Compute the J(a,b) and compare with analytic_answer.
"""
function testjacobian(a, b, analytic_answer, grid)
    isapprox(norm(FourierFlows.jacobian(a, b, grid)), norm(analytic_answer); 
             atol=g.nx*g.ny*1e-14)
end

function testcreatearrays(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  FourierFlows.@createarrays T dims a2 b2
  a1 == a2 && b1 == b2
end


# -----------------------------------------------------------------------------
# Running the tests

# Test on a rectangular grid
nx, ny = 64, 128             # number of points
Lx, Ly = 2π, 3π             # Domain width
g = TwoDGrid(nx, Lx, ny, Ly)

x, y = g.X, g.Y
k0, l0 = 2π/Lx, 2π/Ly
sig = 0.5

# a real-valued gaussian
f1 = exp.(-(x.^2 + y.^2)/(2*sig^2))

# a complex-valued function
f2 = exp.( im*(2k0*x + 3l0*y.^2) ).*( exp.(-(x.^2 + y.^2)/(2sig^2))
                                        + 2im*exp.(-(x.^2 + y.^2)/(5sig^2)) )

# a sine wave
k1, l1 = 2*k0, 6*l0
k2, l2 = 3*k0, -3*l0
s1 = sin.(k1*x + l1*y)
s2 = sin.(k2*x + l2*y)

# the analytical expression for the Jacobian of s1 and s2
Js1s2 = (k1*l2-k2*l1)*cos.(k1*x + l1*y).*cos.(k2*x + l2*y)

@test testparsevalsum(f1, g; realvalued=true)   #real valued f with rfft
@test testparsevalsum(f1, g; realvalued=false)  #real valued f with fft
@test testparsevalsum(f2, g; realvalued=false)  #complex valued f with fft
@test testparsevalsum2(f1, g; realvalued=true)  #real valued f with rfft
@test testparsevalsum2(f1, g; realvalued=false) #real valued f with fft
@test testparsevalsum2(f2, g; realvalued=false) #complex valued f with fft

@test testjacobian(s1, s1, 0*s1, g) # test that J(a,a)=0
@test testjacobian(s1, s2, Js1s2, g) # test J(s1,s2)=Js1s2

@test testcreatearrays() # Test that @createarrays works properly.
