function broadcastmult!(A, x, F; nloops=1000)
  for i = 1:nloops
    @. A = x*F
  end
end

function regularmult!(A, X, F; nloops=1000)
  for i = 1:nloops
    @. A = X*F
  end
end

function broadcastmult!(A::Array{Complex{Float64}, 2}, 
  x::Array{Complex{Float64}, 1}, F::Array{Complex{Float64}, 2}; nloops=1000)
  for i = 1:nloops
    @. A = x*F
  end
end

function broadcastmult!(A::Array{Complex{Float64}, 2}, 
  x::Array{Float64, 1}, F::Array{Complex{Float64}, 2}; nloops=1000)
  for i = 1:nloops
    @. A = x*F
  end
end





  

function testbroadcastspeed(testns)
  srand(123)
  for n in testns

    x = reshape(rand(n), (n, 1))
    X = [ x[i] for i=1:n, j=1:n ]
    F = rand(n, n)
    A = zeros(F)

    xc = Array{Complex{Float64}}(x)
    Fc = Array{Complex{Float64}}(F)
    Ac = Array{Complex{Float64}}(A)

    # Compile
    broadcastmult!(A, x, F; nloops=1)
    regularmult!(A, X, F; nloops=1)

    # Run
    btime = @elapsed broadcastmult!(A, x, F)
    rtime = @elapsed regularmult!(A, X, F)

    @printf "N: %5d^2, broadcast: %.4f, 'regular': %.4f\n" n btime rtime

  end
  nothing
end

function testconvertspeed(testns)
  srand(123)
  for n in testns
    x = reshape(rand(n), (n, 1))
    X = [ x[i] for i=1:n, j=1:n ]
    F = rand(n, n)
    A = zeros(F)

    xc = Array{Complex{Float64}}(x)
    Fc = Array{Complex{Float64}}(F)
    Ac = Array{Complex{Float64}}(A)

    # Compile
    broadcastmult!(Ac, x, Fc; nloops=1)
    broadcastmult!(Ac, xc, Fc; nloops=1)

    cctime = @elapsed broadcastmult!(Ac, xc, Fc)
    rctime = @elapsed broadcastmult!(Ac, x, Fc)
    @printf "N: %5d^2, c-c bcast: %.4f, r-c bcast: %.4f\n" n cctime rctime
  end
  nothing
end

testns = 2.^(6:11)

println("Testing broadcast speed...")
testbroadcastspeed(testns)

println("Testing real-complex convert speed...")
testconvertspeed(testns)
