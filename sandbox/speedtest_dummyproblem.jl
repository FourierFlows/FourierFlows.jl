function dummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=100)
  for i = 1:nloops
    A_mul_B!(a, ifftplan, ah);

    @. b = a*c
    @. c = b*b + a
    @. b = c*c + a
    @. a = c*b + 2*a
     
    A_mul_B!(ah, fftplan, a);
  end
end


# Initialize random number generator
srand(123)
testnthreads = 2.^(0:Int(log2(Sys.CPU_CORES)))
testns = 2.^(5:9)
effort = FFTW.MEASURE
inplace = false

for n in testns
  for nthreads in testnthreads

     a = exp.(im*2π*rand(n, n))
     b = exp.(im*2π*rand(n, n))
     c = exp.(im*2π*rand(n, n))
    ah = ifft(a)

    FFTW.set_num_threads(nthreads)

     fftplan = plan_fft(a; flags=effort);
    ifftplan = plan_ifft(ah; flags=effort);
    # Compile
    dummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=1)
    # Run
    @printf "N: %5d^2, threads: % 2d, %s : " n nthreads "planned"
    @time dummyproblem!(a, b, c, ah, fftplan, ifftplan)
  end
  println()
end
