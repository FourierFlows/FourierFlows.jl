using ParallelAccelerator

function dummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=100)
  for i = 1:nloops
    A_mul_B!(a, ifftplan, ah)

    @. b = a*c
    @. c = b*b + a
    @. b = c*c + a
    @. a = c*b + 2*a
     
    A_mul_B!(ah, fftplan, a)
  end
end


@acc function accdummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=100)
  for i = 1:nloops
    A_mul_B!(a, ifftplan, ah)

    @. b = a*c
    @. c = b*b + a
    @. b = c*c + a
    @. a = c*b + 2*a
     
    A_mul_B!(ah, fftplan, a)
  end
end


# Initialize random number generator
srand(123)
testns = 2.^(5:9)
effort = FFTW.MEASURE

nthreads1 = ENV["OMP_NUM_THREADS"]
nthreads2 = Threads.nthreads()
FFTW.set_num_threads(mininum([nthreads1, nthreads2]))

msg = "Testing ParallelAccelerator with either\n"
msg *= "  Threads.nthreads() = $nthreads1 or\n"
msg *= "  OMP_NUM_THREADS = $nthreads2."
println(msg)
    
for n in testns

   a = exp.(im*2π*rand(n, n))
   b = exp.(im*2π*rand(n, n))
   c = exp.(im*2π*rand(n, n))
  ah = ifft(a)

   fftplan = plan_fft(a; flags=effort);
  ifftplan = plan_ifft(ah; flags=effort);

  # Compile
  dummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=1)
  accdummyproblem!(a, b, c, ah, fftplan, ifftplan; nloops=1)

  # Run
  @printf "N: %5d^2, threads: % 2d, %s : " n nthreads "pure julia"
  @time dummyproblem!(a, b, c, ah, fftplan, ifftplan)

  @printf "N: %5d^2, threads: % 2d, %s : " n nthreads "parallel accelerator"
  @time accdummyproblem!(a, b, c, ah, fftplan, ifftplan)

end
