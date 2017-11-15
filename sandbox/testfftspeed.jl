function fftloop!(a, ah, fftplan, ifftplan; nloops=100)
  for i = 1:nloops
    A_mul_B!(ah, fftplan, a);
    A_mul_B!(a, ifftplan, ah);
  end
end

function fftloop!(a, ah; nloops=100)
  for i = 1:nloops
    fft!(a)
    ifft!(ah)
  end
end


# Initialize random number generator
srand(123)
testnthreads = 1:Sys.CPU_CORES
testns = 2.^(5:11)
effort = FFTW.MEASURE
inplace = false

for n in testns
  for nthreads in testnthreads

     a = exp.(im*2π*rand(n, n))
    ah = ifft(a)

    FFTW.set_num_threads(nthreads)

    if inplace 
      # Compile
      fftloop!(a, ah; nloops=1)
      # Run
      @printf "N: %5d^2, threads: %d, %24s : " n nthreads "in-place"
      @time fftloop!(a, ah)
    else
       fftplan = plan_fft(a; flags=effort);
      ifftplan = plan_ifft(ah; flags=effort);
      # Compile
      fftloop!(a, ah, fftplan, ifftplan; nloops=1)
      # Run
      @printf "N: %5d^2, threads: % 2d, %s : " n nthreads "planned"
      @time fftloop!(a, ah, fftplan, ifftplan)
    end
  end
  println()
end
