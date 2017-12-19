function fftloop!(a, ah,
    fftplan::Base.DFT.FFTW.cFFTWPlan{Complex{Float64}, -1, false, 2},
    ifftplan; nloops=100)

  for i = 1:nloops
    A_mul_B!(ah, fftplan, a)
    A_mul_B!(a, ifftplan, ah)
  end
  nothing
end

function fftloop!(a, ah,
  rfftplan::Base.DFT.FFTW.rFFTWPlan{Float64, -1, false, 2}, irfftplan;
  nloops=100)

  for i = 1:nloops
    A_mul_B!(ah, rfftplan, a)
    A_mul_B!(a, irfftplan, ah)
  end
  nothing
end


function testfftspeed(ns;
  nthreads=nothing, effort=FFTW.MEASURE, nloops=100, ffttype="complex")

  if nthreads == nothing
    nthreads = [1]
    while nthreads[end] < Sys.CPU_CORES
      push!(nthreads, 2*nthreads[end])
    end
  end

  srand(123)
  times = zeros(length(ns), length(nthreads))

  @printf("\ntesting %s fft speed...\n", ffttype)
  for (inth, nth) in enumerate(nthreads)

    # Should we be super conservative about making sure number of
    # threads is right?
    FFTW.set_num_threads(nth)
    ENV["MKL_NUM_THREADS"] = nth
    ENV["JULIA_NUM_THREADS"] = nth

    @printf "(%s fft) nthreads: %d | n=" ffttype nth

    for (i, n) in enumerate(ns)
      @printf "%d... " n

      if ffttype == "complex"
        a = rand(n, n) + im*rand(n, n)
        ah = fft(a)
         fftplan =  plan_fft(Array{Complex{Float64},2}(n, n); flags=effort)
        ifftplan = plan_ifft(Array{Complex{Float64},2}(n, n); flags=effort)
      elseif ffttype == "real"
        a = rand(n, n)
        ah = rfft(a)
        fftplan = plan_rfft(
          Array{Float64,2}(n, n); flags=effort)
        ifftplan = plan_irfft(
          Array{Complex{Float64},2}(Int(n/2+1), n), n; flags=effort)
      end

      # Compile
      fftloop!(a, ah, fftplan, ifftplan; nloops=1)

      # Run
      times[i, inth] = @elapsed fftloop!(a, ah, fftplan, ifftplan;
        nloops=nloops)

    end
    @printf "done.\n"
  end
  ns, nthreads, times
end


"""
Prettily print the results of an fft test across arrays of size ns and
nthreads.
"""
function printresults(nthreads, ns, nloops, times, ffttype)

  # Header
  results = @sprintf("\n*** %s fft results (%d loops) ***\n\n",
    ffttype, nloops)
  results *= @sprintf("threads | n:")
  for n in ns
    results *= @sprintf("% 9d | ", n)
  end

  # Divider
  results *= "\n"
  results *= "----------------------------------------------------------------"
  results *= "---------------\n"

  # Body
  for (ith, nth) in enumerate(nthreads)
    results *= @sprintf("% 7d |   ", nth)
    for (i, n) in enumerate(ns)
      results *= @sprintf("% 9.5f | ", times[i, ith])
    end
    results *= "\n"
  end

  println(results)
  nothing
end


# Default test params
nloops = 200
testns = [64, 128, 256, 512, 1024, 2048]

# Real and complex fft speed tests
ns, nthreads, ctimes = testfftspeed(testns, ffttype="complex", nloops=nloops)
ns, nthreads, rtimes = testfftspeed(testns, ffttype="real", nloops=nloops)

printresults(nthreads, ns, nloops, ctimes, "complex")
printresults(nthreads, ns, nloops, rtimes, "real")
