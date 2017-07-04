include("framework.jl")
include("solver.jl")
include("turbonlysolver.jl")


module Utils

  using Domain, Framework, Solver
  export dealias!, set_q!, set_A!

  function set_q!(q::AbstractArray, v::Vars, p::Params, g::Grid)
      v.q = q
      A_mul_B!(v.qh, g.rfftplan, v.q)
      dealias!(v.qh, g)
  end

  function set_A!(A::AbstractArray, v::Vars, p::Params, g::Grid)
      v.A = Array{Complex128}(A)
      A_mul_B!(v.Ah, g.fftplan, v.A)
      dealias!(v.Ah, g)
  end

end



module PlotUtils

  using Domain, Framework, Solver
  using PyPlot

  function quickplot(nfig::Int, v::Vars, p::Params, g::Grid)

    fig, axs = subplots(nrows=2, ncols=1, sharex=true, sharey=true)
    axes(axs[1])
    imshow(v.q)

    axes(axs[2])
    imshow(sqrt.(v.u.^2.0+v.v.^2.0))

    pause(0.1)

  end

end



module Problems

  using Domain, Framework, TimeSteppers, Utils
  import Solver, TurbOnlySolver

  export test_problem

  function test_turb_problem(nx::Int, dt::Float64)

    Lx    = 2.0*pi              # Domain width
    f0    = 1.0                 # Inertial frequency
    Ro    = 1.0                 # Rossby number
    nuq   = 1e-4                # Vorticity hyperviscosity
    nuqn  = 2                   # Vorticity hyperviscosity order

    # Dummy parameters
    nua   = 1.0                 # Wave hyperviscosity
    nuan  = 2                   # Wave hyperviscosity order
    k0    = 1.0                 # Wavenumber
    sig   = 1.0
    kap   = 1.0

    # Set up problem
    FFTW.set_num_threads(Sys.CPU_CORES)

    g = Grid(nx, Lx)
    p = Params(f0, nuq, nuqn, g)
    v = Vars(p, g)

    qts = ForwardEulerTimeStepper(dt, p.LCq)
    # ats = ETDRK4TimeStepper(dt, p.LCa)

    # Normalize by rms vorticity
    srand(123)
    v.q = rand(g.nx, g.ny)
    v.qh = rfft(v.q)
    TurbOnlySolver.updatevars!(v, p, g)

    return g, p, v, qts
  end



  function test_problem(nx::Int)

    Lx   = 2.0*pi               # Domain width
    f0   = 1.0                  # Inertial frequency
    dt   = 0.20 * 2.0*pi/f0     # Time step

    nuq  = 1e-4                 # Vorticity hyperviscosity
    nuqn = 2                    # Vorticity hyperviscosity order
    nua  = 1e4                  # Wave hyperviscosity
    nuan = 2                    # Wave hyperviscosity order

    # Initial condition
    alpha = 1                   # Frequency parameter
    k0    = 32.0*pi / Lx        # Wavenumber
    Ro    = 0.1                 # Rossby number
    rq    = Lx/20.0             # Radius of Gaussian vortex

    # Computed parameters
    sig  = f0*sqrt(alpha+1)     # Wave frequency
    kap  = k0 / sqrt(alpha)     # Modewise wavenumber

    # Set up problem
    FFTW.set_num_threads(Sys.CPU_CORES)

    g = Grid(nx, Lx)
    p = Params(f0, sig, kap, nuq, nua, nuqn, nuan, g)
    v = Vars(p, g)

    #qts = ETDRK4TimeStepper(dt, p.LCq)
    qts = ETDRK4TimeStepper(dt, p.LCq)
    ats = ETDRK4TimeStepper(dt, p.LCa)

    # Set initial condition
    x0, y0 = g.Lx/8.0, g.Ly/8.0
    x1, y1 = -g.Lx/8.0, -g.Ly/8.0

    A0 = exp.(im*k0*g.X)
    q0 = rand(g.nx, g.ny)

    # Normalize by rms vorticity
    qrms = sqrt(mean(v.q.^2.0))/p.f0
    v.q = f0*Ro * v.q/qrms

    set_q!(q0, v, p, g)
    set_A!(A0, v, p, g)

    A_mul_B!(v.Ah, g.fftplan, v.A)
    A_mul_B!(v.qh, g.rfftplan, v.q)

    return g, p, v, qts, ats
  end



  function two_gaussians(nx::Int)

    Lx   = 1600e3               # Domain width
    f0   = 1e-4                 # Inertial frequency
    dt   = 0.01 * 2.0*pi/f0      # Time step

    nuq  = 1e8                  # Vorticity hyperviscosity
    nuqn = 4                    # Vorticity hyperviscosity order
    nua  = 1e24                 # Wave hyperviscosity
    nuan = 8                    # Wave hyperviscosity order

    # Initial condition
    alpha = 1                   # Frequency parameter
    k0    = 32.0*pi / Lx        # Wavenumber
    Ro    = 0.1                 # Rossby number
    rq    = Lx/20.0             # Radius of Gaussian vortex

    # Computed parameters
    sig  = f0*sqrt(alpha+1)     # Wave frequency
    kap  = k0 / sqrt(alpha)     # Modewise wavenumber

    # Set up problem
    FFTW.set_num_threads(Sys.CPU_CORES)

    g = Grid(nx, Lx)
    p = Params(f0, sig, kap, nuq, nua, nuqn, nuan, g)
    v = Vars(p, g)

    qts = ETDRK4TimeStepper(dt, p.LCq)
    ats = ETDRK4TimeStepper(dt, p.LCa)

    # Set initial condition
    x0, y0 = g.Lx/8.0, g.Ly/8.0
    x1, y1 = -g.Lx/8.0, -g.Ly/8.0

    v.A = exp.(im*k0*g.X)
    v.q = ( f0*Ro * exp.( -((g.X-x0).^2.0 + (g.Y-y0).^2.0) / (2.0*rq^2.0) )
      + f0*Ro * exp.( -((g.X-x1).^2.0 + (g.Y-y1).^2.0) / (2.0*rq^2.0) ) )

    set_q!(v.q, v, p, g)
    set_A!(v.q, v, p, g)

    updatevars!(v, p, g)

    return g, p, v, qts, ats
  end



  function two_gaussians(nx::Int, fdt::Float64)

    Lx   = 1600e3               # Domain width
    f0   = 1e-4                 # Inertial frequency
    dt   = fdt * 2.0*pi/f0

    nuq  = 1e8                  # Vorticity hyperviscosity
    nuqn = 4                    # Vorticity hyperviscosity order
    nua  = 0.0                  # Wave hyperviscosity
    nuan = 4                    # Wave hyperviscosity order

    # Initial condition
    alpha = 1                   # Frequency parameter
    k0    = 32.0*pi / Lx        # Wavenumber
    Ro    = 0.05                # Rossby number
    rq    = Lx/20.0             # Radius of Gaussian vortex

    # Computed parameters
    sig  = f0*sqrt(alpha+1)     # Wave frequency
    kap  = k0 / sqrt(alpha)     # Modewise wavenumber

    # Set up problem
    FFTW.set_num_threads(Sys.CPU_CORES)

    g = Grid(nx, Lx)
    p = Params(f0, sig, kap, nuq, nua, nuqn, nuan, g)
    v = Vars(p, g)

    qts = ForwardEulerTimeStepper(dt, p.LCq)
    ats = ForwardEulerTimeStepper(dt, p.LCa)
    #qts = ETDRK4TimeStepper(dt, p.LCq)
    #ats = ETDRK4TimeStepper(dt, p.LCa)

    # Set initial condition
    x0, y0 = g.Lx/8.0, g.Ly/8.0
    x1, y1 = -g.Lx/8.0, -g.Ly/8.0

    v.A = exp.(im*k0*g.X)
    v.q = ( f0*Ro * exp.( -((g.X-x0).^2.0 + (g.Y-y0).^2.0) / (2.0*rq^2.0) )
      + f0*Ro * exp.( -((g.X-x1).^2.0 + (g.Y-y1).^2.0) / (2.0*rq^2.0) ) )

    set_q!(v.q, v, p, g)
    set_A!(v.A, v, p, g)

    Solver.updatevars!(v, p, g)

    Solver.stepforward!(1, qts, ats, v, p, g)
    TurbOnlySolver.stepforward!(1, qts, v, p, g)

    return g, p, v, qts, ats
  end



end
