using
  FFTW,
  LinearAlgebra,
  Printf,
  JLD2,
  CUDA,
  Test

using
  FourierFlows,
  FourierFlows.Diffusion

using FourierFlows: parsevalsum2

using LinearAlgebra: mul!, ldiv!, norm

# the devices on which tests will run
devices = (CPU(),)
@has_cuda devices = (CPU(), GPU())

const rtol_output = 1e-12
const rtol_utils = 1e-13
const rtol_timesteppers = 1e-12

const steppers = [
  "ForwardEuler",
  "RK4",
  "ETDRK4",
  "AB3",
  "FilteredForwardEuler",
  "FilteredRK4",
  "FilteredETDRK4",
  "FilteredAB3"
]

# Run tests
include("createffttestfunctions.jl")

for dev in devices
  
  const rtol_fft = dev==CPU() ? 1e-12 : 1e-6
  
  @info "testing on " * string(typeof(dev))
  
  @time @testset "Grid tests" begin
    include("test_grid.jl")
    
    nx, Lx =  6, 2π
    ny, Ly =  8, 4π
    nz, Lz = 10, 3.0
    
    # Test 1D grid
    g₁ = OneDGrid(dev, nx, Lx)
    @test testnx(g₁, nx)
    @test testdx(dev, g₁)
    @test testdk(g₁)
    @test testx(g₁)
    @test testk(g₁)
    @test testkr(g₁)
    @test testdealias(g₁)
    @test testmakefilter(dev, g₁)

    # Test 2D rectangular grid
    g₂ = TwoDGrid(dev, nx, Lx, ny, Ly)
    @test testnx(g₂, nx)
    @test testny(g₂, ny)
    @test testdx(dev, g₂)
    @test testdy(dev, g₂)
    @test testdk(g₂)
    @test testdl(g₂)
    @test testx(g₂)
    @test testy(g₂)
    @test testk(g₂)
    @test testkr(g₂)
    @test testl(g₂)
    @test testgridpoints(dev, g₂)
    @test testdealias(g₂)
    @test testmakefilter(dev, g₂)
    
    # Test 3D parallelogram grid
    g₃ = ThreeDGrid(dev, nx, Lx, ny, Ly, nz, Lz)
    @test testnx(g₃, nx)
    @test testny(g₃, ny)
    @test testnz(g₃, nz)
    @test testdx(dev, g₃)
    @test testdy(dev, g₃)
    @test testdz(dev, g₃)
    @test testdk(g₃)
    @test testdl(g₃)
    @test testdm(g₃)
    @test testx(g₃)
    @test testy(g₃)
    @test testz(g₃)
    @test testk(g₃)
    @test testkr(g₃)
    @test testl(g₃)
    @test testm(g₃)
    @test testgridpoints(dev, g₃)
    @test testdealias(g₃)
    @test testmakefilter(dev, g₃)
    @test test_plan_flows_fftrfft(dev)
    
    # Test typed grids
    T = Float32
    @test test_plan_flows_fftrfft(dev, T=T)
    @test testtypedonedgrid(dev, nx, Lx; T=T)
    @test testtypedtwodgrid(dev, nx, Lx, ny, Ly; T=T)
    @test testtypedthreedgrid(dev, nx, Lx, ny, Ly, nz, Lz; T=T)
    
    # Test show() methods
    if dev==CPU()
      FT = "Float64"

      @test repr(g₁) == "OneDimensionalGrid\n  ├─────────── Device: "*FourierFlows.griddevice(g₁)*"\n  ├──────── FloatType: "*FT*"\n  ├────────── size Lx: 6.283185307179586\n  ├──── resolution nx: 6\n  ├── grid spacing dx: 1.0471975511965976\n  └─────────── domain: x ∈ [-3.141592653589793, 2.094395102393195]"
      @test repr(g₂) == "TwoDimensionalGrid\n  ├───────────────── Device: "*FourierFlows.griddevice(g₂)*"\n  ├────────────── FloatType: "*FT*"\n  ├────────── size (Lx, Ly): (6.283185307179586, 12.566370614359172)\n  ├──── resolution (nx, ny): (6, 8)\n  ├── grid spacing (dx, dy): (1.0471975511965976, 1.5707963267948966)\n  └───────────────── domain: x ∈ [-3.141592653589793, 2.094395102393195]\n                             y ∈ [-6.283185307179586, 4.71238898038469]"
      @test repr(g₃) == "ThreeDimensionalGrid\n  ├───────────────────── Device: "*FourierFlows.griddevice(g₃)*"\n  ├────────────────── FloatType: "*FT*"\n  ├────────── size (Lx, Ly, Lz): (6.283185307179586, 12.566370614359172, 12.566370614359172)\n  ├──── resolution (nx, ny, nz): (6, 8, 10)\n  ├── grid spacing (dx, dy, dz): (1.0471975511965976, 1.5707963267948966, 0.3)\n  └────────────────────  domain: x ∈ [-3.141592653589793, 2.094395102393195]\n                                 y ∈ [-6.283185307179586, 4.71238898038469]\n                                 z ∈ [-1.5, 1.2]"
    else
      FT = "Float32"
      
      @test repr(g₁) == "OneDimensionalGrid\n  ├─────────── Device: "*FourierFlows.griddevice(g₁)*"\n  ├──────── FloatType: "*FT*"\n  ├────────── size Lx: 6.2831855\n  ├──── resolution nx: 6\n  ├── grid spacing dx: 1.0471976\n  └─────────── domain: x ∈ [-3.1415927, 2.0943952]"
      @test repr(g₂) == "TwoDimensionalGrid\n  ├───────────────── Device: "*FourierFlows.griddevice(g₂)*"\n  ├────────────── FloatType: "*FT*"\n  ├────────── size (Lx, Ly): (6.2831855f0, 12.566371f0)\n  ├──── resolution (nx, ny): (6, 8)\n  ├── grid spacing (dx, dy): (1.0471976f0, 1.5707964f0)\n  └───────────────── domain: x ∈ [-3.1415927, 2.0943952]\n                             y ∈ [-6.2831855, 4.712389]"
      @test repr(g₃) == "ThreeDimensionalGrid\n  ├───────────────────── Device: "*FourierFlows.griddevice(g₃)*"\n  ├────────────────── FloatType: "*FT*"\n  ├────────── size (Lx, Ly, Lz): (6.2831855f0, 12.566371f0, 12.566371f0)\n  ├──── resolution (nx, ny, nz): (6, 8, 10)\n  ├── grid spacing (dx, dy, dz): (1.0471976f0, 1.5707964f0, 0.3f0)\n  └────────────────────  domain: x ∈ [-3.1415927, 2.0943952]\n                                 y ∈ [-6.2831855, 4.712389]\n                                 z ∈ [-1.5, 1.2]"
    end
  end

  @time @testset "FFT tests" begin
    include("test_fft.jl")

    # Test 1D grid
    nx = 32             # number of points
    Lx = 2π             # Domain width
    g1 = OneDGrid(dev, nx, Lx, T=Float32)
    
    @test test_fft_cosmx(g1)
    @test test_rfft_cosmx(g1)
    @test test_rfft_mul_cosmx(g1)

    # Test 2D rectangular grid
    nx, ny = 32, 64     # number of points
    Lx, Ly = 2π, 3π     # Domain width
    g2 = TwoDGrid(dev, nx, Lx, ny, Ly, T=Float32)

    @test test_fft_cosmxcosny(g2)
    @test test_rfft_cosmxcosny(g2)
    @test test_rfft_mul_cosmxcosny(g2)
    @test test_fft_sinmxny(g2)
    @test test_rfft_sinmxny(g2)
    @test test_rfft_mul_sinmxny(g2)

    # Test 3D parallelogram grid
    nx, ny, nz = 32, 30, 16     # number of points
    Lx, Ly, Lz = 2π, 3π, 4.0    # Domain width
    g3 = ThreeDGrid(dev, nx, Lx, ny, Ly, nz, Lz)

    @test test_fft_cosmxcosny(g3)
    @test test_rfft_cosmxcosny(g3)
    @test test_rfft_mul_cosmxcosny(g3)
    @test test_fft_sinmxny(g3)
    @test test_rfft_sinmxny(g3)
    @test test_rfft_mul_sinmxny(g3)
  end

  @time @testset "IFFT tests" begin
    include("test_ifft.jl")

    # Test 1D grid
    nx = 32             # number of points
    Lx = 2π             # Domain width
    g1 = OneDGrid(dev, nx, Lx)

    @test test_ifft_cosmx(g1)
    @test test_irfft_cosmx(g1)
    @test test_irfft_mul_cosmx(g1)

    # Test 2D rectangular grid
    nx, ny = 32, 64     # number of points
    Lx, Ly = 2π, 3π     # Domain width
    g2 = TwoDGrid(dev, nx, Lx, ny, Ly)

    @test test_ifft_cosmxcosny(g2)
    @test test_irfft_cosmxcosny(g2)
    @test test_irfft_mul_cosmxcosny(g2)
    @test test_ifft_sinmxny(g2)
    @test test_irfft_sinmxny(g2)
    @test test_irfft_mul_sinmxny(g2)

    # Test 3D parallelogram grid
    nx, ny, nz = 32, 30, 16     # number of points
    Lx, Ly, Lz = 2π, 3π, 4.0    # Domain width
    g3 = ThreeDGrid(dev, nx, Lx, ny, Ly, nz, Lz)

    @test test_ifft_cosmxcosny(g3)
    @test test_irfft_cosmxcosny(g3)
    @test test_irfft_mul_cosmxcosny(g3)
    @test test_ifft_sinmxny(g3)
    @test test_irfft_sinmxny(g3)
    @test test_irfft_mul_sinmxny(g3)
  end

  @time @testset "Timestepper tests" begin
    include("test_timesteppers.jl")
    for stepper in steppers
      @test constantdiffusiontest_stepforward(stepper, dev=dev)
      
      @test varyingdiffusiontest_stepforward(stepper, dev=dev)
      
      if FourierFlows.isexplicit(stepper)
        @test constantdiffusiontest_step_until(stepper, dev=dev)
      else
        @test_throws Exception constantdiffusiontest_step_until(stepper, dev=dev)
      end
    end
  end

  @time @testset "Utils tests" begin
    include("test_utils.jl")

    @test test_fltype()
    @test test_cxtype()
    @test test_innereltype()
    @test test_superzeros()
    @test test_supertuplezeros()
    @test test_supersize()
    @test test_arraytype(dev)
    @test test_arraytypeTdim(dev, Float32, 2)

    # Test on a rectangular grid
    nx, ny = 64, 128   # number of points
    Lx, Ly = 2π, 3π    # Domain width
    g = TwoDGrid(dev, nx, Lx, ny, Ly)
    x, y = gridpoints(g)
    k0, l0 = 2π/Lx, 2π/Ly

    # Real and complex-valued functions
    σ = 0.5
    f1 = @. exp(-(x^2 + y^2)/(2σ^2))
    f2 = @. (cos(2k0*x + 3l0*y^2) + im*sin(2k0*x + 3l0*y^2)) * (exp(-(x^2 + y^2)/(2σ^2)) + 2im*exp(-(x^2 + y^2)/(5σ^2)))

    # Sine/Exp waves
    k1, l1 = 2*k0,  6*l0
    k2, l2 = 3*k0, -3*l0
    sinkl1 = @. sin(k1*x + l1*y)
    sinkl2 = @. sin(k2*x + l2*y)
    expkl1 = @. cos(k1*x + l1*y) + im*sin(k1*x + l1*y)
    expkl2 = @. cos(k2*x + l2*y) + im*sin(k2*x + l2*y)

    # Analytical expression for the Jacobian of sin1 and sin2 and of exp1 and exp2
    Jsinkl1sinkl2 = @. (k1*l2-k2*l1)*cos(k1*x + l1*y)*cos(k2*x + l2*y)
    Jexpkl1expkl2 = @. (k2*l1-k1*l2)*(cos((k1+k2)*x + (l1+l2)*y)+im*sin((k1+k2)*x + (l1+l2)*y))

    @test test_parsevalsum(f1, g; realvalued=true)   # Real valued f with rfft
    @test test_parsevalsum(f1, g; realvalued=false)  # Real valued f with fft
    @test test_parsevalsum(f2, g; realvalued=false)  # Complex valued f with fft
    @test test_parsevalsum2(f1, g; realvalued=true)  # Real valued f with rfft
    @test test_parsevalsum2(f1, g; realvalued=false) # Real valued f with fft
    @test test_parsevalsum2(f2, g; realvalued=false) # Complex valued f with fft

    @test test_jacobian(sinkl1, sinkl1, 0*sinkl1, g)      # Test J(a, a) = 0
    @test test_jacobian(sinkl1, sinkl2, Jsinkl1sinkl2, g) # Test J(sin1, sin2) = Jsin1sin2
    @test test_jacobian(expkl1, expkl2, Jexpkl1expkl2, g) # Test J(exp1, exp2) = Jexp1exps2

    @test test_zeros()
    
    g = OneDGrid(dev, nx, Lx)
    σ = 0.5
    f1 = @. exp(-g.x^2/(2σ^2))
    @test test_parsevalsum2(f1, g; realvalued=true)  # Real valued f with rfft
    @test test_parsevalsum2(f1, g; realvalued=false) # Real valued f with fft

    
    # Radial spectrum tests. Note that ahρ = ∫ ah ρ dθ.
    n = 128; δ = n/10                                       # Parameters
    ahkl_isotropic(k, l) = exp(-(k^2 + l^2) / 2δ^2)                     # a  = exp(-ρ²/2δ²)
    ahρ_isotropic(ρ) = 2π * ρ * exp(-ρ^2 / 2δ^2)                        # aᵣ = 2π ρ exp(-ρ²/2δ²)
    @test test_radialspectrum(dev, n, ahkl_isotropic, ahρ_isotropic)
    @test test_radialspectrum(dev, n, ahkl_isotropic, ahρ_isotropic; rfft=true)

    ahkl_anisotropic(k, l) = exp(-(k^2 + l^2) / 2δ^2) * k^2 / (k^2+l^2) # a  = exp(-ρ²/2δ²) cos²θ
    ahρ_anisotropic(ρ) = π * ρ * exp(-ρ^2/2δ^2)                         # aᵣ = π ρ exp(-ρ²/2δ²)
    @test test_radialspectrum(dev, n, ahkl_anisotropic, ahρ_anisotropic)
    @test test_radialspectrum(dev, n, ahkl_anisotropic, ahρ_anisotropic; rfft=true)
    @test test_ongrid(dev)
  end

  @time @testset "Problem instantiation tests" begin
    include("test_instantiate_problem.jl")  

    for stepper in steppers
      @test instantiate_problem(dev, stepper)
      if occursin("Filtered", stepper)
        @test instantiate_problem_with_filter_kwargs(dev, stepper)
      end
    end
  end

  @time @testset "Diagnostics tests" begin
    include("test_diagnostics.jl")

    @test test_diagnosticsteps(dev, freq=1)
    @test test_diagnosticsteps(dev, freq=2)
    @test test_diagnosticsteps(dev, nsteps=100, freq=9, ndata=20)
    @test test_basicdiagnostics(dev)
    @test test_scalardiagnostics(dev, freq=1)
    @test test_scalardiagnostics(dev, freq=2)
    @test test_incrementdiagnostic(dev)
    @test test_extenddiagnostic(dev)
    @test test_getindex(dev)
  end

  @time @testset "Output tests" begin
    include("test_output.jl")
    
    @test test_withoutjld2()
    @test test_uniquepath()
    @test test_outputconstructor(dev)
    @test test_getindex(dev)
    @test test_saveproblem_saveoutput(dev)
    @test test_saveproblemTwoDGrid(dev)
    @test test_savediagnostic(dev)
    
    file = test_savefields(dev; parameter=2.2)
    @test file["params"]["parameter"] == 2.2
    @test_throws KeyError file["params"]["func"] == func
    @test_throws KeyError file["params"]["plan"] == plan
    close(file)
  end

  @time @testset "show() methods tests" begin
    prob = Diffusion.Problem(dev=dev)
    
    struct Params1 <: AbstractParams
        κ1 :: Float64
        κ2 :: Float64
      func :: Function
    end
    
    struct Params2 <: AbstractParams
        κ1 :: Float64
      func :: Function
        κ2 :: Float64
    end
    
    func(x) = sin(x^2)
    
    params1 = Params1(1.0, 2.0, func)
    params2 = Params2(1.0, func, 2.0)
    
    prob1 = FourierFlows.Problem(prob.eqn, "RK4", prob.clock.dt, prob.grid, prob.vars, params1, dev)
    prob2 = FourierFlows.Problem(prob.eqn, "RK4", prob.clock.dt, prob.grid, prob.vars, params2, dev)
    prob = prob1 # prob2 is only useful for testing params show() method
    
    @test repr(prob1.params) == "Parameters\n  ├───── parameter: κ1 -> Float64\n  ├───── parameter: κ2 -> Float64\n  └───── parameter: func -> Function\n"
    @test repr(prob2.params) == "Parameters\n  ├───── parameter: κ1 -> Float64\n  ├───── parameter: func -> Function\n  └───── parameter: κ2 -> Float64\n"
    if VERSION < v"1.6"
      @test repr(prob.vars) == "Variables\n  ├───── variable: c -> 128-element " * string(ArrayType(dev)) * "{Float64,1}\n  ├───── variable: cx -> 128-element " * string(ArrayType(dev)) * "{Float64,1}\n  ├───── variable: ch -> 65-element " * string(ArrayType(dev)) * "{Complex{Float64},1}\n  └───── variable: cxh -> 65-element " * string(ArrayType(dev)) * "{Complex{Float64},1}\n"
      @test repr(prob.eqn) == "Equation\n  ├──────── linear coefficients: L\n  │                              ├───type: Int64\n  │                              └───size: (65,)\n  ├───────────── nonlinear term: calcN!()\n  └─── type of state vector sol: Complex{Float64}"
    else
      if dev == CPU()
        @test repr(prob.vars) == "Variables\n  ├───── variable: c -> 128-element Vector{Float64}\n  ├───── variable: cx -> 128-element Vector{Float64}\n  ├───── variable: ch -> 65-element Vector{ComplexF64}\n  └───── variable: cxh -> 65-element Vector{ComplexF64}\n"
      else
        @test repr(prob.vars) == "Variables\n  ├───── variable: c -> 128-element " * string(ArrayType(dev)) * "{Float64, 1}\n  ├───── variable: cx -> 128-element " * string(ArrayType(dev)) * "{Float64, 1}\n  ├───── variable: ch -> 65-element " * string(ArrayType(dev)) * "{ComplexF64, 1}\n  └───── variable: cxh -> 65-element " * string(ArrayType(dev)) * "{ComplexF64, 1}\n"
      end
      @test repr(prob.eqn) == "Equation\n  ├──────── linear coefficients: L\n  │                              ├───type: Int64\n  │                              └───size: (65,)\n  ├───────────── nonlinear term: calcN!()\n  └─── type of state vector sol: ComplexF64"
    end
    @test repr(prob.clock) == "Clock\n  ├─── timestep dt: 0.01\n  ├────────── step: 0\n  └──────── time t: 0.0"
    @test repr(prob) == "Problem\n  ├─────────── grid: grid (on " * FourierFlows.griddevice(prob.grid) * ")\n  ├───── parameters: params\n  ├────── variables: vars\n  ├─── state vector: sol\n  ├─────── equation: eqn\n  ├────────── clock: clock\n  └──── timestepper: RK4TimeStepper"
  end

end # end loop over devices
