"""
# Exports
$(EXPORTS)
"""
module XYAdvection
   
   export set_2Dgrid, set_Params, set_Vars, set_Problem, set_initialConditions!, 
   run_and_save, get_data, updatevars!, stepforward!, saveoutput

   using DocStringExtensions

   using FourierFlows, JLD2
   using LinearAlgebra: mul!, ldiv!

   """
       struct Params{T} <: AbstractParams

   The parameters for XYAdvection problem:
   $(FIELDS)
   """
   struct Params{T} <: AbstractParams
      "Diffusion coefficient"
      J  :: T
      "Advection coefficient"
      A  :: T
      "LG coefficient"
      Γ  :: T
      "Average order parameter"
      α  :: T
   end

   """
       struct Vars2D{Aphys, Atrans} <: AbstractVars

   The variables for XYAdvection problem:

   $(FIELDS)
   """
   struct Vars2D{Aphys, Atrans} <: AbstractVars
      "px and py field and its x and y derivative pxx, pxy, pyx, pyy"
      px  :: Aphys
      py  :: Aphys

      pxx :: Aphys
      pxy :: Aphys
      pyx :: Aphys
      pyy :: Aphys  
      
      "px_hat, py_hat; FFT of px, py and its derivative fields"
      pxh  :: Atrans
      pyh  :: Atrans
      pxxh :: Atrans
      pxyh :: Atrans
      pyxh :: Atrans
      pyyh :: Atrans
   end

   """
       Vars2D(grid)

   Return the variables `vars` for a XYAdvection problem on `grid`.
   
   $(TYPEDFIELDS)
   """
   function Vars2D(grid::TwoDGrid{T}) where T
      Dev = typeof(grid.device)
   
      @devzeros Dev T (grid.nx, grid.ny) px py pxx pxy pyx pyy
      @devzeros Dev Complex{T} (grid.nkr, grid.nl) pxh pyh pxxh pxyh pyxh pyyh
      
      return Vars2D(px, py, pxx, pxy, pyx, pyy, pxh, pyh, pxxh, pxyh, pyxh, pyyh)
   end

   """
       calcN!(N, sol, t, clock, vars, params, grid)

   Calculate the nonlinear term for the XYAdvection equation.
   """
   function calcN!(N, sol, t, clock, vars, params, grid)
      # multiply p__h with ik to get derivatives
      @. vars.pxxh = im * grid.kr .* sol[:,:,1]
      @. vars.pxyh = im * grid.l  .* sol[:,:,1]

      @. vars.pyxh = im * grid.kr .* sol[:,:,2]
      @. vars.pyyh = im * grid.l  .* sol[:,:,2]

      # get ik*p__h in physical space
      ldiv!(vars.pxx, grid.rfftplan, vars.pxxh) # destroys vars.pxxh when using fftw
      ldiv!(vars.pxy, grid.rfftplan, vars.pxyh) # destroys vars.pxyh when using fftw
      ldiv!(vars.pyx, grid.rfftplan, vars.pyxh) # destroys vars.pyxh when using fftw
      ldiv!(vars.pyy, grid.rfftplan, vars.pyyh) # destroys vars.pyyh when using fftw
      
      # non-linear term
      @. vars.pxx = params.A * ((vars.px * vars.pxx) + (vars.py * vars.pxy)) + params.Γ * (params.α - (vars.px^2 + vars.py^2))*vars.px
      @. vars.pyx = params.A * ((vars.px * vars.pyx) + (vars.py * vars.pyy)) + params.Γ * (params.α - (vars.px^2 + vars.py^2))*vars.py

      # go to fourier space and define N
      mul!(vars.pxxh, grid.rfftplan, vars.pxx)
      mul!(vars.pyxh, grid.rfftplan, vars.pyx)
      N[:,:, 1] = vars.pxxh
      N[:,:, 2] = vars.pyxh

      dealias!(N[:,:,1], grid)
      dealias!(N[:,:,2], grid)
      return nothing
   end

   """
       Equation(params,grid)
   """
   function Equation(params::Params, grid::TwoDGrid)
      dev = grid.device

      # Linear operator
      L = zeros(dev, eltype(grid), (grid.nkr, grid.nl,2))
      @. L[:,:,1] = - params.J * grid.kr^2 - params.J * grid.l^2
      @. L[:,:,2] = - params.J * grid.kr^2 - params.J * grid.l^2
      
      # full equation
      return FourierFlows.Equation(L, calcN!, grid)
   end

   #################################### Exported Functions #################################
   """
   To be called from main
   
       set_2Dgrid(nx::Int64,Lx; dev::Device=CPU(),aliased_fraction = 0,kwargs...)
   
   setup 2D grid given the parameters
   """
   function set_2Dgrid(nx::Int64,Lx; dev::Device=CPU(),aliased_fraction = 0,kwargs...)
      grid = TwoDGrid(dev; nx, Lx, aliased_fraction,kwargs...)
      return grid
   end

   """
   To be called from main
   
       set_Params(;J::T,A::T,Γ::T,α::T)
   
   setup parameter values for the problem
   """
   function set_Params(;J::T,A::T,Γ::T,α::T) where {T<:Float64}
      params = Params(J, A, Γ, α)
      return params
   end

   """
   To be called from main
   
       set_Vars(grid::TwoDGrid)
   
   setup variables for the system
   """
   function set_Vars(grid::TwoDGrid)
      vars = Vars2D(grid)
      return vars
   end

   
   """
   To be called from main

       set_Problem(grid::TwoDGrid,params::Params,vars::Vars2D,dt::Float64=0.02,stepper = "ForwardEuler";stepperkwargs...)

   setup the FourierFlows.Problem 
   """
   function set_Problem(grid::TwoDGrid,params::Params,vars::Vars2D,dt::Float64=0.02,stepper = "ForwardEuler";stepperkwargs...)

      equation = Equation(params,grid)

      prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params; stepperkwargs...)
      return prob
   end

   """
       updatevars!(prob)
   """
   function updatevars!(prob)
      vars, grid, sol = prob.vars, prob.grid, prob.sol
   
      @. vars.pxh  = sol[:,:, 1]
      @. vars.pyh  = sol[:,:, 2]
   
      ldiv!(vars.px,  grid.rfftplan, deepcopy(sol[:,:, 1])) # use deepcopy() because irfft destroys its input
      ldiv!(vars.py,  grid.rfftplan, deepcopy(sol[:,:, 2])) # use deepcopy() because irfft destroys its input 
      return nothing
   end

   """
       set_initialConditions!(prob, px, py)

   Set the solution `sol` as the transform of `px` and `py` and update `vars`.
   """
   function set_initialConditions!(prob, u, v)
      vars, grid, sol = prob.vars, prob.grid, prob.sol
   
      cast_type = typeof(vars.px) # determine the type of vars.px

      # below, e.g., A(px0) converts px0 to the same type as vars expects
      # (useful when px0 is a CPU array but grid.device is GPU)
      mul!(vars.pxh, grid.rfftplan, cast_type(u))
      mul!(vars.pyh, grid.rfftplan, cast_type(v))
   
      @. sol[:,:,1] = vars.pxh
      @. sol[:,:,2] = vars.pyh
   
      updatevars!(prob)
   
      return nothing
   end

   """
   To be called from main

       run_and_save(prob,nsteps::T=10^5,printFreq::T=10^3;filepath::T2 = ".",filename::T2 = "XYAdvection_data.jld2") where {T<:Integer,T2<:String}

   Run the problem and save output to file.
   """
   function run_and_save(prob,nsteps::T=10^5,printFreq::T=10^3;filepath::T2 = ".",filename::T2 = "XYAdvection_data.jld2") where {T<:Integer,T2<:String}
      # assert nsteps/printFreq is an integer
      @assert isinteger(nsteps/printFreq) "requires nsteps/printFreq == Integer"

      fname = joinpath(filepath, filename)
      get_sol(prob) = prob.sol
      out = Output(prob, fname, (:sol, get_sol))
      saveproblem(out)
      for i in 0:Int(nsteps)
         if i % Int(printFreq)==0
            saveoutput(out)
         end
         stepforward!(prob)
         updatevars!(prob)
      end
      return out
   end

   """
   To be called from main

       get_data(titr::Integer,filepath::String,nx::Integer)

   Read output data from file and convert to physical space for analysis.
   """
   function get_data(titr::Integer,filepath::String,nx::Integer)
      file = jldopen(filepath)
      px = irfft(file[string("snapshots/sol/", titr)][:,:, 1], nx)
      py = irfft(file[string("snapshots/sol/", titr)][:,:, 2], nx)
      close(file)
      return px,py
   end
   ############################################################################################################
end #end module