module HeatEquation

export
  Problem,
  updatevars!,
  set_c!,

using 
  FourierFlows, 
  FFTW

using FourierFlows: varsexpression
using LinearAlgebra: mul!, ldiv!

# --
# Problems
# --

noflow(x, y) = 0 # used as defaults for u, v functions in Problem()

function flowargs(u)
  umethods = methods(u)
  length(umethods.ms) == 1 || error("The functions u and v may only have one method.")
  umethods.ms[1].nargs-1
end

"""
    Problem(; parameters...)

Construct a constant diffusivity problem with steady or time-varying flow.
"""
function TwoDProblem(;
            nx = 128,
            Lx = 2π,
            ny = nx,
            Ly = Lx,
           kap = 0,
           eta = kap,
          kaph = 0,
         nkaph = 0,
             u = noflow,
             v = noflow,
            dt = 0.01,
       stepper = "RK4",
  unsteadyflow = false,
             T = Float64,
          grid = TwoDGrid(nx, Lx, ny, Ly; T=T)
  )

  flowargs(u) == flowargs(v) || error("Functions u and v must have same number of arguments.")
  U = flowargs(u) > 2 ? u : u.(grid.x, grid.y)
  V = flowargs(v) > 2 ? v : v.(grid.x, grid.y)

   p = Params(eta, kap, kaph, nkaph, U, V)
   v = Vars(grid)
  eq = Equation(p, grid)

  FourierFlows.Problem(eq, stepper, dt, grid, vars, params)
end

function OneDProblem(;
            nx = 128,
            Lx = 2π,
            dt = 0.01,
       stepper = "RK4",
             T = Float64
  )

   g = OneDGrid(nx, Lx)
   p = EmptyParams()
   v = Vars(g)
  eq = Equation(p, g)

  FourierFlows.Problem(eq, stepper, dt, grid, vars, params)
end

# --
# Params
# --

"""
    Params(eta, kap, kaph, nkaph, u, v)

Returns the params for constant diffusivity problem with time-varying flow.
"""
struct Params{Tη,Tκ,Th,Tu} <: AbstractParams
  eta::Tη         # Horizontal diffusivity
  kap::Tκ         # Vertical diffusivity
  kaph::Th        # Constant isotropic hyperdiffusivity
  nkaph::Int      # Constant isotropic hyperdiffusivity order
  u::Tu           # Advecting x-velocity
  v::Tu           # Advecting y-velocity
end

Params(eta, kap, u, v) = Params(eta, kap, 0eta, 0, u, v)

# --
# Equations
# --

"""
    Equation(p, g)

Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(p, g)
  L = @. -p.eta*g.kr^2 - p.kap*g.l^2 - p.kaph*g.Krsq^p.nkaph
  FourierFlows.Equation(L, calcN!, g)
end

# --
# Vars
# --

# Construct Vars types
const physicalvars = [:c, :cx, :cy]
const transformvars = [:ch, :cxh, :cyh]

eval(varsexpression(:Vars, physicalvars, transformvars))

"""
    Vars(g)

Returns the vars for constant diffusivity problem on grid g.
"""
function Vars(g::FourierFlows.AbstractGrid{T}) where T
  @createarrays T (g.nx, g.ny) c cx cy
  @createarrays Complex{T} (g.nkr, g.nl) ch cxh cyh
  Vars(c, cx, cy, ch, cxh, cyh)
end

# --
# Solvers
# --

"Calculate advection of c by steady flow."
function calcadvection!(Np, u, v, cx, cy, x, y, t)
  @. Np = u*cx + v*cy
  nothing
end

"Calculate advection of c by time-varying flow."
function calcadvection!(Np, u::Function, v::Function, cx, cy, x, y, t)
  @. Np = u(x, y, t)*cx + v(x, y, t)*cy
  nothing
end

"""
    calcN!(N, sol, t, s, v, p, g)

Calculate the advective terms for a tracer equation with constant diffusivity and time-varying flow.
"""
function calcN!(N, sol, t, cl, v, p, g)
  @. v.cxh = im * g.kr * sol
  @. v.cyh = im * g.l  * sol

  ldiv!(v.cx, g.rfftplan, v.cxh) # destroys v.cxh when using fftw
  ldiv!(v.cy, g.rfftplan, v.cyh) # destroys v.cyh when using fftw

  calcadvection!(v.cx, p.u, p.v, v.cx, v.cy, g.x, g.y, t)

  mul!(N, g.rfftplan, v.cx)
  nothing
end

# --
# Helper functions
# --

"""
    updatevars!(v, g, sol)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(v, g, sol)
  @. v.ch = sol
  ldiv!(v.c, g.rfftplan, deepcopy(v.ch))
  nothing
end

updatevars!(prob) = updatevars!(prob.vars, prob.grid, prob.sol)

"""
    set_c!(prob, c)

Set the solution as the transform of `c`.
"""
function set_c!(prob, c)
  mul!(prob.sol, g.rfftplan, c)
  updatevars(prob)
end

set_c!(prob, c::Function) = set_c!(prob, c.(prob.grid.x, prob.grid.y))

end # module
