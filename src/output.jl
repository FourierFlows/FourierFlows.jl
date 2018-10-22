const twodimgridfieldstosave = [:nx, :ny, :Lx, :Ly, :x, :y]

"""
    Output(prob, filename, fieldtuples...)

Define output for the Problem prob with fields and functions that calculate
the output in (fldname, func) tuples.
"""
mutable struct Output
  prob::Problem
  filename::String
  fields::Dict{Symbol,Function}
  init::Bool
end

function Output(prob::Problem, filename::String, fields::Dict{Symbol,Function})
  if filename[end-4:end] == ".jld2" # Initialize output: remove trailing ".jld2" to form 'basefilename'
    basefilename = filename[1:end-5]
  else
    basefilename = filename
  end

  filename = basefilename*".jld2"
  n = 0
  while isfile(filename) # append numbers until unique name is found
    n += 1
    filename = basefilename * "_$n.jld2"
  end

  saveproblem(prob, filename)
  Output(prob, filename, fields, true)
end

function Output(prob::Problem, filename::String, fieldtuples...)
  Output(prob, filename, Dict{Symbol,Function}(
    [ (symfld[1], symfld[2]) for symfld in fieldtuples ]
  ))
end

getindex(out::Output, key) = out.fields[key](out.prob)

"""
    saveoutput(out)

Save current output fields for file in out.filename.
"""
function saveoutput(out::Output)
  groupname = "snapshots"
  jldopen(out.filename, "a+") do file
    file["$groupname/t/$(out.prob.step)"] = out.prob.t
    for fieldname in keys(out.fields)
      file["$groupname/$fieldname/$(out.prob.step)"] = out[fieldname]
    end
  end
  nothing
end

function savefields(file, grid::TwoDGrid)
  for field in [:nx, :ny, :Lx, :Ly, :x, :y]
    file["grid/$field"] = getfield(grid, field)
  end
  nothing
end

function savefields(file, grid::OneDGrid)
  for field in [:nx, :Lx, :x]
    file["grid/$field"] = getfield(grid, field)
  end
  nothing
end

function savefields(file, params::AbstractParams)
  for name in fieldnames(typeof(params))
    field = getfield(params, name)
    if !(typeof(field) <: Function)
      file["params/$name"] = field
    end
  end
  nothing
end

function savefields(file, clock::Clock)
  file["clock/dt"] = clock.dt   # Timestepper
  nothing
end

function savefields(file, eqn::Equation)
  file["eqn/L"] = eqn.L       
  file["eqn/dims"] = eqn.dims 
  file["eqn/T"] = eqn.T 
  nothing
end

"""
    saveproblem(prob, filename)

Save certain aspects of a problem timestepper, grid, and params. Functions
that are fields in params are not saved.
"""
function saveproblem(prob::Problem, filename::String)
  jldopen(filename, "a+") do file
    for field in [:eqn, :clock, :grid, :params]
      savefields(file, getfield(prob, field))
    end
  end
  nothing
end

saveproblem(out::Output) = saveproblem(out.prob, out.filename)

"""
    savediagnostic(diag, diagname)

Save diagnostics to file, labeled by the string diagname.
"""
function savediagnostic(diag::Diagnostic, diagname::String, filename::String)
  jldopen(filename, "a+") do file
    file["diags/$diagname/steps"] = diag.steps
    file["diags/$diagname/t"] = diag.t
    file["diags/$diagname/data"] = diag.data
  end
  nothing
end
