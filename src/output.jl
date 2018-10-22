"""
    Output(prob, filename, fieldtuples...)

Define output for `prob` with fields and functions that calculate
the output in the list of tuples `fieldtuples = (fldname, func)...`.
"""
struct Output
  prob::Problem
  path::String
  fields::Dict{Symbol,Function}
end

withoutjld2(path) = path[end-4:end] == ".jld2" ? path[1:end-5] : path

"""
    uniquepath(path)

Returns `path` with a number appended if `isfile(path)`, incremented until `path` does not exist.
"""
function uniquepath(path)
  n = 0
  while isfile(path)
    n += 1
    path = withoutjld2(path) * "_$n.jld2"
  end
  path
end
    
function Output(prob, path, fields::Dict{Symbol,Function})
  truepath = uniquepath(withoutjld2(path)*".jld2") # ensure path ends in ".jld2"
  saveproblem(prob, truepath)
  Output(prob, filename, fields)
end

Output(prob, path, fields...) = Output(prob, path, Dict{Symbol,Function}([fields...]))

getindex(out::Output, key) = out.fields[key](out.prob)

"""
    saveoutput(out)

Save the fields in `out.fields` to `out.path`.
"""
function saveoutput(out)
  groupname = "snapshots"
  jldopen(out.filename, "a+") do file
    file["$groupname/t/$(out.prob.step)"] = out.prob.t
    for fieldname in keys(out.fields)
      file["$groupname/$fieldname/$(out.prob.step)"] = out[fieldname]
    end
  end
  nothing
end

"""
    savefields(file, field)

Saves some parameters of `prob.field`.
"""
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

Save certain aspects of a problem.
"""
function saveproblem(prob, filename)
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

Save diagnostics in `diag` to file, labeled by `diagname`.
"""
function savediagnostic(diag, diagname, filename)
  jldopen(filename, "a+") do file
    file["diags/$diagname/steps"] = diag.steps
    file["diags/$diagname/t"] = diag.t
    file["diags/$diagname/data"] = diag.data
  end
  nothing
end
