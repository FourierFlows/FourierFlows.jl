using JLD2
import Base: getindex, setindex!, push!, append!, fieldnames
export Output, saveoutput, saveproblem, groupsize, savediagnostic

gridfieldstosave = [:nx, :ny, :Lx, :Ly, :X, :Y]

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

  # Initialize output: first remove trailing ".jld2" to form 'basefilename'
  if filename[end-4:end] == ".jld2"
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
  groupname = "timeseries"
  jldopen(out.filename, "a+") do file
    file["$groupname/t/$(out.prob.step)"] = out.prob.t
    for fieldname in keys(out.fields)
      file["$groupname/$fieldname/$(out.prob.step)"] = out[fieldname]
    end
  end
  nothing
end


""" 
    saveproblem(prob, filename)

Save certain aspects of a problem timestepper, grid, and params. Functions
that are fields in params are not saved.
"""
function saveproblem(prob::AbstractProblem, filename::String)

  jldopen(filename, "a+") do file
      file["timestepper/dt"] = prob.ts.dt   # Timestepper

      for field in gridfieldstosave         # Grid
        file["grid/$field"] = getfield(prob.grid, field)
      end

      for name in fieldnames(prob.params)   # Params
        field = getfield(prob.params, name)
        if !(typeof(field) <: Function)
          file["params/$name"] = field
        end
      end

  end

  nothing
end

saveproblem(out::Output) = saveproblem(out.prob, out.filename)


"""
    savediagnostic(diag, diagname)

Save diagnostics to file, labeled by the string diagname.
"""
function savediagnostic(diag::AbstractDiagnostic, diagname::String,
  filename::String) 

  jldopen(filename, "a+") do file
    file["diags/$diagname/time"] = diag.time
    file["diags/$diagname/data"] = diag.data
  end

  nothing
end
