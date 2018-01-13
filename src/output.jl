using JLD2
import Base: getindex, setindex!, push!, append!, fieldnames
export Output, saveoutput, saveproblem, groupsize, savediagnostic

gridfieldstosave = [:nx, :ny, :Lx, :Ly, :X, :Y]

""" Output type for FourierFlows problems. """
mutable struct Output
  prob::Problem
  filename::String
  fields::Dict{Symbol, Function}
  init::Bool
end

function Output(prob::Problem, filename::String,
  fields::Dict{Symbol, Function})
  saveproblem(prob, filename)
  Output(prob, filename, fields, true)
end

""" Constructor for Outputs with no fields. """
function Output(prob::Problem, filename::String)
  fields = Dict{Symbol, Function}()
  Output(prob, filename, fields)
end

""" Constructor for Outputs in which the name, field pairs are passed as
tupled arguments."""
function Output(prob::Problem, filename::String, fieldtuples...)
  Output(prob, filename,
    Dict{Symbol, Function}([(symfld[1], symfld[2]) for symfld in fieldtuples]))
end

""" Get the current output field. """
function getindex(out::Output, key)
  out.fields[key](out.prob)
end

function setindex!(out::Output, calcfield::Function, fieldname::Symbol)
  out.fields[fieldname] = calcfield
end

""" Add output name, calculator pairs when supplied as tupled arguments. """
function push!(out::Output, newfields...)
  for i = length(newfields)
    out.fields[newfields[i][1]] = newfields[i][2]
  end
end

""" Append a dictionary of name, calculator pairs to the dictionary of
output fields. """
function append!(out::Output, newfields::Dict{Symbol, Function})
  for key in keys(newfields)
    push!(out, (key, newfields[key]))
  end
end

function fieldnames(out::Output)
  fieldnames(out.fields)
end

""" Save the current output fields. """
function saveoutput(out::Output)
  step = out.prob.step
  groupname = "timeseries"

  jldopen(out.filename, "a+") do file
    file["$groupname/t/$step"] = out.prob.t
    for fieldname in keys(out.fields)
      file["$groupname/$fieldname/$step"] = out[fieldname]
    end
  end

  nothing
end

""" Save attributes of the Problem associated with the given Output. """
function saveproblem(out::Output)
  saveproblem(out.prob, out.filename)
end

""" Find the number of elements in a JLD2 group. """
function groupsize(group::JLD2.Group)
  try
    value = length(group.unwritten_links) + length(group.written_links)
  catch
    value = length(group.unwritten_links)
  end
  return value
end

""" 
Save certain aspects of a Problem. Entire problems cannot be saved
in general, because functions cannot be saved (and functions may use
arbitrary numbers of global variables that cannot be included in a saved 
object). 
"""
function saveproblem(prob::AbstractProblem, filename::String)
  jldopen(filename, "a+") do file
      file["timestepper/dt"] = prob.ts.dt
      for field in gridfieldstosave
        file["grid/$field"] = getfield(prob.grid, field)
      end

      for name in fieldnames(prob.params)
        field = getfield(prob.params, name)
        if !(typeof(field) <: Function)
          file["params/$name"] = field
        end
      end
  end

  nothing
end

"""
Save diagnostics to file.
"""
function savediagnostic(diag::AbstractDiagnostic, diagname::String,
  filename::String) 

  jldopen(filename, "a+") do file
    file["diags/$diagname/time"] = diag.time
    file["diags/$diagname/data"] = diag.data
  end

  nothing
end
