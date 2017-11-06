__precompile__()

using JLD2, HDF5

export Output, saveoutput, saveproblem, groupsize




""" Output type for FourierFlows problems. """
type Output
  name::String
  calc::Function
  prob::Problem
  filename::String
end


""" Save output to file. """
function saveoutput(out::Output)
  step = out.prob.step
  groupname = "timeseries"
  name = out.name

  jldopen(out.filename, "a+") do file
    file["$groupname/$name/$step"] = out.calc(out.prob)
    file["$groupname/t/$step"]     = out.prob.t
  end

  nothing
end


""" Save an array of outputs to file. """
function saveoutput(outs::AbstractArray)

  step = outs[1].prob.step
  groupname = "timeseries"

  jldopen(outs[1].filename, "a+") do file
    file["$groupname/t/$step"] = outs[1].prob.t # save timestamp
    for out in outs # save output data 
      name = out.name
      file["$groupname/$name/$step"] = out.calc(out.prob)
    end
  end

  nothing
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


""" Save certain aspects of a Problem. Entire problems cannot be saved
in general, because functions cannot be saved (and functions may use
arbitrary numbers of global variables that cannot be included in a saved 
object). """
function saveproblem(prob::AbstractProblem, filename::String)

  gridfieldstosave = [:nx, :ny, :Lx, :Ly, :x, :y, :X, :Y, :K, :L, :Kr, :Lr]

  jldopen(filename, "a+") do file
      file["setup/dt"] = prob.ts.dt

      for field in gridfieldstosave
        file["setup/grid/$field"] = getfield(prob.grid, field)
      end

      names = fieldnames(prob.params)
      for name in names
        field = getfield(prob.params, name)
        if !(typeof(field) <: Function)
          file["setup/params/$name"] = field
        end
      end
  end
  nothing
end
