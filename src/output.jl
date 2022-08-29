"""
    struct Output

The composite type for output.

$(TYPEDFIELDS)
"""
struct Output
    "the relevant problem for the output"
      prob :: Problem
    "the path for the output file"
      path :: String
    "the fields to be saved; the relevant problem for this diagnostic"
    fields :: Dict{Symbol, Function}
  
  @doc """
      Output(prob, filename, fieldtuples...)

  Define output for `prob` with fields and functions that calculate
  the output in the list of tuples `fieldtuples = (fldname, func)...`.
  """
  function Output(prob, path, fieldtuples::Dict{Symbol,Function})
    truepath = uniquepath(withoutjld2(path) * ".jld2") # ensure path ends in ".jld2"
    return new(prob, truepath, fieldtuples)
  end
end

Output(prob, path, fields...) = Output(prob, path, Dict{Symbol,Function}([fields...]))
Output(prob, path, field::Tuple{Symbol,T}) where T = Output(prob, path, Dict{Symbol,Function}([field]))

withoutjld2(path) = (length(path)>4 && path[end-4:end] == ".jld2") ? path[1:end-5] : path

"""
    uniquepath(path)

Return `path` with a number appended if `isfile(path) == true`. The number is incremented
until `path` does not exist.
"""
function uniquepath(path)
  n = 1
  
  if isfile(path)
    path = withoutjld2(path) * "_$n.jld2"
  end
  
  while isfile(path)
    n += 1
    path = withoutjld2(path)[1:end-length("_$(n-1)")] * "_$n.jld2"
  end
  
  return path
end

getindex(out::Output, key) = out.fields[key](out.prob)

"""
    saveoutput(out)

Save the fields in `out.fields` to `out.path`.
"""
function saveoutput(out)
  groupname = "snapshots"
  jldopen(out.path, "a+") do path
    path["$groupname/t/$(out.prob.clock.step)"] = out.prob.clock.t
    for fieldname in keys(out.fields)
      path["$groupname/$fieldname/$(out.prob.clock.step)"] = out[fieldname]
    end
  end
  
  return nothing
end

"""
    savefield(file, location, data)

Saves a particular field's `data` to `file`.
"""
savefield(file, location, data) = file[location] = data
savefield(file, location, data::AbstractArray) = file[location] = Array(data)

"""
    savefields(file, grid)

Saves some parameters of problem's `grid` to `file`.
"""
function savefields(file::JLD2.JLDFile{JLD2.MmapIO}, grid::TwoDGrid)
  for field in [:nx, :ny, :Lx, :Ly, :x, :y]
    savefield(file, "grid/$field", getfield(grid, field))
  end
  
  return nothing
end

function savefields(file::JLD2.JLDFile{JLD2.MmapIO}, grid::OneDGrid)
  for field in [:nx, :Lx, :x]
    savefield(file, "grid/$field", getfield(grid, field))
  end
  
  return nothing
end

forbidden_output_types = [Function, FFTW.FFTWPlan, CUDA.CUFFT.CuFFTPlan]

function savefields(file::JLD2.JLDFile{JLD2.MmapIO}, params::AbstractParams)
  for name in fieldnames(typeof(params))
    field = getfield(params, name)
    if !(any([typeof(field) <: forbidden for forbidden in forbidden_output_types]))
      savefield(file, "params/$name", field)
    end
  end
  
  return nothing
end

function savefields(file::JLD2.JLDFile{JLD2.MmapIO}, clock::Clock)
  file["clock/dt"] = clock.dt
  
  return nothing
end

function savefields(file::JLD2.JLDFile{JLD2.MmapIO}, eqn::Equation)
  savefield(file, "eqn/L", eqn.L)
  file["eqn/dims"] = eqn.dims
  file["eqn/T"] = eqn.T
  
  return nothing
end

"""
    saveproblem(prob, filename)

Save certain aspects of a problem `prob` to `filename`.
"""
function saveproblem(prob, filename)
  file = jldopen(filename, "a+")

  for field in [:eqn, :clock, :grid, :params]
    savefields(file, getfield(prob, field))
  end
  close(file)

  return nothing
end

saveproblem(out::Output) = saveproblem(out.prob, out.path)

"""
    savediagnostic(diagnostic, diagname, filename)

Save `diagnostic` to `filename` under name `diagname`. Only the computed diagnostic
is saved, that is, everything up to diagnostic's iteration `diagnostic.i`.
"""
function savediagnostic(diagnostic, diagname, filename)
  jldopen(filename, "a+") do file
    file["diagnostics/$diagname/steps"] = diagnostic.steps[1:diagnostic.i]
    file["diagnostics/$diagname/t"] = diagnostic.t[1:diagnostic.i]
    file["diagnostics/$diagname/data"] = diagnostic.data[1:diagnostic.i]
  end
  
  return nothing
end

show(io::IO, out::Output) =
     print(io, "Output\n",
               "  ├──── prob: ", summary(out.prob), '\n', 
               "  ├──── path: ", out.path, '\n', 
               "  └── fields: ", out.fields)
