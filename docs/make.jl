using
  Documenter,
  Literate,
  CairoMakie,  # so that Literate.jl does not capture precompilation output
  Glob,
  FourierFlows

# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

 #####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

examples = [
    "OneDShallowWaterGeostrophicAdjustment.jl",
]

for example in examples
  example_filepath = joinpath(EXAMPLES_DIR, example)
  withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
    Literate.markdown(example_filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
  end
end

#####
##### Build and deploy docs
#####

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(
    collapselevel = 2,
       prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://fourierflows.github.io/FourierFlowsDocumentation/dev/",
)

pages = [
    "Home" => "index.md",
    "Installation Instructions" => "installation_instructions.md",
    "Code Basics" => "basics.md",
    "Grids" => "grids.md",
    "Aliasing" => "aliasing.md",
    "Problem" => "problem.md",
    "Diagnostics" => "diagnostics.md",
    "Output" => "output.md",
    "GPU" => "gpu.md",
    "Examples" => [ 
        "literated/OneDShallowWaterGeostrophicAdjustment.md",
        ],
    "Contributor's guide" => "contributing.md",
    "Library" => [ 
        "Contents" => "library/outline.md",
        "Public" => "library/public.md",
        "Private" => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(
   sitename = "FourierFlows.jl",
    authors = "Gregory L. Wagner and Navid C. Constantinou",
    modules = [FourierFlows],
     format = format,
      pages = pages,
    doctest = true,
     strict = :doctest,
      clean = true,
  checkdocs = :exports
)

withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
  deploydocs(
            repo = "github.com/FourierFlows/FourierFlowsDocumentation.git",
        versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
    push_preview = false,
       devbranch = "main"
  )
end


@info "Cleaning up temporary .jld2 and .nc files created by doctests..."
for file in vcat(glob("docs/*.jld2"))
  rm(file)
end
