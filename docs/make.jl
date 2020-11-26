using
  Documenter,
  Literate,
  Plots,  # so that Literate.jl does not capture precompilation output
  FourierFlows
  
# Gotta set this environment variable when using the GR run-time on Travis CI.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "OneDShallowWaterGeostrophicAdjustment.jl",
]

for example in examples
  example_filepath = joinpath(EXAMPLES_DIR, example)
  withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
    Literate.notebook(example_filepath, OUTPUT_DIR, documenter=true)
    Literate.script(example_filepath, OUTPUT_DIR, documenter=true)
  end
end

#####
##### Build and deploy docs
#####

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid Travis CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/FourierFlowsDocumentation/dev/",
     # mathengine = Documenter.MathJax()
)

pages = [
    "Home" => "index.md",
    "Installation Instructions" => "installation_instructions.md",
    "Code Basics" => "basics.md",
    "Grids" => "grids.md",
    "Problem" => "problem.md",
    "GPU" => "gpu.md",
    "Examples" => [ 
        "generated/OneDShallowWaterGeostrophicAdjustment.md",
        ],
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
     # strict = true,
      clean = true,
  checkdocs = :exports
)

withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
  deploydocs(        repo = "github.com/FourierFlows/FourierFlowsDocumentation.git",
                 versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
             push_preview = true
  )
end
