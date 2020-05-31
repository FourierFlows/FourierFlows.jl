using
  Documenter,
  Literate,
  Plots,  # to not capture precompilation output
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
  Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true, nbviewer_root_url="github.com/FourierFlows/FourierFlowsDocumentation")
  Literate.notebook(example_filepath, OUTPUT_DIR, documenter=true, nbviewer_root_url="github.com/FourierFlows/FourierFlowsDocumentation")
  Literate.script(example_filepath, OUTPUT_DIR, documenter=true)
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
      canonical = "https://fourierflows.github.io/FourierFlowsDocumentation/dev/"
)

makedocs(
    modules = [FourierFlows],
    doctest = false,
      clean = true,
   checkdocs = :all,
     format = format,
    authors = "Gregory L. Wagner and Navid C. Constantinou",
   sitename = "FourierFlows.jl",
  
      pages = Any[
              "Home" => "index.md",
              "Code Basics" => "basics.md",
              "Forcing" => "forcing.md",
              "Examples" => [ 
                  "generated/OneDShallowWaterGeostrophicAdjustment.md",
               ],
              "DocStrings" => Any[
                  "man/types.md",
                  "man/functions.md"]
                 ]
)

withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
  deploydocs(        repo = "github.com/FourierFlows/FourierFlowsDocumentation.git",
                 versions = ["stable" => "v^", "v#.#.#"],
             push_preview = true,
  )
end
