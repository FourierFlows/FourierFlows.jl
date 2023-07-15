using
  Documenter,
  DocumenterCitations,
  Literate,
  CairoMakie,  # so that Literate.jl does not capture precompilation output
  FourierFlows

#####
##### Generate examples
#####

@show bib_filepath = joinpath(@__DIR__, "src/references.bib")
bib = CitationBibliography(bib_filepath)

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

examples = [
  "OneDShallowWaterGeostrophicAdjustment.jl",
  "XYAdvectionSetupRun.jl",
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

format = Documenter.HTML(
    collapselevel = 2,
       prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://fourierflows.github.io/FourierFlowsDocumentation/stable",
)

pages = [
    "Home" => "index.md",
    "Installation Instructions" => "installation_instructions.md",
    "Code Basics" => "basics.md",
    "Grids" => "grids.md",
    "Aliasing" => "aliasing.md",
    "Problem" => "problem.md",
    "Time stepping" => "timestepping.md",
    "Diagnostics" => "diagnostics.md",
    "Output" => "output.md",
    "GPU" => "gpu.md",
    "Examples" => [ 
      "literated/OneDShallowWaterGeostrophicAdjustment.md",
      "literated/XYAdvectionSetupRun.md"
      ],
    "Contributor's guide" => "contributing.md",
    "Library" => [ 
      "Contents" => "library/outline.md",
      "Public" => "library/public.md",
      "Private" => "library/internals.md",
      "Function index" => "library/function_index.md",
      ],
    "References" => "references.md",
]

makedocs(bib, sitename = "FourierFlows.jl",
               authors = "Gregory L. Wagner and Navid C. Constantinou and contributors",
               modules = [FourierFlows],
                format = format,
                 pages = pages,
               doctest = false,
                strict = :doctest,
                 clean = true,
             checkdocs = :exports)

@info "Cleaning up temporary .jld2 and .nc files created by doctests or literated examples..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
  global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
  rm(file)
end

withenv("GITHUB_REPOSITORY" => "FourierFlows/FourierFlowsDocumentation") do
  deploydocs(
             repo = "github.com/FourierFlows/FourierFlowsDocumentation.git",
         versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
        forcepush = true,
     push_preview = true,
        devbranch = "main"
  )
end
