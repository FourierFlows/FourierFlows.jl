# Contributors Guide

This is a short guide for potential FourierFlows.jl contributors.

Please feel free to ask us questions and chat, either by raising an [issue](https://github.com/FourierFlows/FourierFlows.jl/issues) or starting a [discussion](https://github.com/FourierFlows/FourierFlows.jl/discussions).

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative practices. 
New contributors should make sure to read that guide.

## What can I do?

* Tackle an existing [issue](https://github.com/FourierFlows/FourierFlows.jl/issues).

* Improve documentation, docstrings, or comments if you found something is hard to use.

* Implement a new feature (e.g., a new diagnostic into a module).

If you're interested in working on something, let us know by commenting on
existing issues or by opening a new issue. This is to make sure no one else
is working on the same issue and so we can help and guide you in case there
is anything you need to know beforehand.

## Ground Rules

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.
* Do not commit changes to files that are irrelevant to your feature or bugfix
  (e.g., `.gitignore`).
* Be willing to accept criticism and work on improving your code; we don't want
  to break other users' code, so care must be taken not to introduce bugs. We
  discuss pull requests and keep working on them until we believe we've done a
  good job.
* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug

The easiest way to get involved is to report issues you encounter when using FourierFlows.jl 
or by requesting something you think is missing.

* Head over to the [issues](https://github.com/FourierFlows/FourierFlows.jl/issues) page.
* Search to see if your issue already exists or has even been solved previously.
* If you indeed have a new issue or request, click the "New Issue" button.
* Please be as specific as possible. Include the version of the code you were using, as
  well as what operating system you are running. The output of Julia's `versioninfo()`
  and `] status` is helpful to include. If possible, include complete, minimal example
  code that reproduces the problem.

## Setting up your development environment

* Install [Julia](https://julialang.org/) on your system.
* Install git on your system if it is not already there (install XCode command line tools on
  a Mac or git bash on Windows).
* Login to your GitHub account and make a fork of the
  [FourierFlows.jl repository](https://github.com/FourierFlows/FourierFlows.jl) by
  clicking the "Fork" button.
* Clone your fork of the FourierFlows.jl repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/your-user-name/FourierFlows.jl.git
  ```
* Navigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.
* Connect your repository to the upstream (main project).
  ```
  git remote add FourierFlows https://github.com/FourierFlows/FourierFlows.jl.git
  ```
* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the `Project.toml`
  file.
* You can test to make sure FourierFlows.jl works by typing in `] test`. This will run all
  the tests (this can take a while). In an ideal world you should run the tests on a machine
  with a GPU capability but if that's not a possibility that is available to you then don't 
  worry -- simply comment in a PR that you didn't test on GPU.

Your development environment is now ready!

## Pull Requests

Changes and contributions should be made via GitHub pull requests against the `master` branch.

When you're done making changes, commit the changes you made. Chris Beams has written 
a [guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages.

When you think your changes are ready to be merged into the main repository,
push to your fork and [submit a pull request](https://github.com/FourierFlows/FourierFlows.jl/compare/).

**Working on your first Pull Request?** You can learn how from this _free_ video series
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github), Aaron Meurer's [tutorial on the git workflow](https://www.asmeurer.com/git-workflow/), 
or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).

## Documentation

All PRs that introduces new features or new modules should be accompanied with appropriate 
docstrings and documentation. Writing documentation strings is really important to make sure 
others use your functionality properly. Didn't write new functions? That's fine, but be sure 
that the documentation for the code you touched is still in great shape. It is not uncommon 
to find some strange wording or clarification that you can take care of while you are here.

We encourage using [unicode](https://docs.julialang.org/en/v1/manual/unicode-input/) characters 
when writing docstrings, e.g., use `α` instead of `\alpha`. This makes the rendering of the 
docstrings in the Documentation and in the Julia REPL's `help?>` mode as similar as possible.

You can preview how the Documentation will look like after merging by building the documentation 
locally. To do that, from the main directory of your local repository call

```
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```
 
and then open `docs/build/index.html` in your favorite browser.

## Credits

This contributor's guide is heavily based on the [MetPy contributor's guide](https://github.com/Unidata/MetPy/blob/master/CONTRIBUTING.md) 
and on its "cover" made by [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/contributing/).
