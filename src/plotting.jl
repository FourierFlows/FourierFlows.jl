__precompile__()

import PyPlot, PyCall

using PyPlot


abstract type AbstractPlot end
abstract type AbstractTwoComponentPlot <: AbstractPlot end



type PlotComponent
  getfield
  title
  limz
  colors
end

type PlotMessage
  getmessage
  axs
  fontsize
end

type ProblemPlot
  prob
  fig
  axs
  components
  messages
  name
  dpi
end

# Key for figure sizes. Index is the number of plot components.
figsizekey = [(6, 6), (10, 5), (12, 4)]

function ProblemPlot(prob; components=nothing, messages=nothing, 
  name="problemplot", dpi=240)

  ncomps = length(components)

  fig, axs = subplots(nrows=1, ncols=ncomps, sharex=true, sharey=true,
    figsize=figsizekey[ncomps])

  ProblemPlot(prob, fig, axs, components, messages, name, dpi)
end

# One component plots --------------------------------------------------------- 
""" A type for a one-component plot. """
type OneComponentPlot <: AbstractPlot
  fig::PyPlot.Figure
  axs::PyCall.PyObject
  figsize::Tuple{Int64, Int64}
  g::TwoDGrid
  v::AbstractVars
  p::AbstractParams
  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String
  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String
  message::Function
  plotname::String
  iplot::Int
  dpi::Int
end

function OneComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, title1::String, limz1::AbstractArray, colors1::String, 
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  xlabel::String, ylabel::String, message::Function, plotname::String,
  dpi::Int; figsize=(10, 5))

  fig, axs = subplots(nrows=1, ncols=1, sharex=true, sharey=true,
    figsize=figsize)

  OneComponentPlot(fig, axs, figsize, g, v, p, 
    comp1, title1, limz1, colors1,
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, 0, dpi)
end

function OneComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1::Function, limz1::AbstractArray, colors1::String, 
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm,
  message::Function, plotname::String, dpi::Int; figsize=(10, 5))

  OneComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end

function OneComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  message::Function, plotname::String, dpi::Int; figsize=(12, 4))

  xlimz = [minimum(g.x), maximum(g.x)]
  ylimz = [minimum(g.y), maximum(g.y)]
  Lnorm = 1.0

  OneComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end


function makeplot!(pl::OneComponentPlot; save=false, show=false)

  c1 = pl.comp1(pl.v, pl.p, pl.g)

  figure(pl.fig[:number]) 
  axes(pl.axs)
  cla()

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, c1,
    cmap=pl.colors1, vmin=pl.limz1[1], vmax=pl.limz1[2])

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  xlabel(pl.xlabel)
  ylabel(pl.xlabel)

  message = pl.message(pl.v, pl.p, pl.g)
  text(0.02, 0.94, message, transform=pl.axs[:transAxes], fontsize=16)

  if pl.xlabel == "" && pl.ylabel == ""
    pl.axs[:xaxis][:set_visible](false)
    pl.axs[:yaxis][:set_visible](false)
  end

  if (pl.title1 == "" && 
      pl.xlabel == "" &&
      pl.ylabel == "" )
    tight_layout()
  end

  if save
    filename = @sprintf("%s_%03d.png", pl.plotname, pl.iplot)
    savefig(filename, dpi=pl.dpi)
    pl.iplot += 1
  end

  if show
    pause(0.01)
  end

end













# Two component plots --------------------------------------------------------- 
type TwoComponentProblemPlot <: AbstractTwoComponentPlot
  fig::PyPlot.Figure
  axs::Array{PyCall.PyObject, 1}
  figsize::Tuple{Int64, Int64}

  prob::AbstractProblem

  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String

  comp2::Function
  title2::String
  limz2::Array{Float64, 1}
  colors2::String

  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String

  message::Function
  plotname::String
  iplot::Int
  dpi::Int
end

function TwoComponentProblemPlot(prob::AbstractProblem,
  comp1, title1::String, limz1::AbstractArray, colors1::String, 
  comp2, title2::String, limz2::AbstractArray, colors2::String,
  xlimz, ylimz, Lnorm,
  xlabel::String, ylabel::String, message::Function, plotname::String;
  dpi=240, figsize=(10, 5)
)

  fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
    figsize=figsize)

  TwoComponentProblemPlot(fig, axs, figsize, prob, 
    comp1, title1, limz1, colors1, 
    comp2, title2, limz2, colors2, 
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, 0, dpi)
end

""" Constrcut a TwoComponentProblemPlot object with no scaling and 
covering the size of the domain. """
function TwoComponentProblemPlot(prob::AbstractProblem,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  xlabel::String, ylabel::String, message::Function, plotname::String;
  dpi=240, figsize=(10, 5))

  xlimz = [minimum(prob.grid.x), maximum(prob.grid.x)]
  ylimz = [minimum(prob.grid.y), maximum(prob.grid.y)]
  Lnorm = 1.0

  TwoComponentProblemPlot(prob, 
    comp1, "", limz1, colors1, comp2, "", limz2, colors2, xlimz, ylimz, Lnorm,
    xlabel, ylabel, message, plotname; dpi=dpi, figsize=figsize)
end

""" Constrcut a TwoComponentProblemPlot object with no scaling and 
covering the size of the domain. """
function TwoComponentProblemPlot(prob::AbstractProblem,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  message::Function, plotname::String;
  dpi=240, figsize=(10, 5))

  TwoComponentProblemPlot(prob, 
    comp1, limz1, colors1, comp2, limz2, colors2,
    "", "", message, plotname; dpi=dpi, figsize=figsize)
end

function get_components_and_message(pl::TwoComponentProblemPlot)
  pl.comp1(pl.prob), pl.comp2(pl.prob), pl.message(pl.prob)
end










""" A type for a two-component plot. """
type TwoComponentPlot <: AbstractTwoComponentPlot

  fig::PyPlot.Figure
  axs::Array{PyCall.PyObject, 1}
  figsize::Tuple{Int64, Int64}

  g::TwoDGrid
  v::AbstractVars
  p::AbstractParams

  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String

  comp2::Function
  title2::String
  limz2::Array{Float64, 1}
  colors2::String

  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String

  message::Function
  plotname::String
  iplot::Int
  dpi::Int
end

function TwoComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, title1::String, limz1::AbstractArray, colors1::String, 
  comp2, title2::String, limz2::AbstractArray, colors2::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  xlabel::String, ylabel::String, message::Function, plotname::String,
  dpi::Int; figsize=(10, 5))

  fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
    figsize=figsize)

  TwoComponentPlot(fig, axs, figsize, g, v, p, 
    comp1, title1, limz1, colors1,
    comp2, title2, limz2, colors2,
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, 0, dpi)
end

function TwoComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  xlabel::String, ylabel::String, message::Function, plotname::String,
  dpi::Int; figsize=(10, 5))

  TwoComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, dpi;
    figsize=figsize)
end

function TwoComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  message::Function, plotname::String, dpi::Int; figsize=(10, 5))

  TwoComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end

function TwoComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  message::Function, plotname::String, dpi::Int; figsize=(10, 5))

  xlimz = [minimum(g.x), maximum(g.x)]
  ylimz = [minimum(g.y), maximum(g.y)]
  Lnorm = 1.0

  TwoComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end


function get_components_and_message(pl::TwoComponentPlot)
  (pl.comp1(pl.v, pl.p, pl.g), pl.comp1(pl.v, pl.p, pl.g), 
    pl.message(pl.v, pl.p, pl.g))
end







function makeplot!(pl::AbstractTwoComponentPlot; save=false, show=false)

  c1, c2, message = get_components_and_message(pl)

  g = try
    pl.prob.grid
  catch
    pl.g
  end

  figure(pl.fig[:number]) 
  axes(pl.axs[1])
  cla()

  pcolormesh(g.X/pl.Lnorm, g.Y/pl.Lnorm, c1, 
    cmap=pl.colors1, vmin=pl.limz1[1], vmax=pl.limz1[2])

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  xlabel(pl.xlabel)
  ylabel(pl.xlabel)


  axes(pl.axs[2])
  cla()

  pcolormesh(g.X/pl.Lnorm, g.Y/pl.Lnorm, c2,
    cmap=pl.colors2, vmin=pl.limz2[1], vmax=pl.limz2[2])

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)
  xlabel(pl.xlabel)


  text(0.02, 0.94, message, transform=pl.axs[1][:transAxes], fontsize=16)

  if pl.xlabel == "" && pl.ylabel == ""
    pl.axs[1][:xaxis][:set_visible](false)
    pl.axs[1][:yaxis][:set_visible](false)
    pl.axs[2][:xaxis][:set_visible](false)
    pl.axs[2][:yaxis][:set_visible](false)
  end

  if (pl.title1 == "" && 
      pl.title2 == "" &&
      pl.xlabel == "" &&
      pl.ylabel == "" )
    tight_layout()
  end

  if save
    filename = @sprintf("%s_%06d.png", pl.plotname, pl.iplot)
    savefig(filename, dpi=pl.dpi)
    pl.iplot += 1
  end

  if show
    pause(0.01)
  end

end











# Three component plots ------------------------------------------------------- 

type ThreeComponentPlot <: AbstractPlot

  fig::PyPlot.Figure
  axs::Array{PyCall.PyObject, 1}
  figsize::Tuple{Int64, Int64}

  g::TwoDGrid
  v::AbstractVars
  p::AbstractParams

  comp1::Function
  title1::String
  limz1::Array{Float64, 1}
  colors1::String

  comp2::Function
  title2::String
  limz2::Array{Float64, 1}
  colors2::String

  comp3::Function
  title3::String
  limz3::Array{Float64, 1}
  colors3::String

  xlimz::Array{Float64, 1}
  ylimz::Array{Float64, 1}
  Lnorm::Float64
  xlabel::String
  ylabel::String

  message::Function
  plotname::String
  iplot::Int
  dpi::Int
end

function ThreeComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, title1::String, limz1::AbstractArray, colors1::String, 
  comp2, title2::String, limz2::AbstractArray, colors2::String,
  comp3, title3::String, limz3::AbstractArray, colors3::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  xlabel::String, ylabel::String, message::Function, plotname::String, 
  dpi::Int; figsize=(12, 4))

  fig, axs = subplots(nrows=1, ncols=3, sharex=true, sharey=true,
    figsize=figsize)

  ThreeComponentPlot(fig, axs, figsize, g, v, p, 
    comp1, title1, limz1, colors1,
    comp2, title2, limz2, colors2,
    comp3, title3, limz3, colors3,
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, 0, dpi)
end

function ThreeComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  comp3, limz3::AbstractArray, colors3::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  xlabel::String, ylabel::String, message::Function, plotname::String, 
  dpi::Int; figsize=(12, 4))

  ThreeComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    comp3, "", limz3, colors3,
    xlimz, ylimz, Lnorm, xlabel, ylabel, message, plotname, dpi;
    figsize=figsize)
end

function ThreeComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  comp3, limz3::AbstractArray, colors3::String,
  xlimz::AbstractArray, ylimz::AbstractArray, Lnorm, 
  message::Function, plotname::String, dpi::Int; figsize=(12, 4))

  ThreeComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    comp3, "", limz3, colors3,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end

function ThreeComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  comp3, limz3::AbstractArray, colors3::String,
  message::Function, plotname::String; figsize=(12, 4))

  xlimz = [minimum(g.x), maximum(g.x)]
  ylimz = [minimum(g.y), maximum(g.y)]
  Lnorm = 1.0

  ThreeComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    comp3, "", limz3, colors3,
    xlimz, ylimz, Lnorm, "", "", message, plotname, 240;
    figsize=figsize)
end



function ThreeComponentPlot(
  g::TwoDGrid, v::AbstractVars, p::AbstractParams,
  comp1, limz1::AbstractArray, colors1::String, 
  comp2, limz2::AbstractArray, colors2::String,
  comp3, limz3::AbstractArray, colors3::String,
  message::Function, plotname::String, dpi::Int; figsize=(12, 4))

  xlimz = [minimum(g.x), maximum(g.x)]
  ylimz = [minimum(g.y), maximum(g.y)]
  Lnorm = 1.0

  ThreeComponentPlot(g, v, p, 
    comp1, "", limz1, colors1,
    comp2, "", limz2, colors2,
    comp3, "", limz3, colors3,
    xlimz, ylimz, Lnorm, "", "", message, plotname, dpi;
    figsize=figsize)
end





function makeplot!(pl::ThreeComponentPlot; save=false, show=false)

  figure(pl.fig[:number]) 

  c1 = pl.comp1(pl.v, pl.p, pl.g)
  c2 = pl.comp2(pl.v, pl.p, pl.g)
  c3 = pl.comp3(pl.v, pl.p, pl.g)


  axes(pl.axs[1])
  cla()

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c1, 
    cmap=pl.colors1,
    vmin=pl.limz1[1], vmax=pl.limz1[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)

  xlabel(pl.xlabel)
  ylabel(pl.xlabel)


  axes(pl.axs[2])
  cla()

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c2, 
    cmap=pl.colors2,
    vmin=pl.limz2[1], vmax=pl.limz2[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)

  xlabel(pl.xlabel)


  axes(pl.axs[3])
  cla()

  pcolormesh(pl.g.X/pl.Lnorm, pl.g.Y/pl.Lnorm, 
    c3, 
    cmap=pl.colors3,
    vmin=pl.limz3[1], vmax=pl.limz3[2]
  )

  xlim(pl.xlimz[1]/pl.Lnorm, pl.xlimz[2]/pl.Lnorm)
  ylim(pl.ylimz[1]/pl.Lnorm, pl.ylimz[2]/pl.Lnorm)

  xlabel(pl.xlabel)


  message = pl.message(pl.v, pl.p, pl.g)
  text(0.02, 1.03, message, transform=pl.axs[1][:transAxes], fontsize=16)

  if pl.xlabel == "" && pl.ylabel == ""
    pl.axs[1][:xaxis][:set_visible](false)
    pl.axs[1][:yaxis][:set_visible](false)
    pl.axs[2][:xaxis][:set_visible](false)
    pl.axs[2][:yaxis][:set_visible](false)
    pl.axs[3][:xaxis][:set_visible](false)
    pl.axs[3][:yaxis][:set_visible](false)
  end

  tight_layout(rect=(0.05, 0.05, 0.95, 0.95))

  if show
    pause(0.01)
  end

  if save
    filename = @sprintf("%s_%03d.png", pl.plotname, pl.iplot)
    savefig(filename, dpi=pl.dpi)
    pl.iplot += 1
  end


end
