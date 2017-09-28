__precompile__()

import SpecialFunctions, 
       PyPlot,
       PyCall,
       JLD

using PyPlot


"Return the fftwavenumber vector with length n and domain size L."
fftwavenums(n::Int; L=1.0) = 2.0*pi/L*cat(1, 0:n/2, -n/2+1:-1)

""" Return the root-mean-square of an array. """
rms(q) = sqrt(mean(q.^2))


"Fast loop-based dealiasing method for complex, spectral-space vars."
function dealias!(a::Array{Complex{Float64}, 2}, g::AbstractGrid)
  if size(a)[1] == g.nk       # Transform of a complex var
    for j in g.lderange
      @simd for i in g.kderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end
  else                        # Transform of a real var
    for j in g.lderange
      @simd for i in g.krderange
        @inbounds a[i, j] = 0.0 + 0.0*im
      end
    end
  end
end


"Dealiasing method for 3-arrays with broadcasting for third dimension."
function dealias!(a::Array{Complex{Float64}, 3}, g::AbstractGrid)
  if size(a)[1] == g.nk       # Transform of a complex var
    for j in g.lderange
      @simd for i in g.kderange
        @inbounds @views @. a[i, j, :] = Complex{Float64}(0.0)
      end
    end
  else                        # Transform of a real var
    for j in g.lderange
      @simd for i in g.krderange
        @inbounds @views @. a[i, j, :] = Complex{Float64}(0.0)
      end
    end
  end
end


""" 
Generate a real and random two-dimensional distribution phi(x, y) with
a Fourier spectrum peaked around a central non-dimensional wavenumber kpeak.
The spectrum is normalized either by setting the root-mean-square value of phi
with the keyword 'rms', or the maximum value of phi with the keyword 'maxval'.
"""
function peaked_isotropic_spectrum(nkl::Tuple{Int, Int}, kpeak::Real; 
  ord=4.0, rms=1.0, maxval=0.0)

  # Non-dimensional wavenumbers
  nk, nl = nkl
  k, l   = fftwavenums(nk), fftwavenums(nl)

  K = zeros(Float64, nk, nl)
  for j = 1:nl, i = 1:nk
    K[i, j] = sqrt(k[i]^2.0 + l[j]^2.0)
  end
  
  # Generate random spectrum and then normalize 
  phih = exp.(2.0*im*pi*rand(nk, nl)) ./ (1.0 .+ K./kpeak).^ord

  # Normalize by maximum value if specified
  if maxval > 0
    phi = real.(ifft(phih))
    phi = maxval * phi / maximum(abs.(phi))
  else
    phih .*= rms ./ sqrt.(sum(abs.(phih).^2.0)) 
    phi = real.(ifft(phih))
  end

  return phi
end

"Alternative input form for peaked_isotropic_spectrum for square grids."
function peaked_isotropic_spectrum(nx::Int, npeak::Real; ord=4.0, rms=1.0, 
  maxval=0.0)
  peaked_isotropic_spectrum((nx, nx), npeak; ord=ord, rms=rms, maxval=maxval)
end




""" Return a 2D vorticity field corresponding to the Lamb Dipole with
strength Ue, radius R, and wavenumber k, and centered around 
(xc, yc)=center. The default value of 'center' is the middle of the grid."""
function lambdipole(Ue::Real, R::Real, g::TwoDGrid; center=(nothing, nothing))

  if center == (nothing, nothing)
    xc = mean(g.x)
    yc = mean(g.y)
  else
    xc = center[1]
    yc = center[2]
  end

  # Wavenumber corresponding to radius R and the first bessel func zero.
  k = 3.8317 / R 
  q0 = -2*Ue*k/SpecialFunctions.besselj(0, k*R)

  r = sqrt.((g.X-xc).^2.0 + (g.Y-yc).^2.0)
  q = q0 * SpecialFunctions.besselj.(1, k*r) .* (g.Y-yc)./r

  q[r .== 0.0] = 0.0 # just in case.
  q[r .> R] = 0.0

  return q
end





""" Return a vorticity field with magnitude q0, radius R, and center at
center[1], center[2] on a TwoDGrid g corresponding to a 'Gaussian vortex' with 
Gaussian streamfunction. """
function gaussianvortex(q0::Real, R::Real, g::TwoDGrid; 
  center=(nothing, nothing))

  if center == (nothing, nothing)
    xc = mean(g.x)
    yc = mean(g.y)
  else
    xc = center[1]
    yc = center[2]
  end

  ( q0/R^2.0 * ( (g.X-xc).^2.0 + (g.Y-yc).^2.0 - 2*R^2.0 )
        .* exp.( -((g.X-xc).^2.0 + (g.Y-yc).^2.0) / (2.0*R^2.0)) )
end



""" Return an array of random numbers on a TwoDGrid normalized to have a 
specifed rms value. """
function rmsrand(g::TwoDGrid, rmsval::Real)
  q = rand(g.nx, g.ny)
  q .*= rmsval / rms(q) 
  return q
end




""" Integrate the square of a variables's Fourier transform on a 2D grid
using Parseval's theorem, taking into account for FFT normalization and
testing whether the coefficients are the product of a real or complex
Fourier transform. """
function parsint(uh, g::TwoDGrid)

  # Weird normalization (hopefully holds for both Julia and MKL FFT)
  norm = 2*g.Lx*g.Ly/(g.nx^2*g.ny^2)

  nk, nl = size(uh) 

  # Different summing techniques for complex or real transform
  if nk == g.nkr
    U = sum(abs2.(uh[1, :]))
    U += 2*sum(abs2.(uh[2:end, :]))
  else
    U = sum(abs2.(uh))
  end

  return U*norm
end




# ----------------------------------------------------------------------------- 
# Plot utilities
abstract type AbstractPlot end

""" A type for a two-component plot. """
type TwoComponentPlot <: AbstractPlot

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






function makeplot!(pl::TwoComponentPlot; save=false, show=false)

  figure(pl.fig[:number]) 

  c1 = pl.comp1(pl.v, pl.p, pl.g)
  c2 = pl.comp2(pl.v, pl.p, pl.g)

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


  message = pl.message(pl.v, pl.p, pl.g)
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
    filename = @sprintf("%s_%03d.png", pl.plotname, pl.iplot)
    savefig(filename, dpi=pl.dpi)
    pl.iplot += 1
  end

  if show
    pause(0.01)
  end

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
  text(0.02, 0.94, message, transform=pl.axs[1][:transAxes], fontsize=16)

  if pl.xlabel == "" && pl.ylabel == ""
    pl.axs[1][:xaxis][:set_visible](false)
    pl.axs[1][:yaxis][:set_visible](false)
    pl.axs[2][:xaxis][:set_visible](false)
    pl.axs[2][:yaxis][:set_visible](false)
    pl.axs[3][:xaxis][:set_visible](false)
    pl.axs[3][:yaxis][:set_visible](false)
  end

  if (pl.title1 == "" && 
      pl.title2 == "" &&
      pl.xlabel == "" &&
      pl.ylabel == "" )
    tight_layout()
  end

  if show
    pause(0.01)
  end

  if save
    filename = @sprintf("%s_%03d.png", pl.plotname, pl.iplot)
    savefig(filename, dpi=pl.dpi)
    pl.iplot += 1
  end


end




function saveproblem(name::String, v::AbstractVars, 
  p::AbstractParams, g::AbstractGrid, ts::AbstractTimeStepper, 
  eq::AbstractEquation)

  savename = @sprintf("%s_problem.jld", name)

  JLD.save(savename, "v", v, "ts", ts, "p", p, "g", g, "eq", eq)

  nothing
end




function savesnapshot(name::String, v::AbstractVars, 
  ts::AbstractTimeStepper) 

  savename = @sprintf("%s_%06d.jld", name, ts.step)

  if symbol("c") in fieldnames(ts) # Problem has real and complex parts
    JLD.save(savename, "t", v.t, "step", ts.step, 
      "solc", v.solc, "solr", v.solr)
  else
    JLD.save(savename, "t", v.t, "step", ts.step, "sol", v.sol)
  end

  nothing
end




function savesnapshot(name::String, pb::Problem)

  savename = @sprintf("%s_%06d.jld", name, ts.step)

  if symbol("c") in fieldnames(pb.ts) # Problem has real and complex parts
    JLD.save(savename, "t", pb.v.t, "step", pb.ts.step, 
      "solc", pb.v.solc, "solr", pb.v.solr)
  else
    JLD.save(savename, "t", pb.v.t, "step", pb.ts.step, "sol", pb.v.sol)
  end

  nothing
end
