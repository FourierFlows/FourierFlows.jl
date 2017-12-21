include("../src/FourierFlows.jl")

using FourierFlows, PyPlot

nx = 128
Lx = 2π
Ly = Lx
ny = nx

g  = TwoDGrid(nx, Lx, ny, Ly)


# High-wavenumber filter for K, L and Kr, Lr wavenumber grids
cphi=0.65π
filterfac=23.6

wv = sqrt.((g.K*g.dx).^2 + (g.L*g.dy).^2)
wvr = sqrt.((g.Kr*g.dx).^2 + (g.Lr*g.dy).^2)

filter = exp.(-filterfac*(wv-cphi).^4);
filterr = exp.(-filterfac*(wvr-cphi).^4);

for i = 1:g.nk, j = 1:g.nl
	if real(wv[i, j]) < cphi
		filter[i, j] = 1
	end
end
for i = 1:g.nkr, j=1:g.nl
	if real(wvr[i, j]) < cphi
		filterr[i, j] = 1
	end
end

filterorder=4.0
innerfilterK=0.65
outerfilterK=0.95




function makefilter(g::TwoDGrid; order=4.0, innerK=0.65, outerK=1.0,
  realvars=true)

  # Get decay rate for filter
  decay = 15.0*log(10.0) / (outerK-innerK)^order
  decay = 23.6*π^order

  # Non-dimensional square wavenumbers
  if realvars
    KK = sqrt.( (g.Kr*g.dx/π).^2 + (g.Lr*g.dy/π).^2 )
  else
    KK = sqrt.( (g.K*g.dx/π).^2  + (g.L*g.dy/π).^2  )
  end

  filt = exp.( -decay*(KK-innerK).^order )
  filt[ real.(KK).<innerK ] = 1

  return filt
end

filter_Greg = makefilter(g; order=filterorder, innerK=innerfilterK,
  outerK=outerfilterK, realvars=true)


figure(1)
pcolormesh(abs.(filterr))

figure(2)
pcolormesh(abs.(filter_Greg))

figure(3)
pcolormesh(abs.(filter_Greg)-abs.(filterr))
colorbar()
