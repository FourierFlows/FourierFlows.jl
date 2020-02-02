sizeofphysicalgrid(g::OneDGrid) = (g.nx,)
sizeofphysicalgrid(g::TwoDGrid) = (g.nx, g.ny)
sizeofphysicalgrid(g::ThreeDGrid) = (g.nx, g.ny, g.nz)

sizeofcoefficientgrid(g::OneDGrid; realvalued=true) = realvalued ? (g.nkr,) : (g.nk,)
sizeofcoefficientgrid(g::TwoDGrid; realvalued=true) = realvalued ? (g.nkr, g.nl) : (g.nk, g.nl)
sizeofcoefficientgrid(g::ThreeDGrid; realvalued=true) = realvalued ? (g.nkr, g.nl, g.nm) : (g.nk, g.nl, g.nm)

function test_field_from_grid_values(gr::AbstractGrid{T}; realvalued=true, dev=CPU()) where T
  Ta = realvalued ? T : Complex{T}
   a = ArrayType(dev)(rand(Ta, sizeofphysicalgrid(gr)))
  ah = realvalued ? rfft(a) : fft(a)
  
   F = field_from_grid_values(a, gr)
  isapprox(F.grid_values, a) && isapprox(F.coefficient_values, ah)
end

function test_field_from_coefficient_values(gr::AbstractGrid{T}; realvalued=true, dev=CPU()) where T
  ah = ArrayType(dev)(rand(Complex{T}, sizeofcoefficientgrid(gr; realvalued=realvalued)))
  a = realvalued ? irfft(ah, gr.nx) : ifft(ah)
  
  F = field_from_coefficient_values(ah, gr; realvalued=realvalued)
  isapprox(F.grid_values, a) && isapprox(F.coefficient_values, ah)
end
