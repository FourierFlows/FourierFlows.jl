include("../src/domain.jl")

using Domain
# using PyPlot


# # test on square grids
nx = 8                  # number of points
Lx = 2.0*pi             # Domain width
g = Grid(nx, Lx)

# test on rectangular grids
# nx, ny = 8, 12                  # number of points
# Lx, Ly = 2.0, 3.4             # Domain widths
# g = Grid(nx, ny, Lx, Ly);


println(" ")
################################################################################
# Test that the X & Y grids are actually created

if maximum(abs.(g.X))<1e-10
  info("X grid is ALL zeros!")
else
  println("X  grid seems OK.")
end

if maximum(abs.(g.Y))<1e-10
  info("Y grid is ALL zeros!")
else
  println("Y  grid seems OK.")
end


################################################################################
# Test that K, L, Kr, Lr grids are created in the correct order

test = 0;
for i in g.krange
  test += abs(maximum(abs.(g.K[i,:]))-minimum(abs.(g.K[i,:])));
  #the inner abs's are needed because Greg defined K as Complex128
end
if test>1e-12
  info("K grid is not correcty oriented!")
else
  println("K  grid seems OK.")
end

test = 0;
for i in g.lrange
  test += abs(maximum(abs.(g.L[:,i]))-minimum(abs.(g.L[:,i])));
  #the inner abs's are needed because Greg defined L as Complex128
end
if test>1e-12
  info("L grid is not correcty oriented!")
else
  println("L  grid seems OK.")
end


test = 0;
for i in g.krrange
  test += abs(maximum(abs.(g.Kr[i,:]))-minimum(abs.(g.Kr[i,:])));
  #the inner abs's are needed because Greg defined Kr as Complex128
end
if test>1e-12
  info("Kr grid is not correcty oriented!")
else
  println("Kr grid seems OK.")
end

println(" ")
################################################################################
# Test fft's are giving what you expect them to give

if nx/2>=4

  m = 1;
  n = 2;

  f1   = cos.(m*real(g.k[2])*g.X).*cos.(n*real(g.l[2])*g.Y);
  f2   = sin.(m*real(g.k[2])*g.X + n*real(g.l[2])*g.Y);
  f1h  = fft(f1);
  f2h  = fft(f2);
  f1hr = rfft(f1);
  f2hr = rfft(f2);

  f2hr_mul = Array{Complex128}(g.nkr, g.nl)
  A_mul_B!( f2hr_mul, g.rfftplan, f2 )

  f1h_th   = zeros(size(f1h));  f2h_th   = zeros(size(f2h));
  f1hr_th  = zeros(size(f1hr)); f2hr_th  = zeros(size(f2hr));

  for i in g.krange, j in g.lrange
    if ( abs(real(g.K[i,j])) == m*real(g.k[2]) && abs(real(g.L[i,j])) == n*real(g.l[2]) )
      f1h_th[i,j] = - g.nx*g.ny/4;
    end
    if ( real(g.K[i,j]) == m*real(g.k[2]) && real(g.L[i,j]) == n*real(g.l[2]) )
      f2h_th[i,j] = -g.nx*g.ny/2;
    elseif ( real(g.K[i,j]) == -m*real(g.k[2]) && real(g.L[i,j]) == -n*real(g.l[2]) )
      f2h_th[i,j] = g.nx*g.ny/2;
    end
  end
  f2h_th = -im*f2h_th;


  for i in g.krrange, j in g.lrange
    if ( abs(g.Kr[i,j])==m*real(g.k[2]) && abs(g.L[i,j])==n*real(g.l[2]) )
      f1hr_th[i,j] = - g.nx*g.ny/4;
    end
    if ( real(g.Kr[i,j]) == m*real(g.k[2]) && real(g.L[i,j]) == n*real(g.l[2]) )
      f2hr_th[i,j] = -g.nx*g.ny/2;
    elseif ( real(g.Kr[i,j]) == -m*real(g.k[2]) && real(g.L[i,j]) == -n*real(g.l[2]) )
      f2hr_th[i,j] = g.nx*g.ny/2;
    end
  end
  f2hr_th = -im*f2hr_th;

  if norm(f1h-f1h_th)/norm(f1h_th)>1e-12
    info("  fft for cos(mx)cos(ny) is not correcty calculated!")
  else
    println("  fft for cos(mx)cos(ny) seems OK.")
  end

  if norm(f1hr-f1hr_th)/norm(f1hr_th)>1e-12
    info(" rfft for cos(mx)cos(ny) is not correcty calculated!")
  else
    println(" rfft for cos(mx)cos(ny) seems OK.")
  end

  if norm(f2h-f2h_th)/norm(f2h_th)>1e-12
    info("  fft for sin(mx+ny) is not correcty calculated!")
  else
    println("  fft for sin(mx+ny) seems OK.")
  end

  if norm(f2hr-f2hr_th)/norm(f2hr_th)>1e-12
    info(" rfft for sin(mx+ny) is not correcty calculated!")
  else
    println(" rfft for sin(mx+ny) seems OK.")
  end

  if norm(f2hr_mul-f2hr_th)/norm(f2hr_th)>1e-12
    info(" rfft with A_mul_B for sin(mx+ny) is not correcty calculated!")
  else
    println(" rfft with A_mul_B for sin(mx+ny) seems OK.")
  end

end

println(" ")
################################################################################
# Test ifft's taking you back where you started

f1b = real(ifft(f1h));

if norm(f1-f1b)/norm(f1) > 1e-12
  info(" ifft for cos(mx)cos(ny) is not correcty calculated!")
else
  println(" ifft for cos(mx)cos(ny) seems OK.")
end

f1b = irfft(f1hr, nx);

if norm(f1-f1b)/norm(f1) > 1e-12
  info("irfft for cos(mx)cos(ny) is not correcty calculated!")
else
  println("irfft for cos(mx)cos(ny) seems OK.")
end

a = deepcopy(f2hr);

f2b = Array{Float64}(g.nx, g.ny)
A_mul_B!( f2b, g.irfftplan, f2hr )
if norm(f2-f2b)/norm(f2) > 1e-12
  info("irfft with A_mul_B for sin(mx+ny) is not correcty calculated!")
else
  println("irfft with A_mul_B for sin(mx+ny) seems OK.")
end

b = deepcopy(f2hr);

f2b3 = irfft(f2hr, g.nx);
if norm(f2-f2b3)/norm(f2) > 1e-12
  info("irfft for sin(mx+ny) is not correcty calculated!")
else
  println("irfft for sin(mx+ny) seems OK.")
end

c = deepcopy(f2hr);

println("why is a different from b?   norm(a-b)=", norm(a-b))


f2b4 = real(ifft(f2h));

if norm(f2-f2b4)/norm(f2) > 1e-12
  info(" ifft for sin(mx+ny) is not correcty calculated!")
else
  println(" ifft for sin(mx+ny) seems OK.")
end
