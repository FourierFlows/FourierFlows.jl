include("../src/domain.jl")

using Domain
# using PyPlot


# # test on square grids
# nx = 8;
# Lx = 2.0;
# g = Grid(nx, Lx);

# test on rectangular grids
nx = 32;
ny = 64;
Lx = 2.0;
Ly = 3.4;
g = Grid(nx, ny, Lx, Ly);





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

################################################################################
# Test fft's are giving what you expect them to give

if nx/2>=4

  m = 2;
  n = 1;

  f   = cos.(m*real(g.k[2])*g.X).*cos.(n*real(g.l[2])*g.Y);
  fh  = fft(f);
  fhr = rfft(f);

  fh_th  = zeros(size(fh));
  for i in g.krange, j in g.lrange
    if ( abs(g.K[i,j])==m*real(g.k[2]) && abs(g.L[i,j])==n*real(g.l[2]) )
      fh_th[i,j] = - g.nx*g.ny/4;
    end
  end

  fhr_th  = zeros(size(fhr));
  for i in g.krrange, j in g.lrange
    if ( abs(g.Kr[i,j])==m*real(g.k[2]) && abs(g.L[i,j])==n*real(g.l[2]) )
      fhr_th[i,j] = - g.nx*g.ny/4;
    end
  end

  if norm(fh-fh_th)>1e-12
    info("  fft for sin(mx)sin(ny) is not correcty calculated!")
  else
    println("  fft for sin(mx)sin(ny) seems OK.")
  end

  if norm(fhr-fhr_th)>1e-12
    info(" rfft for sin(mx)sin(ny) is not correcty calculated!")
  else
    println(" rfft for sin(mx)sin(ny) seems OK.")
  end

end


################################################################################
# Test ifft's taking you back where you started

f2 = real(ifft(fh));

if sqrt(mean((f-f2).^2)) > 1e-12
  info(" ifft for sin(mx)sin(ny) is not correcty calculated!")
else
  println(" ifft for sin(mx)sin(ny) seems OK.")
end

f2 = irfft(fhr,nx);


if sqrt(mean((f-f2).^2)) > 1e-12
  info("irfft for sin(mx)sin(ny) is not correcty calculated!")
else
  println("irfft for sin(mx)sin(ny) seems OK.")
end
