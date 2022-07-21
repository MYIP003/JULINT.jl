# Write .int file from input amplitude and phase (.int format = float 32)
# Written by Michelle YIP
# Edited 31/05/2021

using LinearAlgebra

function julint_w(ph::Matrix{Float32}, A::Matrix{Float32})
  
   INT_r = A.*cos.(ph);
   INT_i = A.*sin.(ph);
   INT = Complex.(INT_r,INT_i);
   
   write("test.int",INT)
   return INT
end

function julint_w_fast(ph::Matrix{Float32}, A::Matrix{Float32})
   
   INT_r = INT_i = similar(ph);
   mul!(INT_r,A,cos.(ph));
   mul!(INT_i,A,sin.(ph));

   lines,width = size(ph);
   INT = Array{Float32}(undef,lines,Int(width*2));
   INT[:,1:2:end] = INT_r;
   INT[:,2:2:end] = INT_i;

   return transpose(INT)
end
