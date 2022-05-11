# Extract phase and amplitude from complex floating matrix retrieved from julint_r
# Written by Michelle YIP on 31/05/2021

# Get phase and amplitude  from INT
# ph, A = julint_i2pha(INT);
function julint_i2pha(INT::Matrix{ComplexF32})
   INT_r = real(INT);
   INT_i = imag(INT);
   ph = atan.(INT_i,INT_r);                 # Phase
   A = sqrt.(abs2.(INT_r) + abs2.(INT_i));  # Amplitude
   return ph, A
end

