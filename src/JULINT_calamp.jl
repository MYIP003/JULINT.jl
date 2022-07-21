########################################################
# Calculate Amplitude calibration constant for SLC files
# (i.e. calculate mean Amplitude values from each SLC)
# Author: Michelle Yip
# Last Edited: 20-09-2021
########################################################
#using Mmap
#using LightXML
#using LinLogQuantization
#using HDF5
"""
Calculate mean amplitude values from each uncompressed SLC file
```
julint_calamp_uncompressed(infilename::String,width::Int64,outfilename::String)
```
# Arguments
- `infilename::String` : filename of calamp.in that contains names of the list of SLCs.
- `width::Int64` : width of SLC.
- `outfilename::String` : calamp.out that contain input SLCs and the mean amplitude values.
# Examples
```julia-repl
julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp.in";
julia> outfilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp.out";
julia> julint_calamp_uncompressed(infilename,25483,outfilename)
```
"""
function julint_calamp_uncompressed(infilename::String,width::Int64,outfilename::String)
  outfileios = open(outfilename,"a")
  for i in 1:length(readlines(infilename))
   ## Initialize
    sumamp = 0.0;
    A = 0.0;
    nof_pixels = 0;
    nof_zero_pixels = 0;
    calib = 0.0;
    infile = readlines(infilename)[i];

    ## Read from SLC
     s = open(infile);
     println("Opening ",infile," ...")
     lines = Int64(filesize(infile)/(width*sizeof(ComplexF32)));
     data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));

     for j in 1:lines
       for k in 1:width
         A = abs(data[j,k]);
         if A>0.001
          sumamp += A;
          nof_pixels += 1;
         else
          nof_zero_pixels += 1;
         end
       end
     end

    close(s);

    calib = round(sumamp/nof_pixels,digits=4);
    println("Mean amplitude = ",calib)
    println("Number of pixels with zero amplitude = ",nof_zero_pixels)
    println("Number of pixels with amplitude different than zero = ",nof_pixels)  

    outfile = string(infile," ",calib,"\n");
    write(outfileios,outfile);

  end
  close(outfileios);
end
#########################################################################
"""
Calculate mean amplitude values from each uncompressed/compressed amplitude files (CFloat32 to UInt8/UInt16)
```
julint_calamp(infilename::String,width::Int64,outfilename::String,nbits::Int64)
```
# Arguments
- `infilename::String` : filename of calamp.in that contains names of the compressed amplitude files.
- `width::Int64` : width of SLC.
- `outfilename::String` : calamp.out that contain names of compressed SLCs and their mean amplitude values.
- `nbits::Int64` : no. of bits of the input SLCs.
# Examples
```julia-repl
julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_8.in";
julia> outfilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_8.out";
julia> julint_calamp(infilename,width,outfilename,8)
```
"""
function julint_calamp(infilename::String,width::Int64,outfilename::String,nbits::Int64)
  outfileios = open(outfilename,"a")
  for i in 1:length(readlines(infilename))
   ## Initialize
    sumamp = 0.0;
    A = 0.0;
    nof_pixels = 0;
    nof_zero_pixels = 0;
    calib = 0.0;
    infile = readlines(infilename)[i];

    ## Read from SLC
     if nbits == 32
       s = open(infile);
       println("Opening ",infile," ...")
       lines = Int64(filesize(infile)/(width*sizeof(ComplexF32)));
       data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));

       for j in 1:lines
         for k in 1:width
           A = abs(data[j,k]);
           if A>0.001
            sumamp += A;
            nof_pixels += 1;
           else
            nof_zero_pixels += 1;
           end
         end
       end

     elseif nbits == 16 || nbits == 8 
       println("Opening ",infile," ...")
       T = typeof(h5read(infile,"A")).parameters[1];
       N = typeof(h5read(infile,"A")).parameters[2];
       println("De-quantizing amplitudes"," from ",T," to"," Float32"," ...")
       data = LogQuantArray{T,N}(h5read(infile,"A"), h5read(infile,"min"), h5read(infile,"max"));
       data = Array(data); 
     

       for j in 1:size(data)[1]
         for k in 1:size(data)[2]
           A = data[j,k];
           if A>0.001
            sumamp += A;
            nof_pixels += 1;
           else
            nof_zero_pixels += 1;
           end
         end
       end
    end # end if nbits

    calib = round(sumamp/nof_pixels,digits=4);
    println("Mean amplitude = ",calib)
    println("Number of pixels with zero amplitude = ",nof_zero_pixels)
    println("Number of pixels with amplitude different than zero = ",nof_pixels)

    outfile = string(infile," ",calib,"\n");
    write(outfileios,outfile);

  end
  close(outfileios);
end

#########################################################################
## Test the script
## Define the variables

#=
dir = "/Users/michelleymw/Desktop/Interferogram2/SLC/"
outfilename = string(dir,"calamp.out"); # Output filename
infilename = string(dir,"calamp.in");   # Input filename
linebytes = sizeof(ComplexF32)*width; #Define no. of bytes in a line
outfile = 0;
=#
#########################################################################
########Calculate Amplitude dispersion without using for loop############
#########################################################################

#=
function julint_calamp_without_forloop(infilename::String,width::Int,outfilename::String)
  for i in 1:length(readlines(infilename))
    ## Initialize
    sumamp = 0.0;
    nof_pixels = 0;
    nof_zero_pixels = 0;
    calib_factor = 0;
   
    infile = readlines(infilename)[i];
    s = open(infile);
    println("Opening ",infile," ...")
    lines = Int64(filesize(infile)/(width*sizeof(ComplexF32)))
    data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
    A = abs.(data);
    sumamp = sum(A[findall(>(0.001),A)]);
    nof_pixels = length(findall(>(0.001),A));
    nof_zero_pixels = lines*width-nof_pixels;
    calib = sumamp/nof_pixels;

    println("Mean amplitude = ",calib)
    println("Number of pixels with zero amplitude = ",nof_zero_pixels)
    println("Number of pixels with amplitude different than zero = ",nof_pixels)     
    outfileios = open(outfilename,"a")
    outfile = string(infile," ",calib,"\n");
    write(outfileios,outfile);
  end

end
=#
