########################################################
# Compress SLC files from CFloat32 to UInt16/UInt8
# using Linear/Logarithmic Quantization
# Author: Michelle Yip
# Last Edited: 23-09-2021
########################################################
# using Mmap
# using LinLogQuantization
# using HDF5

"""
Compress SLC files from CFloat32 to UInt16/UInt8 using Linear/Logarithmic Quantization
# Arguments
- `infilename::String` : name of list of SLCs to be compressed.
- `width::Int64` : width of SLC.
- `linlogbit::Int64` : number of quantization bit, 16 or 8.
# Examples
```julia-repl
julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp.in";
julia> julint_compress(infilename,25483,16)
``` 
"""
function julint_compress(infilename::String,width::Int64,linlogbit::Int64)
   for i in 1:length(readlines(infilename))
     infile = readlines(infilename)[i];
     s = open(infile);
     println("Opening ",infile," ...")
     lines = Int64(filesize(infile)/(width*sizeof(ComplexF32)));
     data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
    # ph = angle.(data);
    # A = abs.(data);
     
     if linlogbit == 16
       println("Compressing to UInt16 ...")
       log16A = LogQuant16Array(abs.(data));
       lin16ph = LinQuant16Array(angle.(data));
       # Save compressed phase and amplitude to outout file
       println("Compressed amplitude saving as ",infile,"_","A","_",linlogbit,".h5")
       println("Compressed phase saving as ",infile,"_","ph","_",linlogbit,".h5")
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"A",log16A.A);
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"max",log16A.max);
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"min",log16A.min);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"A",lin16ph.A);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"max",lin16ph.max);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"min",lin16ph.min);
  
     elseif linlogbit == 8
       println("Compressing to UInt8 ...")
       log8A = LogQuant8Array(abs.(data));
       lin8ph = LinQuant8Array(angle.(data));
       # Save compressed phase and amplitude to outout file
       println("Compressed amplitude saving as ",infile,"_","A","_",linlogbit,".h5")
       println("Compressed phase saving as ",infile,"_","ph","_",linlogbit,".h5")
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"A",log8A.A);
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"max",log8A.max);
       h5write(string(infile,"_","A","_",linlogbit,".h5"),"min",log8A.min);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"A",lin8ph.A);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"max",lin8ph.max);
       h5write(string(infile,"_","ph","_",linlogbit,".h5"),"min",lin8ph.min);
     
     else
       println("Error, number of quantization bit can only be 16 or 8")
     end
   end
end
