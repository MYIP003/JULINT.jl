########################################################
# Calculate interferogram from SLC images
# -e=a*conj(b)
# Author: Michelle Yip
# 30/03/2022 Added imagemath function for SB option
########################################################
#using Mmap
#using LightXML
#using LinLogQuantization
#using HDF5
"""
Calculate interferogram from compressed SLC image pairs (for generating Small Baselines interferograms)
# Arguments
- `infilename::String` : list of primary and secondary image (i.e. calamp.in), first line is the primary image.
- `width::Int64` : width of each SLC image.
- `nbits::Int64` : number of bits of compressed SLC images (8 or 16 or 32).
# Examples
```julia-repl
julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp.in";
julia> width = 25483;
julia> nbits = 8;
julia> julint_imagemath_sb(infilename,width,nbits)
```
"""
function julint_imagemath_sb(infilename::String,width::Int64,nbits::Int64)
 if nbits == 32

  for i in [1:2:length(readlines(infilename));]
  
    # Read primary image
    pri_infile = readlines(infilename)[1];
    s = open(pri_infile);
    println("Opening ",pri_infile," ...")
    lines = Int64(filesize(pri_infile)/(width*sizeof(ComplexF32)));
    pri_complex = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));

    # Read Secondary image
    sec_infile = readlines(infilename)[i+1];
    s = open(sec_infile);
    println("Opening ",sec_infile," ...")
    lines = Int64(filesize(sec_infile)/(width*sizeof(ComplexF32)));
    secondary_data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));

    # Calculate interferogram
    int_data = primary_data.*conj(secondary_data);
    int_data = transpose(int_data); # transpose to save as same dimension as .int files
    # Save interferogram
    write(intname,int_data);
    println("Saved as ",intname)
  end # end for loop

 else 

   for i in [1:2:length(readlines(infilename));] 

    # Read Primary image (Read line # 1,3,5,7...)
    pri_infile_A = readlines(infilename)[i];
    if nbits == 16
      pri_infile_ph = string(chop(readlines(infilename)[1],tail=length("_A_16.h5")),"_ph_",nbits,".h5");
    elseif nbits == 8
      pri_infile_ph = string(chop(readlines(infilename)[1],tail=length("_A_8.h5")),"_ph_",nbits,".h5");
    end
    println("Opening ",pri_infile_A, " ...")
    println("Opening ",pri_infile_ph, " ...")
    T = typeof(h5read(pri_infile_A,"A")).parameters[1];
    N = typeof(h5read(pri_infile_A,"A")).parameters[2];

    # De-quantize primary amplitude and phase linearly/logarithmically
    println("De-quantizing primary amplitude and phase"," from ",T," to"," Float32"," ...")
    pri_A = LogQuantArray{T,N}(h5read(pri_infile_A,"A"), h5read(pri_infile_A,"min"), h5read(pri_infile_A,"max"));
    pri_A = Array(pri_A);
    pri_ph = LinQuantArray{T,N}(h5read(pri_infile_ph,"A"), h5read(pri_infile_ph,"min"), h5read(pri_infile_ph,"max"));
    pri_ph = Array(pri_ph);

    # Reconstruct complex matrix (CFloat32) using amplitude and phase
    println("Reconstruct complex matrix using amplitude and phase")
    pri_r = pri_A.*cos.(pri_ph);
    pri_i = pri_A.*sin.(pri_ph);
    pri_complex = Complex.(pri_r,pri_i);


    # Read Secondary image (Read line # 2,4,6,8...)
    sec_infile_A = readlines(infilename)[i+1];
    if nbits == 16
      sec_infile_ph = string(chop(readlines(infilename)[i+1],tail=length("_A_16.h5")),"_ph_",nbits,".h5");
      sec_infile_position = chop(readlines(infilename)[i+1],tail=length("salve.slc_A_16.h5"));  # Remove last 8 characters("slave.slc_A_16.h5") from string
    elseif nbits == 8
      sec_infile_ph = string(chop(readlines(infilename)[i+1],tail=length("_A_8.h5")),"_ph_",nbits,".h5");
      sec_infile_position = chop(readlines(infilename)[i+1],tail=length("salve.slc_A_8.h5"));
    end 

    println("Opening ",sec_infile_A, " ...")
    println("Opening ",sec_infile_ph, " ...")
        
    # De-quantize primary amplitude and phase linearly/logarithmically
    println("De-quantizing secondary amplitude and phase"," from ",T," to"," Float32"," ...")

    sec_A = LogQuantArray{T,N}(h5read(sec_infile_A,"A"), h5read(sec_infile_A,"min"), h5read(sec_infile_A,"max"));
    sec_A = Array(sec_A);
    sec_ph = LinQuantArray{T,N}(h5read(sec_infile_ph,"A"), h5read(sec_infile_ph,"min"), h5read(sec_infile_ph,"max"));     
    sec_ph = Array(sec_ph);

    # Reconstruct complex matrix (CFloat32) using amplitude and phase
    println("Reconstruct complex matrix using amplitude and phase")
    sec_r = sec_A.*cos.(sec_ph);
    sec_i = sec_A.*sin.(sec_ph);
    sec_complex = Complex.(sec_r,sec_i); 

    # Calculate interferogram
    println("Calculating interferogram ... ")
    int_data = pri_complex.*conj(sec_complex);
    int_data = transpose(int_data); # transpose to save as same dimension as .int files
    
    write(string(sec_infile_position,"isce_minrefdem_",nbits,".int"),int_data);
    println(string("Interferogram saved as ",sec_infile_position,"isce_minrefdem_",nbits,".int"))

  end # and for loop
 end # end if nbits
end

"""
Calculate interferogram from compressed SLC image pairs (for generating Single Reference interferograms)
# Arguments
- `infilename::String` : list of primary and secondary image (i.e. calamp.in), first line is the primary image.
- `width::Int64` : width of each SLC image.
- `nbits::Int64` : number of bits of compressed SLC images (8 or 16 or 32).
# Examples
```julia-repl
julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp.in";
julia> width = 25483;
julia> nbits = 8;
julia> julint_imagemath(infilename,width,nbits)
```
"""
function julint_imagemath_ps(infilename::String,width::Int64,nbits::Int64) 
 if nbits == 32

  # Read primary image
  pri_infile = readlines(infilename)[1];
  s = open(pri_infile);
  println("Opening ",pri_infile," ...")
  lines = Int64(filesize(pri_infile)/(width*sizeof(ComplexF32)));
  pri_complex = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));

  for i in 1:length(readlines(infilename))-1
   
    # Read Secondary image
    sec_infile = readlines(infilename)[i+1];
    s = open(sec_infile);
    println("Opening ",sec_infile," ...")
    lines = Int64(filesize(sec_infile)/(width*sizeof(ComplexF32)));
    secondary_data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
 
    # Calculate interferogram
    int_data = primary_data.*conj(secondary_data);
    int_data = transpose(int_data); # transpose to save as same dimension as .int files
    # Save interferogram
    write(intname,int_data);
    println("Saved as ",intname)
  end # end for loop

 else 
  # Read Primary image
  pri_infile_A = readlines(infilename)[1];
  if nbits == 16
    pri_infile_ph = string(chop(readlines(infilename)[1],tail=length("_A_16.h5")),"_ph_",nbits,".h5");
  elseif nbits == 8
    pri_infile_ph = string(chop(readlines(infilename)[1],tail=length("_A_8.h5")),"_ph_",nbits,".h5");
  end
  println("Opening ",pri_infile_A, " ...")
  println("Opening ",pri_infile_ph, " ...")
  T = typeof(h5read(pri_infile_A,"A")).parameters[1];
  N = typeof(h5read(pri_infile_A,"A")).parameters[2];

  # De-quantize primary amplitude and phase linearly/logarithmically
  println("De-quantizing primary amplitude and phase"," from ",T," to"," Float32"," ...")
  pri_A = LogQuantArray{T,N}(h5read(pri_infile_A,"A"), h5read(pri_infile_A,"min"), h5read(pri_infile_A,"max"));
  pri_A = Array(pri_A);
  pri_ph = LinQuantArray{T,N}(h5read(pri_infile_ph,"A"), h5read(pri_infile_ph,"min"), h5read(pri_infile_ph,"max"));
  pri_ph = Array(pri_ph);

  # Reconstruct complex matrix (CFloat32) using amplitude and phase
  println("Reconstruct complex matrix using amplitude and phase")
  pri_r = pri_A.*cos.(pri_ph);
  pri_i = pri_A.*sin.(pri_ph);
  pri_complex = Complex.(pri_r,pri_i);  

  for i in 1:length(readlines(infilename))-1

    # Read Secondary image
    sec_infile_A = readlines(infilename)[i+1];
    if nbits == 16
      sec_infile_ph = string(chop(readlines(infilename)[i+1],tail=length("_A_16.h5")),"_ph_",nbits,".h5");
      sec_infile_position = chop(readlines(infilename)[i+1],tail=length("salve.slc_A_16.h5"));  # Remove last 8 characters("slave.slc_A_16.h5") from string
    elseif nbits == 8
      sec_infile_ph = string(chop(readlines(infilename)[i+1],tail=length("_A_8.h5")),"_ph_",nbits,".h5");
      sec_infile_position = chop(readlines(infilename)[i+1],tail=length("salve.slc_A_8.h5"));
    end 

    println("Opening ",sec_infile_A, " ...")
    println("Opening ",sec_infile_ph, " ...")
    #T = typeof(h5read(sec_infile_A,"A")).parameters[1];
    #N = typeof(h5read(sec)infile_A,"A")).parameters[2];
        
    # De-quantize primary amplitude and phase linearly/logarithmically
    println("De-quantizing secondary amplitude and phase"," from ",T," to"," Float32"," ...")

    sec_A = LogQuantArray{T,N}(h5read(sec_infile_A,"A"), h5read(sec_infile_A,"min"), h5read(sec_infile_A,"max"));
    sec_A = Array(sec_A);
    sec_ph = LinQuantArray{T,N}(h5read(sec_infile_ph,"A"), h5read(sec_infile_ph,"min"), h5read(sec_infile_ph,"max"));     
    sec_ph = Array(sec_ph);

    # Reconstruct complex matrix (CFloat32) using amplitude and phase
    println("Reconstruct complex matrix using amplitude and phase")
    sec_r = sec_A.*cos.(sec_ph);
    sec_i = sec_A.*sin.(sec_ph);
    sec_complex = Complex.(sec_r,sec_i); 

    # Calculate interferogram
    println("Calculating interferogram ... ")
    int_data = pri_complex.*conj(sec_complex);
    int_data = transpose(int_data); # transpose to save as same dimension as .int files
    #ph = angle.(int_data);        # Extract Phase
    #heatmap(ph[3000:5000,10000:15000],c = :cyclic_mygbm_30_95_c78_n256_s25)  
    # Save the complex interferogram (CFloat32)
    write(string(sec_infile_position,"isce_minrefdem_",nbits,".int"),int_data);
    println(string("Interferogram saved as ",sec_infile_position,"isce_minrefdem_",nbits,".int"))

    # Read saved interferogram
    # read_int = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
  end # end for loop
 end # end if nbits
end # end function
