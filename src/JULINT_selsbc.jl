######################################################## 
# Select SB candidates for compressed SLC files
# Author: Michelle Yip
# Last Edited: 31-03-2022
########################################################
#using Mmap
#using LightXML
"""
Make `selsbc.in` or `selsbc_compressed_nbits.in` file contain amplitude dispersion threshold, width and calamp.out
```
julint_make_selsbc(da_thresh::Float64,width::Int64,calamp_out_name::String,compress_byte::Int64)
```
# Arguments
- `da_thresh::Float64` : amplitude dispersion threshold.
- `width::Int64` : width of SLC files.
- `calamp_out_name::String` : calamp.out file.
- `compress_byte::Int64` : No. of bytes of the input compressed SLCs
# Outputs
- `selpsc.in` : input file for julint_selpsc.

# Examples
```julia-repl
julia> calamp_out_name = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_16.out";
julia> julint_make_selsbc(0.3,25483,calamp_out_name,16)
```
"""
function julint_make_selsbc(da_thresh::Float64,width::Int64,calamp_out_name::String,compress_byte::Int64)

   if compress_byte == 32
     outfileios = open(string("selsbc.in"),"a")
   elseif compress_byte == 16 || compress_byte == 8
     outfileios = open(string("selsbc_compressed_",compress_byte,".in"),"a")
   end

   outfile = string(da_thresh,"\n");
   write(outfileios,outfile);
   write(outfileios,string(width,"\n"));
   write(outfileios,read(calamp_out_name,String));
   
   close(outfileios);

end
"""
Select small baselines candidates
```
julint_selsbc(infilename::String,infile_patch::String,flag_compress::Bool)
```
# Arguments
- `infilename::String` : selsbc.in file.
- `infile_patch::String` : patch.in file.
- `flag_compress::Bool` : false for uncompressed SLC, true for compressed SLC
# Outputs
- `pscands.1.ij` : PS candidate locations.
- `pscands.1.da` : PS candidate amplitude dispersion.
- `mean_amp.flt` : mean amplitude of every pixel.
# Examples
```julia-repl
Julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/selsbc_compressed_16.in";
julia> infile_patch = "/Users/michelleymw/Desktop/Interferogram2/SLC/Test/PATCH_1/patch.in";
julia> julint_selsbc(infilename,infile_patch,flag_compress);
```
"""
function julint_selsbc(infilename::String,infile_patch::String,flag_compress::Bool)

   # Flag_compress == true or false, true for compressed file, false for uncompressed file
   ampfilename = split(readlines(infilename)[3])[1]; # Strat reading file from 3rd line
   if flag_compress == true
        if occursin(".h5",ampfilename)
          #println("Processing compressed file ...")
        else
          println("The input files do not contain compressed SLCs, please change flag to false")
          return(0)
        end
   else # flag_compress == false
        if occursin(".h5",ampfilename)
          println("The input files is compressed SLCs, please change flag to true")
          return(0)
        else
          #println("Processing uncompressed file ...")
        end
   end

   # Define output file names
   ijname="pscands.1.ij";
   daoutname="pscands.1.da"
   meanoutname="mean_amp.flt"
 
   D_thresh = parse(Float64,readlines(infilename)[1]);
   #D_thresh_sq = D_thresh*D_thresh;
   width = parse(Int64,readlines(infilename)[2]);   
   println("dispersion threshold = ",D_thresh)
   println("width = ",width)
   num_files = countlines(infilename)-2; # minus the first two lines that only contain D_thresh and width
   # Initialize
   calib_factor = zeros(Float32,num_files);
   lines =  parse(Int64,readlines(joinpath("../len.txt"))[1]);
   abs_master_amp = abs(1);   


   # Initialize
   rg_start = parse(Int64,readlines(infile_patch)[1]);
   rg_end = parse(Int64,readlines(infile_patch)[2]);
   az_start = parse(Int64,readlines(infile_patch)[3]);
   az_end = parse(Int64,readlines(infile_patch)[4]);

   # Load in patch to the above varables
   # Determine size of a patch
   patch_lines = az_end - az_start + 1;
   patch_width = rg_end - rg_start + 1;

   linebytes = width*sizeof(ComplexF32); 		     # Bytes per line in SLC files
   patch_linebytes = patch_width*sizeof(ComplexF32);    # Bytes per line in each patch
   #patch_amp_linebytes =  patch_width*sizeofelement; # Only amplitude 
   
   println(string("number of lines per file = ",lines))
   println(string("patch lines = ", patch_lines))
   println(string("patch width = ", patch_width))
   
   # Initialize
   pscid = 0;  # PS candidate ID number
   pix_start = (az_start-1)*width+(rg_start-1); #define pixel number of start of 1st line of patch
   #pos_start = pix_start*sizeofelement*2;       #define position of start of 1st line of patch on SLC file
   
   sumamp= zeros(Float32,patch_lines,patch_width);
   sumampdiffsq = zeros(Float32,patch_lines,patch_width);

   # Open output files
   meanoutios = open(meanoutname,"w");
   ijnameios = open(ijname,"w");
   daoutnameios = open(daoutname,"w");
   
   #Fill calib_factor array
   for i in 1:num_files
     calib_factor[i] = parse(Float32,split(readlines(infilename)[i+2])[2]);
   end

   for i in 1:Int64(num_files/2)
     ampfilename1 = split(readlines(infilename)[i*2+1])[1]; # Start reading file from 3rd line
     ampfilename2 = split(readlines(infilename)[i*2+2])[1];
 
     println(string("Opening ",ampfilename1," [file ",i*2-1,"]","..."))
     println(string("Opening",ampfilename2," [file ",i*2,"]","..."))
     if flag_compress == true
       T = typeof(h5read(ampfilename1,"A")).parameters[1];
       N = typeof(h5read(ampfilename1,"A")).parameters[2];

       println("De-quantizing amplitudes"," from ",T," to"," Float32"," ...")
       # Read primary SLC
       data = LogQuantArray{T,N}(h5read(ampfilename1,"A"), h5read(ampfilename1,"min"), h5read(ampfilename1,"max"));
       ampfile1= Array(data);
       # Read secondary SLC
       data = LogQuantArray{T,N}(h5read(ampfilename2,"A"), h5read(ampfilename2,"min"), h5read(ampfilename2,"max"));
       ampfile2= Array(data);

     else # Advoid processing compressed file with flag_compress indicated as false
        # Read primary SLC 
        s = open(ampfilename1);
        ampfile1 = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines))); 
        # Read secondary SLC
        s = open(ampfilename2);
        ampfile2 = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
     end  # end of flag_compress


     for y in 1: patch_lines
       for x in 1: patch_width

         camp1 = ampfile1[az_start-1+y,rg_start-1+x];    # buffer(i*patch_width+x) ??
         camp2 = ampfile2[az_start-1+y,rg_start-1+x];

         amp1 = abs(camp1)/calib_factor[i*2-1];
         amp2 = abs(camp2)/calib_factor[i*2];

         sumamp[y,x] += amp1;
         sumamp[y,x] += amp2;
         sumampdiffsq[y,x] += (amp1-amp2)*(amp1-amp2);

        end # end for x loop
      end # end for y loop

     println(string(i, " SLC processed"))
     
    end # end for i loop

    # Write to mean_amp.flt
    write(meanoutios,transpose(sumamp));

    #LineInBuffer = 10;
    lines_per_file = lines;
    #BufferSize = num_files*linebytes*LineInBuffer;
    println(string("number of amplitude files = ",num_files))
    println(string("number of lines per file = ",lines_per_file))
    #println(string(num_files," files, ",linebytes," lines bytes, ",LineInBuffer," lines in the buffer")
    #println(string(Buffer size = ", BufferSize, " bytes"))


    for y in 1:patch_lines
      for x in 1:patch_width

	 if sumamp[y,x] > 0
           D_a = sqrt(sumampdiffsq[y,x]/(num_files/2))/(sumamp[y,x]/num_files); 

	   if D_a < D_thresh
             pscid += 1;
             
             # Write to pscands.1.ij
             write(ijnameios,string(pscid," ",(az_start-1)+y-1," ",(rg_start-1)+x-1,"\n"));
             
             # Write to pscands.1.da
             write(daoutnameios,string(D_a, "\n")) 

	   end # end for if D_a
          end # end for if sumamp
       end # end for x loop
     
     if y/100 == round(y/100)
       #println(string(y, " lines processed. Buffer contains "," lines"))
       println(string(y, " lines processed."))
     end

   end  # end for y for loop
 
 
   close(meanoutios);
   close(ijnameios);
   close(daoutnameios);

end # end of function
   

