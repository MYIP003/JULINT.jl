######################################################## 
# Select PS candidates for compressed SLC files
# Author: Michelle Yip
# Last Edited: 31-03-2022
########################################################
#using Mmap
#using LightXML
"""
Make `selpsc.in` file contain amplitude dispersion threshold, width and calamp.out
```
julint_make_selpsc(da_thresh::Float64,width::Int64,calamp_out_name::String,compress_byte::Int64)
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
julia> julint_make_selpsc(0.3,25483,calamp_out_name,16)
```
"""
function julint_make_selpsc(da_thresh::Float64,width::Int64,calamp_out_name::String,compress_byte::Int64)

   if compress_byte == 32
     outfileios = open(string("selpsc.in"),"a")
   elseif compress_byte == 16 || compress_byte == 8
     outfileios = open(string("selpsc_compressed_",compress_byte,".in"),"a")
   end
    # T = typeof(h5read(infile,"A")).parameters[1];
    # N = typeof(h5read(infile,"A")).parameters[2];

   outfile = string(da_thresh,"\n");
   write(outfileios,outfile);
   write(outfileios,string(width,"\n"));
   write(outfileios,read(calamp_out_name,String));
   
   close(outfileios);

end
"""
Select persistent scatterer candidates
```
julint_selpsc(infilename::String,infile_patch::String,flag_compress::Bool)
```
# Arguments
- `infilename::String` : selpsc.in file.
- `infile_patch::String` : patch.in file.
- `flag_compress::Bool` : false for uncompressed SLC, true for compressed SLC
# Outputs
- `pscands.1.ij` : PS candidate locations.
- `pscands.1.da` : PS candidate amplitude dispersion.
- `mean_amp.flt` : mean amplitude of every pixel.
- `pscands.1.ij0`: location of pixel with zero amplitude.
# Examples
```julia-repl
Julia> infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/selpsc_compressed_16.in";
julia> infile_patch = "/Users/michelleymw/Desktop/Interferogram2/SLC/Test/PATCH_1/patch.in";
julia> julint_selpsc(infilename,infile_patch,flag_compress);
```
"""
function julint_selpsc(infilename::String,infile_patch::String,flag_compress::Bool)


   # Flag_compress == true of false, true for compressed file, false for uncompressed file
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
          println("The inpit files is compressed SLCs, please change flag to true")
          return(0)
        else
          #println("Processing uncompressed file ...")
        end
   end

   # Define output file names
   ijname="pscands.1.ij";
   #jiname="pscands.1.ij.int";
   daoutname="pscands.1.da"
   meanoutname="mean_amp.flt"
   jiname = string(ijname,".int");
   ijname0 = string(ijname,"0");
   println("File name for zero amplitude PS: ",ijname0);
   
   #
   D_thresh = parse(Float64,readlines(infilename)[1]);
   D_thresh_sq = D_thresh*D_thresh;
   width = parse(Int64,readlines(infilename)[2]);   
   println("dispersion threshold = ",D_thresh)
   println("width = ",width)
   num_files = countlines(infilename)-2; # minus the first two lines that only contain D_thresh and width
   # Initialize
   calib_factor = zeros(Float32,num_files);
   lines =  parse(Int64,readlines(joinpath("../len.txt"))[1]);
   abs_master_amp = abs(1);   


######### This section use up memory space when there is more SLC to be processed
# deleted
###########

   # Initialize
   rg_start = parse(Int64,readlines(infile_patch)[1]);
   rg_end = parse(Int64,readlines(infile_patch)[2]);
   az_start = parse(Int64,readlines(infile_patch)[3]);
   az_end = parse(Int64,readlines(infile_patch)[4]);

   # Load in patch to the above varables
   # Determine size of a patch
   patch_lines = az_end - az_start + 1;
   patch_width = rg_end - rg_start + 1;

   sizeoffloat = 4; # sizeof(Float32) = 4
   sizeofelement = sizeoffloat;
   linebytes = width*sizeofelement*2; 		     # Bytes per line in SLC files
   patch_linebytes = patch_width*sizeofelement*2;    # Bytes per line in each patch
   patch_amp_linebytes =  patch_width*sizeofelement; # Only amplitude 
   #ampfilesize = 5533;
   numline = 5533;
   println(string("number of lines per file = ",numline))
   println(string("patch lines = ", patch_lines))
   println(string("patch width = ", patch_width))
   
   #buffer = Array{ComplexF32}(undef,num_files*patch_linebytes);
   #masterline[1:25483] = 1;

   # Initialize
   pscid = 0;  # PS candidate ID number
   pix_start = (az_start-1)*width+(rg_start-1); #define pixel number of start of 1st line of patch
   pos_start = pix_start*sizeofelement*2;       #define position of start of 1st line of patch on SLC file
   
   sumamp= zeros(Float32,patch_lines,patch_width);
   sumampsq = zeros(Float32,patch_lines,patch_width);
   amp_0 = zeros(Float32,patch_lines,patch_width);


   meanoutios = open(meanoutname,"w");
   ijnameios = open(ijname,"w");
   daoutnameios = open(daoutname,"w");
   ijname0ios = open(ijname0,"w");
   #jinameios = open(jiname,"w");
   
   for i in 1:num_files
     ampfilename = split(readlines(infilename)[i+2])[1]; # Strat reading file from 3rd line
     calib_factor[i] = parse(Float32,split(readlines(infilename)[i+2])[2]);
     println(string("Opening",ampfilename, "..."))
     if flag_compress == true
       T = typeof(h5read(ampfilename,"A")).parameters[1];
       N = typeof(h5read(ampfilename,"A")).parameters[2];

       println("De-quantizing amplitudes"," from ",T," to"," Float32"," ...")
       data = LogQuantArray{T,N}(h5read(ampfilename,"A"), h5read(ampfilename,"min"), h5read(ampfilename,"max"));
       ampfile= Array(data);
     else # Advoid processing compressed file with flag_compress indicated as false 
        s = open(ampfilename);
        ampfile = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines))); 
     end  # end of flag_compress


     for y in 1: patch_lines
       for x in 1: patch_width

         camp = ampfile[az_start-1+y,rg_start-1+x];    # buffer(i*patch_width+x) ?? # Read which part?
         amp = abs(camp)/calib_factor[i]/abs_master_amp;  # What is master_amp? =1?
	 if amp <= 0.00005
           amp_0[y,x] = 1;
	   sumamp[y,x] = 0;
         else
	   sumamp[y,x] += amp;
	   sumampsq[y,x] += amp*amp;
	 end # end for if amp
        end # end for x loop
      end # end for y loop

     println(string(i, " SLC processed"))
     
    end # end for i loop

    # Write to mean_amp.flt
    write(meanoutios,transpose(sumamp));

    for y in 1:patch_lines
      for x in 1:patch_width

	 if sumamp[y,x] > 0
	   D_sq=num_files*sumampsq[y,x]/(sumamp[y,x]*sumamp[y,x]) - 1; 

	   if D_sq < D_thresh_sq
	     if amp_0[y,x] != 1
	       
	       pscid += 1;
	       #if y > 2501 && y < 2505 && x > 15001 && x <  15005
	       #println(string("y = ",y," and x = ",x," pscid = ",pscid))
	       #end
	       # Write to pscands.1.ij
               write(ijnameios,string(pscid," ",(az_start-1)+y-1," ",(rg_start-1)+x-1,"\n")); # minus one to be consistant with c code:selpsc.c
	       ##### save ji file (skip here, for gamma only)
	       #J = (rg_strat-1)+x-1;
	       #I = (az_start-1)+y-1;               

	       D_a = sqrt(D_sq);
	       # Write to pscands.1.da
	       write(daoutnameios,string(D_a, "\n"));

	     else
		# Write to pscands.1.ij0
	        write(ijname0ios,string(pscid," ",(az_start-1)+y-1," ",(rg_start-1)+x-1,"\n"));
	     end # end for amp_0 !=1
	   end # end for if D_sq
          end # end for if sumamp
       end # end for x loop

     if y/100 == round(y/100)
     println(string(y, " lines processed for pscands.1.ij, pscands.1.da, pscands.1.ij0"))
     end
   
#=
     # Read in next line in each amp file
     for i in 1:num_files
        ampfile[i,1,:];  #####################
     end # end for i for loop

     if y/100.0 == Int64(y/100.0)
       println(y," lines processed");
       println(pscid," selected pixels");
       println(D_thresh_sq," D_thresh_sq")
     end
=#
   #close(meanoutios);

   end  # end for y for loop
 
 
   close(meanoutios);
   close(ijnameios);
   close(daoutnameios);
   close(ijname0ios);
   #close(jinameios);

end # end of function
   

