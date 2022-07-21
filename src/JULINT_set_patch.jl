#"#######################################################
# Setup patches
# Author: Michelle Yip
# Last Edited: 24-09-2021
########################################################

"""
Make `patch.in` and `patch_noover.in` that contain range and azimuth boundary of each patch
```
julint_set_patch(prg::Int64,paz::Int64,overlap_rg::Int64,overlap_az::Int64,width::Int64,line::Int64)
```
# Arguments
- `prg::Int64` : number of patches in range.
- `paz::Int64`   : number of patches in azimuth.
- `overlap_rg::Int64` : overlapping pixels between patches in range.
- `overlap_az::Int64` : overlapping pixels between patches in azimuth.
- `width::Int64` : width of SLC files.
- `length::Int64` : length of SLC files.
# Outputs
- `patch.in` : input file for julint_selpsc.
- `patch_noover.in` : input file for julint_selpsc(?).

# Examples
```julia-repl
julia> julint_set_patch(2,2,50,200,25483,5533)
```
"""
function julint_set_patch(prg::Int64,paz::Int64,overlap_rg::Int64,overlap_az::Int64,width::Int64,line::Int64)
   # Initialize
   irg = 0; # Index of range patch
   iaz = 0; # Index of azimuth patch
   ip = 0;  # Index of overall patch
   width_p = floor(width/prg);
   length_p = floor(line/paz);
   WORKDIR = pwd();   

   while irg < prg
     irg = irg+1;
     while iaz < paz
       iaz = iaz+1;
       ip = ip+1;
       #println(ip)
       start_rg1 = width_p*(irg-1)+1;
       start_rg = start_rg1-overlap_rg;
       if start_rg < 1
         start_rg = 1;
       end
       end_rg1 = width_p*irg;
       end_rg = end_rg1+overlap_rg;
       if end_rg > width
         end_rg = width;
       end
       start_az1 = length_p*(iaz-1)+1;
       start_az = start_az1-overlap_az;
       if start_az < 1
	 start_az = 1;
       end
       end_az1 = length_p*iaz;
       end_az = end_az1+overlap_az;
       if end_az > line
	 end_az = line;
       end
       
       # Create PATCH_X folder and patch.in, patch_noover.in
       mkdir(string(WORKDIR,"/PATCH_",ip));
       cd(string(WORKDIR,"/PATCH_",ip));

       # Location of element of start and end of patch with overlap region (patch.in)
       outfileios = open(string(WORKDIR,"/PATCH_",ip,"/","patch.in"),"a")
       write(outfileios, string(Int(start_rg),"\n"));
       write(outfileios, string(Int(end_rg),"\n"));
       write(outfileios, string(Int(start_az),"\n"));
       write(outfileios, string(Int(end_az)));
       close(outfileios);
       
       # Location of element of start and end of patch without overlap region (patch_noover.in)
       outfileios = open(string(WORKDIR,"/PATCH_",ip,"/","patch_noover.in"),"a")
       write(outfileios, string(Int(start_rg1),"\n"));
       write(outfileios, string(Int(end_rg1),"\n"));
       write(outfileios, string(Int(start_az1),"\n"));
       write(outfileios, string(Int(end_az1)));
       close(outfileios);
     end
     iaz = 0; 
   end 
end

