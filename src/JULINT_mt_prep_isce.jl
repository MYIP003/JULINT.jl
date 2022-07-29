######################################################## 
# Select PS candidates for compressed SLC files
# Author: Michelle Yip
# 09/03/2022 completed mt_prep_isce for PS(single primary image)
# 25/03/2022 Added small baselines option: 25-03-2022
########################################################
using DelimitedFiles

"""
```
julint_mt_prep_isce(da_thresh::Float64,rg_patches::Int64,az_patches::Int64,rg_overlap::Int64,az_overlap::Int64,nbits::Int64)
```
# Arguments
- `da_thresh::Float64` : amplitude dispersion threshold.
- `rg_patches::Int64` : number of patches in range.
- `az_patches::Int64` : number of patches in azimuth.
- `rg_overlap::Int64` : overlapping pixels between patches in range.
- `az_overlap::Int64` : overlapping pixels between patches in azimuth.
- `nbits::Int64` : number of bits of compressed data (32/16/8)
# Outputs
- `pscphase.in` : input file for pscphase.
- `psclonlat.in` : input file for psclonlat.
- `pscdem.in` : input file for pscdem.
- `selpsc.in` : input file for julint_selpsc.

# Examples
```julia-repl
julia> julint_mt_prep_isce(0.2,2,2,50,100,8)
```
"""
function julint_mt_prep_isce(da_thresh::Float64,rg_patches::Int64,az_patches::Int64,rg_overlap::Int64,az_overlap::Int64,nbits::Int64)

   println(string("Amplitude Dispersion Threshold:",da_thresh))
   println(string("Processing ",rg_patches," patch(es) in range and ",az_patches," in azimuth"))
   WORKDIR = pwd();
   dirname  = SubString(pwd(),length(pwd())-length("SMALL_BASELINES")+1,length(pwd())); # Both INSAR_20220312 and SMALLBASELINES have length=14, this will extract directory name of INSAR_20220312 for PS option and SMALL_BASELINES for SB option
 
   # Get INSAR directory
   if dirname == "SMALL_BASELINES"
     INSARDIR = SubString(pwd(),1:length(pwd())-length("SMALL_BASELINES")-1);
   else
     INSARDIR = WORKDIR;
   end

   # Retrieve width and length
   width = parse(Int64,readlines(joinpath(WORKDIR,"width.txt"))[1]); 
   line = parse(Int64,readlines(joinpath(WORKDIR,"len.txt"))[1]);
 
 # Create calamp.in
   if dirname == "SMALL_BASELINES"
     println("Small Baseline Processing")
     if nbits == 32
       calampinname = string(WORKDIR,"/calamp.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/*.slc"`));

     elseif nbits == 16
       calampinname = string(WORKDIR,"/calamp_compressed_16.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/*.slc_A_16.h5"`));

     elseif nbits == 8
       calampinname = string(WORKDIR,"/calamp_compressed_8.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/*.slc_A_8.h5"`));
     end
   
   else # if current directory is INSAR directory
     if nbits == 32
       calampinname = string(WORKDIR,"/calamp.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/master/master.slc"`));
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/slave.slc"`));

     elseif nbits == 16
       calampinname = string(WORKDIR,"/calamp_compressed_16.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/master/master.slc_A_16.h5"`));
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/slave.slc_A_16.h5"`));

     elseif nbits == 8
       calampinname = string(WORKDIR,"/calamp_compressed_8.in");
       calampinios = open(calampinname,"w");
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/master/master.slc_A_8.h5"`));
       writedlm(calampinios,readlines(`sh -c "ls $WORKDIR/*/slave.slc_A_8.h5"`));
     end
 
   end # if the current directory is SMALL_BASELINES
   close(calampinios);

   # Calculate interferograms from compressed SLC images pairs
   if dirname == "SMALL_BASELINES"
     JULINT.julint_imagemath_sb(calampinname,width,nbits)
   else
     JULINT.julint_imagemath_ps(calampinname,width,nbits)
   end


   # Calibrate amplitudes
   if nbits == 32
     calampoutname = "calamp.out";
   else  
     calampoutname = string("calamp_compressed_",nbits,".out");
   end

   JULINT.julint_calamp(calampinname,width,calampoutname,nbits)
   
   # Set up patches
   JULINT.julint_set_patch(rg_patches,az_patches,rg_overlap,az_overlap,width,line)
   cd("../") # Directory was in PATCH_*, go back to INSAR directory

   # generating the patch list
   patchlistname = string(WORKDIR,"/patch.list");
   patchlistios = open(patchlistname,"w");
   writedlm(patchlistios,readlines(`sh -c "ls -d PATCH_*"`))
   close(patchlistios);

#=
   # Dumping the interferogram (Create pscphase.in)
   # 1st line: width
   # From 2nd line : /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/*/isce_minrefdem.int
   pscphasename = string(WORKDIR,"/pscphase.in");
   pscphaseios = open(pscphasename,"w");
   write(pscphaseios,string(width,"\n"));
   if nbits == 32
     writedlm(pscphaseios,readlines(`sh -c "ls -d $WORKDIR/*/isce_minrefdem.int"`));
   elseif nbits == 16
     writedlm(pscphaseios,readlines(`sh -c "ls -d $WORKDIR/*/isce_minrefdem_16.int"`));
   elseif nbits == 8
   writedlm(pscphaseios,readlines(`sh -c "ls -d $WORKDIR/*/isce_minrefdem_8.int"`));
   end
   close(pscphaseios);

   # Dumping the geocoordinates (Create psclonlat.in)
   # 1st line: width
   # 2nd line: /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/lon.raw
   # 3rd line: /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/lat.raw
   psclonlatname = string(WORKDIR,"/psclonlat.in");
   psclonlatios = open(psclonlatname,"w");
   write(psclonlatios,string(width,"\n"));
   writedlm(psclonlatios,readlines(`sh -c "ls -d $INSARDIR/lon.raw"`));
   writedlm(psclonlatios,readlines(`sh -c "ls -d $INSARDIR/lat.raw"`));
   close(psclonlatios);   


   # Dumping the radar-coded DEM(Create pscdem.in)
   # 1st line: width
   # 2nd line: /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/dem.raw
   pscdemname = string(WORKDIR,"/pscdem.in");
   pscdemios = open(pscdemname,"w");
   write(pscdemios,string(width,"\n"));
   writedlm(pscdemios,readlines(`sh -c "ls -d $INSARDIR/dem.raw"`));
   close(pscdemios);
=#

   # mt_extract_cands
   if dirname == "SMALL_BASELINES"
     # Make selsbc.in
     JULINT.julint_make_selsbc(da_thresh,width,calampoutname,nbits)
     if nbits == 32
       selsbcinname = string(WORKDIR,"/selsbc.in");
     elseif nbits == 16 || nbits == 8
       selsbcinname = string(WORKDIR,"/selsbc_compressed_",nbits,".in");
     end

     # Read each PATCH in patch.list
     for i in 1:countlines("patch.list")
       cur_patch = readlines("patch.list")[i];
       println("Patch: ",cur_patch)
       cd(cur_patch)

       if nbits == 32
         JULINT.julint_selsbc(selsbcinname,"patch.in",false)
       else
         JULINT.julint_selsbc(selsbcinname,"patch.in",true)
       end
#=
       # Retrieve lon/lat for PS candidates (c script from StaMPS/src)
       run(`sh -c "psclonlat $WORKDIR/psclonlat.in pscands.1.ij pscands.1.ll"`);
       # Retrieve hgt for PS candidates
       run(`sh -c "pscdem $WORKDIR/pscdem.in pscands.1.ij pscands.1.hgt"`);
       # Retrieve phase for PS candidates
       run(`sh -c "pscphase $WORKDIR/pscphase.in pscands.1.ij pscands.1.ph"`);
=#
       cd("../") # go back to INSAR directory
     end # for i
   else
     # Make selpsc.in
     JULINT.julint_make_selpsc(da_thresh,width,calampoutname,nbits)
     if nbits == 32
       selpscinname = string(WORKDIR,"/selpsc.in");
     elseif nbits == 16 || nbits == 8
       selpscinname = string(WORKDIR,"/selpsc_compressed_",nbits,".in");
     end

     # Read each PATCH in patch.list
     for i in 1:countlines("patch.list")
       cur_patch = readlines("patch.list")[i];
       println("Patch: ",cur_patch)
       cd(cur_patch)

       if nbits == 32
         JULINT.julint_selpsc(selpscinname,"patch.in",false)
       else
         JULINT.julint_selpsc(selpscinname,"patch.in",true)
       end
#=       
       # Retrieve lon/lat for PS candidates (c script from StaMPS/src)
       run(`sh -c "psclonlat $WORKDIR/psclonlat.in pscands.1.ij pscands.1.ll"`);
       # Retrieve hgt for PS candidates
       run(`sh -c "pscdem $WORKDIR/pscdem.in pscands.1.ij pscands.1.hgt"`);
       # Retrieve phase for PS candidates
       run(`sh -c "pscphase $WORKDIR/pscphase.in pscands.1.ij pscands.1.ph"`);
=#
       cd("../") # go back to INSAR directory
     end # for i 
   end # if dirname == "SMALL_BASELINES"

    
end #end function

#=

infilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_8.in";
outfilename = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_8.out";
julint_calamp_compressed(infilename,outfilename)



julint_set_patch(2,2,50,200,width,line)

infilename=
infile_patch=

# Define variables
# path = "/media/michelle/Michelle/HK_ASC_2016_2021/StaMPS/INSAR_20181220"
# path = pwd();
path = "/Users/michelleymw/Desktop/Interferogram2/INSAR20220202/"



# Generate patch.list
julint_set_patch(2,2,50,200,25483,5533)
cd PATCH_*
### Generate pscands.1.ij pscands.1.da mean_amp.flt pscands.1.ij0
julint_selpsc(infilename,infile_patch,true)



# Generate pscands.1.hgt pscands.1.ll pscands.1.ph pscands.1.ij.int

   if ($dolonlat == 1) then
        # Retrieve lon/lat for PS candidates
        echo ""
        echo "psclonlat $WORKDIR/psclonlat.in pscands.1.ij pscands.1.ll"
        psclonlat $WORKDIR/psclonlat.in pscands.1.ij pscands.1.ll
    endif

psclonlat /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/psclonlat.in pscands.1.ij pscands.1.ll
opening pscands.1.ij...
width = 48592
opening /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/lon.raw...
opening /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/lat.raw...
100000 PS candidates processed
200000 PS candidates processed
300000 PS candidates processed
400000 PS candidates processed
500000 PS candidates processed

    if ($dodem == 1) then
        # Retrieve hgt for PS candidates
        echo ""
        echo "pscdem $WORKDIR/pscdem.in pscands.1.ij pscands.1.hgt"
        pscdem $WORKDIR/pscdem.in pscands.1.ij pscands.1.hgt
    endif
println(string("opening pscands.1.ij..."))
println(string("width = ",width))
println(string("opening "))
println("input file specified as single precision")

pscdem /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/pscdem.in pscands.1.ij pscands.1.hgt
opening pscands.1.ij...
pscdem: width = 48592
pscdem: opening /media/michelle/Michelle/HK_ASC_2015_2021/StaMPS/INSAR_20181220/dem.raw...
pscdem: input file specified as single precision


    if ($dophase == 1) then
        # Retrieve phase for PS candidates
        echo ""
        echo "pscphase $WORKDIR/pscphase.in pscands.1.ij pscands.1.ph"
        pscphase $WORKDIR/pscphase.in pscands.1.ij pscands.1.ph
    endif

=#

