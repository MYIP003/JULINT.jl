######################################################## 
# Select lon and lat from lon.raw and lat.raw given the selected 
#PS candidates from pscands.1.ji
# Author: Michelle Yip
# Last Edited: 07-01-2022
########################################################
"""
Make selpsc.in file the contain amplitude dispersion threshold
# Arguments
- `parmfile::String` : psclonlat.in
- `pscands.1.ij::String : location of permanent scatterer candidiates.
# Outputs
- `pscands.1.ll::String` : lon/lat of permanent scatterer candidiates.

# Examples
```julia-repl
julia> calamp_out_name = "/Users/michelleymw/Desktop/Interferogram2/SLC/calamp_compressed_16.out";
julia> julint_make_selpsc_compressed(0.3,25483,calamp_out_name,16)
```
"""

function julint_psclonlat(parmfile::String,pscands_ij_name::String,pscands_ll_name::String)

  width = parse(Float64,readlines(parmfile)[1]);
  lonfile = parse(Float64,readlines(parmfile)[2]);
  latfile = parse(Float64,readlines(parmfile)[3]);
  #pscands_ll_name = "pscands.1.ll";

  s = open(pscands_ij_name);
ampfile = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
camp= ampfile[y,:]; # Read line by line


end

function julint_pscphase(parmfile::String,pscands_ij_name::String,pscands_ll_name::String)

  width = parse(Float64,readlines(parmfile)[1]);
  lonfile = parse(Float64,readlines(parmfile)[2]);
  latfile = parse(Float64,readlines(parmfile)[3]);
  #pscands_ll_name = "pscands.1.ll";

  s = open(pscands_ij_name);
ampfile = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
camp= ampfile[y,:]; # Read line by line


end
