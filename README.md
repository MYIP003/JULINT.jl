# JULINT
A Julia package for compressing coregistered Sentinel-1 single look complex (SLC) images from the InSAR Scientific Computing Environment (ISCE) software (Rosen et al., 2012), computing interferogram from the compressed images and select initial pixel candidates based on amplitude dispersion. The output files are saved in format that are readable by the Stanford Method for Persistent Scatterers (StaMPS) software (Hooper et al., 2004; Hooper, 2006) that extract time-series ground displacement of selected PS candidates from the interferograms. We also thank M. Klöwer for development of the LinLogQuantization.jl library package (https://github.com/milankl/LinLogQuantization.jl), which allowed us to explore rapidly the compression and decompression of the SAR data (Klöwer et al., 2021).

This package requires pre-installation of the following Julia package:
1. LinLogQuantization
2. Mmap
3. LightXML
4. HDF5
5. DelimitedFiles
