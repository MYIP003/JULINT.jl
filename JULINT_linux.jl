module JULINT
  using Mmap
  using LightXML
  using LinLogQuantization
  using HDF5
  using DelimitedFiles

  include("/media/michelle/Michelle/JULINT/src/JULINT_r.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_i2pha.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_w.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_calamp.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_compress.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_selpsc.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_selsbc.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_set_patch.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_selpsc_line_by_line.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_imagemath.jl")
  include("/media/michelle/Michelle/JULINT/src/JULINT_mt_prep_isce.jl")

end

