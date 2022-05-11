module JULINT
  using Mmap
  using LightXML
  using LinLogQuantization
  using HDF5
  using DelimitedFiles

  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_r.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_i2pha.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_w.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_calamp.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_compress.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_selpsc.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_selsbc.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_set_patch.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_imagemath.jl")
  include("/Users/michelleymw/Desktop/PhD/JULINT/src/JULINT_mt_prep_isce.jl")

end

