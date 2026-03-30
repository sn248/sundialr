file(REMOVE_RECURSE
  "libsundials_sunmatrixband.a"
  "libsundials_sunmatrixband.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunmatrixband_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
