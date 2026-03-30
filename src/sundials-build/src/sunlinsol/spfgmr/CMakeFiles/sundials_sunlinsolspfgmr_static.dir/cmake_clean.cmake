file(REMOVE_RECURSE
  "libsundials_sunlinsolspfgmr.a"
  "libsundials_sunlinsolspfgmr.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunlinsolspfgmr_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
