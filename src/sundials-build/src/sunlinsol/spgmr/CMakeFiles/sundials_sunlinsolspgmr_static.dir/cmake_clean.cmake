file(REMOVE_RECURSE
  "libsundials_sunlinsolspgmr.a"
  "libsundials_sunlinsolspgmr.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunlinsolspgmr_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
