file(REMOVE_RECURSE
  "libsundials_sunlinsolspbcgs.a"
  "libsundials_sunlinsolspbcgs.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunlinsolspbcgs_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
