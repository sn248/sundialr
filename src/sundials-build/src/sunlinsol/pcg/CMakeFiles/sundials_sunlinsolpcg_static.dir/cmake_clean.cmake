file(REMOVE_RECURSE
  "libsundials_sunlinsolpcg.a"
  "libsundials_sunlinsolpcg.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunlinsolpcg_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
