file(REMOVE_RECURSE
  "libsundials_arkode.a"
  "libsundials_arkode.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_arkode_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
