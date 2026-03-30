file(REMOVE_RECURSE
  "libsundials_cvodes.a"
  "libsundials_cvodes.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_cvodes_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
