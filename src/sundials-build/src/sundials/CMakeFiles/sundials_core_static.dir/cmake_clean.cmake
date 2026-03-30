file(REMOVE_RECURSE
  "libsundials_core.a"
  "libsundials_core.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_core_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
