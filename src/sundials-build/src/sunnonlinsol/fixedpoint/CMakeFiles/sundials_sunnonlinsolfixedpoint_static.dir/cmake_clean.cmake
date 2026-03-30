file(REMOVE_RECURSE
  "libsundials_sunnonlinsolfixedpoint.a"
  "libsundials_sunnonlinsolfixedpoint.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunnonlinsolfixedpoint_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
