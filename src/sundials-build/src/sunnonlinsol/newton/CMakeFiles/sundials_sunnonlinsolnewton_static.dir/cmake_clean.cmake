file(REMOVE_RECURSE
  "libsundials_sunnonlinsolnewton.a"
  "libsundials_sunnonlinsolnewton.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_sunnonlinsolnewton_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
