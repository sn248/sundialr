file(REMOVE_RECURSE
  "libsundials_cvode.a"
  "libsundials_cvode.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_cvode_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
