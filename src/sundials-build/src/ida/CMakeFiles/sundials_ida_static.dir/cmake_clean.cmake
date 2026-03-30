file(REMOVE_RECURSE
  "libsundials_ida.a"
  "libsundials_ida.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/sundials_ida_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
