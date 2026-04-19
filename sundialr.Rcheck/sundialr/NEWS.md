sundialr v0.1.7
===============
* Changes in `/tools/cmake_call.sh` file to remove warnings with respect to stdout/stderr/fprintf etc in upstream C library
* Updated the underlying SUNDIALS library to v7.6.0 (Jan 2026)

sundialr v0.1.6
================
Updated build system to use cmake

sundialr v0.1.5
================
Bug fixes:
* Fixed the bug in assigning absolute tolerances

sundialr v0.1.4.1
================
Bug fixes:
* Fixed the linking problem due to multiple defined symbols

sundialr v0.1.4
================

New features:
* added the `cvsolve` function to solve ODEs with multiple discontinuities

Bug fixes:
* None

API change:
* NONE is existing functions
* A new function `cvsolve` is added

Internal changes:
* Now uses SUNDIALS version 5.2.0 (released March 2020) at the back end
