cmake_minimum_required(VERSION 3.8)

add_library(mexbind0x INTERFACE)
target_include_directories(mexbind0x INTERFACE .)
target_compile_definitions(mexbind0x INTERFACE MATLAB_MEX_FILE)
add_library(mexbind0x::mexbind0x ALIAS mexbind0x)
