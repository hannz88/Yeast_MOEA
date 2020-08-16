cmake_minimum_required(VERSION 3.16)
project(yeastfitness)

set(CMAKE_CXX_STANDARD 17)
include_directories("/usr/lib/R/site-library/Rcpp/include/")
add_library(yeastfitness yeastfitness.cpp)