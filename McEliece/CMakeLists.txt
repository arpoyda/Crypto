cmake_minimum_required(VERSION 3.10)
project(McElieceCrypto)

set(CMAKE_CXX_STANDARD 20)

add_compile_options(-fPIC)

add_library(McElieceCryptoLib SHARED GaloisField.cpp GaloisField.h McEliceCryptoC.h)

add_executable(McElieceCrypto main.cpp)

target_link_libraries(McElieceCrypto McElieceCryptoLib -lpthread)
