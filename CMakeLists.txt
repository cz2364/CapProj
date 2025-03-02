cmake_minimum_required(VERSION 3.14...3.22)
project(CapProj
        VERSION 1.0
        LANGUAGES CXX)

# Add Google Test as a submodule (or use find_package)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()


# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Create library
add_library(CapProj src/defs.cpp)
target_include_directories(CapProj PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


# Create test executable
add_executable(CapProjTest test/defs_test.cpp)
target_link_libraries(CapProjTest CapProj GTest::gtest_main)

# Enable testing with CTest
include(GoogleTest)
enable_testing()
gtest_discover_tests(CapProjTest)
