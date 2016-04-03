cmake_minimum_required(VERSION 2.8)
 
project(MSMNoduleDetection)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

include_directories(../../Common)

add_executable(MSMNoduleDetection MSMNoduleDetection.cxx)
target_link_libraries(MSMNoduleDetection  ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()