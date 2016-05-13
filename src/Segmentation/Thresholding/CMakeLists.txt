cmake_minimum_required(VERSION 2.8)
 
project(ThresholdingNoduleDetection)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

include_directories(../../Common)

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

add_executable(thresholding ThresholdingNoduleDetection.cxx)
target_link_libraries(thresholding  ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
