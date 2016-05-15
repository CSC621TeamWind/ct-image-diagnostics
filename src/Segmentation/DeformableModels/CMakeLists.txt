cmake_minimum_required(VERSION 2.8)
 
project(MSMNoduleDetection)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

include_directories(../../Common)

add_executable(MSMNoduleDetection MSMNoduleDetection.cxx)
target_link_libraries(MSMNoduleDetection ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
