# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 3.14 )

# discover the compilers from PATH (when system install is overridden)
find_program(gcc_loc gcc)
find_program(gpp_loc g++)

set(CMAKE_C_COMPILER ${gcc_loc})
set(CMAKE_CXX_COMPILER ${gpp_loc})

# declare the project name
PROJECT(GBL)

## ========= Versioning ==========
set( ${PROJECT_NAME}_VERSION_MAJOR 4 )
set( ${PROJECT_NAME}_VERSION_MINOR 0 )
set( ${PROJECT_NAME}_VERSION_PATCH 2 )
set( ${PROJECT_NAME}_SOVERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}" )
set( ${PROJECT_NAME}_VERSION   "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_VERSION_PATCH}" )



## ========= Options ==========
option( SUPPORT_ROOT "Support ROOT for user I/O" False )
option(JNA_DEBUG "Debug printouts for helping debug JNA memory issues" OFF)
option(JNA_MEMORY_MONITOR 
  "Only keep track of structures allocated and destructed within JNA wrappers" OFF)


## ========= Compiler flags and install prefix ==========
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -Wno-long-long -Wuninitialized -Wl,-no-undefined  -Wno-unused-function" )

# install prefix - default to a local install 
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX  PROPERTY VALUE "${PROJECT_BINARY_DIR}/GBLInstall" )
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)


## ========= Dependencies =====   

# obtain ROOT 
if( SUPPORT_ROOT )
  find_package(ROOT REQUIRED CONFIG)
  add_definitions(-DGBL_EIGEN_SUPPORT_ROOT)
endif()

# Eigen3 
find_package(Eigen3 REQUIRED CONFIG)

# Mille 
find_package(Mille QUIET )
if (${Mille_FOUND})
    message ("Found existing Mille installation")
else()
    message ("No existing Mille installation found - will clone and build the recommended version")
    include ("cmake/getMille.cmake")
endif()

# ======== Configure the actual build ===== 

add_library(${PROJECT_NAME} SHARED)
target_sources(${PROJECT_NAME} 
              PRIVATE 
                src/BorderedBandMatrix.cpp  src/GblData.cpp  src/GblMeasurement.cpp  src/GblPoint.cpp  src/GblTrajectory.cpp  src/GblUtilities.cpp  src/JavaNativeAccessWrappers.cpp  src/VMatrix.cpp
              PUBLIC 
              FILE_SET HEADERS 
                BASE_DIRS include
                FILES include/BorderedBandMatrix.h  include/GblData.h  include/GblMeasurement.h  include/GblPoint.h  include/GblTrajectory.h  include/GblUtilities.h  include/VMatrix.h
              )

install( TARGETS ${PROJECT_NAME} DESTINATION lib )

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen Mille2)

if( SUPPORT_ROOT )
  # Use imported targets for correct transitive dependency handling.
  # ROOT::Matrix provides TMatrixD/TMatrixDSym/TMatrixTBase, ROOT::MathCore provides TVectorD.
  target_link_libraries(${PROJECT_NAME} PUBLIC ROOT::Core ROOT::Matrix ROOT::MathCore )
endif() 

if( JNA_DEBUG )
  target_compile_definitions(${PROJECT_NAME} PUBLIC JNA_DEBUG)
endif()
# JNA_MEMORY_MONITOR is a subset of JNA_DEBUG,
#   basically, it turns off the runtime printing and only keeps track of counts
#   printing the final tally when the program is finished running
#   this is very helpful when developing a java program using GBL via JNA
#   to make sure your java program is cleaning up the memory allocated by GBL
if( JNA_MEMORY_MONITOR )
  target_compile_definitions(${PROJECT_NAME} PUBLIC JNA_MEMORY_MONITOR)
endif()

# install GBL
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION lib
  FILE_SET HEADERS 
  DESTINATION include
)

# the examples
add_subdirectory( examples/ )
  
# ======== Configure documentation ===== 

# add a target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

else()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generatred!" )

endif(DOXYGEN_FOUND)

# ======== Add a setup script for the user ==== 


configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/setup.sh.cmake.in
    ${PROJECT_BINARY_DIR}/gblsetup.sh @ONLY)

install(FILES ${PROJECT_BINARY_DIR}/gblsetup.sh
         DESTINATION . )


# ======== Create and install package config and version files ===== 
configure_file(
${PROJECT_SOURCE_DIR}/cmake/GBLConfig.cmake.in
${PROJECT_BINARY_DIR}/cmake/GBLConfig.cmake @ONLY)

configure_file(
${PROJECT_SOURCE_DIR}/cmake/GBLConfig-version.cmake.in
${PROJECT_BINARY_DIR}/cmake/GBLConfig-version.cmake @ONLY)

install(FILES ${PROJECT_BINARY_DIR}/cmake/GBLConfig.cmake
        ${PROJECT_BINARY_DIR}/cmake/GBLConfig-version.cmake
        DESTINATION lib/cmake/${PROJECT_NAME}/)
