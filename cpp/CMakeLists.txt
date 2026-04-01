# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 3.10 )

# discover the compilers from PATH (when system install is overridden)
find_program(gcc_loc gcc)
find_program(gpp_loc g++)

set(CMAKE_C_COMPILER ${gcc_loc})
set(CMAKE_CXX_COMPILER ${gpp_loc})

# declare the project name
PROJECT(GBL)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 4 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )

# make life easier and simply use the ilcsoft default settings
# load default ilcsoft settings (install prefix, build type, rpath, etc.)
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
INCLUDE( ilcsoft_default_settings )

# enable ROOT support, check required CXX standard (see below)
OPTION( SUPPORT_ROOT "Support ROOT for user I/O" False )
IF( SUPPORT_ROOT )
  FIND_PACKAGE(ROOT REQUIRED CONFIG)
  ADD_DEFINITIONS(-DGBL_EIGEN_SUPPORT_ROOT)
ENDIF()
FIND_PACKAGE(Eigen3 REQUIRED CONFIG)

find_package(Mille QUIET )
if (${Mille_FOUND})
    message ("Found existing Mille installation")
else()
    message ("No existing Mille installation found - will clone and build the recommended version")
    include ("cmake/getMille.cmake")
endif()

# require proper c++
#ADD_DEFINITIONS( "-Wall -g -ansi -pedantic -Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL( DIRECTORY ./include DESTINATION . PATTERN ".svn" EXCLUDE )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

IF( SUPPORT_ROOT )
  # Use imported targets for correct transitive dependency handling.
  # ROOT::Matrix provides TMatrixD/TMatrixDSym/TMatrixTBase, ROOT::MathCore provides TVectorD.
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} Eigen3::Eigen ROOT::Core ROOT::Matrix ROOT::MathCore  Mille2)
ELSE()
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} Eigen3::Eigen Mille2)
ENDIF()

OPTION(JNA_DEBUG "Debug printouts for helping debug JNA memory issues" OFF)
IF( JNA_DEBUG )
  TARGET_COMPILE_DEFINITIONS(${PROJECT_NAME} PUBLIC JNA_DEBUG)
ENDIF()
# JNA_MEMORY_MONITOR is a subset of JNA_DEBUG,
#   basically, it turns off the runtime printing and only keeps track of counts
#   printing the final tally when the program is finished running
#   this is very helpful when developing a java program using GBL via JNA
#   to make sure your java program is cleaning up the memory allocated by GBL
OPTION(JNA_MEMORY_MONITOR 
  "Only keep track of structures allocated and destructed within JNA wrappers" OFF)
IF( JNA_MEMORY_MONITOR )
  TARGET_COMPILE_DEFINITIONS(${PROJECT_NAME} PUBLIC JNA_MEMORY_MONITOR)
ENDIF()

# the examples, build with "make examples"
ADD_SUBDIRECTORY( examples/ )
  
# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)

    INCLUDE(GlobalSVNRevision)
    
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generatred!" )

ENDIF(DOXYGEN_FOUND)

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( GBLConfig.cmake GBLConfigVersion.cmake )
