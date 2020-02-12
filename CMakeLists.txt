# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 3.1 )

# Set default build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# declare the project name
PROJECT(GBL VERSION "2.2.1" LANGUAGES CXX)

###############################
# Setup the build environment #
###############################

# Additional packages to be searched for by cmake
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Configure the installation prefix to allow both local as system-wide installation
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# Set up the RPATH so executables find the libraries even when installed in non-default location
SET(CMAKE_MACOSX_RPATH TRUE)
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Add the automatically determined parts of the RPATH which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" IsSystemDir)
IF("${IsSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${IsSystemDir}" STREQUAL "-1")

######################
# Define build flags #
######################

# Set standard build flags
SET(COMPILER_FLAGS -pedantic -Wall -Wextra -Wcast-align -Wcast-qual   -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op   -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel  -Wswitch-default -Wundef -Wshadow -Wformat-security -Wdeprecated -fdiagnostics-color=auto )

# Set linker flags
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  SET ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-undefined,error")
ELSEIF(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
  SET ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
ENDIF()

# Require a C++11 compiler but do allow extensions
SET(CMAKE_CXX_STANDARD "11" CACHE STRING "Required C++ Standard")
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

INCLUDE("cmake/compiler-flag-checks.cmake")

############################
# Check for Prerequisistes #
############################

FIND_PACKAGE(Eigen3 REQUIRED)

# enable ROOT support
OPTION( SUPPORT_ROOT "Support ROOT for user I/O" True )
IF( SUPPORT_ROOT )
  FIND_PACKAGE(ROOT REQUIRED)
  ADD_DEFINITIONS(-DGBL_EIGEN_SUPPORT_ROOT)
  INCLUDE_DIRECTORIES(SYSTEM ${ROOT_INCLUDE_DIRS})
ENDIF()

INSTALL( DIRECTORY ./include DESTINATION . )

INCLUDE_DIRECTORIES(include)

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_LIBRARY(${PROJECT_NAME} SHARED ${library_sources} )
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} )

# the examples, build with "make examples"
ADD_SUBDIRECTORY( examples/ )

# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
ELSE()
    MESSAGE( STATUS "Doxygen not installed or found, no API documentation can be generated" )
ENDIF()

# generate and install following configuration files
INCLUDE("cmake/MacroGeneratePackageConfigFiles.cmake")
GENERATE_PACKAGE_CONFIGURATION_FILES( GBLConfig.cmake GBLConfigVersion.cmake )
