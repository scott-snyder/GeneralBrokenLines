# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 3.1 )

# declare the project name
PROJECT(GBL VERSION "2.1.3" LANGUAGES CXX)

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

# Require a C++11 compiler but do allow extensions
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

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

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} Eigen3::Eigen)

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

# generate and install following configuration files
INCLUDE("cmake/MacroGeneratePackageConfigFiles.cmake")
GENERATE_PACKAGE_CONFIGURATION_FILES( GBLConfig.cmake GBLConfigVersion.cmake )
