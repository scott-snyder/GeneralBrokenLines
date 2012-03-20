# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )

# declare the project name
PROJECT(GBL)

# project version: *nix style MAJOR.MINOR.PATCH
SET( ${PROJECT_NAME}_SOVERSION "0.0.0" )
SET( ${PROJECT_NAME}_VERSION   "0.0.0" )

# general cmake settings
SET( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" CACHE PATH "Path to CMakeModules" FORCE)
# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# add install path to the rpath list (apple)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH ) 
### end general settings ####


FIND_PACKAGE(ROOT REQUIRED)

# set default install prefix to project root directory
# instead of the cmake default /usr/local
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()
# write this variable to cache
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )


# set cmake output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# project options -- build a shared library
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )

# declare the include directory for resolving the header locations
INCLUDE_DIRECTORIES(include ${ROOT_INCLUDE_DIRS})

# declare the stupid ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# add install path to the rpath list (apple)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# require proper c++
ADD_DEFINITIONS( "-Wall -g -ggdb -O3 -ansi -pedantic -Wno-long-long" )

AUX_SOURCE_DIRECTORY( src/ cc_srcs )

# build the shared library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${cc_srcs})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBS})

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/GBLConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/GBLConfig.cmake" @ONLY
)

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/GBLConfig.cmake" DESTINATION . )


# the examples, build with "make examples"
ADD_SUBDIRECTORY( examples/ )

# change lib_target properties
SET_TARGET_PROPERTIES( GBL PROPERTIES
        # create *nix style library versions + symbolic links
        VERSION ${${PROJECT_NAME}_VERSION}
        SOVERSION ${${PROJECT_NAME}_SOVERSION}
)


# install library
INSTALL( TARGETS GBL DESTINATION lib PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "GBLBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/GBLBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "GBLLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/GBLLibDeps.cmake" DESTINATION lib/cmake )
