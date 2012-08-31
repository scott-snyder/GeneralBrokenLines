# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

# declare the project name
PROJECT(GBL)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# make life easier and simply use the ilcsoft default settings
# load default ilcsoft settings (install prefix, build type, rpath, etc.)
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# include helper macros by ilcsoft/Jan Engels, DESY
INCLUDE( MacroDisplayStandardVariables )
INCLUDE( MacroGeneratePackageConfigFiles )

# set default install prefix to project root directory
# instead of the cmake default /usr/local
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

# write this variable to cache
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )


FIND_PACKAGE(ROOT REQUIRED)

# require proper c++
ADD_DEFINITIONS( "-Wall -g -O0 -ansi -pedantic -Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( BEFORE ./include ${ROOT_INCLUDE_DIRS})
INSTALL( DIRECTORY ./include DESTINATION . PATTERN ".svn" EXCLUDE )

# declare the stupid ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${library_sources} )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES
        # create *nix style library versions + symbolic links
        VERSION ${${PROJECT_NAME}_VERSION}
        SOVERSION ${${PROJECT_NAME}_SOVERSION}
    )

INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBS})

# the examples, build with "make examples"
ADD_SUBDIRECTORY( examples/ )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( GBLConfig.cmake GBLConfigVersion.cmake )
