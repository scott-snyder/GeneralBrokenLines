# set default install prefix to project root directory
# instead of the cmake default /usr/local
IF( ${CMAKE_VERSION} VERSION_LESS "3.7.1" )
    # in versions of cmake less than 3.7.1, the only
    #   way to check if the install prefix was defaulted was to
    #   check if it was set to the default path
    IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
        SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
    ENDIF()
ELSE()
    # in versions of cmake later than 3.7.1, cmake added a variable
    #   specifically to check if the install location was set to the
    #   default path
    # this allows the user to specifically change the path to happen to
    #   correspond to the cmake default path
    IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        SET(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/GBLInstall" )
    ENDIF()
ENDIF()

# write this variable to cache
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )

