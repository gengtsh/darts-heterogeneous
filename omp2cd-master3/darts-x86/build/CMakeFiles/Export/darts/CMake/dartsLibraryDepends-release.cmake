#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "amm" for configuration "RELEASE"
set_property(TARGET amm APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(amm PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/lib64/libhwloc.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libamm.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS amm )
list(APPEND _IMPORT_CHECK_FILES_FOR_amm "${_IMPORT_PREFIX}/darts/lib/libamm.a" )

# Import target "common" for configuration "RELEASE"
set_property(TARGET common APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(common PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "-lpthread"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libcommon.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS common )
list(APPEND _IMPORT_CHECK_FILES_FOR_common "${_IMPORT_PREFIX}/darts/lib/libcommon.a" )

# Import target "rsmanager" for configuration "RELEASE"
set_property(TARGET rsmanager APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(rsmanager PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "common"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/librsmanager.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS rsmanager )
list(APPEND _IMPORT_CHECK_FILES_FOR_rsmanager "${_IMPORT_PREFIX}/darts/lib/librsmanager.a" )

# Import target "codelet" for configuration "RELEASE"
set_property(TARGET codelet APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(codelet PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libcodelet.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS codelet )
list(APPEND _IMPORT_CHECK_FILES_FOR_codelet "${_IMPORT_PREFIX}/darts/lib/libcodelet.a" )

# Import target "threadlocal" for configuration "RELEASE"
set_property(TARGET threadlocal APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(threadlocal PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libthreadlocal.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS threadlocal )
list(APPEND _IMPORT_CHECK_FILES_FOR_threadlocal "${_IMPORT_PREFIX}/darts/lib/libthreadlocal.a" )

# Import target "darts" for configuration "RELEASE"
set_property(TARGET darts APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(darts PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "common;amm;rsmanager;codelet;scheduler;threadlocal;rt;/lib64/libtbb.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libdarts.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS darts )
list(APPEND _IMPORT_CHECK_FILES_FOR_darts "${_IMPORT_PREFIX}/darts/lib/libdarts.a" )

# Import target "scheduler" for configuration "RELEASE"
set_property(TARGET scheduler APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(scheduler PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/darts/lib/libscheduler.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS scheduler )
list(APPEND _IMPORT_CHECK_FILES_FOR_scheduler "${_IMPORT_PREFIX}/darts/lib/libscheduler.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
