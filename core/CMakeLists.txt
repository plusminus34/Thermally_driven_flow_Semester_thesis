SET(folder_package_depend )
IF(WIN32)
SET(folder_link_libs)
ELSE() 
SET(folder_link_libs ${VTK_LIBRARIES} ${NETCDF_LIBRARIES})
ENDIF()
SET(folder_include_dirs ${VTK_USE_FILE})
SET(folder_lib_depend)

SET(current_dir_lists ${CMAKE_CURRENT_LIST_FILE})

IF(WIN32)
  SET(folder_lib_type STATIC)
ENDIF()

INCLUDE(${CMAKE_HOME_DIRECTORY}/cmake/BuildSrcFolderGenericLists.txt NO_POLICY_SCOPE)
