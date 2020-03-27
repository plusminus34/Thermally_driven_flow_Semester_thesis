SET(folder_package_depend )
IF(WIN32)
	SET(folder_link_libs ${NETCDF_LIBRARIES} ${VTK_LIBRARIES})
ELSE()
	SET(folder_link_libs ${NETCDF_LIBRARIES} ${VTK_LIBRARIES} -lstdc++fs)
ENDIF()

SET(folder_include_dirs ${EIGEN3_INCLUDE_DIR} ${NETCDF_INCLUDES} ${VTK_USE_FILE})
SET(folder_lib_depend core )

SET(current_dir_lists ${CMAKE_CURRENT_LIST_FILE})

IF(WIN32)
  SET(folder_lib_type STATIC)
ENDIF()

INCLUDE(${CMAKE_HOME_DIRECTORY}/cmake/BuildSrcFolderGenericLists.txt NO_POLICY_SCOPE)

if(CMAKE_CL_64)
set(CMAKE_MSVCIDE_RUN_PATH ${CMAKE_SOURCE_DIR}/extern/x64)
else()
set(CMAKE_MSVCIDE_RUN_PATH ${CMAKE_SOURCE_DIR}/extern/x86)
endif()
 