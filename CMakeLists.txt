PROJECT (mincbeast)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
 cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

SET(BEAST_SRC 
  alloc.c beast_lib.c distance_patch.c moments.c nlm_utils.c nlmfilter.c nlmsegFuzzy.c label.c ParseArgv.c
)

OPTION(USE_MINC      "Use minc1/minc2 library" ON)

IF(USE_MINC)
  SET(BEAST_SRC ${BEAST_SRC} mincio.c)
  SET(MINC_ROOT "/usr/local/bic" CACHE PATH "minc library base")
  ADD_DEFINITIONS( -DHAVE_MINC )
  OPTION(USE_MINC2      "Use minc 2 library" ON)
  LINK_DIRECTORIES(
    "${MINC_ROOT}/lib"
  )
  INCLUDE_DIRECTORIES( 
    "${MINC_ROOT}/include"
  )
  IF(USE_MINC2)
    ADD_DEFINITIONS( -DMINC2 )
    LINK_LIBRARIES(volume_io2 minc2 netcdf hdf5 z m)
  ELSE(USE_MINC2)
    LINK_LIBRARIES(volume_io minc netcdf m)
  ENDIF(USE_MINC2)
ENDIF(USE_MINC)

OPTION(USE_NIFTI      "Use niftiio library" ON)
OPTION(USE_DMALLOC    "Use dmalloc library to debug memory leaks" OFF)

IF(USE_NIFTI)
  ADD_DEFINITIONS( -DHAVE_NIFTI )
  LINK_LIBRARIES(niftiio znz z m)
  SET(NIFTI_ROOT "/usr/local/bic" CACHE PATH "niftiio library base")
  SET(BEAST_SRC ${BEAST_SRC} niftiio.c)
  LINK_DIRECTORIES(
    "${NIFTI_ROOT}/lib"
  )
  INCLUDE_DIRECTORIES( 
    "${NIFTI_ROOT}/include/nifti"
  )
ENDIF(USE_NIFTI)


IF(USE_DMALLOC)
  ADD_DEFINITIONS( -DUSE_DMALLOC -DDMALLOC_FUNC_CHECK  )
  LINK_LIBRARIES(dmalloc)
ENDIF(USE_DMALLOC)

ADD_LIBRARY(beast 
  ${BEAST_SRC}
)

ADD_EXECUTABLE(mincbeast
  mincbeast.c
)

TARGET_LINK_LIBRARIES(mincbeast beast)

IF(USE_MINC)
IF(USE_NIFTI)
  ADD_EXECUTABLE(minc2nifti
    minc2nifti.c
    )
  TARGET_LINK_LIBRARIES(minc2nifti beast)
ENDIF(USE_NIFTI)
ENDIF(USE_MINC)


INSTALL(TARGETS 
        mincbeast
        RUNTIME DESTINATION bin )