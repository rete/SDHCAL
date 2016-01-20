#################################################
# cmake file for building SDHCAL utils
# @author Eté Rémi ,IPNL
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# project name
PROJECT( SDHCAL )

SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# --------------- ilcsoft dependencies -----------
LIST( APPEND CMAKE_MODULE_PATH $ILCSOFT )

# load default settings from ILCSOFT_CMAKE_MODULES
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
INCLUDE( ilcsoft_default_settings )

# -------------- packages -------------
FIND_PACKAGE( Marlin REQUIRED )
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( ROOT 5.26.00 REQUIRED COMPONENTS Core Eve Geom RGL EG ) # minimum required ROOT version

FOREACH( pkg Marlin ROOT LCIO CLHEP )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

# include directories
INCLUDE_DIRECTORIES( include )
INCLUDE_DIRECTORIES( include/processors )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")

AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )

INSTALL (
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION lib
)

MACRO ( ADD_MARLIN_PROCESSOR _name _processor_src_file )

	ADD_SHARED_LIBRARY ( ${_name} ${_processor_src_file} )
	
	INSTALL (
		TARGETS ${_name}
		LIBRARY DESTINATION lib
	)
	
	TARGET_LINK_LIBRARIES( ${_name} ${PROJECT_NAME} )
	
ENDMACRO ( ADD_MARLIN_PROCESSOR )

# build processors
ADD_MARLIN_PROCESSOR( CutProcessor src/processors/CutProcessor.cc )
ADD_MARLIN_PROCESSOR( OverlayEstimatorProcessor src/processors/OverlayEstimatorProcessor.cc )
ADD_MARLIN_PROCESSOR( SDHCALProcessor src/processors/SDHCALProcessor.cc )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
