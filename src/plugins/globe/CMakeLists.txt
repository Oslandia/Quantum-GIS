# set path to additional CMake modules
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH})

FIND_PACKAGE(OSG REQUIRED)
FIND_PACKAGE(OSGEARTH REQUIRED)
FIND_PACKAGE(Threads)

IF(HAVE_OSGEARTH_CHILD_SPACING)
  ADD_DEFINITIONS(-DHAVE_OSGEARTH_CHILD_SPACING)
ENDIF(HAVE_OSGEARTH_CHILD_SPACING)

########################################################
# Files

SET (globe_plugin_SRCS
     globe_plugin.cpp
     qgsosgearthtilesource.cpp
     globe_plugin_dialog.cpp
     osgEarthQt/ViewerWidget.cpp
     osgEarthUtil/Controls.cpp
)

SET (globe_plugin_UIS
     globe_plugin_dialog_guibase.ui
     )

SET (globe_plugin_MOC_HDRS
     globe_plugin.h
     globe_plugin_dialog.h
)

SET (globe_plugin_RCCS  globe_plugin.qrc)

########################################################
# Build

QT4_WRAP_UI (globe_plugin_UIS_H  ${globe_plugin_UIS})

QT4_WRAP_CPP (globe_plugin_MOC_SRCS  ${globe_plugin_MOC_HDRS})

QT4_ADD_RESOURCES(globe_plugin_RCC_SRCS ${globe_plugin_RCCS})

ADD_LIBRARY (globeplugin MODULE ${globe_plugin_SRCS} ${globe_plugin_MOC_SRCS} ${globe_plugin_RCC_SRCS} ${globe_plugin_UIS_H})

INCLUDE_DIRECTORIES(
     ${CMAKE_CURRENT_BINARY_DIR}
     ${OSGEARTH_INCLUDE_DIR}
     ${GEOS_INCLUDE_DIR}
     ../../core ../../core/raster ../../core/renderer ../../core/symbology ../../core/symbology-ng
     ../../gui
     ..
)

TARGET_LINK_LIBRARIES(globeplugin
  qgis_core
  qgis_gui
  ${QT_QTOPENGL_LIBRARY}
  ${OSGDB_LIBRARY}
  ${OSGGA_LIBRARY}
  ${OSGUTIL_LIBRARY}
  ${OSG_LIBRARY}
  ${OSGQT_LIBRARY}
  ${OSGVIEWER_LIBRARY}
  ${OSGEARTH_LIBRARY}
  ${OSGEARTHFEATURES_LIBRARY}
  ${OSGEARTHSYMBOLOGY_LIBRARY}
  ${OSGEARTHUTIL_LIBRARY}
  ${OPENTHREADS_LIBRARY}
)

###################################################
# osgEarth plugin
###################################################

SET(OSGEARTH_PLUGIN_PREFIX "")
IF(UNIX AND NOT WIN32 AND NOT APPLE)
  IF(CMAKE_SIZEOF_VOID_P MATCHES "8")
    SET(LIB_POSTFIX "64" CACHE STRING "suffix for 32/64 dir placement")
    MARK_AS_ADVANCED(LIB_POSTFIX)
  ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
ENDIF(UNIX AND NOT WIN32 AND NOT APPLE)
IF(NOT DEFINED LIB_POSTFIX)
  SET(LIB_POSTFIX "")
ENDIF(NOT DEFINED LIB_POSTFIX)
SET( DYNAMIC_OSGEARTH "ON" )
SET(OSGEARTH_USER_DEFINED_DYNAMIC_OR_STATIC "SHARED")

INCLUDE( OsgEarthMacroUtils )

SET(TARGET_SRC
    QGISFeatureSource.cpp
)

SET(TARGET_H
    QGISFeatureOptions
)

SET(TARGET_COMMON_LIBRARIES ${TARGET_COMMON_LIBRARIES}
  qgis_core
  qgis_gui
#  ${QT_QTOPENGL_LIBRARY}
  ${OSGDB_LIBRARY}
#  ${OSGGA_LIBRARY}
#  ${OSGUTIL_LIBRARY}
  ${OSG_LIBRARY}
#  ${OSGQT_LIBRARY}
#  ${OSGVIEWER_LIBRARY}
  ${OSGEARTH_LIBRARY}
  ${OSGEARTHFEATURES_LIBRARY}
  ${OSGEARTHSYMBOLOGY_LIBRARY}
  ${OSGEARTHUTIL_LIBRARY}
  ${OPENTHREADS_LIBRARY}
)
SETUP_PLUGIN(osgearth_feature_qgis)

# to install public driver includes:
SET(LIB_NAME feature_qgis)
SET(LIB_PUBLIC_HEADERS ${TARGET_H})
INCLUDE(ModuleInstallOsgEarthDriverIncludes OPTIONAL)


########################################################
# Install

INSTALL(TARGETS globeplugin
  RUNTIME DESTINATION ${QGIS_PLUGIN_DIR}
  LIBRARY DESTINATION ${QGIS_PLUGIN_DIR})

INSTALL (FILES images/world.tif
         DESTINATION ${QGIS_DATA_DIR}/globe)

INSTALL (DIRECTORY images/gui
         DESTINATION ${QGIS_DATA_DIR}/globe)
