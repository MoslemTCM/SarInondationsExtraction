set(OTBAppSarInondationsExtraction_LINK_LIBS
  ${OTBTextures_LIBRARIES}
  ${OTBImageBase_LIBRARIES}
  ${OTBObjectList_LIBRARIES}
  ${OTBApplicationEngine_LIBRARIES}
  ${OTBImageManipulation_LIBRARIES}
  ${OTBGDAL_LIBRARIES}
)

set(OTBlib_c
   bimage.c
   path_queue.cxx
   path_support.c
   pathopen.cxx
  )


OTB_CREATE_APPLICATION(NAME SarInondationsExtraction
                       SOURCES SarInondationsExtraction.cxx ${OTBlib_c}
					   LINK_LIBRARIES ${${otb-module}_LIBRARIES})
					   
                      
#add_library(otbapp_SarRiverDetection ${OTBlib_c})