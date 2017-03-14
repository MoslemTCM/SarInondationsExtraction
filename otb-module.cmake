set(DOCUMENTATION "Extraction des zones inondées en utilisant les images radars.")

otb_module(OTBAppSarInondationsExtraction

  DEPENDS
    OTBITK
    OTBApplicationEngine
	OTBGdalAdapters
    OTBApplicationEngine
    OTBImageBase
    OTBCommon
    OTBImageManipulation
	OTBTextures
	OTBEdge
	
	
  TEST_DEPENDS
    OTBTestKernel
    OTBCommandLine

  DESCRIPTION
    "${DOCUMENTATION}"
)
