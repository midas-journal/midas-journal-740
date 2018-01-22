set(DOCUMENTATION "This module implements cuberille implicit surface
polygonization for ITK. This method operates by diving the surface into a
number of small cubes called cuberilles. Each cuberille is centered at a
pixel lying on the iso-surface and then quadrilaterals are generated for each
face. The original approach is improved by projecting the vertices of each
cuberille onto the implicit surface, smoothing the typical block-like
resultant mesh.

A more detailed description can be found in the Insight Journal article:

  Mueller, D. \"Cuberille Implicit Surface Polygonization for ITK\"
  http://hdl.handle.net/10380/3186
  http://www.insight-journal.org/browse/publication/740
  July 20, 2010.
")

itk_module(Cuberille
  DEPENDS
    ITKCommon
    ITKImageFunction
    ITKImageGradient
    ITKMesh
  TEST_DEPENDS
    ITKTestKernel
    ITKQuadEdgeMesh
    ITKQuadEdgeMeshFiltering
    ITKThresholding
    ITKIOImageBase
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
)
