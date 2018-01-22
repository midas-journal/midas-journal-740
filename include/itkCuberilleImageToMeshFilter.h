/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkCuberilleImageToMeshFilter_h
#define itkCuberilleImageToMeshFilter_h

#define DEBUG_PRINT 0
#define USE_GRADIENT_RECURSIVE_GAUSSIAN 0
#define USE_ADVANCED_PROJECTION 0
#define USE_LINESEARCH_PROJECTION 0

#include "itkMacro.h"
#include "itkMesh.h"
#include "itkImageToMeshFilter.h"
#include "itkCellInterface.h"
#include "itkTriangleCell.h"
#include "itkQuadrilateralCell.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace itk
{

/** \class CuberilleImageToMeshFilter
 *
 * This filter uses the 'cuberille' method to convert an implicit surface
 * (image) to a mesh.
 *
 * The 'cuberille' model was proposed over 30 years ago [1,2].
 * A basic summary of the algorithm is as follows:
 * step over all pixels, for each pixel determine if it lies on the surface,
 * center a cube on the surface pixel, create quadrilateral faces aligned with
 * the cube (taking care of neighbouring pixels also on the surface), use a
 * gradient descent based method to project cube vertices to iso-surface [3].
 *
 * \par Parmeters
 * (Required) Input: specifies the input image containing the implicit surface
 * for polygonization. Currently, the input must NOT have iso-surface pixels
 * on the edge of the image. If this is the case, you MUST use
 * itkConstantPadImageFilter to pad the edges by at least 1 pixel.
 *
 * (Required) IsoSurfaceValue: specifies the value of the iso-surface for which
 * to generate the mesh. Pixels equal to or greater than this value are considered
 * to lie on the surface or inside the resultant mesh.
 *
 * (Optional) GenerateTriangleFaces: specifies whether triangle or
 * quadrilateral faces should be generated. The default is to generate
 * triangle faces.
 *
 * (Optional) ProjectVerticesToIsoSurface: specifies whether the vertices
 * should be projected onto the iso-surface. If projection is disabled, the
 * resultant mesh exhibits the traditional blocky features. Projection takes
 * roughly half of the algorithm time. The default is to project the vertices.
 *
 * (Optional) ProjectVertexSurfaceDistanceThreshold: specifies the threshold
 * for the 'distance' from iso-surface during vertex projection. Note that the
 * distance is actually measured in pixel value units (not space).
 * The smaller this value, the closer the vertices will be to the iso-surface.
 * Small values result in longer convergence time (i.e. slower).
 * Values are clamped to the range [0.0, max pixel value].
 * The default value is 0.5.
 *
 * (Optional) ProjectVertexStepLength: specifies the threshold for the step
 * length during vertex projection.
 * The smaller this value, the more likely the vertices will end up closer to
 * the surface. Small values cause the projection to take longer to converge.
 * Values are clamped to the range [0.0, large].
 * The default value is max spacing * 0.25 (expressed in physical space).
 *
 * (Optional) ProjectVertexStepLengthRelaxationFactor: specifies the step
 * length relaxation factor during vertex projection. The step length is
 * multiplied by this factor each iteration to allow convergence.
 * Values are clamped to the range [0.0, 1.0].The default value is 0.95.
 *
 * (Optional) ProjectVertexMaximumNumberOfSteps: specifies the maximum number
 * of steps used during vertex projection. The default value is 50.
 *
 * \par References
 * [1] G. Herman and H. Liu, "Three-dimensional Display of Human organs
 *     from Computed Tomograms", Computer Graphics and Images Processing,
 *     Volume 9, Issue 1, Pages 1-21, 1979.
 * [2] D. Gordon and J.K. Udupa, "Fast surface tracking in Three-dimensional
       Binary Images", Computer Vision, Graphics and Image Processing,
       Volume 45, Pages 196-214, 1989.
 * [3] http://www2.imm.dtu.dk/~jab/gallery/polygonization.html
 *
 * This implementation was taken from the Insight Journal:
 *     http://hdl.handle.net/10380/3186
 *
 * \author Dan Mueller, Philips Healthcare, dan dot muel at gmail dot com
 *
 * \ingroup Cuberille
 *
 */
template < class TInputImage, class TOutputMesh, class TInterpolator=itk::LinearInterpolateImageFunction<TInputImage> >
class CuberilleImageToMeshFilter : public ImageToMeshFilter< TInputImage, TOutputMesh >
{
public:
  /** Standard "Self" typedef. */
  typedef CuberilleImageToMeshFilter                    Self;
  typedef ImageToMeshFilter< TInputImage, TOutputMesh > Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CuberilleImageToMeshFilter, ImageToMeshFilter);

  /** Some convenient typedefs. */
  typedef TOutputMesh                                     OutputMeshType;
  typedef typename OutputMeshType::Pointer                OutputMeshPointer;
  typedef typename OutputMeshType::MeshTraits             OutputMeshTraits;
  typedef typename OutputMeshType::PointType              OutputPointType;
  typedef typename OutputMeshTraits::PixelType            OutputPixelType;
  typedef typename OutputMeshType::CellTraits             CellTraits;
  typedef typename OutputMeshType::PointsContainerPointer PointsContainerPointer;
  typedef typename OutputMeshType::PointsContainer        PointsContainer;
  typedef typename OutputMeshType::CellsContainerPointer  CellsContainerPointer;
  typedef typename OutputMeshType::CellsContainer         CellsContainer;
  typedef typename OutputMeshType::PointIdentifier        PointIdentifier;
  typedef typename OutputMeshType::CellIdentifier         CellIdentifier;
  typedef CellInterface<OutputPixelType, CellTraits>      CellInterfaceType;
  typedef TriangleCell<CellInterfaceType>                 TriangleCellType;
  typedef typename TriangleCellType::SelfAutoPointer      TriangleAutoPointer;
  typedef typename TriangleCellType::CellAutoPointer      TriangleCellAutoPointer;
  typedef QuadrilateralCell<CellInterfaceType>            QuadrilateralCellType;
  typedef typename QuadrilateralCellType::SelfAutoPointer QuadrilateralAutoPointer;
  typedef typename QuadrilateralCellType::CellAutoPointer QuadrilateralCellAutoPointer;

  typedef TInputImage                               InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename InputImageType::PixelType        InputPixelType;
  typedef typename InputImageType::SizeType         SizeType;
  typedef typename InputImageType::SpacingType      SpacingType;
  typedef typename InputImageType::SpacingValueType SpacingValueType;
  typedef typename InputImageType::IndexType        IndexType;
  typedef typename OutputMeshType::PointType        PointType;

  typedef TInterpolator                         InterpolatorType;
  typedef typename InterpolatorType::Pointer    InterpolatorPointer;
  typedef typename InterpolatorType::OutputType InterpolatorOutputType;

  /** Other convenient typedefs. */
  typedef ConstShapedNeighborhoodIterator< InputImageType >              InputImageIteratorType;
#if USE_GRADIENT_RECURSIVE_GAUSSIAN
  typedef GradientRecursiveGaussianImageFilter< InputImageType >         GradientFilterType;
#else
  typedef GradientImageFilter< InputImageType >                          GradientFilterType;
#endif
  typedef typename GradientFilterType::Pointer                           GradientFilterPointer;
  typedef typename GradientFilterType::OutputImageType                   GradientImageType;
  typedef typename GradientImageType::Pointer                            GradientImagePointer;
  typedef typename GradientFilterType::OutputPixelType                   GradientPixelType;
  typedef itk::VectorLinearInterpolateImageFunction< GradientImageType > GradientInterpolatorType;
  typedef typename GradientInterpolatorType::Pointer                     GradientInterpolatorPointer;

  /** Get/Set the iso-surface value.
    * This parameter specifies the value of the iso-surface for which to
    * generate the mesh. Pixels equal to or less than this value are
    * considered on the surface or inside the resultant mesh.
    */
  itkGetMacro( IsoSurfaceValue, InputPixelType );
  itkSetMacro( IsoSurfaceValue, InputPixelType );

  /** Accept the input image. */
  using Superclass::SetInput;
  virtual void SetInput( const InputImageType * inputImage );

  /** Get/Set interpolate function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get/Set whether triangle or quadrilateral faces should be generated.
    * True = triangle faces, False = quadrilateral faces.
      Default = true (triangle faces). */
  itkGetMacro( GenerateTriangleFaces, bool );
  itkSetMacro( GenerateTriangleFaces, bool );
  itkBooleanMacro( GenerateTriangleFaces );

  /** Get/Set whether the vertices should be project to the iso-surface.
      Default = true. */
  itkGetMacro( ProjectVerticesToIsoSurface, bool );
  itkSetMacro( ProjectVerticesToIsoSurface, bool );
  itkBooleanMacro( ProjectVerticesToIsoSurface );

  /** Get/Set the threshold for the "distance" from iso-surface during vertex projection.
      Note that the distance is actually measured in pixel value units (not space).
      The smaller this value, the closer the vertices will be to the iso-surface.
      Small values result in longer convergence time (i.e. slower).
      Values are clamped to the range [0.0, max pixel value].
      Default = 0.5. */
  itkGetMacro( ProjectVertexSurfaceDistanceThreshold, double );
  itkSetClampMacro( ProjectVertexSurfaceDistanceThreshold, double, 0.0, NumericTraits<InputPixelType>::max() );

  /** Get/Set the the initial step length for vertex projection.
      Values are clamped to the range [0.0, large].
      Default = max spacing * 0.25 (expressed in physical space). */
  itkGetMacro( ProjectVertexStepLength, double );
  itkSetClampMacro( ProjectVertexStepLength, double, 0.0, 100000.0 );

  /** Get/Set the step length relaxation factor during vertex projection.
      The step length is multiplied by this factor each iteration to allow convergence.
      Values are clamped to the range [0.0, 1.0].
      Default = 0.95. */
  itkGetMacro( ProjectVertexStepLengthRelaxationFactor, double );
  itkSetClampMacro( ProjectVertexStepLengthRelaxationFactor, double, 0.0, 1.0 );

  /** Get/Set the maximum number of steps used during vertex projection.
      Default = 50. */
  itkGetMacro( ProjectVertexMaximumNumberOfSteps, unsigned int );
  itkSetMacro( ProjectVertexMaximumNumberOfSteps, unsigned int );

protected:
  CuberilleImageToMeshFilter();
  ~CuberilleImageToMeshFilter() override;
  void PrintSelf(std::ostream& os, Indent indent) const override;

  void GenerateData() override;
  void GenerateOutputInformation() override{ }; // do nothing

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CuberilleImageToMeshFilter);

  /** \class VertexLookupNode A private class containing lookup details for vertices.
   *  \ingroup Cuberille */
  class VertexLookupNode
    {
    public:
      /** Convenient typedefs */
      typedef VertexLookupNode Self;

      /** Constructors */
      VertexLookupNode() : m_X( 0 ), m_Y( 0 ) { }
      VertexLookupNode(unsigned long x, unsigned long y) : m_X( x ), m_Y( y ) { }

      /** Parameters */
      unsigned long GetX() { return m_X; }
      unsigned long GetY() { return m_Y; }

      /** Comparison operators for sorting */
      bool operator> ( const Self& node ) const
        { return (m_Y > node.m_Y) || ( (m_Y == node.m_Y) && (m_X > node.m_X) ); }
      bool operator>= ( const Self& node ) const
        { return (m_Y >= node.m_Y) || ( (m_Y == node.m_Y) && (m_X >= node.m_X) ); }
      bool operator< ( const Self& node ) const
        { return (m_Y < node.m_Y) || ( (m_Y == node.m_Y) && (m_X < node.m_X) ); }
      bool operator<= ( const Self& node ) const
        { return (m_Y <= node.m_Y) || ( (m_Y == node.m_Y) && (m_X <= node.m_X) ); }

    private:
      unsigned long m_X;
      unsigned long m_Y;
  };

  /** \class VertexLookupMap A private class providing vertex lookup functionality.
   *  \ingroup Cuberille */
  template< class TMeshType >
  class VertexLookupMap
    {
    public:
      /** Convenient typedefs */
      typedef VertexLookupMap                               Self;
      typedef typename TMeshType::PointIdentifier           PointIdentifier;
      typedef std::map< VertexLookupNode, PointIdentifier > MapType;

      /** Constructors */
      VertexLookupMap() { }

      /** Clear the lookup map. */
      void Clear() { m_Map.clear(); }

      /** Add the given vertex identifer to the given [x,y] position. */
      void AddVertex( unsigned int x, unsigned int y, PointIdentifier id )
        {
          VertexLookupNode node( x, y );
          m_Map.insert( typename MapType::value_type(node, id) );
        }

      /** Get the vertex identifer for the given [x,y] position.
        * Returns true if the vertex exists and id contains the identifer.
        * Returns false if the vertex does not exist and id is undefined. */
      bool GetVertex( unsigned int x, unsigned int y, PointIdentifier &id )
        {
        bool result = false;
        VertexLookupNode node( x, y );
        typename MapType::iterator it = m_Map.find( node );
        if ( it != m_Map.end() )
          {
          result = true;
          id = it->second;
          }
        return result;
        }

    private:
      MapType m_Map;
  };

  /** Some convenient typedefs. */
  typedef VertexLookupMap< OutputMeshType > VertexLookupMapType;

  /** Private functions to implement the algorithm. */

  /** Compute gradient image. */
  inline void ComputeGradientImage();

  /** Set a flag activating each vertex for the given face. */
  inline void SetVerticesFromFace( unsigned int face, bool *vertexHasQuad );

  /** Get the vertex lookup index from the given index and vertex number. */
  inline IndexType GetVertexLookupIndex( unsigned int vertex, IndexType index );

  /** Project vertex to the iso-surface by stepping along normal. */
  inline void ProjectVertexToIsoSurface( PointType &vertex );

  /** Add a vertex to the given mesh. Increments point identifier. */
  inline void AddVertex( PointIdentifier &id, IndexType index, const InputImageType* image, OutputMeshType* mesh );

  /** Add quadrilateral face to the given mesh. Increments cell identifier. */
  inline void AddQuadFace( CellIdentifier &id, PointIdentifier f[4], OutputMeshType* mesh );

  InputPixelType              m_IsoSurfaceValue;
  InterpolatorPointer         m_Interpolator;
  GradientInterpolatorPointer m_GradientInterpolator;
  SpacingValueType            m_MaxSpacing;
  bool                        m_GenerateTriangleFaces;
  bool                        m_ProjectVerticesToIsoSurface;
  double                      m_ProjectVertexSurfaceDistanceThreshold;
  double                      m_ProjectVertexStepLength;
  double                      m_ProjectVertexStepLengthRelaxationFactor;
  unsigned int                m_ProjectVertexMaximumNumberOfSteps;
#if DEBUG_PRINT
  unsigned int m_ProjectVertexTerminate[3];
#endif
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCuberilleImageToMeshFilter.hxx"
#endif

#endif
