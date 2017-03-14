/*========================================================================= */

// CRIM : Centre de recherche en informatique de Montréal

// Équipe Télédetection pour les catastrophes majeures (TCM).

// Programme : Extraction des zones inondées en se  basant des images Radar.

// Auteur : Moslem Ouled Sghaier

// Version : 1

/*========================================================================= */

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "otbSFSTexturesImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include <itkGrayscaleMorphologicalOpeningImageFilter.h>
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>
#include "itkVector.h"
#include "itkBresenhamLine.h"
#include <vector>
#include <math.h>
#include <ctime>

#include "itkImageDuplicator.h"

#include "itkFlipImageFilter.h"
#include "itkFixedArray.h"

#include <cstdlib>

// Seuillage
#include "itkBinaryThresholdImageFilter.h"

//Bresenham
#include "itkBresenhamLine.h"

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkImageRegistrationMethod.h"

// Image  to LabelMap and LabelImage
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include <itkLabelMapToBinaryImageFilter.h>
#include "otbWrapperApplication.h" 
#include "otbWrapperApplicationRegistry.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperTags.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "otbPersistentVectorizationImageFilter.h"
#include "otbVectorDataProjectionFilter.h"

//Utils
#include "itksys/SystemTools.hxx"
#include "itkListSample.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"

// Image thining
#include "itkBinaryThinningImageFilter.h"


// ceci sera utile pour les données vectorielles
#include "otbVectorData.h" 
#include "otbVectorDataFileReader.h" 
#include "otbVectorDataFileWriter.h"

#include "itkPreOrderTreeIterator.h" 
#include "otbObjectList.h" 
#include "otbPolygon.h"
#include <otbVectorDataTransformFilter.h>

#include "itkIntensityWindowingImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"

//using namespace std;

extern "C" {
	#include "pde_toolbox_bimage.h"
	#include "pde_toolbox_defs.h"
//	#include "pde_toolbox_LSTB.h"
//	#include "ImageMagickIO.h"
}

#include "pathopenclose.h"

#include "itkChangeInformationImageFilter.h"
#include "itkVersor.h"

namespace otb
{
namespace Wrapper
{

class SarInondationsExtraction : public Application
{

public:
	typedef SarInondationsExtraction Self;
    typedef Application              Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

#define PI 3.14159265


/** Standard macro */
    itkNewMacro(Self);
    itkTypeMacro(SarInondationsExtraction, otb::Application);

typedef unsigned char CharPixelType; // IO
typedef otb::Image<CharPixelType, 2> CharImageType;


private:

void DoInit()
{

SetName("SarInondationsExtraction"); // Nécessaire
SetDocName("SarInondationsExtraction");
SetDocLongDescription("Un module pour l'extraction des zones inondées en se basant sur des images radar");
SetDocLimitations("Les autres paramètres seront ajoutés plus tard");
SetDocAuthors("Moslem Ouled Sghaier");


AddParameter(ParameterType_InputImage,"in1", "Input Image 1");
SetParameterDescription("in1", "The input image 1");
AddParameter(ParameterType_InputImage,"in2", "Input Image 2");
SetParameterDescription("in2", "The input image 2");
AddParameter(ParameterType_OutputVectorData,"out1", "Output Vector 1");
SetParameterDescription("out1","The output vector 1");
AddParameter(ParameterType_OutputVectorData,"out2", "Output Vector 2");
SetParameterDescription("out2","The output vector 2");
AddParameter(ParameterType_OutputVectorData,"out3", "Output Vector 3");
SetParameterDescription("out3","The output vector 3");

}

void DoUpdateParameters()
{
	// Nothing to do here : all parameters are independent
}

void DoExecute()
{
  
clock_t start, stop;
start = clock();

   /////////////////////////////////////////////////////////////////////// déclaration des structures globales

typedef unsigned char CharPixelType; // IO
const unsigned int Dimension = 2;
typedef otb::Image<CharPixelType, Dimension> CharImageType;

// Conversion des deux images

CharImageType::Pointer rescale1  = Convert (GetParameterUInt32Image("in1"));
CharImageType::Pointer rescale2  = Convert (GetParameterUInt32Image("in2"));

// Application du SFS-SD sur les deux images

CharImageType::Pointer SFSFilter1 = SFSTextures(rescale1);
CharImageType::Pointer SFSFilter2 = SFSTextures(rescale2);

typedef otb::ImageFileWriter<CharImageType> WriterType;
WriterType::Pointer writer20 = WriterType::New();
WriterType::Pointer writer10 = WriterType::New();
writer20->SetFileName("C:/Users/ouledsmo/Desktop/richard10.tif");
writer10->SetFileName("C:/Users/ouledsmo/Desktop/richard20.tif");
writer20->SetInput(SFSFilter1);
writer10->SetInput(SFSFilter2);
writer20->Update();
writer10->Update();

// Filtrer les bords des deux images

SFSFilter1 = border  (SFSFilter1);
SFSFilter2 = border  (SFSFilter2);

// Application de l'ouverture morhologique en utilisant un élément structurant sous la forme d'un chemin, images 1 et 2

CharImageType::Pointer openingFilter1 = PathOpening(SFSFilter1,100);
CharImageType::Pointer openingFilter2 = PathOpening(SFSFilter2,50);

// Seuillage pour afficher les zones homogènes

CharImageType::Pointer seuil1 = Seuillage(openingFilter1);
CharImageType::Pointer seuil2 = Seuillage(openingFilter2);

// Application de l'ouverture morphologique (closing) pour les deux images

CharImageType::Pointer closingFilter1 = Closing(seuil1);
CharImageType::Pointer closingFilter2 = Closing(seuil2);

typedef otb::ImageFileWriter<CharImageType> WriterType;
WriterType::Pointer writer = WriterType::New();
WriterType::Pointer writer1 = WriterType::New();
writer->SetFileName("C:/Users/ouledsmo/Desktop/richard.tif");
writer1->SetFileName("C:/Users/ouledsmo/Desktop/richard1.tif");
writer->SetInput(closingFilter1);
writer1->SetInput(closingFilter2);
writer->Update();
writer1->Update();

// extraire l'image de différence 

CharImageType::Pointer difference = diff(closingFilter1, closingFilter2);

stop = clock();
std::cout << "CPU time elapsed:" << ((double)stop-start)/CLOCKS_PER_SEC << std::endl;

//Label the objects in a binary image

VectorDataType::Pointer vectorDataProjection =  Affichage(closingFilter1,1);
SetParameterOutputVectorData("out1",vectorDataProjection); 

VectorDataType::Pointer vectorDataProjection1 =  Affichage(closingFilter2,2);
SetParameterOutputVectorData("out2",vectorDataProjection1);

VectorDataType::Pointer vectorDataProjection2 =  Affichage(difference,2);
SetParameterOutputVectorData("out3",vectorDataProjection2); 
}

  // Les differents fonctions utilisées dans le code

  CharImageType::Pointer Convert (otb::Wrapper::UInt32ImageType* input)
  {
     typedef itk::RescaleIntensityImageFilter<UInt32ImageType,CharImageType> RescaleFilter0;
     RescaleFilter0::Pointer rescale1 = RescaleFilter0::New();
     rescale1->SetInput(input);
     rescale1->UpdateLargestPossibleRegion();
	 
	 return (rescale1->GetOutput());
  }

 CharImageType::Pointer SFSTextures (CharImageType::Pointer image)
  {
    typedef otb::SFSTexturesImageFilter<CharImageType, CharImageType> SFSFilterType;
    SFSFilterType::Pointer SFSFilter1 = SFSFilterType::New();
    SFSFilter1->SetSpectralThreshold(8);
    SFSFilter1->SetSpatialThreshold(100); 
    SFSFilter1->SetNumberOfDirections(20); 
    SFSFilter1->SetRatioMaxConsiderationNumber(5); 
    SFSFilter1->SetAlpha(1.00);
    SFSFilter1->SetInput(image);
    SFSFilter1->Update();

	return (SFSFilter1->GetOutput());
  } 

 CharImageType::Pointer border (CharImageType::Pointer SFSFilter1)

  {
  unsigned int x= SFSFilter1->GetLargestPossibleRegion().GetSize()[0];
  unsigned int y= SFSFilter1->GetLargestPossibleRegion().GetSize()[1];

  for (unsigned int i=0 ; i< 2;i++ )
  for (unsigned int j=0 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  SFSFilter1->SetPixel(pixelIndex,0);}}

  for (unsigned int i=0 ; i< x;i++ )
  for (unsigned int j=0 ; j< 2;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  SFSFilter1->SetPixel(pixelIndex,0);}}

  for (unsigned int i=x-2 ; i< x;i++ )
  for (unsigned int j=0 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  SFSFilter1->SetPixel(pixelIndex,0);}}

  for (unsigned int i=0 ; i< x;i++ )
  for (unsigned int j=y-2 ; j< y;j++)
  {{
  CharImageType::IndexType pixelIndex; 
  pixelIndex[0] = i;   // x position 
  pixelIndex[1] = j;   // y position
  SFSFilter1->SetPixel(pixelIndex,0);}}

  return (SFSFilter1);
  }

 CharImageType::Pointer opening (CharImageType::Pointer input)
  
 {
  typedef itk::FlatStructuringElement<2>  StructuringElementType;
  typedef itk::GrayscaleMorphologicalOpeningImageFilter<CharImageType,CharImageType,StructuringElementType>  GrayscaleDilateImageFilterType;
  GrayscaleDilateImageFilterType::Pointer openingFilter1 = GrayscaleDilateImageFilterType::New();

  openingFilter1->SetInput(input);
  StructuringElementType::RadiusType elementRadius;
  elementRadius.Fill(1);
  StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);
  openingFilter1->SetKernel(structuringElement);
  openingFilter1->GetSafeBorder();
  openingFilter1->Update();
 
  return (openingFilter1->GetOutput());
 }

 CharImageType::Pointer diff (CharImageType::Pointer avant, CharImageType::Pointer apres)
 {
  CharImageType::SizeType size, size1, size2;

  CharImageType::Pointer resultat = CharImageType::New();

  size=avant->GetLargestPossibleRegion().GetSize();
  size1=apres->GetLargestPossibleRegion().GetSize();

  size2[0]=std::min(size[0],size1[0]);
  size2[1]=std::min(size[1],size1[1]);
  CharImageType::RegionType region; 
  CharImageType::IndexType start; 
  start[0] =   0;  // first index on X 
  start[1] =   0;  // first index on Y

  region.SetIndex(start);
  region.SetSize(size2);

  resultat->SetRegions(region);
  resultat->Allocate();
  resultat->FillBuffer(0);
  resultat->Update();

  for (int i = 0; i < size2[0]; i++)
  {
    for (int j = 0 ; j < size2[1]; j++ )
	{
       CharImageType::IndexType pixelIndex; 
	   pixelIndex[0] = i;   // x position 
       pixelIndex[1] = j;   // y position
	   CharImageType::PixelType pixelValue = avant->GetPixel(pixelIndex);
	   CharImageType::PixelType pixelValue1 = apres->GetPixel(pixelIndex);
	  if ( (pixelValue == 0) && (pixelValue1 == 255) )
	  resultat->SetPixel(pixelIndex,255);
	}
  
  }
 resultat->Update();

 typedef otb::ImageFileWriter<CharImageType> WriterType;
WriterType::Pointer writer123 = WriterType::New();
writer123->SetFileName("C:/Users/ouledsmo/Desktop/richard3.tif");
writer123->SetInput(resultat);
writer123->Update();

 return resultat;
 }

 CharImageType::Pointer PathOpening (CharImageType::Pointer openingFilter, int L)
 {
    int   i, K;

    K=1;

  unsigned int nx = openingFilter->GetLargestPossibleRegion().GetSize()[0];//input_bimage->dim->buf[0];
  unsigned int ny = openingFilter->GetLargestPossibleRegion().GetSize()[1];;//input_bimage->dim->buf[1];
  unsigned int num_pixels = nx*ny;

  std::cout << num_pixels << std::endl;

  PATHOPEN_PIX_TYPE * input_image = new PATHOPEN_PIX_TYPE[nx * ny];
  PATHOPEN_PIX_TYPE * output_image = new PATHOPEN_PIX_TYPE[nx * ny];
  
  // Convert intermediate float to PATHOPEN_PIX_TYPE (unsigned char)
	for (i = 0; i < num_pixels; ++i) {
		CharImageType::IndexType index = {i % nx, i / nx};
		input_image[i] = openingFilter->GetPixel(index);//static_cast<PATHOPEN_PIX_TYPE>(input_bimage->buf[i]);
	}

	std::cout << "Calling pathopen ()" << std::endl;
        
	pathopen(
            input_image, // The input image //
            nx, ny,	 // Image dimensions //
            L,		 // The threshold line length //
            K,		 // The maximum number of gaps in the path //
            output_image // Output image //
            ); 

	for (i = 0; i < num_pixels; ++i) {
		CharImageType::IndexType index = {i % nx, i / nx};
		openingFilter->SetPixel(index,  static_cast<CharImageType::PixelType>(output_image[i]));
	}

	return (openingFilter);
 }

 CharImageType::Pointer Seuillage (CharImageType::Pointer openingFilter)
 {
  typedef itk::BinaryThresholdImageFilter<CharImageType, CharImageType>  FilterType1;
  FilterType1::Pointer filter1 = FilterType1::New();
  filter1->SetInput(openingFilter);
  filter1->SetLowerThreshold(14);
  filter1->Update();

  return (filter1->GetOutput());
 }

 CharImageType::Pointer Closing (CharImageType::Pointer input)
  
 {
  typedef itk::FlatStructuringElement<2>  StructuringElementType;
  typedef itk::BinaryMorphologicalClosingImageFilter<CharImageType,CharImageType,StructuringElementType>  ClosingImageFilterType;
  ClosingImageFilterType::Pointer closingFilter = ClosingImageFilterType::New();

  closingFilter->SetInput(input);
  StructuringElementType::RadiusType elementRadius;
  elementRadius.Fill(2);
  StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);
  closingFilter->SetKernel(structuringElement);
  closingFilter->GetSafeBorder();
  closingFilter->Update();
  return (closingFilter->GetOutput());
 }

 VectorDataType::Pointer Affichage(CharImageType::Pointer seuil, int proj)
  {
  typedef otb::Polygon<double>             PolygonType;
  typedef PolygonType::Pointer             PolygonPointerType;
  typedef PolygonType::ContinuousIndexType PolygonIndexType;
  typedef otb::ObjectList<PolygonType>     PolygonListType;
  typedef PolygonListType::Pointer         PolygonListPointerType;
  typedef unsigned long LabelPixelType;
  typedef otb::Image<LabelPixelType, 2>  LabeledImageType;

  typedef itk::ConnectedComponentImageFilter<CharImageType,LabeledImageType> ConnectedFilterType;
  typedef otb::PersistentVectorizationImageFilter<LabeledImageType,PolygonType> PersistentVectorizationFilterType;

  ConnectedFilterType::Pointer connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput(seuil);

  //Perform vectorization in a persistent way
    PersistentVectorizationFilterType::Pointer persistentVectorization = PersistentVectorizationFilterType::New();
    persistentVectorization->Reset();
    persistentVectorization->SetInput(connectedFilter->GetOutput());
    try
      {
      persistentVectorization->Update();
      }
    catch (itk::ExceptionObject& err)
      {
      std::cout << "\nExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      }

	PolygonListPointerType OutputPolyList = persistentVectorization->GetPathList();
	//Display results
    std::cout << "nb objects found = " << OutputPolyList->Size() << std::endl;

	VectorDataType::Pointer outVectorData = VectorDataType::New(); 
    typedef VectorDataType::DataNodeType DataNodeType;
	typedef VectorDataType::DataTreeType            DataTreeType; 
    typedef itk::PreOrderTreeIterator<DataTreeType> TreeIteratorType;

	DataNodeType::Pointer document = DataNodeType::New(); 
    document->SetNodeType(otb::DOCUMENT); 
    document->SetNodeId("polygon"); 
    DataNodeType::Pointer folder = DataNodeType::New(); 
    folder->SetNodeType(otb::FOLDER); 
    DataNodeType::Pointer multiPolygon = DataNodeType::New(); 
    multiPolygon->SetNodeType(otb::FEATURE_MULTIPOLYGON);

	DataTreeType::Pointer tree = outVectorData->GetDataTree(); 
    DataNodeType::Pointer root = tree->GetRoot()->Get(); 
 
    tree->Add(document, root); 
    tree->Add(folder, document); 
    tree->Add(multiPolygon, folder);

	for (PolygonListType::Iterator pit = OutputPolyList->Begin(); 
       pit != OutputPolyList->End(); ++pit) 
    { 
    DataNodeType::Pointer newPolygon = DataNodeType::New(); 
    newPolygon->SetPolygonExteriorRing(pit.Get()); 
    tree->Add(newPolygon, multiPolygon); 
    }

	typedef itk::AffineTransform<VectorDataType::PrecisionType, 2> TransformType;
	typedef otb::VectorDataTransformFilter <VectorDataType,VectorDataType > VectorDataFilterType;
    VectorDataFilterType:: Pointer vectorDataProjection = VectorDataFilterType:: New ();

	char* projection;
	if (proj == 1)
	{projection = "in1";}
	if (proj == 2)
	{projection = "in2";}

	TransformType::ParametersType params;
    params.SetSize(6);
    params[0] = GetParameterUInt32Image(projection)->GetSpacing()[0];
    params[1] = 0;
    params[2] = 0;
    params[3] = GetParameterUInt32Image(projection)->GetSpacing()[1];
    params[4] = GetParameterUInt32Image(projection)->GetOrigin()[0];
    params[5] = GetParameterUInt32Image(projection)->GetOrigin()[1];
	TransformType::Pointer transform = TransformType::New();
    transform->SetParameters(params);

	vectorDataProjection->SetTransform(transform);
	vectorDataProjection->SetInput(outVectorData);
	vectorDataProjection->Update();

	std::cout << "OK" << std::endl;
	std::cout << "L'origine1:" << vectorDataProjection->GetOutput()->GetOrigin() << GetParameterUInt32Image(projection)->GetOrigin() << std::endl;
    std::cout << "L'espacement1:" << vectorDataProjection->GetOutput()->GetSpacing() << GetParameterUInt32Image(projection)->GetSpacing() << std::endl << std::endl;
	std::cout << "La projection1:"  << vectorDataProjection->GetOutput()->GetProjectionRef() << GetParameterUInt32Image(projection)->GetProjectionRef() << std::endl << std::endl;
	
	return (vectorDataProjection->GetOutput());
  }

     };
}
   
}
OTB_APPLICATION_EXPORT(otb::Wrapper::SarInondationsExtraction);