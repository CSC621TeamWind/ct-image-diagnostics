#ifndef utilFunctionsAlt_hxx
#define utilFunctionsAlt_hxx

#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "QuickView.h"
#include "itkPointSetToImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkPNGImageIO.h" // TODO: Check if this is needed.
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"


typedef float PixelType;
typedef itk::Image< PixelType, 3 > ImageType3D;
typedef itk::ImageSeriesReader< ImageType3D > ReaderType;
typedef itk::Image< PixelType, 2 > ImageType2D;
typedef itk::ExtractImageFilter < ImageType3D, ImageType2D > FilterType2D;
typedef itk::PointSet< PixelType, 3 > PointSetType;
typedef itk::Image< unsigned char, 2 >  BinaryImageType;

namespace utility {

  /*
  Reads a DICOM image series given the directory path containing
  the image series.
  */
  ReaderType::Pointer readDicomImageSeries(std::string dirname) {
    ReaderType::Pointer reader = ReaderType::New();

    typedef itk::GDCMImageIO       ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();

    reader->SetImageIO(dicomIO);
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetDirectory(dirname);

    try {
      typedef std::vector< std::string >    SeriesIdContainer;
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

      while (seriesItr != seriesEnd) {
        ++seriesItr;
      }
      std::string seriesIdentifier;
      seriesIdentifier = seriesUID.begin()->c_str();

      typedef std::vector< std::string >   FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      reader->SetFileNames(fileNames);

      try {
        reader->Update();
        return reader;
      }
      catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
      }
    }
    catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
    }
  }

  /*
  Display the given  image. No return type.
  */
  void display2DImage(ImageType2D::Pointer image) {
    typedef itk::RescaleIntensityImageFilter<ImageType2D, ImageType2D> RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(image);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->Update();

    QuickView viewer;
    viewer.AddImage(image.GetPointer());
    viewer.AddImage(rescaleFilter->GetOutput());
    viewer.Visualize();
  }

  /*
  Extracts a 2D image slice from a 3D image.
  */
  ImageType2D::Pointer extract2DImageSlice(ImageType3D::Pointer reader, int plane, int slice) {
    FilterType2D::Pointer filter = FilterType2D::New();

    ImageType3D::RegionType inputRegion =
      reader->GetLargestPossibleRegion();

    ImageType3D::SizeType size = inputRegion.GetSize();
    size[plane] = 0;  // collapsing the 3rd dimension

    ImageType3D::IndexType start = inputRegion.GetIndex();
    const unsigned int sliceNumber = slice;
    start[plane] = sliceNumber;  // the required index of the 3rd dimension

    ImageType3D::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);

    filter->SetExtractionRegion(desiredRegion);

    // setting the direction of collapse
    filter->SetDirectionCollapseToSubmatrix();
    filter->SetInput(reader);

    ImageType2D::Pointer img = filter->GetOutput();
    img->Update();

    return img;
  }

  void WriteSliceAsPNG(ImageType3D::Pointer myimage, int sliceNumber, std::string filename) {
	  ImageType2D::Pointer testimage = extract2DImageSlice(myimage, 2, sliceNumber);
	  std::cout << "Displaying slice number : " << sliceNumber << std::endl;
	  std::cout << "Displaying the filename : " << filename << std::endl;
	  display2DImage(testimage);
	  const float spacing[2] = { 1.0, 1.0 };
	  testimage->SetSpacing(spacing);
	  // Since png only supports unsigned char and unsigned float, 
	  // casting the image to unsigned char.
	  typedef unsigned char                             OutputPixelType;
	  typedef itk::Image< OutputPixelType, 2 >  OutputImageType;

	  typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageType2D >
		  RescaleType;
	  RescaleType::Pointer rescale = RescaleType::New();
	  rescale->SetInput(testimage);
	  rescale->SetOutputMinimum(0);
	  rescale->SetOutputMaximum(itk::NumericTraits< OutputPixelType >::max());

	  typedef itk::CastImageFilter< ImageType2D, OutputImageType > FilterType;
	  FilterType::Pointer filter = FilterType::New();
	  filter->SetInput(rescale->GetOutput());

	  typedef  itk::ImageFileWriter< OutputImageType > WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetImageIO(itk::PNGImageIO::New());
	  writer->SetFileName(filename);
	  writer->SetInput(filter->GetOutput());
	  writer->Update();
  }
  

  /* Creates a binary png image of size width x height. Draws the points of the pointset of slice 'slice'
  on the image.
  */
  void WriteBinaryLabelImage(PointSetType::Pointer pointSet, int width, int height, float sliceNumber, std::string filename) {
	  typedef unsigned short PixelType2; // Since PNG only uses unsigned char and unsigned short
	  typedef itk::Image< PixelType2, 2 >  ImageTypeUnsigned;
	  ImageTypeUnsigned::Pointer binaryImage = ImageTypeUnsigned::New();

	  ImageTypeUnsigned::IndexType start;
	  start.Fill(0);

	  ImageTypeUnsigned::SizeType size;
	  size[0] = width;
	  size[1] = height;

	  ImageTypeUnsigned::RegionType region;
	  region.SetSize(size);
	  region.SetIndex(start);
	  binaryImage->SetRegions(region);
	  binaryImage->Allocate();
	  binaryImage->FillBuffer(0);

	  /*typedef itk::RescaleIntensityImageFilter< ImageTypeUnsigned, ImageTypeUnsigned > RescaleType;
	  RescaleType::Pointer rescale = RescaleType::New();
	  rescale->SetInput(binaryImage);
	  rescale->SetOutputMinimum(0);
	  rescale->SetOutputMaximum(255);
	  rescale->SetOutputMaximum(itk::NumericTraits< PixelType2 >::max());*/

	  typedef PointSetType::PointsContainer      PointsContainer;
	  PointsContainer::Pointer  points2 = pointSet->GetPoints();
	  typedef PointsContainer::Iterator     PointsIterator;
	  PointsIterator  pointIterator = points2->Begin();
	  PointsIterator end = points2->End();
	  while (pointIterator != end)
	  {
		  itk::Point<PixelType, 3> p = pointIterator.Value();
		  if (sliceNumber == p[2]) {
			  ImageTypeUnsigned::IndexType pixelIndexSpherePoint;
			  pixelIndexSpherePoint[0] = p[0];
			  pixelIndexSpherePoint[1] = p[1];
			  //std::cout << "The index of the point is " << pointIterator.Index() << std::endl;
			  binaryImage->SetPixel(pixelIndexSpherePoint, 65530);  // TODO: Rescale intensity range.
			  //std::cout << "The required points are " << pixelIndexSpherePoint << std::endl;
		  }
		  ++pointIterator;
	  }

	  // Writing out the image.
	  typedef  itk::ImageFileWriter< ImageTypeUnsigned > WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetImageIO(itk::PNGImageIO::New());
	  writer->SetFileName(filename);
	  writer->SetInput(binaryImage);
	  writer->Update();

  }

  /* Create a test black image with two white squares
  */
  void CreateImage(ImageType2D::Pointer image )
  {
	  ImageType2D::IndexType start;
	  start.Fill(0);

	  ImageType2D::SizeType size;
	  size.Fill(100);

	  ImageType2D::RegionType region;
	  region.SetSize(size);
	  region.SetIndex(start);
	  image->SetRegions(region);
	  image->Allocate();
	  image->FillBuffer(0);

	  itk::ImageRegionIterator<ImageType2D> imageIterator(image, image->GetLargestPossibleRegion());

	  // Make two squares
	  while (!imageIterator.IsAtEnd())
	  {
		  if ((imageIterator.GetIndex()[0] > 5 && imageIterator.GetIndex()[0] < 20) &&
			  (imageIterator.GetIndex()[1] > 5 && imageIterator.GetIndex()[1] < 20))
		  {
			  imageIterator.Set(255);
		  }

		  if ((imageIterator.GetIndex()[0] > 50 && imageIterator.GetIndex()[0] < 70) &&
			  (imageIterator.GetIndex()[1] > 50 && imageIterator.GetIndex()[1] < 70))
		  {
			  imageIterator.Set(255);
		  }
		  ++imageIterator;
	  }
  }
}

#endif