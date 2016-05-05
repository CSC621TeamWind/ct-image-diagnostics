#ifndef utilFunctions_hxx
#define utilFunctions_hxx

#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "QuickView.h"
#include "itkPointSetToImageFilter.h"

typedef float PixelType;
typedef itk::Image< PixelType, 3 > ImageType3D;
typedef itk::ImageSeriesReader< ImageType3D > ReaderType;
typedef itk::Image< PixelType, 2 > ImageType2D;
typedef itk::ExtractImageFilter < ImageType3D, ImageType2D > FilterType2D;
typedef itk::PointSet< PixelType, 3 > PointSetType;
typedef itk::Image< PixelType, 2 >  BinaryImageType;

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

  /*
  BinaryImageType::Pointer getImageFromPoints(PointSetType::Pointer pointSet) {
	  typedef itk::PointSetToImageFilter<PointSetType, BinaryImageType> PointsToImageFilterType;
	  PointsToImageFilterType::Pointer pointsToImageFilter = PointsToImageFilterType::New();
	  pointsToImageFilter->SetInput(pointSet);
	  BinaryImageType::SpacingType spacing;
	  spacing.Fill(1.0);
	  BinaryImageType::PointType origin;
	  origin.Fill(0.0);
	  pointsToImageFilter->SetSpacing(spacing);
	  pointsToImageFilter->SetOrigin(origin);
	  pointsToImageFilter->Update();
	  BinaryImageType::Pointer binaryImage = pointsToImageFilter->GetOutput();
	  return binaryImage;
  }
  */

  void CreateImage(ImageType2D::Pointer image, PointSetType::Pointer pointSet)
  {
	  // Create a black image with a white square
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