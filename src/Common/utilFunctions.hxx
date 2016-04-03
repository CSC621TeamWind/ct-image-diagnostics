#ifndef utilFunctions_hxx
#define utilFunctions_hxx

#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

namespace utility {
  typedef signed short PixelType;
  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ImageSeriesReader< ImageType > ReaderType;

  /* 
  Reads a DICOM image series given the directory path containing
  the image series.
  */
  ReaderType::Pointer readDicomImageSeries(std::string dirname) {
    ReaderType::Pointer reader = ReaderType::New();

    typedef itk::GDCMImageIO       ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    
    reader->SetImageIO( dicomIO );
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails( true );
    nameGenerator->AddSeriesRestriction("0008|0021" );
    nameGenerator->SetDirectory(dirname);
    
    try {
      typedef std::vector< std::string >    SeriesIdContainer;
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

      while( seriesItr != seriesEnd ) {
        ++seriesItr;
      }
      std::string seriesIdentifier;
      seriesIdentifier = seriesUID.begin()->c_str();

      typedef std::vector< std::string >   FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames( seriesIdentifier );

      reader->SetFileNames( fileNames );

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
}

#endif