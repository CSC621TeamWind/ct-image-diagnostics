/**
 * URN-NBN-SI-DOC-HAVH299II
 */

#include <string>

// FIXME: use quotes always?
#include "itkImage.h"
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "utilFunctions.hxx" 
#include "itkOtsuThresholdCalculator.h"
#include "SegmentedLungFilter.hxx"
#include "itkImageToImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"

//#include "MultipleThresholding.hxx"

#include "util.hxx"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

typedef itk::PointSet< PixelType, 3 > PointSetType;

using itk::Indent;
using itk::Object;
using itk::BinaryThresholdImageFilter;
using itk::SliceBySliceImageFilter;
using itk::BinaryFunctorImageFilter;
using itk::ExceptionObject;

//typedef typename TImage::ConstPointer ConstInputImagePointer;

int usage(char * prog) {
    cerr << "Usage: " << prog << " DICOM_DIR" << endl;
    return -1;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        return usage(argv[0]);
    }

    bool verbose = false;

    for (int i = 1; i < argc - 1; i++) {
        if (string("-v") == argv[i]) {
            verbose = true;
        } else if (string("--verbose") == argv[i]) {
            verbose = true;
        } else {
            return usage(argv[0]);
        }
    }

    string dicom_path = argv[argc - 1];

    int MultiThreshold[10] = {-600, -550, -500, -450, -400, -350, -300, -250, -200, -150};
    

    int max_iterations =  10;

    unsigned int Dimension = 2;
    unsigned int radius = 1;

    typedef unsigned char InputPixelType;
    typedef unsigned char OutputPixelType;

    typedef itk::Image< InputPixelType,  3 >   InputImageType;
    typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

    typename itk::Image< InputImageType, 3 >::Pointer TempImageType;

    typedef itk::BinaryThresholdImageFilter< InputImageType, InputImageType > ThresholdFilterType;
    ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
    
    typedef std::vector<ThresholdFilterType::Pointer> ThresholdVector;
    ThresholdVector  Tvector;
 
    //For Erosion
    typedef itk::FlatStructuringElement<2> StructuringElementType;
    StructuringElementType::RadiusType elementRadius;
    elementRadius.Fill(radius);
 
    StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);
 
    typedef itk::BinaryErodeImageFilter <InputImageType, InputImageType, StructuringElementType>
    BinaryErodeImageFilterType;
     
     for (int k = 0; k < 10; k++)
        Tvector[k] = ThresholdFilterType::New();
    
    const PixelType OutsideValue = static_cast<PixelType>( 0 );
    const PixelType InsideValue = static_cast<PixelType>( 255 );

    // read the DICOM image series and print out the number of slices.
    try 
    {
        ReaderType::Pointer reader = ReaderType::New();
        reader = utility::readDicomImageSeries(dicom_path);

	    typedef signed short DICOMPixelType;
	    typename itk::Image<DICOMPixelType, 3>::Pointer image = reader->GetOutput();

	    typedef SegmentedLungFilter<itk::Image<DICOMPixelType, 3> > LungFilter;
	    typename LungFilter::Pointer lungs = LungFilter::New();
	    typedef typename LungFilter::OutputImageType SegmentedLungImage;
	    lungs->SetInput(image);
	    lungs->Update();
  
        /* Attemp erosion

        BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
        erodeFilter->SetInput(lungs->Update);
        erodeFilter->SetKernel(structuringElement);
        erodeFilter->SetErodeValue(255);*/

        thresholder->SetInput( lungs->GetOutput() );
  
    for (int i = 0; i < max_iterations; i++)
    {
      const PixelType LowerThreshold = static_cast<PixelType>(MultiThreshold[i]);
      
      thresholder->SetLowerThreshold( LowerThreshold );
   
      thresholder->SetOutsideValue( OutsideValue );
      thresholder->SetInsideValue( InsideValue );
      
      thresholder->Update();
      Tvector[i] = thresholder;
    }
    
	//displaySlice<SegmentedLungImage>(lungs->GetOutput(), 2, 63);

    for (int j = 0; j < 10; j++)
        displaySlice<SegmentedLungImage>(Tvector[j]->GetOutput(), 2, 63);

    /*
    /Image Iterator to extract data from will go here 
    */
    
    //InputImageType::IndexType pixelIndex = {{0,0,0}}; //Initial Position of {x,y,z}
    //InputImageType::PixelType pixelValue;

    
   /* for (int l = 0; l < 10; l++)
    {
       pixelValue = Tvector[l]->GetPixel(pixelIndex);

    }*/

    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writerThreshold = WriterType::New();
    writerThreshold->SetFileName("~/Desktop/test.jpeg");
    writerThreshold->SetInput(lungs->GetOutput());
    writerThreshold->Update();
    
    }
    catch (ExceptionObject &ex) {
        cerr << "The program encountered an exception: " << ex << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
