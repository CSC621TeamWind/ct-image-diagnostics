/**
 * URN-NBN-SI-DOC-HAVH299II
 */

#include <string>

// FIXME: use quotes always?
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "utilFunctions.hxx" 
#include "itkOtsuThresholdCalculator.h"
#include "SegmentedLungFilter.hxx"

#include "util.hxx"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

typedef itk::PointSet< PixelType, 3 > PointSetType;

using itk::Indent;
using itk::Object;
using itk::BinaryThresholdImageFilter;
using itk::SliceBySliceImageFilter;
using itk::BinaryFunctorImageFilter;
using itk::ExceptionObject;

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

    // read the DICOM image series and print out the number of slices.
    try {
        ReaderType::Pointer reader = ReaderType::New();
        reader = utility::readDicomImageSeries(dicom_path);

	typedef signed short DICOMPixelType;
	typename itk::Image<DICOMPixelType, 3>::Pointer image = reader->GetOutput();

	typedef SegmentedLungFilter<itk::Image<DICOMPixelType, 3> > LungFilter;
	typename LungFilter::Pointer lungs = LungFilter::New();
	typedef typename LungFilter::OutputImageType SegmentedLungImage;
	lungs->SetInput(image);
	lungs->Update();

	typedef itk::Image<unsigned char, 2> SegmentedLungSlice;
	typename SegmentedLungSlice::Pointer image2d = extract2DImageSlice<SegmentedLungImage, SegmentedLungSlice >(lungs->GetOutput(), 2, 63);

	typedef itk::RescaleIntensityImageFilter<SegmentedLungSlice, ImageType2D> RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image2d);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	QuickView viewer;
	//viewer.AddImage(image2d);
	viewer.AddImage(rescaleFilter->GetOutput());
	viewer.Visualize();

    }
    catch (ExceptionObject &ex) {
        cerr << "The program encountered an exception: " << ex << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
