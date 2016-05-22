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
#include "itkMaskImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkConnectedComponentImageFilter.h"

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

    typedef signed int InputPixelType;
    typedef unsigned char OutputPixelType;

    typedef signed short DICOMPixelType;
    typedef itk::Image<DICOMPixelType, 3> DICOMImageType;
    typedef DICOMImageType   InputImageType;
    typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

    typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > ThresholdFilterType;

    typedef std::vector<ThresholdFilterType::Pointer> ThresholdVector;
    ThresholdVector  Tvector;

    //For Erosion
    typedef itk::FlatStructuringElement<2> StructuringElementType;
    StructuringElementType::RadiusType elementRadius;
    elementRadius.Fill(radius);

    StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);

    typedef itk::BinaryErodeImageFilter <InputImageType, InputImageType, StructuringElementType>
        BinaryErodeImageFilterType;

    const PixelType OutsideValue = static_cast<PixelType>( 0 );
    const PixelType InsideValue = static_cast<PixelType>( 255 );

    // read the DICOM image series and print out the number of slices.
    try 
    {
        ReaderType::Pointer reader = ReaderType::New();
        reader = utility::readDicomImageSeries(dicom_path);

        DICOMImageType::Pointer image = reader->GetOutput();

	typename DICOMImageType::SpacingType spacing = image->GetSpacing();

	// Calculate minimum nodule size in integer coordinates
	int minNoduleSize[3]; 
	for (int i = 0; i < 3; i++) {
		minNoduleSize[i] = 4 /* mm */ / spacing[i];
	}

	int maxNoduleSize[3]; 
	for (int i = 0; i < 3; i++) {
		maxNoduleSize[i] = 30 /* mm */ / spacing[i];
	}

        typedef SegmentedLungFilter<DICOMImageType> LungFilter;
        typename LungFilter::Pointer lungs = LungFilter::New();
        typedef typename LungFilter::OutputImageType SegmentedLungImage;
        lungs->SetInput(image);
        lungs->Update();

        /* Attemp erosion

           BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
           erodeFilter->SetInput(lungs->Update);
           erodeFilter->SetKernel(structuringElement);
           erodeFilter->SetErodeValue(255);*/

        typedef itk::MaskImageFilter<DICOMImageType, typename LungFilter::OutputImageType, DICOMImageType > MaskedLungImage;
        typedef typename MaskedLungImage::Pointer MaskedLungImagePointer;

        MaskedLungImagePointer maskedImage = MaskedLungImage::New();
        maskedImage->SetInput(image);
        maskedImage->SetOutsideValue(itk::NumericTraits<DICOMPixelType>::min());
        maskedImage->SetMaskImage(lungs->GetOutput());
        maskedImage->Update();

        typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructuringElementType;
        StructuringElementType openingElement;
        openingElement.SetRadius(1);
        openingElement.CreateStructuringElement();



        for (int i = 0; i < max_iterations; i++)
        {
            const PixelType LowerThreshold = static_cast<PixelType>(MultiThreshold[i]);

            ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
            thresholder->SetInput( maskedImage->GetOutput() );

            thresholder->SetLowerThreshold( LowerThreshold );

            thresholder->SetOutsideValue( OutsideValue );
            thresholder->SetInsideValue( InsideValue );

            thresholder->Update();

            typedef itk::BinaryMorphologicalOpeningImageFilter<typename ThresholdFilterType::OutputImageType, typename ThresholdFilterType::OutputImageType, StructuringElementType> OpeningFilter;
            typedef typename OpeningFilter::Pointer OpeningFilterPointer;

            OpeningFilterPointer filter = OpeningFilter::New();
            filter->SetInput(thresholder->GetOutput());
            filter->SetKernel(openingElement);
            filter->Update();

            typedef itk::BinaryImageToLabelMapFilter<typename OpeningFilter::OutputImageType> LabelFilter;
            typedef typename LabelFilter::Pointer LabelFilterPointer;

            typedef itk::LabelObject<long unsigned int, OpeningFilter::OutputImageType::ImageDimension > LabelObjectType;
            typedef itk::LabelMap<LabelObjectType> LabelMap;
            typedef typename LabelMap::Pointer LabelMapPointer;

            LabelFilterPointer labelFilter = LabelFilter::New();
            labelFilter->SetInput(filter->GetOutput());
            labelFilter->Update();

            LabelMapPointer map = labelFilter->GetOutput();

            // FIXME: i shadows other i
            // XXX: i = 0 is the background
            for (int i = 1; i < map->GetNumberOfLabelObjects(); i++) {
                auto label = map->GetLabelObject(i);

                // Manually calculate bounds
                // FIXME: 99% certain itk already includes functionality for this
                int min[3] = {999999, 999999, 999999};
                int max[3] = {-999999, -999999, -999999};
                for(unsigned int j = 0; j < label->Size(); j++) {
                    auto index = label->GetIndex(j);

                    for (int k = 0; k < 3; k++) {
                        if (index[k] < min[k]) {
                            min[k] = index[k];
                        }
                        if (index[k] > max[k]) {
                            max[k] = index[k];
                        }
                    }
                }

                int mid[3];
                int size[3];

                for (int j = 0; j < 3; j++) {
                    size[j] = max[j] - min[j];
                    mid[j] = (min[j] + max[j]) / 2;
                }

		bool tooSmall = false; // XXX: dumb way to break out of outer loop
		for (int j = 0; j < 3; j++) {
			if (size[j] < minNoduleSize[j]) {
				tooSmall = true;
			}
		}
		if (tooSmall) {
			continue;
		}

		bool tooBig = false;
		for (int j = 0; j < 3; j++) {
			if (size[j] > maxNoduleSize[j]) {
				tooBig = true;
			}
		}
		if (tooBig) {
			continue;
		}

                // Output point to stdout
                std::cout << mid[0] << ',' << mid[1] << ',' << mid[2] << ' ';
                std::cout << size[0] << ", " << size[1] << ", " << size[2] << std::endl;
            }

#           ifdef DISPLAY_SEGMENTED_IMAGES
            displaySlice<SegmentedLungImage>(filter->GetOutput(), 2, 269 - 51);
#           endif
        }
    }
    catch (ExceptionObject &ex) {
        cerr << "The program encountered an exception: " << ex << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
