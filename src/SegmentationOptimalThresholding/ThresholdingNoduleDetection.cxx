/**
 * URN-NBN-SI-DOC-HAVH299II
 */

#include <string>

#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "utilFunctions.hxx" 
#include <itkHistogramAlgorithmBase.h>
#include <itkHistogramThresholdCalculator.h>
#include <itkImageToHistogramFilter.h>
#include <itkOtsuThresholdCalculator.h>
#include <itkBinaryThresholdImageFilter.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

//typedef unsigned short PixelType;
typedef itk::PointSet< PixelType, 3 > PointSetType;

using itk::HistogramAlgorithmBase;
using itk::SmartPointer;
using itk::NumericTraits;
using itk::Indent;
using itk::Object;
using itk::Statistics::ImageToHistogramFilter;
using itk::BinaryThresholdImageFilter;

template <typename THistogram>
class OptimalThresholdCalculator : public HistogramAlgorithmBase<THistogram>
{
public:
    typedef OptimalThresholdCalculator      Self;
    typedef Object                          Superclass;
    typedef SmartPointer<Self>              Pointer;
    typedef SmartPointer<const Self>        ConstPointer;

    typedef typename THistogram::MeasurementType       MeasurementType;
    typedef typename THistogram::AbsoluteFrequencyType FrequencyType;
    typedef size_t SizeValueType;

    typedef typename NumericTraits< MeasurementType >::RealType MeanType;

    itkNewMacro(Self);

    itkTypeMacro(OptimalThresholdCalculator, Object);

    typedef THistogram  HistogramType;

    /** Typedef for the thresholds output */
    typedef std::vector< MeasurementType > OutputType;

    /** Returns the thresholds vector */
    const OutputType & GetOutput() { return m_Output; }

    const MeasurementType GetThreshold() const { return threshold; }

    void Compute(void) ITK_OVERRIDE {
        if ( this->GetInputHistogram()->GetSize().Size() != 1 ) {
            itkExceptionMacro(<< "Histogram must be 1-dimensional.");
        }

        // FIXME: other types of divergence
        bool converged = false;

        for (int i = 0; i < max_iterations; i++) {
            MeasurementType next_threshold = next();

            MeasurementType delta = next_threshold - threshold;

            threshold = next_threshold;

            if ((-max_error) <= delta && delta <= max_error) {
                converged = true;
                break;
            }
        }

        m_Output.resize(1);
        m_Output[0] = threshold;
    }

    MeasurementType next() const {
        typename THistogram::ConstPointer histogram = this->GetInputHistogram();
        typename THistogram::ConstIterator iter = histogram->Begin();
        typename THistogram::ConstIterator end = histogram->End();

        MeasurementType foreground = NumericTraits<MeasurementType>::ZeroValue(); 
        FrequencyType nforeground = NumericTraits<FrequencyType>::ZeroValue();

        MeasurementType background = NumericTraits<MeasurementType>::ZeroValue(); 
        FrequencyType nbackground = NumericTraits<FrequencyType>::ZeroValue();

        while (iter != end) {
            MeasurementType value = iter.GetMeasurementVector()[0];

            if (value <= threshold) {
                background += value * iter.GetFrequency();
                nbackground += iter.getFrequency();
            } else {
                foreground += value * iter.GetFrequency();
                nforeground += iter.getFrequency();
            }
        }

        MeanType mean_foreground = static_cast<MeanType>(foreground) / static_cast<MeanType>(nforeground);
        MeanType mean_background = static_cast<MeanType>(background) / static_cast<MeanType>(nbackground);

        return static_cast<MeasurementType> ((mean_foreground + mean_background) / 2);
    }

    OptimalThresholdCalculator(MeasurementType initial_threshold) {
        threshold = initial_threshold;
        max_error = static_cast<MeasurementType>(0.001);
        max_iterations = NumericTraits<SizeValueType>::max() - 1;
    }
protected:
    virtual ~OptimalThresholdCalculator() {}

    void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE {
        Superclass::PrintSelf(os, indent);

        os << indent << "Threshold: " << threshold << endl;
    }

    private:
    MeasurementType threshold;
    MeasurementType max_error;
    OutputType m_Output;
    SizeValueType max_iterations;
    OptimalThresholdCalculator(const Self&) ITK_DELETE_FUNCTION;
    void operator=(const Self&) ITK_DELETE_FUNCTION;
};

template <typename TImage>
class ThresholdingNoduleDetection {
public:
    typedef typename TImage::Pointer ImagePointer;
    typedef ImageToHistogramFilter<TImage> HistogramFilterType;
    typedef typename HistogramFilterType::HistogramType HistogramType;
    typedef typename HistogramFilterType::Pointer HistogramFilterPointer;
    typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;
    typedef typename HistogramFilterType::HistogramMeasurementType HistogramMeasurementType;
    typedef OptimalThresholdCalculator<HistogramType> ThresholdCalculatorType;

    typedef unsigned char SegmentedImagePixelType;
    typedef itk::Image< SegmentedImagePixelType, 3 > SegmentedImageType;
    typedef SegmentedImageType::Pointer SegmentedImagePointer;
    typedef BinaryThresholdImageFilter<TImage, SegmentedImageType> BinaryImageFilter;
    typedef typename BinaryImageFilter::Pointer BinaryImageFilterPointer;

    ThresholdingNoduleDetection(ImagePointer image) {
        this->image = image;
    }
    void initial_lung_segmentation() {
        HistogramFilterPointer filter = HistogramFilterType::New();

        HistogramSizeType size(1);
        size[0] = 4096; // FIXME: should be more flexible
        filter->SetInput(image);
        filter->SetAutoMinimumMaxium(true);
        filter->setHistogramSize(size);

        filter->SetMarginalScale(10); // FIXME: parameter
        filter->Update();

        const HistogramType * histogram = filter->GetOutput();

        // 1000 HU = Air (recommended initial guess from paper)
        ThresholdCalculatorType thresholdCalculator(1000);

        //TODO: ImageThresholdFilter

        thresholdCalculator->SetInput(histogram); // FIXME
        thresholdCalculator->Compute();
        HistogramMeasurementType threshold = thresholdCalculator->GetThreshold();

        BinaryImageFilterPointer air = BinaryImageFilter::New();
        air->SetInput(image);
        air->SetLowerThreshold(threshold);


    }
protected:
    ImagePointer image;
    //Histogram histogram;
};


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


    }
    catch (itk::ExceptionObject &ex) {
        std::cout << "The program encountered an exception: " << ex << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
