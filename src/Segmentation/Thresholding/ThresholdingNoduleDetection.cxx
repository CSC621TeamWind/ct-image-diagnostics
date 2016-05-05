/**
 * URN-NBN-SI-DOC-HAVH299II
 */

#include <string>

// FIXME: use quotes always?
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "utilFunctions.hxx" 
#include <itkHistogramAlgorithmBase.h>
#include <itkHistogramThresholdCalculator.h>
#include <itkImageToHistogramFilter.h>
#include <itkOtsuThresholdCalculator.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkImageSliceConstIteratorWithIndex.h"
#include <itkBinaryFillholeImageFilter.h>
#include "itkSliceBySliceImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

typedef itk::PointSet< PixelType, 3 > PointSetType;

using itk::HistogramAlgorithmBase;
using itk::SmartPointer;
using itk::NumericTraits;
using itk::Indent;
using itk::Object;
using itk::Statistics::ImageToHistogramFilter;
using itk::BinaryThresholdImageFilter;
using itk::SliceBySliceImageFilter;
using itk::BinaryFunctorImageFilter;
using itk::ExceptionObject;

/**
 * Simple boolean operator to combine the mask of air with the torso mask
 */
template <typename TPixel>
class LungTorsoSegment
{
public:
    ~LungTorsoSegment() {}

    bool operator!=(const LungTorsoSegment &) const { return false; }

    bool operator==(const LungTorsoSegment & other) const { return !( *this != other ); }

    inline TPixel operator()(const TPixel & A, const TPixel & B) const
    { return static_cast<TPixel>((!A) && B); }
};

/**
 * Calculate the optimal image threshold based on a histogram
 */
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
class SegmentedLungFilter {
public:
    typedef typename TImage::Pointer InputImagePointer;

    typedef unsigned char SegmentedImagePixelType;
    typedef itk::Image< SegmentedImagePixelType, 3 > SegmentedImageType;
    typedef SegmentedImageType::Pointer SegmentedImagePointer;

    typedef BinaryThresholdImageFilter<TImage, SegmentedImageType> BinaryImageFilter;
    typedef BinaryImageFilter AirFleshSegmentedImage;
    typedef typename AirFleshSegmentedImage::Pointer AirFleshSegmentedImagePointer;

    typedef itk::Image< SegmentedImagePixelType, 2 > SegmentedSliceImage;

    SegmentedLungFilter(InputImagePointer image): image(image) { }

    typedef ImageToHistogramFilter<TImage> HistogramFilterType;
    typedef typename HistogramFilterType::HistogramType HistogramType;
    typedef typename HistogramType::Pointer HistogramPointer;
    typedef typename HistogramFilterType::Pointer HistogramFilterPointer;
    typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;
    typedef typename HistogramFilterType::HistogramMeasurementType HistogramMeasurementType;
    /**
     * Get histogram of input image
     */
    HistogramPointer getHistogram() {
        HistogramFilterPointer filter = HistogramFilterType::New();

        HistogramSizeType size(1);
        size[0] = 4096; // FIXME: should be more flexible
        filter->SetInput(image);
        filter->SetAutoMinimumMaxium(true);
        filter->setHistogramSize(size);

        filter->SetMarginalScale(10); // FIXME: parameter
        filter->Update();

        return filter->GetOutput();
    }

    typedef OptimalThresholdCalculator<HistogramType> ThresholdCalculatorType;
    /**
     * Segment air from flesh using an optimal threshold algorithm
     */
    AirFleshSegmentedImagePointer segmentAirFromFlesh() {
        HistogramPointer histogram = getHistogram();

        ThresholdCalculatorType thresholdCalculator(1000);
        thresholdCalculator->SetInput(histogram);
        thresholdCalculator->Compute();

        HistogramMeasurementType threshold = thresholdCalculator->GetThreshold();

        AirFleshSegmentedImagePointer flesh = BinaryImageFilter::New();
        flesh->SetInput(image);
        flesh->SetLowerThreshold(threshold);

        return flesh;
    }

    // FIXME: Triple-check this
    typedef itk::BinaryFillholeImageFilter<SegmentedSliceImage> HoleFillingFilter;
    typedef SliceBySliceImageFilter<BinaryImageFilter, SegmentedImageType, HoleFillingFilter> SegmentedSliceFilter;
    typedef typename SegmentedSliceFilter::Pointer SegmentedSliceFilterPointer;
    typedef SegmentedImageType TorsoSegmentedImage;
    typedef typename TorsoSegmentedImage::Pointer TorsoSegmentedImagePointer;

    /**
     * Uses a hole filling algorithm to segment the entire torso
     */
    TorsoSegmentedImagePointer segmentTorso(AirFleshSegmentedImagePointer flesh) {
        SegmentedSliceFilterPointer slices = SegmentedSliceFilter::New();

        TorsoSegmentedImagePointer torso = HoleFillingFilter::New();
        slices->SetFilter(torso);
        slices->SetInput(flesh);

        // FIXME: 3d-connected labeling to isolate torso (may not be necessary)

        return torso;
    }

    typedef BinaryFunctorImageFilter<AirFleshSegmentedImage, TorsoSegmentedImage, LungTorsoSegment<SegmentedImagePixelType>, SegmentedImageType> InitialLungsSegmentedImage;
    typedef typename InitialLungsSegmentedImage::Pointer InitialLungsSegmentedImagePointer;

    /**
     * Segment lungs by taking the intersection between air and the torso
     */
    InitialLungsSegmentedImagePointer segmentLungsInitial(AirFleshSegmentedImagePointer air, TorsoSegmentedImagePointer torso) {
        InitialLungsSegmentedImagePointer initialLungs = InitialLungsSegmentedImage::New();

        initialLungs->SetInput1(air);
        initialLungs->SetInput2(torso);
        initialLungs->Update(); // FIXME: don't need to update... right?

        return initialLungs;
    }

    typedef HoleFillingFilter FinalLungsSegmentedImage;
    typedef typename FinalLungsSegmentedImage::Pointer FinalLungsSegmentedImagePointer;

    /**
     * Calculate the final lung mask from the initial mask
     */
    FinalLungsSegmentedImagePointer segmentLungs(InitialLungsSegmentedImagePointer initialLungs) {
        FinalLungsSegmentedImagePointer lungs = FinalLungsSegmentedImage::New();

        SegmentedSliceFilterPointer slices = SegmentedSliceFilter::New();

        slices->SetFilter(lungs);
        slices->SetInput(initialLungs);

        return lungs;
    }

    FinalLungsSegmentedImagePointer segmentLungs() {
        // p.16
        // Segment flesh, this mask is true where there is flesh (!M_i in paper)
        AirFleshSegmentedImagePointer flesh = segmentAirFromFlesh();

        // Fill holes in the flesh mask to obtain the torso mask
        // body mask in paper (M_b)
        TorsoSegmentedImagePointer torso = segmentTorso(flesh);

        // Get the intersection between the torso and the air masks to obtain initial lung mask
        // Called secondary lung mask in paper (M_s)
        InitialLungsSegmentedImagePointer initialLungs = segmentLungsInitial(flesh, torso);

        // Fill holes in the initial lung mask to obtain the final lung mask
        FinalLungsSegmentedImagePointer lungs = segmentLungs(initialLungs); 

        return lungs;
    }
protected:
    InputImagePointer image;
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
    catch (ExceptionObject &ex) {
        cerr << "The program encountered an exception: " << ex << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
