#pragma once

#ifndef CSC621WIND_SEGMENTEDLUNGFILTER_H_
#define CSC621WIND_SEGMENTEDLUNGFILTER_H_

#include "itkImage.h"
#include "itkImageToHistogramFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkImageToImageFilter.h"

#include "OptimalThresholdCalculator.hxx"

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
 *
 */
template <typename TImage, typename TOutputImage = itk::Image<unsigned char, 3> >
class SegmentedLungFilter: public itk::ImageToImageFilter<TImage, TOutputImage> {
public:
    typedef SegmentedLungFilter           Self;
    typedef itk::Object                   Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef TOutputImage OutputImageType;

    itkNewMacro(Self);

    itkTypeMacro(SegmentedLungFilter, itk::ImageToImageFilter);
    typedef typename TImage::Pointer InputImagePointer;
    typedef typename TImage::ConstPointer ConstInputImagePointer;

    //typedef unsigned char SegmentedImagePixelType;
    //typedef itk::Image< SegmentedImagePixelType, 3 > SegmentedImageType;
    typedef TOutputImage SegmentedImage;
    typedef typename SegmentedImage::Pointer SegmentedImagePointer;
    typedef typename TOutputImage::PixelType SegmentedImagePixel;

    typedef itk::BinaryThresholdImageFilter<TImage, SegmentedImage> BinaryImageFilter;
    //typedef BinaryImageFilter AirFleshSegmentFilter;
    typedef BinaryImageFilter AirFleshSegmentedImage;
    typedef typename AirFleshSegmentedImage::Pointer AirFleshSegmentedImagePointer;

    typedef itk::Image<SegmentedImagePixel, 2 > SegmentedSliceImage;

    //SegmentedLungFilter(InputImagePointer image): image(image) { }

    typedef itk::Statistics::ImageToHistogramFilter<TImage> HistogramFilterType;
    typedef typename HistogramFilterType::HistogramType HistogramType;
    typedef typename HistogramType::Pointer HistogramPointer;
    typedef typename HistogramFilterType::Pointer HistogramFilterPointer;
    typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;
    typedef typename HistogramFilterType::HistogramMeasurementType HistogramMeasurementType;
    /**
     * Get histogram of input image
     */
    HistogramPointer getHistogram(ConstInputImagePointer image) {
        HistogramFilterPointer filter = HistogramFilterType::New();

        HistogramSizeType size(1);
        size[0] = 4096; // FIXME: should be more flexible
        filter->SetInput(image);
        //filter->SetAutoMinimumMaxium(true);
        //filter->setHistogramSize(size);

        filter->SetMarginalScale(10); // FIXME: parameter
        filter->Update();

        return filter->GetOutput();
    }

    typedef OptimalThresholdCalculator<HistogramType> ThresholdCalculatorType;
    /**
     * Segment air from flesh using an optimal threshold algorithm
     */
    AirFleshSegmentedImagePointer segmentAirFromFlesh(ConstInputImagePointer image) {
        HistogramPointer histogram = getHistogram(image);

	typename ThresholdCalculatorType::Pointer thresholdCalculator = ThresholdCalculatorType::New();
        thresholdCalculator->SetInputHistogram(histogram);
        thresholdCalculator->Compute();

        HistogramMeasurementType threshold = thresholdCalculator->GetThreshold();

	AirFleshSegmentedImagePointer flesh = AirFleshSegmentedImage::New();
        flesh->SetInput(image);
        flesh->SetLowerThreshold(threshold);

        return flesh;
    }

    // FIXME: Triple-check this
    typedef itk::BinaryFillholeImageFilter<SegmentedSliceImage> HoleFillingFilter;
    typedef itk::SliceBySliceImageFilter<typename BinaryImageFilter::OutputImageType, SegmentedImage, HoleFillingFilter> SegmentedSliceFilter;
    typedef typename SegmentedSliceFilter::Pointer SegmentedSliceFilterPointer;

    typedef itk::SliceBySliceImageFilter<typename AirFleshSegmentedImage::OutputImageType, SegmentedImage, HoleFillingFilter> TorsoHoleFillingFilter;
    typedef TorsoHoleFillingFilter TorsoSegmentedImage;
    typedef typename TorsoSegmentedImage::Pointer TorsoSegmentedImagePointer;

    /**
     * Uses a hole filling algorithm to segment the entire torso
     */
    TorsoSegmentedImagePointer segmentTorso(AirFleshSegmentedImagePointer flesh) {
        TorsoSegmentedImagePointer slices = TorsoSegmentedImage::New();

	// Update flesh before using it
	flesh->Update();

	typename HoleFillingFilter::Pointer filter = HoleFillingFilter::New();
        slices->SetFilter(filter);
        slices->SetInput(flesh->GetOutput());

        // FIXME: 3d-connected labeling to isolate torso (may not be necessary)

        return slices;
    }

    typedef itk::BinaryFunctorImageFilter<SegmentedImage, SegmentedImage, SegmentedImage, LungTorsoSegment<SegmentedImagePixel> > InitialLungsSegmentedImage;
    typedef typename InitialLungsSegmentedImage::Pointer InitialLungsSegmentedImagePointer;

    /**
     * Segment lungs by taking the intersection between air and the torso
     */
    InitialLungsSegmentedImagePointer segmentLungsInitial(AirFleshSegmentedImagePointer air, TorsoSegmentedImagePointer torso) {
        InitialLungsSegmentedImagePointer initialLungs = InitialLungsSegmentedImage::New();

	//Update air and torso before use
	air->Update();
	torso->Update();

        initialLungs->SetInput1(air->GetOutput());
        initialLungs->SetInput2(torso->GetOutput());

        return initialLungs;
    }

    typedef SegmentedImage FinalLungsSegmentedImage;
    typedef typename FinalLungsSegmentedImage::Pointer FinalLungsSegmentedImagePointer;

    /**
     * Calculate the final lung mask from the initial mask
     */
    SegmentedSliceFilterPointer segmentLungsFinal(InitialLungsSegmentedImagePointer initialLungs) {
        FinalLungsSegmentedImagePointer lungs = FinalLungsSegmentedImage::New();

	// Update initial lungs before use
	initialLungs->Update();

        SegmentedSliceFilterPointer slices = SegmentedSliceFilter::New();
	
	typename HoleFillingFilter::Pointer filter = HoleFillingFilter::New();
        slices->SetFilter(filter);
        slices->SetInput(initialLungs->GetOutput());

        return slices;
    }

    void GenerateData() {
        ConstInputImagePointer image = this->GetInput();

        // p.16
        // Segment flesh, this mask is true where there is flesh (!M_i in paper)
	AirFleshSegmentedImagePointer flesh = segmentAirFromFlesh(image);

        // p.16
        // Fill holes in the flesh mask to obtain the torso mask
        // body mask in paper (M_b)
        TorsoSegmentedImagePointer torso = segmentTorso(flesh);

        // Get the intersection between the torso and the air masks to obtain initial lung mask
        // Called secondary lung mask in paper (M_s)
        InitialLungsSegmentedImagePointer initialLungs = segmentLungsInitial(flesh, torso);

        // Fill holes in the initial lung mask to obtain the final lung mask
        SegmentedSliceFilterPointer lungs = segmentLungsFinal(initialLungs);
	lungs->Update();
	this->GraftOutput(lungs->GetOutput());
    }
};

#endif /* CSC621WIND_SEGMENTEDLUNGFILTER_H_ */
