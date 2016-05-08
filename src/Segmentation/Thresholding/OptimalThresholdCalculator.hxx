#pragma once

#ifndef CSC621WIND_OPTIMALTHRESHOLDCALCULATOR_HXX_
#define CSC621WIND_OPTIMALTHRESHOLDCALCULATOR_HXX_

#include <iostream>
#include <ostream>
#include "itkHistogramAlgorithmBase.h"
#include "itkHistogramThresholdCalculator.h"

/**
 * Calculate the optimal image threshold based on a histogram
 */
template <typename THistogram>
class OptimalThresholdCalculator : public itk::HistogramAlgorithmBase<THistogram>
{
public:
    typedef OptimalThresholdCalculator      Self;
    typedef itk::Object                     Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    typedef typename THistogram::MeasurementType       MeasurementType;
    typedef typename THistogram::AbsoluteFrequencyType FrequencyType;
    typedef size_t SizeValueType;

    typedef typename itk::NumericTraits< MeasurementType >::RealType MeanType;

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

        MeasurementType foreground = itk::NumericTraits<MeasurementType>::ZeroValue(); 
        FrequencyType nforeground = itk::NumericTraits<FrequencyType>::ZeroValue();

        MeasurementType background = itk::NumericTraits<MeasurementType>::ZeroValue(); 
        FrequencyType nbackground = itk::NumericTraits<FrequencyType>::ZeroValue();

        while (iter != end) {
            MeasurementType value = iter.GetMeasurementVector()[0];

            if (value <= threshold) {
                background += value * iter.GetFrequency();
                nbackground += iter.GetFrequency();
            } else {
                foreground += value * iter.GetFrequency();
                nforeground += iter.GetFrequency();
            }
	    ++iter;
        }

        MeanType mean_foreground = static_cast<MeanType>(foreground) / static_cast<MeanType>(nforeground);
        MeanType mean_background = static_cast<MeanType>(background) / static_cast<MeanType>(nbackground);

        return static_cast<MeasurementType> ((mean_foreground + mean_background) / 2);
    }

protected:
    OptimalThresholdCalculator(/*MeasurementType initial_threshold*/) {
	// FIXME: hardcoded
        //threshold = initial_threshold;
        threshold = 1000;
        max_error = static_cast<MeasurementType>(2);
        max_iterations = itk::NumericTraits<SizeValueType>::max() - 1;
    }
    //OptimalThresholdCalculator() {}
    virtual ~OptimalThresholdCalculator() {}

    void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE {
        Superclass::PrintSelf(os, indent);

        os << indent << "Threshold: " << threshold << std::endl;
    }

private:
    MeasurementType threshold;
    MeasurementType max_error;
    OutputType m_Output;
    SizeValueType max_iterations;
    OptimalThresholdCalculator(const Self&) ITK_DELETE_FUNCTION;
    void operator=(const Self&) ITK_DELETE_FUNCTION;
};


#endif /* CSC621WIND_OPTIMALTHRESHOLDCALCULATOR_HXX_ */
