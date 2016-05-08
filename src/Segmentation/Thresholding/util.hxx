#pragma once

#ifndef CSC621WIND_UTIL_H_
#define CSC621WIND_UTIL_H_

#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "QuickView.h"

template <typename TInputImage, typename TOutputImage = itk::Image<typename TInputImage::PixelType, TInputImage::ImageDimension> >
typename TOutputImage::Pointer extract2DImageSlice(typename TInputImage::Pointer input, int plane, int slice) {
	typedef itk::ExtractImageFilter < TInputImage, TOutputImage > FilterType;
	typename FilterType::Pointer filter = FilterType::New();

	typename TInputImage::RegionType inputRegion = input->GetLargestPossibleRegion();

	typename TInputImage::SizeType size = inputRegion.GetSize();
	size[plane] = 0;  // collapsing the 3rd dimension

	typename TInputImage::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;  // the required index of the 3rd dimension

	typename TInputImage::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	// setting the direction of collapse
	filter->SetDirectionCollapseToSubmatrix();
	filter->SetInput(input);

	typename TOutputImage::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}

template <typename TViewer, typename TImage, typename TImageDisplay = itk::Image<unsigned char, 2> >
void addImage(TViewer & viewer, typename TImage::Pointer image) {
	typedef itk::RescaleIntensityImageFilter<TImage, TImageDisplay> RescaleFilterType;
	typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	viewer.AddImage(rescaleFilter->GetOutput());
}

template <typename TImage>
void displaySlice(typename TImage::Pointer image, int plane, int sliceIndex ) {
	typedef itk::Image<typename TImage::PixelType, 2> SliceType;
	typename SliceType::Pointer slice = extract2DImageSlice<TImage, SliceType>(image, plane, sliceIndex);

	QuickView viewer;
	addImage<QuickView, SliceType>(viewer, slice);
	viewer.Visualize();
}

#endif /* CSC621WIND_UTIL_H_ */
