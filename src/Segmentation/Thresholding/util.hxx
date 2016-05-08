#pragma once

#ifndef CSC621WIND_UTIL_H_
#define CSC621WIND_UTIL_H_

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


#endif /* CSC621WIND_UTIL_H_ */
