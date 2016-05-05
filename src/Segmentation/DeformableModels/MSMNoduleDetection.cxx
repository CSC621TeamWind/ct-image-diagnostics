
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientImageFilter.h"
#include "utilFunctions.hxx" 
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "math.h"
/*
Detection of nodules using stable mass-spring models.
To be integrated into the segmentation framework to use with different datasets.
*/

// Temporarily using this dataset with a hardcoded seedpoint and radius.
const std::string TEMP_DICOM_DATASET_DIR = "../../../../../datasets/Cornell/SS0016/SS0016-20000101/SS0016/20000101-094600-0-2";
const int TEMP_NODULE_RADUIS = 40;  // the temporary radius of the nodule in question in pixels.
const int TEMP_SEEDPOINT[] = { 112, 229, 83 };  // the center point of the nodule detected through preprocessing.
const int TEMP_NODULE_SLICES = 22;  // number of image slices over which the nodule is present.
const int TEMP_n = 20;  // number of mass points per image slice.

float standard_deviation(float data[], int n)
{
	float mean = 0.0, sum_deviation = 0.0;
	int i;
	for (i = 0; i<n; ++i)
	{
		mean += data[i];
	}
	mean = mean / n;
	for (i = 0; i<n; ++i)
		sum_deviation += (data[i] - mean)*(data[i] - mean);
	return sqrt(sum_deviation / n);
}

int main(int, char *[]) {


  // read the DICOM image series and print out the number of slices.
  std::cout << "Reading the DICOM image directory : " << TEMP_DICOM_DATASET_DIR << std::endl;
  try {
    ReaderType::Pointer reader = ReaderType::New();
    reader = utility::readDicomImageSeries(TEMP_DICOM_DATASET_DIR);
    typedef std::vector< std::string > FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = reader->GetFileNames();
    std::cout << "The total number of slices are " << fileNames.size() << std::endl;
    std::cout << "Using the file:" << fileNames[fileNames.size()-TEMP_SEEDPOINT[2]] << std::endl;  // this seems to read the images backwards fileName[0] has the largest instead

    // Preparing the image for displaying the pointset.
    ImageType3D::Pointer image = reader->GetOutput();
    const ImageType3D::IndexType pixelIndex = { {TEMP_SEEDPOINT[0], TEMP_SEEDPOINT[1], TEMP_SEEDPOINT[2]} };
    ImageType3D::PixelType pixelValue = image->GetPixel(pixelIndex);

    // Generating a PointSet for the initial spherical MSM.
    PointSetType::Pointer  spherePointSet = PointSetType::New();
    typedef PointSetType::PointType PointType;
   
    int radius = TEMP_NODULE_RADUIS;    
    int M = TEMP_NODULE_SLICES;   // M is the number of slices. (latitude)
    int N = TEMP_n;   // N is the number of points per slice. (longitude)

    PointType seedPoint;
    seedPoint[0] = TEMP_SEEDPOINT[0];
    seedPoint[1] = TEMP_SEEDPOINT[1];
    seedPoint[2] = TEMP_SEEDPOINT[2];

	// Finding the size of the image.
	// Extracting the 2d slice of interest from the 3d image and finding its dimensions. 
	ImageType2D::Pointer image2d = utility::extract2DImageSlice(image, 2, fileNames.size() - seedPoint[2]);
	ImageType2D::SizeType imageSize = image2d->GetLargestPossibleRegion().GetSize();
	int imageWidth = imageSize[0];
	int imageHeight = imageSize[1];
	
    // Set the initial point ID
    int spherePointSetID = 0;
    int zFactor = -(M / 2);
    // Generate the spherical mass points model.
    for (int m = 1; m <= M; m++) {
      std::cout << "Processing slice no: " << m << std::endl;
      for (int n = 1; n <= N; n++) {
          // (x, y, z) = (sin(Pi * m/M) cos(2Pi * n/N), sin(Pi * m/M) sin(2Pi * n/N), cos(Pi * m/M))
        PointType pSphere;

        pSphere[0] = seedPoint[0] + ( radius * (sin(itk::Math::pi * (float)m / M) * cos(2 * itk::Math::pi * (float)n / N)));  // x
        pSphere[1] = seedPoint[1] + ( radius * (sin(itk::Math::pi * (float)m / M) * sin(2 * itk::Math::pi * (float)n / N)));; // y 
        pSphere[2] = seedPoint[2] + zFactor;  // radius * cos(itk::Math::pi * (float)m / M); // z

        // Set the point location
        spherePointSet->SetPoint(++spherePointSetID, pSphere);

        // Set the point data
        spherePointSet->SetPointData(spherePointSetID, 255);  // All points are white for now.

        // Printing out the point generated.
        std::cout << "ID: "<< spherePointSetID << "(" << pSphere[0] << ", " << pSphere[1] << ", " << pSphere[2] << ")" << std::endl;

		// Draw the seed point.
		ImageType3D::IndexType pixelIndexSeed;
		pixelIndexSeed[0] = seedPoint[0];
		pixelIndexSeed[1] = imageHeight - seedPoint[1];  // The seedpoint given seems to be the point from bottom left
		pixelIndexSeed[2] = fileNames.size() - pSphere[2]; 
		image->SetPixel(pixelIndexSeed, 255);

		// Draw the sphere point.
		ImageType3D::IndexType pixelIndexSpherePoint;
		pixelIndexSpherePoint[0] = pSphere[0];
		pixelIndexSpherePoint[1] = imageHeight - pSphere[1];
		pixelIndexSpherePoint[2] = fileNames.size() - pSphere[2];

        // Setting the pixel value.
        //image->SetPixel(pixelIndexSpherePoint, 255);
      }
      zFactor++; // Z refers to the index of the axial slice within the dataset.
      std::cout << std::endl;

    }

	// Extracting the slice with the largest radius containing the mass points.
	image2d = utility::extract2DImageSlice(image, 2, fileNames.size() - seedPoint[2]);

    // displaying the 2D slice.
    std::cout << "Displaying the file:" << fileNames[fileNames.size() - TEMP_SEEDPOINT[2]] << std::endl;
    utility::display2DImage(image2d);
	
    // Calculating the external forces.
    // Finding the mangitude of the gradient of the image.
    typedef itk::Image< unsigned char, 2 >  UnsignedCharImageType;
    typedef itk::Image< float, 2 >   FloatImageType;
    typedef itk::GradientMagnitudeImageFilter<ImageType2D, ImageType2D >  filterType;
    filterType::Pointer gradientFilter = filterType::New();
    gradientFilter->SetInput(image2d);
    gradientFilter->Update();
    std::cout << "Finding the magnitude of the image gradient:" << std::endl;
    utility::display2DImage(gradientFilter->GetOutput());  // Displaying the gradient magnitude image

	// Test printing the values of a line after calculating its gradient magnitude.
	/*for (int j = 0; j < imageWidth; j++) {
		const ImageType2D::IndexType pixelIndex2D = { { j, 20 } };
		ImageType2D::PixelType pixelValue2D = gradientFilter->GetOutput()->GetPixel(pixelIndex2D);
		std::cout << pixelValue2D << "  ";
	}*/



	//for extracting a scalar from the vector image
	typedef float       OutputPixelTypeImage;
	typedef float       ComponentType;
	typedef  itk::CovariantVector<ComponentType, 2> OutputPixelType;
	typedef  itk::Image <OutputPixelType, 2> NewOutputImageType;
	typedef  itk::VectorIndexSelectionCastImageFilter<NewOutputImageType, ImageType2D> SelectionFilterType; // < intputType , outputType

		// Testing the gradient filter.
	typedef float OperatorType;
	typedef itk::GradientImageFilter<ImageType2D, OperatorType, OperatorType, NewOutputImageType> gradientImageFilterType;

	gradientImageFilterType::Pointer gradientPointer = gradientImageFilterType::New();
	gradientPointer->SetInput(image2d);
	gradientPointer->SetUseImageSpacingOff();
	gradientPointer->Update();

	SelectionFilterType::Pointer componentExtractor_x = SelectionFilterType::New();
	SelectionFilterType::Pointer componentExtractor_y = SelectionFilterType::New();

	componentExtractor_x->SetIndex(0);// x component of the gradient
	componentExtractor_y->SetIndex(1);// y component of the gradient

	componentExtractor_x->SetInput(gradientPointer->GetOutput());
	componentExtractor_y->SetInput(gradientPointer->GetOutput());

	componentExtractor_x->Update();
	componentExtractor_y->Update();

	std::cout << " The gradient is ......" << std::endl;
	// Test printing the values of a line after calculating its gradient magnitude.
	/*for (int j = 0; j < imageWidth; j++) {
		const ImageType2D::IndexType pixelIndex2D = { { j, 20 } };
		ImageType2D::PixelType pixelValue2D = componentExtractor_y->GetOutput()->GetPixel(pixelIndex2D);
		std::cout << pixelValue2D << "  ";
	}*/

	// Potential energy. 

	// Functional energy. 

	// printing out the origin, spacing and direction of the image.
	const ImageType3D::SpacingType& sp = image->GetSpacing();
	std::cout << "Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << "," << sp[2] << std::endl;

	const ImageType3D::PointType& opt = image->GetOrigin();
	std::cout << "Origin = ";
	std::cout << opt[0] << ", " << opt[1] << ", " << opt[2] << std::endl;

	const ImageType3D::DirectionType& direct = image->GetDirection();
	std::cout << "Direction = ";
	std::cout << direct << std::endl;

	// Getting the coordinates of the points at slice 83. (12th slice within the sphere)
	int offset = 12 -1;
	float Egrad[TEMP_n];
	std::cout << "The gradients at the selected points are "<< std::endl;
	for (int j = 1; j < TEMP_n; j++) { // TEMP_n points on each slice.
		const ImageType2D::IndexType pixelIndex2D = { { spherePointSet->GetPoint((offset*TEMP_n) + j)[0] , spherePointSet->GetPoint((offset*TEMP_n) + j)[1] } };
		std::cout << "Printing the values at ID: "<< (offset*TEMP_n) + j << "   "  << spherePointSet->GetPoint((offset*TEMP_n) + j)[0] << " and " << spherePointSet->GetPoint((offset*TEMP_n) + j)[1] << std::endl;
		ImageType2D::PixelType pixelValue2D = componentExtractor_x->GetOutput()->GetPixel(pixelIndex2D) + componentExtractor_y->GetOutput()->GetPixel(pixelIndex2D);
		std::cout << pixelValue2D << "  ";
		Egrad[j - 1] = pixelValue2D;
	}

	// Calculating the potential energy.
	float Epot[TEMP_n];
	offset = 12 - 1;
	std::cout << "The potential energy at the points are " << std::endl;
	for (int j = 1; j < TEMP_n; j++) { // TEMP_n points on each slice.
		const ImageType2D::IndexType pixelIndex2D = { { spherePointSet->GetPoint((offset*TEMP_n) + j)[0] , spherePointSet->GetPoint((offset*TEMP_n) + j)[1] } };
		std::cout << "Printing the values at ID: " << (offset*TEMP_n) + j << "   " << spherePointSet->GetPoint((offset*TEMP_n) + j)[0] << " and " << spherePointSet->GetPoint((offset*TEMP_n) + j)[1] << std::endl;
		ImageType2D::PixelType pixelValue2D = image2d->GetPixel(pixelIndex2D);
		std::cout << pixelValue2D << "  ";

	}
	std::cout << std::endl << "The second potential energy values are " << std::endl;
	for (int j = 1; j < TEMP_n; j++) { // TEMP_n points on each slice.

		const ImageType3D::IndexType pixelIndex3D = { { spherePointSet->GetPoint((offset*TEMP_n) + j)[0], imageHeight - spherePointSet->GetPoint((offset*TEMP_n) + j)[1], fileNames.size() - spherePointSet->GetPoint((offset*TEMP_n) + j)[2] } };
		ImageType3D::PixelType pixelValue3D = image->GetPixel(pixelIndex3D);
		std::cout << pixelValue3D << " ";
		Epot[j - 1] = pixelValue3D;
	}

	// Elastic energy.
	float Eelastic[TEMP_n];
	for (int j = 2; j < TEMP_n; j++) {
		// Finding the distance between two points. 
		itk::Point<float, 3> p0 = spherePointSet->GetPoint((offset*TEMP_n) + j -1);
		itk::Point<float, 3> p1 = spherePointSet->GetPoint((offset*TEMP_n) + j);
		double dist = p0.EuclideanDistanceTo(p1);
		std::cout << " The distance between both the points is " << dist << std::endl;
		Eelastic[j] = dist;
	}
	

	// Bending energy.
	float Ebend[TEMP_n];
	for (int j = 2; j < TEMP_n-1; j++) {
		// Finding the X and Y calculation for position.  
		itk::Point<float, 3> p0 = spherePointSet->GetPoint((offset*TEMP_n) + j - 1);
		itk::Point<float, 3> p1 = spherePointSet->GetPoint((offset*TEMP_n) + j);
		itk::Point<float, 3> p2 = spherePointSet->GetPoint((offset*TEMP_n) + j + 1);
		// Calculating the bending energy for X
		float EbendX = p0[0] - (2 * p1[0]) + p2[0];

		// Calculating the bending energy for Y
		float EbendY = p0[1] - (2 * p1[1]) + p2[1];

		std::cout << "The Bending energy X = " << EbendX << " The Bending Y = " << EbendY << std::endl;
		Ebend[j - 1] = EbendX + EbendY;
	}

	// Calculating the attraction energy. 
	// Calculate the distance from the points to the center of the slice.
	float seedPointDist[TEMP_n];
	float d[TEMP_n];
	float avgDist = 0;
	for (int j = 1; j < TEMP_n; j++) {
		// Finding the distance between two points. 
		itk::Point<float, 3> p0 = spherePointSet->GetPoint((offset*TEMP_n) + j - 1);
		itk::Point<float, 3> pCenter;
		pCenter[0] = seedPoint[0];
		pCenter[1] = seedPoint[1];
		pCenter[2] = offset;
		seedPointDist[j-1] = p0.EuclideanDistanceTo(pCenter);
		avgDist += seedPointDist[j - 1];
		std::cout << " The distance between both the point and seedpoint is " << seedPointDist[j - 1] << std::endl;
	}
	avgDist = avgDist / (TEMP_n - 1);
	std::cout << "The average is " << avgDist;
	std::cout << "The standard deviation is " << standard_deviation(seedPointDist, TEMP_n);

	// Computing the energy functional.
	float Efunctional[TEMP_n];
	float Eattr[TEMP_n];
	for (int i = 1; i < TEMP_n; i++) {
		Eattr[i] = seedPointDist[i] / avgDist;
		Efunctional[i] = Eelastic[i] + Ebend[i] + Eattr[i] + Egrad[i] + Epot[i];
		std::cout << "The functional energy is " << Efunctional[i];
	}

	ImageType2D::Pointer binaryImage = ImageType2D::New();
	utility::CreateImage(binaryImage);
	utility::display2DImage(binaryImage);
  }
  catch (itk::ExceptionObject &ex) {
    std::cout << "The program encountered an exception: " << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
