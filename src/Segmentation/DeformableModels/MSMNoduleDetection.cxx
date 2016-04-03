
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "utilFunctions.hxx" 

/*
Detection of nodules using stable mass-spring models.
To be integrated into the segmentation framework to use with different datasets.
*/

// Temporarily using this dataset with a hardcoded seedpoint and radius.
const std::string TEMP_DICOM_DATASET_DIR = "../../../../../datasets/Cornell/SS0016/SS0016-20000101/SS0016/20000101-094600-0-2";
const int TEMP_NODULE_RADUIS = 40;  // the temporary radius of the nodule in question in pixels.
const int TEMP_SEEDPOINT[] = { 112, 229, 83};  // the center point of the nodule detected through preprocessing.
const int TEMP_NODULE_SLICES = 22;  // number of image slices over which the nodule is present.
const int TEMP_n = 10;  // number of mass points per image slice.

int main(int, char *[]) {

  typedef unsigned short PixelType;
  typedef itk::PointSet< PixelType, 3 > PointSetType;

  // read the DICOM image series and print out the number of slices.
  std::cout << "Reading the DICOM image directory : " << TEMP_DICOM_DATASET_DIR << std::endl;
  try {
    utility::ReaderType::Pointer reader = utility::ReaderType::New();
    reader = utility::readDicomImageSeries(TEMP_DICOM_DATASET_DIR);
    typedef std::vector< std::string > FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = reader->GetFileNames();
    std::cout << "The total number of slices are " << fileNames.size() << std::endl;
    std::cout << "Using the file:" << fileNames[fileNames.size()-TEMP_SEEDPOINT[2]] << std::endl;  // this seems to read the images backwards fileName[0] has the largest instead

    // Preparing the image for displaying the pointset.
    utility::ImageType::Pointer image = reader->GetOutput();
    const utility::ImageType::IndexType pixelIndex = { {TEMP_SEEDPOINT[0], TEMP_SEEDPOINT[1], TEMP_SEEDPOINT[2]} };
    utility::ImageType::PixelType pixelValue = image->GetPixel(pixelIndex);

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

    // Set the initial point ID
    int spherePointSetID = 0;
    int zFactor = -(M / 2);
    // Generate the spherical mass points model.
    for (int m = 1; m <= M; m++) {
  	  std::cout << "Processing slice no: " << m + 1 << std::endl;
  	  for (int n = 1; n <= N; n++) {
            // (x, y, z) = (sin(Pi * m/M) cos(2Pi * n/N), sin(Pi * m/M) sin(2Pi * n/N), cos(Pi * m/M))
  		  PointType pSphere;

  		  pSphere[0] = seedPoint[0] + ( radius * (sin(itk::Math::pi * (float)m / M) * cos(2 * itk::Math::pi * (float)n / N)));  // x
  		  pSphere[1] = seedPoint[1] + ( radius * (sin(itk::Math::pi * (float)m / M) * sin(2 * itk::Math::pi * (float)n / N)));; // y 
  		  pSphere[2] = seedPoint[2] + zFactor;  // radius * cos(itk::Math::pi * (float)m / M); // z

  		  // Set the point location
  		  spherePointSet->SetPoint(++spherePointSetID, pSphere);
  		  
  		  // Set the point data
  		  spherePointSet->SetPointData(spherePointSetID, 255);

  		  // Printing out the point generated.
  		  std::cout << "(" << pSphere[0] << ", " << pSphere[1] << ", " << pSphere[2] << ")" << std::endl;

        // Setting the pixel value.
        image->SetPixel(pixelIndex, 255);
  	  }
      zFactor++; // Z refers to the index of the axial slice within the dataset.
      std::cout << std::endl;
    }
  }
  catch (itk::ExceptionObject &ex) {
    std::cout << "The program encountered an exception: " << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
