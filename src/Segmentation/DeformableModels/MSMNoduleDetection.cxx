
#include "itkPointSet.h"
#include "itkPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientImageFilter.h"
#include "utilFunctionsAlt.hxx" 
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"

#include "math.h"
#include <map>
#include <array>

/*
Detection of nodules using stable mass-spring models.
To be integrated into the segmentation framework to use with different datasets.
*/

// Debug data - Temporarily using this dataset with a hardcoded seedpoint and radius.
const std::string TEMP_DICOM_DATASET_DIR = "../../../../../datasets/Cornell/SS0016/SS0016-20000101/SS0016/20000101-094600-0-2";
const std::string TEMP_CANDIDATE_NODULE_FILE = "../../CandidateNodules.txt";
const int TEMP_NODULE_RADUIS = 40;  // the temporary radius of the nodule in question in pixels.
const int TEMP_SEEDPOINT[] = { 112, 229, 83 };  // the center point of the nodule detected through preprocessing.
const int TEMP_NODULE_SLICES = 22;  // number of image slices over which the nodule is present.
const int TEMP_n = 20;  // number of mass points per image slice.

// Globals
ImageType3D::Pointer image;
ImageType3D::Pointer gradientImage;
int height;
int width;
int depth;

const float k = 0.01;   // Spring constant.
const float alpha = 0.2;
const float beta = 0.2;
const float gamma = 0.2;
const float delta = 0.2;
const float epsilon = 0.2;
const int nNeighborhood = 8;

// forward declaration
template <typename TPoint>
class MSMVertexNeighborhood;

/* Represents a single mass spring point.
*/
template <typename TPoint>
class MSMVertex {
public:
	TPoint point;
	int id;
	float eFunctional;		// Total functional energy.
	float eGradient;
	float ePotential;
	float eElastic;
	float eBending;
	float eAttraction;

	MSMVertex() {
	}
	~MSMVertex() {}
};

/* Represents the 8-neighborhood of the mass-sping point.
*/
template<typename TPoint>
class MSMVertexNeighborhood {
public:
	//MSMVertex<TPoint> pixel;
	int id;	// Same as the point identifier in the pointset.
	MSMVertex<TPoint> vertices[8];  // To calculate functional energy in the 8 neighborhood.
	PointSetType::Pointer nPointset[8]; // Represents the 8 neighbors of a pixel.

};

/* Represents a single slice of a candidate nodule.
*/
template <typename TPoint>
class CandidateNoduleSlice {
public:
	CandidateNoduleSlice() {
	}
	~CandidateNoduleSlice() {}
	TPoint midPoint; // Slice mid-point
	int sliceNumber;
	int minId;  // lowest ID for the slice.
	int maxId;  // highest ID for the slice.
	bool isFirst = false; // Is this the first slice ?
	bool isLast = false;  // Is this the last slice ?
	MSMVertex<TPoint> vertices[TEMP_n];
	MSMVertexNeighborhood<TPoint> neighborhood[TEMP_n];  // Each vertex has an 8 neighborhood
	PointSetType::Pointer nodulePointset;
	float eFunctional = 0;

	//typedef itk::SobelEdgeDetectionImageFilter< ImageType2D, ImageType2D > FilterType;
	//FilterType::Pointer gradientFilter = FilterType::New();

	void process() {
		//std::cout << "Find the total functional energy for slice: " << sliceNumber;
		computeFunctionalEnergy();
	}

	// Computes the functional energy of each point in the pointset.
	void computeFunctionalEnergy() {

		computeElasticEnergy();
		computeBendingEnergy();
		computeAttractionEnergy();
		computeGradientEnergy();
		computePotentialEnergy();
		
		float totalFunctional = 0;
		for (int i = 0; i < TEMP_n; i++) {
			// Calculates the weighted sum of all the energies.
			vertices[i].eFunctional = (alpha * vertices[i].eElastic) +
				(beta * vertices[i].eBending) +
				(gamma * vertices[i].eAttraction) +
				(delta * vertices[i].eGradient) +
				(epsilon * vertices[i].ePotential);
			//std::cout << vertices[i].eElastic << " %  " << vertices[i].eBending << " % " << vertices[i].eAttraction << " %  " << vertices[i].eGradient << " % " << vertices[i].ePotential << std::endl;
			//std::cout << "eFunction of Vertex: " << i << "   " << vertices[i].id << " = " << vertices[i].eFunctional << std::endl;
			totalFunctional += vertices[i].eFunctional;
			if (vertices[i].eFunctional < -5000)
				exit(1);
			// Calculate the functional energy of each pixel in the neighborhood.
			for (int j = 0; j < nNeighborhood; j++) {
				neighborhood[i].vertices[j].eFunctional = (alpha * neighborhood[i].vertices[j].eElastic) +
					(beta * neighborhood[i].vertices[j].eBending) +
					(gamma * neighborhood[i].vertices[j].eAttraction) +
					(delta * neighborhood[i].vertices[j].eGradient) +
					(epsilon * neighborhood[i].vertices[j].ePotential);
				//std::cout << neighborhood[i].vertices[j].eElastic << " %  " << neighborhood[i].vertices[j].eBending << " % " << neighborhood[i].vertices[j].eAttraction << " %  " << neighborhood[i].vertices[j].eGradient << " % " << neighborhood[i].vertices[j].ePotential << std::endl;
				//std::cout << "eFunction of Neighborhood Vertex: " << neighborhood[i].vertices[j].id << " = " << neighborhood[i].vertices[j].eFunctional << std::endl;
				if (neighborhood[i].vertices[j].eFunctional < -5000)
					exit(1);
			}
		}
		eFunctional = totalFunctional;
	}

	// Compares the functional energies of the point in question and its 8 neighbors.
	// Evolves the slice by selecting the point with the minimum functional energy.
	void evolveSlice(PointSetType::Pointer myPointset) {
		float candidatePointEnergies[9];
		std::array<int, 2> offset;
		for (int i = 0, pointId = minId; i < TEMP_n, pointId <= maxId; i++, pointId++) {
			for (int j = 0; j < nNeighborhood; j++) {
				candidatePointEnergies[j] = neighborhood[i].vertices[j].eFunctional;
			}
			candidatePointEnergies[nNeighborhood] = vertices[i].eFunctional;
			// Finding the index of the array with the minimum value.
			int index = 0;
			for (int k = 0; k <= nNeighborhood; k++)
			{
				//std::cout << "The Value is " << candidatePointEnergies[k] << std::endl;
				if (candidatePointEnergies[k] < candidatePointEnergies[index])
					index = k;
			}
			offset = getXYOffset(index);

			std::cout << "Selecting value " << candidatePointEnergies[index] << " with index: " << index << std::endl;
			std::cout << "Using transformation " << offset[0] << ", " << offset[1] << std::endl;
			
			TPoint p1 = nodulePointset->GetPoint(pointId);
			//std::cout << "Old PointId is " << pointId << "The point is " << p1[0] << "," << p1[1] << ", " << p1[2] << std::endl;
			p1[0] = p1[0] + offset[0];
			p1[1] = p1[1] + offset[1];
			myPointset->SetPoint(pointId, p1);
			//std::cout << "PointId is " << pointId << "The point is " << p1[0] << "," << p1[1] << ", " << p1[2] << std::endl;
		}
	}

	void computeElasticEnergy() {
		// Compute elastic energy for each point
		std::array<int, 2> offset;  // neighborhood offsets
		for (int i = minId; i <= maxId; i++) {
			// Finding the distance between the point and its neighbors.
			itk::Point<float, 3> p1 = nodulePointset->GetPoint(i);  

			float eSpringLeft = 0.0;
			itk::Point<float, 3> p0;
			if ((i - 1) < minId) { // use the last point
				p0 = nodulePointset->GetPoint(maxId);
			}
			else {
				p0 = nodulePointset->GetPoint(i - 1);
			}
			double dist = p0.EuclideanDistanceTo(p1);
			eSpringLeft = 0.5 * k * dist * dist;

			float eSpringRight = 0.0;
			if ((i + 1) > maxId) {
				p0 = nodulePointset->GetPoint(minId);
			}
			else {
				p0 = nodulePointset->GetPoint(i + 1);
			}
			dist = p0.EuclideanDistanceTo(p1);
			eSpringRight = 0.5 * k * dist * dist;

			float eSpringUp = 0.0;
			if (!isFirst) {
				p0 = nodulePointset->GetPoint(i - TEMP_n);
				dist = p0.EuclideanDistanceTo(p1);
				eSpringUp = 0.5 * k * dist * dist;
			}
			
			float eSpringDown = 0.0;
			if (!isLast) {
				p0 = nodulePointset->GetPoint(i + TEMP_n);
				dist = p0.EuclideanDistanceTo(p1);
				eSpringDown = 0.5 * k * dist * dist;
			}
			vertices[i - minId].eElastic = eSpringLeft + eSpringRight + eSpringUp + eSpringDown;
			//std::cout << "eElastic of ID " << i << "is " << vertices[i-minId].eElastic << std::endl;

			// Finding the elastic energy of the neighbours
			for (int j = 0; j < nNeighborhood; j++) {
				offset = getXYOffset(j);
				itk::Point<float, 3> pNeighbor = p1;
				pNeighbor[0] = pNeighbor[0] + offset[0];
				pNeighbor[1] = pNeighbor[1] + offset[1];


				eSpringLeft = 0.0;
				if ((i - 1) < minId) { // use the last point
					p0 = nodulePointset->GetPoint(maxId);
				}
				else {
					p0 = nodulePointset->GetPoint(i - 1);
				}

				dist = p0.EuclideanDistanceTo(pNeighbor);
				eSpringLeft = 0.5 * k * dist * dist;

				eSpringRight = 0.0;
				if ((i + 1) > maxId) {
					p0 = nodulePointset->GetPoint(minId);
				}
				else {
					p0 = nodulePointset->GetPoint(i + 1);
				}
				dist = p0.EuclideanDistanceTo(pNeighbor);
				eSpringRight = 0.5 * k * dist * dist;

				eSpringUp = 0.0;
				if (!isFirst) {
					p0 = nodulePointset->GetPoint(i - TEMP_n);
					dist = p0.EuclideanDistanceTo(pNeighbor);
					eSpringUp = 0.5 * k * dist * dist;
				}

				eSpringDown = 0.0;
				if (!isLast) {
					p0 = nodulePointset->GetPoint(i + TEMP_n);
					dist = p0.EuclideanDistanceTo(pNeighbor);
					eSpringDown = 0.5 * k * dist * dist;
				}
				neighborhood[i-minId].vertices[j].eElastic = eSpringLeft + eSpringRight + eSpringUp + eSpringDown;
				//std::cout << "Elastic Value " << neighborhood[i - minId].vertices[j].eElastic << std::endl;
			}
		}
	}

	void computeGradientEnergy() {
		// Compute gradient energy for each point (External force)
		std::array<int, 2> offset;
		// Compute the gradient of the image slice.
		typedef itk::SobelEdgeDetectionImageFilter< ImageType2D, ImageType2D > SobelFilterType;
		SobelFilterType::Pointer gradientFilter = SobelFilterType::New();
		ImageType2D::Pointer image2d = ImageType2D::New();
		image2d = utility::extract2DImageSlice(image, 2, depth - midPoint[2]);
		gradientFilter->SetInput(image2d);
		ImageType2D::Pointer gradient2d = ImageType2D::New();
		gradient2d = gradientFilter->GetOutput();
		typedef itk::RescaleIntensityImageFilter<ImageType2D, ImageType2D> RescaleFilterType;
		RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetOutputMinimum(0);
		rescaleFilter->SetOutputMaximum(255);
		rescaleFilter->SetInput(gradient2d);
		rescaleFilter->Update();
		//utility::display2DImage(gradient2d);
		// Storing the intensity value at each pixel.
		for (int i = minId; i <= maxId; i++) {
			  ImageType2D::IndexType pixelIndex2D = { {
					nodulePointset->GetPoint(i)[0],
					height - nodulePointset->GetPoint(i)[1],
				} };
			ImageType2D::PixelType pixelValue2D = rescaleFilter->GetOutput()->GetPixel(pixelIndex2D);
			vertices[i - minId].eGradient = -pixelValue2D;
			//std::cout << "The gradient energy is " << vertices[i - minId].eGradient << std::endl;
			// Compute the value of the corresponding energy in the neighborhood . 
			for (int j = 0; j < nNeighborhood; j++) {
				offset = getXYOffset(j);
				const ImageType2D::IndexType npixelIndex2D = { {
						(nodulePointset->GetPoint(i)[0]) + offset[0],
						(height - nodulePointset->GetPoint(i)[1]) + offset[1],
					} };
				ImageType2D::PixelType nValue = rescaleFilter->GetOutput()->GetPixel(npixelIndex2D);
				neighborhood[i - minId].vertices[j].eGradient = -nValue;
				//std::cout << "Grad Value is " << neighborhood[i - minId].vertices[j].eGradient << std::endl;
			}
		}
	}

	void computePotentialEnergy() {
		std::array<int, 2> offset;
		for (int i = minId; i <= maxId; i++) { 
			const ImageType3D::IndexType pixelIndex3D = { {
				nodulePointset->GetPoint(i)[0],
				height - nodulePointset->GetPoint(i)[1],
				depth - nodulePointset->GetPoint(i)[2] 
			} };
			ImageType3D::PixelType pixelValue3D = image->GetPixel(pixelIndex3D);
			//std::cout << pixelValue3D << " ";
			vertices[i - minId].ePotential = pixelValue3D;
			//std::cout << "The potential energy is " << vertices[i - minId].ePotential << std::endl;
			
			// Compute the value of the neighborhood. 
			for (int j = 0; j < nNeighborhood; j++) {
				offset = getXYOffset(j);
				const ImageType3D::IndexType npixelIndex3D = { {
						(nodulePointset->GetPoint(i)[0]) + offset[0],
						(height - nodulePointset->GetPoint(i)[1]) + offset[1],
						depth - nodulePointset->GetPoint(i)[2]
					} };
				ImageType3D::PixelType nValue = image->GetPixel(npixelIndex3D);
				neighborhood[i - minId].vertices[j].ePotential = nValue;
				//std::cout << "Pot Value is " << nValue << std::endl;
			}
		}
	}

	std::array<int, 2> getXYOffset(int i) {
		std::array<int, 2> offset;
		switch (i) {
		case 0:
			offset[0] = -1;
			offset[1] = -1;
			break;
		case 1:
			offset[0] = 0;
			offset[1] = -1;
			break;
		case 2:
			offset[0] = 1;
			offset[1] = -1;
			break;
		case 3:
			offset[0] = -1;
			offset[1] = 0;
			break;
		case 4:
			offset[0] = 1;
			offset[1] = 0;
			break;
		case 5:
			offset[0] = -1;
			offset[1] = 1;
			break;
		case 6:
			offset[0] = 0;
			offset[1] = 1;
			break;
		case 7:
			offset[0] = 1;
			offset[1] = 1;
			break;
		default:
			offset[0] = 0;
			offset[1] = 0;
			break;
		}
		
		return offset;
	}

	void computeBendingEnergy() {
		std::array<int, 2> offset;
		for (int i = minId; i <= maxId; i++) {
			// Finding the X and Y calculation for position. 
			itk::Point<float, 3> p0;
			itk::Point<float, 3> p1;
			itk::Point<float, 3> p2;
			p1 = nodulePointset->GetPoint(i);

			if (i - 1 < minId) {
				p0 = nodulePointset->GetPoint(maxId); // wrapping around
			}
			else {
				p0 = nodulePointset->GetPoint(i - 1);
			}

			if (i + 1 > maxId) {
				p2 = nodulePointset->GetPoint(minId);
			}
			else {
				p2 = nodulePointset->GetPoint(i + 1);
			}
			// Calculating the bending energy for X
			float EbendX = p0[0] - (2 * p1[0]) + p2[0];
			// Calculating the bending energy for Y
			float EbendY = p0[1] - (2 * p1[1]) + p2[1];
			vertices[i - minId].eBending = EbendX + EbendY;
			//std::cout << "The bending energy is " << vertices[i - minId].eBending << std::endl;

			// Compute the bending energy of the neighbors. 
			// Finding the elastic energy of the neighbours
			for (int j = 0; j < nNeighborhood; j++) {
				offset = getXYOffset(j);
				itk::Point<float, 3> pNeighbor = p1;
				pNeighbor[0] = pNeighbor[0] + offset[0];
				pNeighbor[1] = pNeighbor[1] + offset[1];

				if (i - 1 < minId) {
					p0 = nodulePointset->GetPoint(maxId); // wrapping around
				}
				else {
					p0 = nodulePointset->GetPoint(i - 1);
				}

				if (i + 1 > maxId) {					  // wrapping around
					p2 = nodulePointset->GetPoint(minId);
				}
				else {
					p2 = nodulePointset->GetPoint(i + 1);
				}
				// Calculating the bending energy for X
				EbendX = p0[0] - (2 * pNeighbor[0]) + p2[0];
				// Calculating the bending energy for Y
				EbendY = p0[1] - (2 * pNeighbor[1]) + p2[1];
				neighborhood[i-minId].vertices[j].eBending = EbendX + EbendY;

				//std::cout << "The value of " << neighborhood[i - minId].vertices[j].eBending << std::endl;
			}

		}
	}

	void computeAttractionEnergy() {
		std::array<int, 2> offset;
		// Calculate the average distance of each point from the centroid.
		// Calculate the standard deviation of all the distances.
		// if the distance is greater than the average + SD compute the attraction energy
		float seedPointDist[TEMP_n]; // Distance of the point from the midpoint of the slice.
		float avgDist = 0;
		float sd = 0;
		for (int i = minId; i <= maxId; i++) {
			// Finding the distance between two points. 
			itk::Point<float, 3> p0 = nodulePointset->GetPoint(i);
			seedPointDist[i-minId] = p0.EuclideanDistanceTo(midPoint);
			avgDist += seedPointDist[i - 1];
		}
		avgDist = avgDist / TEMP_n;  // Find the average of all the distances of the point from the midpoint.
		//std::cout << "The average is " << avgDist;
		sd = standard_deviation(seedPointDist, TEMP_n);
		std::cout << "The slice number is " << sliceNumber << " min id " << minId << "Max id " << maxId << std::endl;
		for (int i = minId; i <= maxId; i++) {
			// Finding the distance between two points. 
			if (seedPointDist[i - minId] > avgDist + sd) {
				vertices[i - minId].eAttraction = seedPointDist[i - minId] / avgDist;
			}
			else {
				vertices[i - minId].eAttraction = 0;
			}
			//std::cout << "The attraction energy is " << vertices[i - minId].eAttraction << std::endl;
		}

		for (int j = 0; j < nNeighborhood; j++) {
			offset = getXYOffset(j);
			itk::Point<float, 3> pNeighbor;

			avgDist = 0;
			sd = 0;
			itk::Point<float, 3> p0;
			for (int i = minId; i <= maxId; i++) {
				// Finding the distance between two points. 
				p0 = nodulePointset->GetPoint(i);
				pNeighbor = p0;
				pNeighbor[0] = pNeighbor[0] + offset[0];
				pNeighbor[1] = pNeighbor[1] + offset[1];
				seedPointDist[i - minId] = pNeighbor.EuclideanDistanceTo(midPoint);
				avgDist += seedPointDist[i - 1];
			}
			avgDist = avgDist / TEMP_n;  // Find the average of all the distances of the point from the midpoint.
										 //std::cout << "The average is " << avgDist;
			sd = standard_deviation(seedPointDist, TEMP_n);

			for (int i = minId; i <= maxId; i++) {
				// Finding the distance between two points. 
				if (seedPointDist[i - minId] > avgDist + sd) {
					neighborhood[i-minId].vertices[j].eAttraction = seedPointDist[i - minId] / avgDist;
				}
				else {
					neighborhood[i - minId].vertices[j].eAttraction = 0;
				}
				//std::cout << "The attraction Value is " << vertices[i - minId].eAttraction << std::endl;
			}
		}
	}
};

/* Represents a candidate nodule within the segmentation operation.
*/
template <typename TPoint>
class CandidateNodule {
public:
	TPoint seedPoint;
	float maxRadius;
	float minRadius;
	int numSlices;
	std::map<int, CandidateNoduleSlice<TPoint>> sliceMap;
	std::map<int, MSMVertexNeighborhood<TPoint>> neighborhoodMap;
	PointSetType::Pointer nodulePointset;
	PointSetType::Pointer evolvedNodulePointset;
	typedef itk::Size<3> noduleSize;
	typedef itk::Index<3> noduleIndex;
	typedef itk::ImageRegion<3> NoduleRegionType;
	NoduleRegionType nRegion(noduleSize, noduleIndex);
	float prevTotalFunctional = 0;

	CandidateNodule() { }
	CandidateNodule(TPoint point, float maxRadius, float minRadius) {
		this.seedPoint = point;
		this.maxRadius = maxRadius;
		this.minRadius = minRadius;
	}

	void initializePointNeighborhood(int spherePointSetID, TPoint pSphere) {
		// instantiate Neighborhood with point ID. 

	}

	PointSetType::Pointer generateInitialSphereModel(int radius) {
		nodulePointset = PointSetType::New();
		evolvedNodulePointset = PointSetType::New();
		PointSetType::Pointer spherePointSet = PointSetType::New();
		// Set the initial point ID
		int spherePointSetID = 0;
		int zFactor = -(radius / 2);
		int M = numSlices;
		int N = TEMP_n;
		// Generate the spherical mass points model.
		for (int m = 1; m <= M; m++) {
			std::cout << "Processing slice no: " << m << std::endl;
			CandidateNoduleSlice<TPoint> noduleSlice;
			noduleSlice.sliceNumber = m;
			noduleSlice.minId = spherePointSetID + 1;
			noduleSlice.maxId = spherePointSetID + TEMP_n;
			noduleSlice.midPoint[0] = seedPoint[0];
			noduleSlice.midPoint[1] = seedPoint[1];
			noduleSlice.midPoint[2] = seedPoint[2] + zFactor;
			if (m == 1)
				noduleSlice.isFirst = true;
			if (m == M)
				noduleSlice.isLast = true;
			for (int n = 1; n <= N; n++) {
				// (x, y, z) = (sin(Pi * m/M) cos(2Pi * n/N), sin(Pi * m/M) sin(2Pi * n/N), cos(Pi * m/M))
				TPoint pSphere;

				pSphere[0] = seedPoint[0] + (radius * (sin(itk::Math::pi * (float)m / M) * cos(2 * itk::Math::pi * (float)n / N)));  // x
				pSphere[1] = seedPoint[1] + (radius * (sin(itk::Math::pi * (float)m / M) * sin(2 * itk::Math::pi * (float)n / N)));; // y 
				pSphere[2] = seedPoint[2] + zFactor;  // radius * cos(itk::Math::pi * (float)m / M); // z

				// Set the ID and Point 
				spherePointSetID++;
				spherePointSet->SetPoint(spherePointSetID, pSphere);

				// Set the point data
				//spherePointSet->SetPointData(spherePointSetID, 255);  // All points are white for now.

				// Printing out the point generated.
				//std::cout << "ID: " << spherePointSetID << "(" << pSphere[0] << ", " << pSphere[1] << ", " << pSphere[2] << ")" << std::endl;

				// Draw the seed point.
				ImageType3D::IndexType pixelIndexSeed;
				pixelIndexSeed[0] = seedPoint[0];
				pixelIndexSeed[1] = height - seedPoint[1];  // The seedpoint given seems to be the point from bottom left
				pixelIndexSeed[2] = depth - pSphere[2];
				//image->SetPixel(pixelIndexSeed, 255);

				// Draw the sphere point.
				ImageType3D::IndexType pixelIndexSpherePoint;
				pixelIndexSpherePoint[0] = pSphere[0];
				pixelIndexSpherePoint[1] = height - pSphere[1];
				pixelIndexSpherePoint[2] = depth - pSphere[2];

				// Setting the pixel value.
				//image->SetPixel(pixelIndexSpherePoint, 255);

				// Adding the MSMVertex to the slice.
				MSMVertex<TPoint> vertex;
				vertex.id = spherePointSetID;
				vertex.point = pSphere;
				noduleSlice.vertices[n-1] = vertex;
			}
			
			zFactor++; // Z refers to the index of the axial slice within the dataset.
			std::cout << std::endl;
			sliceMap[m] = noduleSlice;  // Add the slice into the nodule set. Associate sliceNumber -> CandidateNodule.
		}
		nodulePointset = spherePointSet;
		return spherePointSet;
	}

	void generateInitialCubeModel(int side) {
		// TODO - good to test with.
	}

	int processNodule() {
		// for each slice in the map, process the vertices. 
		// Check for stopping condition.

		float totalFunctional = 0;
		typedef std::map<int, CandidateNoduleSlice<TPoint>>::iterator itType;
		for (itType iterator = sliceMap.begin(); iterator != sliceMap.end(); iterator++) {
			//if (iterator->first == 11) { // TODO: Figure out which label maps need to be written out.
				CandidateNoduleSlice<TPoint> noduleSlice = iterator->second;
				std::cout << "Processing slice: " << noduleSlice.sliceNumber;
				noduleSlice.nodulePointset = PointSetType::New();
				noduleSlice.nodulePointset = nodulePointset;
				//utility::WriteBinaryLabelImage(nodulePointset, width, height, noduleSlice.midPoint[2], "results/TestLabel.png");
				noduleSlice.process();  // calculates all the energies.
				totalFunctional += noduleSlice.eFunctional;
				noduleSlice.evolveSlice(evolvedNodulePointset);
				// std::string filename = "TestLabel_" + std::to_string(noduleSlice.midPoint[2]) + ".png";
				//utility::WriteBinaryLabelImage(evolvedNodulePointset, width, height, noduleSlice.midPoint[2], "results/TestLabel.png");
			//}
		}
		std::cout << "Total functional is " << totalFunctional << std::endl;
		std::cout << "Previous functional is " << prevTotalFunctional << std::endl;

		if (prevTotalFunctional < totalFunctional) {
			std::cout << "Stopping condition reached.";
			return 1;
		}
		// Replace the nodulePointset with the evolved pointset. 
		// FIXME: Figure if i need a pointer or value replacement. For now pointer seems to be working fine.
		nodulePointset = evolvedNodulePointset;
		prevTotalFunctional = totalFunctional;
		return 0;
	}

	void writeImageSlices() {
		typedef std::map<int, CandidateNoduleSlice<TPoint>>::iterator itType;
		for (itType iterator = sliceMap.begin(); iterator != sliceMap.end(); iterator++) {
			if (iterator->first == 11) { // only write out this for now.
				CandidateNoduleSlice<TPoint> noduleSlice = iterator->second;
				std::cout << "Writing slice: " << noduleSlice.midPoint[2];

				//std::string filename = sprintf("results/Slice_%d.png", noduleSlice.midPoint[2]);
				std::string filename = "results/TestSlice_" + std::to_string((int)noduleSlice.midPoint[2]) + ".png";
				utility::WriteSliceAsPNG(image, noduleSlice.midPoint[2], filename);
			}
		}
	}

	float getTotalFunctionalEnergy() {
		float totalEnergy = 0;
		typedef std::map<int, CandidateNoduleSlice<TPoint>>::iterator itType;
		for (itType iterator = sliceMap.begin(); iterator != sliceMap.end(); iterator++) {
			CandidateNoduleSlice<TPoint> noduleSlice = iterator->second;
			std::cout << "Summing up functional energy of slice: " << noduleSlice.sliceNumber;
			std::cout << noduleSlice.vertices[0].eFunctional << "SOMETHING";
			for (int i = 0; i < TEMP_n; i++) {
				// Calculates the weighted sum of all the energies.
				totalEnergy = totalEnergy + (iterator->second).vertices[i].eFunctional;
				std::cout << (iterator->second).vertices[i].eFunctional << "   ";
			}
		}
		std::cout << "The total functional energy of this shape is " << totalEnergy << std::endl;
		return totalEnergy;
	}
};

/* Represents the segmentation operation using Deformable Mass spring models.
*/
template <typename TPoint>
class DeformableMSMSegmentation {
public :
	int numIterations;
	int noduleCount = 0;
	std::map<int, CandidateNodule<TPoint>> noduleList;

	void setIterations(int num) {
		numIterations = num;
	}
	void addNodule(CandidateNodule<TPoint> nodule) {
		noduleList[noduleCount] = nodule;
		noduleCount++;
	}
	void process() {
		typedef std::map<int, CandidateNodule<TPoint>>::iterator itType;
		for (itType iterator = noduleList.begin(); iterator != noduleList.end(); iterator++) {
			std::cout << "Processing nodule " << iterator->first << std::endl;
			CandidateNodule<TPoint> nodule = iterator->second;
			PointSetType::Pointer  spherePointSet = PointSetType::New();
			spherePointSet = nodule.generateInitialSphereModel(nodule.maxRadius);
			nodule.writeImageSlices();

			// Main processing starts here.
			for (int k = 0; k < numIterations; k++) {
				// process each nodule in the list.
				int result = nodule.processNodule();
				if (result == 1) {
					std::cout << "Reached the final result at iteration " << (k + 1) << std::endl;
					break;
				}
				std::cout << "Completed " << k + 1 << " iterations." << std::endl;
			}
		}
	}
};

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

int validateArguments(std::string dicomPath, std::string candidateListPath) {
	// Check whether the dicom directory exists.
	struct stat info;

	if (stat(dicomPath.c_str(), &info) != 0) {
		std::cout << "Cannot access Dicom directory: " << dicomPath << std::endl;
		return -1;
	}

	if (stat(candidateListPath.c_str(), &info) != 0) {
		std::cout << "Cannot access nodule candidate list file: " << candidateListPath << std::endl;
		return -1;
	}
	return 0;
}

int getFileLineCount(std::string filename) {
	std::ifstream f(filename);
	std::string line;
	int i = 0;
	for (i = 0; std::getline(f, line); ++i);
	return i;
}

class CandidateNoduleFileData {
public:
	int id;
	float x;
	float y;
	float z;
	float xDim;
	float yDim;
	float zDim;
};

int main(int argc, char *argv[]) {


  // Accepting arguments from the command line.
  if (argc < 2) {
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << std::endl;
		std::cerr << " <Dicom_Directory> <Candidate_list_file>" << std::endl;
		return EXIT_FAILURE;
   }

  // read the DICOM image series and print out the number of slices.
  std::string dicomDir = argv[1];
  std::string candidateFilePath = argv[2];
  if (dicomDir == "DEBUG") {
	  dicomDir = TEMP_DICOM_DATASET_DIR;
  }
  if (candidateFilePath == "DEBUG") {
	  candidateFilePath = TEMP_CANDIDATE_NODULE_FILE;
  }

  // Validate the arguments
  if (validateArguments(dicomDir, candidateFilePath) != 0) {
	  std::cerr << "One or more arguments are invalid." << std::endl;
	  return EXIT_FAILURE;
  }

  std::cout << "Reading the DICOM image directory : " << dicomDir << std::endl;
  try {

    ReaderType::Pointer reader = ReaderType::New();
    reader = utility::readDicomImageSeries(dicomDir);
    typedef std::vector< std::string > FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = reader->GetFileNames();
    std::cout << "The total number of slices are " << fileNames.size() << std::endl;
    //std::cout << "Using the file:" << fileNames[fileNames.size()-TEMP_SEEDPOINT[2]] << std::endl;  // this seems to read the images backwards fileName[0] has the largest instead

	int numLines = getFileLineCount(candidateFilePath);
	// Read the nodule file.
	std::cout << "Reading the nodule candidate list file: " << candidateFilePath << std::endl;
	std::ifstream ifs;
	ifs.open(candidateFilePath, std::ifstream::in); 
	std::string fileLine;
	int noduleCount = 0;
	std::map<int, CandidateNoduleFileData> noduleMap;

	while (std::getline(ifs, fileLine)) {
		char* dup = _strdup(fileLine.c_str());
		char* token = std::strtok(dup, ",\t");
		if (fileLine.find(",") != std::string::npos) {
			while (token != NULL) {
				CandidateNoduleFileData cf;
				cf.id = ++noduleCount;
				cf.x = atof(token);
				token = std::strtok(NULL, ",\t");
				cf.y = atof(token);
				token = std::strtok(NULL, ",\t");
				cf.z = atof(token);
				token = std::strtok(NULL, ",\t");
				cf.xDim = atof(token);
				token = std::strtok(NULL, ",\t");
				cf.yDim = atof(token);
				token = std::strtok(NULL, ",\t");
				cf.zDim = atof(token);
				token = std::strtok(NULL, ",\t");

				noduleMap[cf.id] = cf;
			}
		}
		if (noduleCount == numLines-1)
			break;
		free(dup);
	}

	// Assigning global variables. TODO: Convert to static ?
    image = reader->GetOutput();
	ImageType3D::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();
	width = imageSize[0];
	height = imageSize[1];
	depth = imageSize[2];
	std::cout << "Width is " << width << " HEight is " << height << " Depth is " << depth;
	// Set origin and spacing for viewing.
	std::cout << "The origin is " << image->GetOrigin(); 

	// Initialize the set of candidate nodules.
	DeformableMSMSegmentation<itk::Point<PixelType, 3>> seg;
	seg.setIterations(20); 
	// Create an add candidate nodules for processing.
	typedef std::map<int, CandidateNoduleFileData>::iterator itType;
	for (itType iterator = noduleMap.begin(); iterator != noduleMap.end(); iterator++) {
		CandidateNodule<itk::Point<PixelType, 3>> nodule;
		CandidateNoduleFileData cf = iterator->second;
		nodule.seedPoint[0] = cf.x;
		nodule.seedPoint[1] = cf.y;
		nodule.seedPoint[2] = cf.z;
		nodule.numSlices = cf.zDim;
		nodule.maxRadius = cf.xDim;
		seg.addNodule(nodule);
	}

	seg.process();

	/*CandidateNodule<itk::Point<PixelType, 3>> nodule;
	nodule.seedPoint[0] = TEMP_SEEDPOINT[0];
	nodule.seedPoint[1] = TEMP_SEEDPOINT[1];
	nodule.seedPoint[2] = TEMP_SEEDPOINT[2];
	nodule.numSlices = TEMP_NODULE_SLICES;

	// Create the nodule list.

	DeformableMSMSegmentation<itk::Point<PixelType, 3>> seg;
	seg.setIterations(3);
	seg.addNodule(nodule);
	seg.process();*/

  }
  catch (itk::ExceptionObject &ex) {
    std::cout << "The program encountered an exception: " << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
