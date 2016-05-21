#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCenterOfMass.h>
#include <vtkMassProperties.h>
#include <vtkTriangle.h>
#include <string>

int main ( int argc, char *argv[] )
{
  // Ensure a filename was specified
  if(argc != 2)
    {
    std::cerr << "Usage: " << argv[0] << " InputFilename" << std::endl;
    return EXIT_FAILURE;
    }
 
  // Get the filename from the command line
  std::string inputFilename = argv[1];
 
  // Get all data from the file
  vtkSmartPointer<vtkGenericDataObjectReader> reader = 
      vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();
 
  // All of the standard data types can be checked and obtained like this:
  if (reader->IsFilePolyData())
  {
	  std::cout << "output is a polydata" << std::endl;
	  vtkPolyData* output = reader->GetPolyDataOutput();
	  std::cout << "output has " << output->GetNumberOfPoints() << " points." << std::endl;


	  // Compute the center of mass
	  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
		  vtkSmartPointer<vtkCenterOfMass>::New();

	  centerOfMassFilter->SetInputData(output);

	  centerOfMassFilter->SetUseScalarsAsWeights(false);
	  centerOfMassFilter->Update();

	  double center[3];
	  centerOfMassFilter->GetCenter(center);

	  std::cout << "Center of mass is " << center[0] << " " << center[1] << " " << center[2] << std::endl;
	  vtkMassProperties *massProperty = vtkMassProperties::New();

	  //massProperty->SetInputConnection(triangle->GetOutputPort());

	  massProperty->SetInputData(output);

	  massProperty->Update();

	  double vol = massProperty->GetVolume();

	  std::cout << "The volume is " << vol << std::endl;

	  // Calculating the area: 
	  double totalArea = 0;
	  for (vtkIdType i = 0; i < output->GetNumberOfCells(); i++)
	  {
		  vtkCell* cell = output->GetCell(i);

		  vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
		  double p0[3];
		  double p1[3];
		  double p2[3];
		  triangle->GetPoints()->GetPoint(0, p0);
		 // std::cout << "p0: " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
		  triangle->GetPoints()->GetPoint(1, p1);
		 // std::cout << "p1: " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		  triangle->GetPoints()->GetPoint(2, p2);
		 // std::cout << "p2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;

		  double area = vtkTriangle::TriangleArea(p0, p1, p2);
		  //std::cout << "area of triangle " << i << ": " << area << std::endl;
		  totalArea += area;
	  }
	  std::cout << "The total area is " << totalArea << std::endl;
  }
  return EXIT_SUCCESS;
}