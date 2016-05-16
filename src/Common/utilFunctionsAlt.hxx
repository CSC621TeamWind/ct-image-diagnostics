#ifndef utilFunctionsAlt_hxx
#define utilFunctionsAlt_hxx

// ITK
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkPointSetToImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkPNGImageIO.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "itkQuadrilateralCell.h"
#include "itkMeshFileWriter.h"
#include "itkSimplexMesh.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkSimplexMeshVolumeCalculator.h"
#include "itkRegularSphereMeshSource.h"
#include "itkBoundingBox.h"

// VTK
#include "QuickView.h"
#include "vtkVersion.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridWriter.h"

typedef float PixelType;
typedef itk::Image< PixelType, 3 > ImageType3D;
typedef itk::ImageSeriesReader< ImageType3D > ReaderType;
typedef itk::Image< PixelType, 2 > ImageType2D;
typedef itk::ExtractImageFilter < ImageType3D, ImageType2D > FilterType2D;
typedef itk::PointSet< PixelType, 3 > PointSetType;
typedef itk::Image< unsigned char, 2 >  BinaryImageType;
typedef itk::Mesh< float, 3 >   MeshType;
typedef itk::SimplexMesh< float, 3 > SimplexType;
typedef itk::TriangleMeshToSimplexMeshFilter< MeshType, SimplexType > TConvert;
typedef itk::SimplexMeshVolumeCalculator< SimplexType > TVolume;


namespace utility {

  /*
  Reads a DICOM image series given the directory path containing
  the image series.
  */
  ReaderType::Pointer readDicomImageSeries(std::string dirname) {
    ReaderType::Pointer reader = ReaderType::New();

    typedef itk::GDCMImageIO       ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();

    reader->SetImageIO(dicomIO);
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetDirectory(dirname);

    try {
      typedef std::vector< std::string >    SeriesIdContainer;
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

      while (seriesItr != seriesEnd) {
        ++seriesItr;
      }
      std::string seriesIdentifier;
      seriesIdentifier = seriesUID.begin()->c_str();

      typedef std::vector< std::string >   FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      reader->SetFileNames(fileNames);

      try {
        reader->Update();
        return reader;
      }
      catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
      }
    }
    catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
    }
	return NULL;
  }

  /* typedef itk::TriangleCell<MeshType::CellType> TriangleType;
  CellAutoPointer line0;
  line0.TakeOwnership(new TriangleType);
  line0->SetPointId(0, 0); // line between points 0 and 1
  line0->SetPointId(1, 1);
  line0->SetPointId(2, 2);
  mesh->SetCell(0, line0);

  */
  /*
  Display the given  image. No return type.
  */
  void display2DImage(ImageType2D::Pointer image) {
    typedef itk::RescaleIntensityImageFilter<ImageType2D, ImageType2D> RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(image);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->Update();

    QuickView viewer;
    viewer.AddImage(image.GetPointer());
    viewer.AddImage(rescaleFilter->GetOutput());
    viewer.Visualize();
  }

  /*
  Extracts a 2D image slice from a 3D image.
  */
  ImageType2D::Pointer extract2DImageSlice(ImageType3D::Pointer reader, int plane, int slice) {
    FilterType2D::Pointer filter = FilterType2D::New();

    ImageType3D::RegionType inputRegion =
      reader->GetLargestPossibleRegion();

    ImageType3D::SizeType size = inputRegion.GetSize();
    size[plane] = 0;  // collapsing the 3rd dimension

    ImageType3D::IndexType start = inputRegion.GetIndex();
    const unsigned int sliceNumber = slice;
    start[plane] = sliceNumber;  // the required index of the 3rd dimension

    ImageType3D::RegionType desiredRegion;
    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(start);

    filter->SetExtractionRegion(desiredRegion);

    // setting the direction of collapse
    filter->SetDirectionCollapseToSubmatrix();
    filter->SetInput(reader);

    ImageType2D::Pointer img = filter->GetOutput();
    img->Update();

    return img;
  }

  void WriteSliceAsDCM(ImageType3D::Pointer myimage, int sliceNumber, std::string filename) {
	  ImageType2D::Pointer testimage = extract2DImageSlice(myimage, 2, sliceNumber);
	  const float spacing[2] = { 1.0, 1.0 };
	  testimage->SetSpacing(spacing);
	  typedef  itk::ImageFileWriter< ImageType2D > WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetFileName(filename);
	  writer->SetInput(testimage);
	  writer->Update();

  }
  void WriteSliceAsPNG(ImageType3D::Pointer myimage, int sliceNumber, std::string filename) {
	  ImageType2D::Pointer testimage = extract2DImageSlice(myimage, 2, sliceNumber);
	  //WriteSliceAsDCM(myimage, sliceNumber, "Testing.dcm");
	  const float spacing[2] = { 1.0, 1.0 };
	  testimage->SetSpacing(spacing);
	  // Since png only supports unsigned char and unsigned float, 
	  // casting the image to unsigned char.
	  typedef unsigned char                             OutputPixelType;
	  typedef itk::Image< OutputPixelType, 2 >  OutputImageType;

	  typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageType2D >
		  RescaleType;
	  RescaleType::Pointer rescale = RescaleType::New();
	  rescale->SetInput(testimage);
	  rescale->SetOutputMinimum(0);
	  rescale->SetOutputMaximum(itk::NumericTraits< OutputPixelType >::max());

	  typedef itk::CastImageFilter< ImageType2D, OutputImageType > FilterType;
	  FilterType::Pointer filter = FilterType::New();
	  filter->SetInput(rescale->GetOutput());

	  typedef  itk::ImageFileWriter< OutputImageType > WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetImageIO(itk::PNGImageIO::New());
	  writer->SetFileName(filename);
	  writer->SetInput(filter->GetOutput());
	  writer->Update();
  }
  

  /* Creates a binary png image of size width x height. Draws the points of the pointset of slice 'slice'
  on the image.
  */
  void WriteBinaryLabelImage(PointSetType::Pointer pointSet, int width, int height, float sliceNumber, std::string filename) {
	  typedef unsigned short PixelType2; // Since PNG only uses unsigned char and unsigned short
	  typedef itk::Image< PixelType2, 2 >  ImageTypeUnsigned;
	  ImageTypeUnsigned::Pointer binaryImage = ImageTypeUnsigned::New();

	  ImageTypeUnsigned::IndexType start;
	  start.Fill(0);

	  ImageTypeUnsigned::SizeType size;
	  size[0] = width;
	  size[1] = height;

	  ImageTypeUnsigned::RegionType region;
	  region.SetSize(size);
	  region.SetIndex(start);
	  binaryImage->SetRegions(region);
	  binaryImage->Allocate();
	  binaryImage->FillBuffer(0);

	  /*typedef itk::RescaleIntensityImageFilter< ImageTypeUnsigned, ImageTypeUnsigned > RescaleType;
	  RescaleType::Pointer rescale = RescaleType::New();
	  rescale->SetInput(binaryImage);
	  rescale->SetOutputMinimum(0);
	  rescale->SetOutputMaximum(255);
	  rescale->SetOutputMaximum(itk::NumericTraits< PixelType2 >::max());*/

	  typedef PointSetType::PointsContainer      PointsContainer;
	  PointsContainer::Pointer  points2 = pointSet->GetPoints();
	  typedef PointsContainer::Iterator     PointsIterator;
	  PointsIterator  pointIterator = points2->Begin();
	  PointsIterator end = points2->End();
	  while (pointIterator != end)
	  {
		  itk::Point<PixelType, 3> p = pointIterator.Value();
		  if (sliceNumber == p[2]) {
			  ImageTypeUnsigned::IndexType pixelIndexSpherePoint;
			  pixelIndexSpherePoint[0] = p[0];
			  pixelIndexSpherePoint[1] = height - p[1]; // TODO: Figure out the final coordinate system.
			  //std::cout << "The index of the point is " << pointIterator.Index() << std::endl;
			  binaryImage->SetPixel(pixelIndexSpherePoint, 65534);  // TODO: Rescale intensity range.
			  //std::cout << "The required points are " << pixelIndexSpherePoint << std::endl;
		  }
		  ++pointIterator;
	  }

	  // Writing out the image.
	  typedef  itk::ImageFileWriter< ImageTypeUnsigned > WriterType;
	  WriterType::Pointer writer = WriterType::New();
	  writer->SetImageIO(itk::PNGImageIO::New());
	  writer->SetFileName(filename);
	  writer->SetInput(binaryImage);
	  writer->Update();

  }


  class VisitVTKCellsClass
  {
	  vtkCellArray* m_Cells;
	  int* m_LastCell;
	  int* m_TypeArray;
  public:
	  // typedef the itk cells we are interested in
	  typedef itk::CellInterface<
		  MeshType::PixelType,
		  MeshType::CellTraits >  CellInterfaceType;

	  typedef itk::LineCell<CellInterfaceType> floatLineCell;
	  typedef itk::TriangleCell<CellInterfaceType>      floatTriangleCell;
	  typedef itk::QuadrilateralCell<CellInterfaceType> floatQuadrilateralCell;

	  // Set the vtkCellArray that will be constructed
	  void SetCellArray(vtkCellArray* a)
	  {
		  m_Cells = a;
	  }

	  // Set the cell counter pointer
	  void SetCellCounter(int* i)
	  {
		  m_LastCell = i;
	  }

	  // Set the type array for storing the vtk cell types
	  void SetTypeArray(int* i)
	  {
		  m_TypeArray = i;
	  }

	  // Visit a triangle and create the VTK_TRIANGLE cell
	  void Visit(unsigned long, floatTriangleCell* t)
	  {
		  m_Cells->InsertNextCell(3, (vtkIdType*)t->PointIdsBegin());
		  m_TypeArray[*m_LastCell] = VTK_TRIANGLE;
		  (*m_LastCell)++;
	  }

	  // Visit a triangle and create the VTK_QUAD cell
	  void Visit(unsigned long, floatQuadrilateralCell* t)
	  {
		  m_Cells->InsertNextCell(4, (vtkIdType*)t->PointIdsBegin());
		  m_TypeArray[*m_LastCell] = VTK_QUAD;
		  (*m_LastCell)++;
	  }

	  // Visit a line and create the VTK_LINE cell
	  void Visit(unsigned long, floatLineCell* t)
	  {
		  m_Cells->InsertNextCell(2, (vtkIdType*)t->PointIdsBegin());
		  m_TypeArray[*m_LastCell] = VTK_LINE;
		  (*m_LastCell)++;
	  }
  };

  /* Converts a pointset to a mesh with triangle cells. Using TriangleCells since area and volume calculations seem
  to require it.
  */
  MeshType::Pointer ConvertPointsetToTriangleCellMesh2(PointSetType::Pointer pointSet, int n) {

	  MeshType::Pointer mesh = MeshType::New();

	  // Add the points to the mesh
	  int numPoints = pointSet->GetNumberOfPoints();
	  std::cout << "Num points are " << numPoints << std::endl;
	  for (int i = 0; i < numPoints; i++) {
		  MeshType::PointType p = pointSet->GetPoint(i);
		  p[2] = p[2] * 2.5;
		  mesh->SetPoint(i, p);
		  //std::cout << p << std::endl;
		  // FIXME: Make sure we dont have any wrong values.
		  if (p[0] < -5000 || p[1] < -5000 || p[2] < -5000 || p[0] > 5000 || p[1] > 5000 || p[2] > 5000) {
			  std::cout << "Aborting. Potential wrong value at index " << i;
			  exit(1);
		  }
	  }

	  int cellIndex = 0;
	  typedef MeshType::CellType::CellAutoPointer CellAutoPointer;
	  typedef itk::TriangleCell< MeshType::CellType > TriangleType;
	  for (int i = 0; i <numPoints - n; i++) {
		  CellAutoPointer triangle;
		  triangle.TakeOwnership(new TriangleType);
		  triangle->SetPointId(0, i);
		  triangle->SetPointId(1, i + n);
		  int offset = i + n + 1;
		  if (offset%n == 0) {
			  offset = i + 1;
		  }
		  //cout << "Triangle " << i << "  " << (i + n) << "   " << offset;
		  triangle->SetPointId(2, offset);
		  mesh->SetCell(cellIndex, triangle);
		  cellIndex++;

		  // Create the other half of the quad
		  CellAutoPointer triangle2;
		  triangle2.TakeOwnership(new TriangleType);
		  triangle2->SetPointId(0, i);
		  int offset2 = i + 1;
		  if ((i) % n == (n-1)) {
			  offset2 = i - (n - 1);
		  }
		  triangle2->SetPointId(1, offset2);
		  triangle2->SetPointId(2, offset);
		  mesh->SetCell(cellIndex, triangle2);
		  cellIndex++;
		  //cout << "Alt Triangle " << i << "  " << offset2 << "   " << offset << std::endl;
	  }

	  std::cout << "Created a mesh";
	  return mesh;
  }


  /* Converts a pointset to a mesh with triangle cells. Using TriangleCells since area and volume calculations seem
	 to require it.
  */
  MeshType::Pointer ConvertPointsetToTriangleCellMesh(PointSetType::Pointer pointSet, int n) {

	  MeshType::Pointer mesh = MeshType::New();

	  // Add the points to the mesh
	  int numPoints = pointSet->GetNumberOfPoints();
	  std::cout << "Num points are " << numPoints << std::endl;
	  for (int i = 1; i < numPoints; i++) {
		  MeshType::PointType p = pointSet->GetPoint(i);
		  p[2] = p[2] * 2.5;
		  mesh->SetPoint(i, p);
		  //std::cout << p << std::endl;
		  // FIXME: Make sure we dont have any wrong values.
		  if (p[0] < -5000 || p[1] < -5000 || p[2] < -5000 || p[0] > 5000 || p[1] > 5000 || p[2] > 5000) {
			  std::cout << "Aborting. Potential wrong value at index " << i;
			  exit(1);
		  }
	  }

	  int cellIndex = 0;
	  typedef MeshType::CellType::CellAutoPointer CellAutoPointer;
	  typedef itk::TriangleCell< MeshType::CellType > TriangleType;
	  for (int i = 1; i <numPoints-n; i++) {
		  CellAutoPointer triangle;
		  triangle.TakeOwnership(new TriangleType);
		  triangle->SetPointId(0, i);
		  triangle->SetPointId(1, i + n);  
		  int offset = i + n + 1;
		  if (offset%n == 1) {
			  offset = i + 1;
		  }
		  //cout << "Triangle " << i << "  " << (i + n) << "   " << offset;
		  triangle->SetPointId(2, offset);
		  mesh->SetCell(cellIndex, triangle);
		  cellIndex++;

		  // Create the other half of the quad
		  CellAutoPointer triangle2;
		  triangle2.TakeOwnership(new TriangleType);
		  triangle2->SetPointId(0, i);
		  int offset2 = i + 1;
		  if ((i) % n == 0) {
			  offset2 = i - (n-1);
		  }
		  triangle2->SetPointId(1, offset2);
		  triangle2->SetPointId(2, offset);
		  mesh->SetCell(cellIndex, triangle2);
		  cellIndex++;
		  //cout << "Alt Triangle " << i << "  " << offset2 << "   " << offset << std::endl;
	  }

	  std::cout << "Created a mesh";
	  return mesh;
  }

  /* Converts a pointset to a mesh.
  */
  MeshType::Pointer ConvertPointsetToMesh(PointSetType::Pointer pointSet) {

	  MeshType::Pointer mesh = MeshType::New();

	  // Add the points to the mesh
	  int numPoints = pointSet->GetNumberOfPoints();
	  for (int i = 1; i < numPoints; i++) {
		  MeshType::PointType p = pointSet->GetPoint(i);
		  p[2] = p[2] * 2.5;
		  mesh->SetPoint(i, p);
		  //std::cout << p << std::endl;
		  // FIXME: Make sure we dont have any wrong values.
		  if (p[0] < -5000 || p[1] < -5000 || p[2] < -5000 || p[0] > 5000 || p[1] > 5000 || p[2] > 5000) {
			  std::cout << "Aborting. Potential wrong value at index " << i;
			  exit(1);
		  }
	  }

	  // Create the mesh cells - connect all the points

	  int cellIndex = 0;
	  typedef MeshType::CellType::CellAutoPointer CellAutoPointer;
	  typedef itk::LineCell< MeshType::CellType > LineType;
	  for (int i = 1; i < numPoints - 1; i++) {
		  CellAutoPointer line;
		  line.TakeOwnership(new LineType);
		  line->SetPointId(0, i);
		  line->SetPointId(1, i + 1);
		  mesh->SetCell(cellIndex, line);
		  cellIndex++;
		  // connect to longitudinal points.
		  if (i > 20) {
			  CellAutoPointer line2;
			  line2.TakeOwnership(new LineType);
			  line2->SetPointId(0, i);
			  line2->SetPointId(1, i - 20);
			  mesh->SetCell(cellIndex, line2);
			  cellIndex++;
		  }
	  }

	  std::cout << "Created a mesh";
	  return mesh;
  }


  /* Write out an ITK mesh
  */
  void WriteITKMesh(MeshType::Pointer mesh, std::string filename) {
	  typedef itk::MeshFileWriter<MeshType> MeshTypeWriter;
	  MeshTypeWriter::Pointer writer = MeshTypeWriter::New();
	  writer->SetInput(mesh);
	  writer->SetFileName(filename);
	  writer->Update();
	  std::cout << "Written out the itk mesh" << std::endl;
  }

  /* Convert an ITK mesh into VTK unstructured grid.
  */
  void ConvertMeshToUnstructuredGrid(MeshType::Pointer mesh, vtkUnstructuredGrid* unstructuredGrid)
  {
	  // Get the number of points in the mesh
	  int numPoints = mesh->GetNumberOfPoints();
	  if (numPoints == 0)
	  {
		  mesh->Print(std::cerr);
		  std::cerr << "no points in Grid " << std::endl;
		  exit(-1);
	  }

	  // Create the vtkPoints object and set the number of points
	  vtkPoints* vpoints = vtkPoints::New();
	  vpoints->SetNumberOfPoints(numPoints);
	  // Iterate over all the points in the itk mesh filling in
	  // the vtkPoints object as we go
	  MeshType::PointsContainer::Pointer points = mesh->GetPoints();

	  // In ITK the point container is not necessarily a vector, but in VTK it is
	  vtkIdType VTKId = 0;
	  std::map< vtkIdType, int > IndexMap;

	  for (MeshType::PointsContainer::Iterator i = points->Begin();
	  i != points->End(); ++i, VTKId++)
	  {
		  // Get the point index from the point container iterator
		  IndexMap[VTKId] = i->Index();

		  // Set the vtk point at the index with the the coord array from itk
		  // itk returns a const pointer, but vtk is not const correct, so
		  // we have to use a const cast to get rid of the const
		  vpoints->SetPoint(VTKId, const_cast<float*>(i->Value().GetDataPointer()));
	  }

	  // Set the points on the vtk grid
	  unstructuredGrid->SetPoints(vpoints);

	  // Setup some VTK things
	  int vtkCellCount = 0; // running counter for current cell being inserted into vtk
	  int numCells = mesh->GetNumberOfCells();
	  int *types = new int[numCells]; // type array for vtk
									  // create vtk cells and estimate the size
	  vtkCellArray* cells = vtkCellArray::New();
	  cells->EstimateSize(numCells, 4);

	  // Setup the line visitor
	  typedef itk::CellInterfaceVisitorImplementation<
		  float, MeshType::CellTraits,
		  itk::LineCell< itk::CellInterface<MeshType::PixelType, MeshType::CellTraits > >,
		  VisitVTKCellsClass> LineVisitor;
	  LineVisitor::Pointer lv = LineVisitor::New();
	  lv->SetTypeArray(types);
	  lv->SetCellCounter(&vtkCellCount);
	  lv->SetCellArray(cells);

	  // Setup the triangle visitor
	  typedef itk::CellInterfaceVisitorImplementation<
		  float, MeshType::CellTraits,
		  itk::TriangleCell< itk::CellInterface<MeshType::PixelType, MeshType::CellTraits > >,
		  VisitVTKCellsClass> TriangleVisitor;
	  TriangleVisitor::Pointer tv = TriangleVisitor::New();
	  tv->SetTypeArray(types);
	  tv->SetCellCounter(&vtkCellCount);
	  tv->SetCellArray(cells);

	  // Setup the quadrilateral visitor
	  typedef itk::CellInterfaceVisitorImplementation<
		  float, MeshType::CellTraits,
		  itk::QuadrilateralCell< itk::CellInterface<MeshType::PixelType, MeshType::CellTraits > >,
		  VisitVTKCellsClass> QuadrilateralVisitor;
	  QuadrilateralVisitor::Pointer qv = QuadrilateralVisitor::New();
	  qv->SetTypeArray(types);
	  qv->SetCellCounter(&vtkCellCount);
	  qv->SetCellArray(cells);

	  // Add the visitors to a multivisitor

	  MeshType::CellType::MultiVisitor::Pointer mv =
		  MeshType::CellType::MultiVisitor::New();

	  mv->AddVisitor(tv);
	  mv->AddVisitor(qv);
	  mv->AddVisitor(lv);

	  // Now ask the mesh to accept the multivisitor which
	  // will Call Visit for each cell in the mesh that matches the
	  // cell types of the visitors added to the MultiVisitor
	  mesh->Accept(mv);

	  // Now set the cells on the vtk grid with the type array and cell array
	  unstructuredGrid->SetCells(types, cells);
	  std::cout << "Unstructured grid has " << unstructuredGrid->GetNumberOfCells() << " cells." << std::endl;

	  // Clean up vtk objects
	  cells->Delete();
	  vpoints->Delete();

  }

  /* Write out a VTK unstructured grid (.vtu file)
  */
  void WriteVTKUnstructuredGrid(MeshType::Pointer mesh, std::string filename) {
	  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	  utility::ConvertMeshToUnstructuredGrid(mesh, unstructuredGrid);

	  // Write file
	  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	  writer->SetFileName(filename.c_str());
		#if VTK_MAJOR_VERSION <= 5
		writer->SetInputConnection(unstructuredGrid->GetProducerPort());
		#else
		writer->SetInputData(unstructuredGrid);
		#endif
	  writer->Write();
	  std::cout << "Finished writing out the mesh" << std::endl;
  }

  void CalculateBoundingBox(PointSetType::Pointer pointset) {
	  typedef itk::BoundingBox<itk::IdentifierType, 3, PixelType> BoundingBoxType;
	  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
	  boundingBox->SetPoints(pointset->GetPoints());
	  boundingBox->ComputeBoundingBox();
	  itk::Point<float, 3 > p; 
	  p[0] = 0;
	  p[1] = 0;
	  p[2] = 0;

	  for (int i = 0; i < pointset->GetNumberOfPoints(); i++) {
		  cout << pointset->GetPoint(i) << "  ";
	  }

	  std::cout << "bounds: " << boundingBox->GetBounds() << std::endl;
	  std::cout << "center: " << boundingBox->GetCenter() << std::endl;
	  std::cout << "diagonal length squared: " << boundingBox->GetDiagonalLength2() << std::endl;
  }

  void CalculateAreaAndVolume(MeshType::Pointer triangleMesh) {

	  // Ensure that all cells of the mesh are triangles.
	  for (MeshType::CellsContainerIterator it = triangleMesh->GetCells()->Begin();
	  it != triangleMesh->GetCells()->End();
		  ++it)
	  {
		  MeshType::CellAutoPointer cell;
		  triangleMesh->GetCell(it->Index(), cell);
		  if (3 != cell->GetNumberOfPoints())
		  {
			  std::cerr << "ERROR: All cells must be trianglar." << std::endl;
			  exit(1);
		  }
	  } 
	  /*MeshType::Pointer mesh = MeshType::New();
	  // Create points
	  MeshType::PointType p0, p1, p2, p3;

	  p0[0] = -1.0; p0[1] = -1.0; p0[2] = 0.0; // first  point ( -1, -1, 0 )
	  p1[0] = 1.0; p1[1] = -1.0; p1[2] = 0.0; // second point (  1, -1, 0 )
	  p2[0] = 1.0; p2[1] = 1.0; p2[2] = 0.0; // third  point (  1,  1, 0 )

	  mesh->SetPoint(0, p0);
	  mesh->SetPoint(1, p1);
	  mesh->SetPoint(2, p2);

	  typedef itk::TriangleCell<MeshType::CellType> TriangleType;
	  typedef MeshType::CellType::CellAutoPointer CellAutoPointer;
	  CellAutoPointer line0;
	  line0.TakeOwnership(new TriangleType);
	  line0->SetPointId(0, 0); // line between points 0 and 1
	  line0->SetPointId(1, 1);
	  line0->SetPointId(2, 2);
	  mesh->SetCell(0, line0);*/


	  // Convert the triangle mesh to a simplex mesh.
	  TConvert::Pointer convert = TConvert::New();
	  convert->SetInput(triangleMesh);
	  try
	  {
		  convert->Update();
	  }
	  catch (itk::ExceptionObject &ex) {
		  std::cout << "The program encountered an exception: " << ex << std::endl;
	  }
	  
	  // Calculate the volume and area of the simplex mesh.
	 /* TVolume::Pointer volume = TVolume::New();
	  volume->SetSimplexMesh(convert->GetOutput());
	  volume->Compute();

	  // Compare with the volume and area of an ideal sphere.
	  //std::cout << "Ideal Volume: " << 4.0 / 3.0*M_PI*pow(5.0, 3) << std::endl;
	  std::cout << "Mesh Volume: " << volume->GetVolume() << std::endl;
	  //std::cout << "Ideal Surface Area: " << 4.0*M_PI*pow(5.0, 2) << std::endl;
	  std::cout << "Mesh Surface Area: " << volume->GetArea() << std::endl;*/

  }
  /* DEBUG : Create a test black image with two white squares
  */
  void CreateImage(ImageType2D::Pointer image )
  {
	  ImageType2D::IndexType start;
	  start.Fill(0);

	  ImageType2D::SizeType size;
	  size.Fill(100);

	  ImageType2D::RegionType region;
	  region.SetSize(size);
	  region.SetIndex(start);
	  image->SetRegions(region);
	  image->Allocate();
	  image->FillBuffer(0);

	  itk::ImageRegionIterator<ImageType2D> imageIterator(image, image->GetLargestPossibleRegion());

	  // Make two squares
	  while (!imageIterator.IsAtEnd())
	  {
		  if ((imageIterator.GetIndex()[0] > 5 && imageIterator.GetIndex()[0] < 20) &&
			  (imageIterator.GetIndex()[1] > 5 && imageIterator.GetIndex()[1] < 20))
		  {
			  imageIterator.Set(255);
		  }

		  if ((imageIterator.GetIndex()[0] > 50 && imageIterator.GetIndex()[0] < 70) &&
			  (imageIterator.GetIndex()[1] > 50 && imageIterator.GetIndex()[1] < 70))
		  {
			  imageIterator.Set(255);
		  }
		  ++imageIterator;
	  }
  }
}

#endif