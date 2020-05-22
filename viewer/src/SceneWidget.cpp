#include "SceneWidget.h"
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkOutlineFilter.h>
#include <vtkDelaunay2D.h>

#include "core/Math.hpp"
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vector>
#include "core/ParticleTracer.hpp"
#include "core/ISampleField.hpp"
#include "core/NetCDF.hpp"

#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

//copypasted(mostly) from Exercise 4
// creates an actor that draws a line in a given color
vtkSmartPointer<vtkActor> createLineActor(const std::vector<Vec3f>& line, const Vec3f& color)
{
	// create polyline vertices
	vtkNew<vtkPoints> points;
	for (const Vec3f& vertex : line)
		points->InsertNextPoint(vertex[0], vertex[1], vertex[2]);

	// create polyline indices
	vtkNew<vtkPolyLine> polyLine;
	polyLine->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
	for (unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
		polyLine->GetPointIds()->SetId(i, i);

	// Create a cell array to store the lines in and add the lines to it
	vtkNew<vtkCellArray> cells;
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);	// Add the points to the dataset
	polyData->SetLines(cells);			// Add the lines to the dataset

	// Setup actor and mapper
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputData(polyData);
	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color[0], color[1], color[2]);
	actor->GetProperty()->SetLineWidth(2);
	return actor;
}

SceneWidget::SceneWidget()
{
	mRenderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
	SetRenderWindow(mRenderWindow);

	CreateTestScene();
}

const double PI = 3.14159265358979;

class ABCFlow : public ISampleField<Vec4f,4> {
public:
	virtual Vec4f Sample(const Vec<double, 4>& coord) const override {
		return Vec4f(
			((1 - exp(-0.1*coord[3]))*sin(2.0*PI*coord[3]) + sqrt(3))*sin(coord[2]) + cos(coord[1]),
			((1 - exp(-0.1*coord[3]))*sin(2.0*PI*coord[3]) + sqrt(3))*cos(coord[2]) + sqrt(2)*sin(coord[0]),
			sin(coord[1]) + sqrt(2)*cos(coord[0]),
			1
		);
	}
};

class CenterField : public ISampleField<Vec2f, 2> {
public:
	virtual Vec2f Sample(const Vec<double, 2>& coord) const override {
		return Vec2f( -coord[1], coord[0] );
	}
};

void SceneWidget::CreateTestScene()
{
	// Settings

	const int sampled_field_resolution = 30;
	const double t0 = 0.0;
	const double t1 = 100000.0;
	const double dt = 5;

	const double bb_size = 7.0;

	// display surface
	std::string constantsfile = "../../../lfff00000000c.nc";
	RegScalarField2f* hsurf = NetCDF::ImportScalarField2f(constantsfile, "HSURF", "rlat", "rlon");
	vtkSmartPointer<vtkPoints> hfield_points = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < hsurf->GetData().size(); ++i) {
		Vec2i gridCoord = hsurf->GetGridCoord(i);
		Vec2d coord = hsurf->GetCoordAt(gridCoord);
		hfield_points->InsertNextPoint(coord[1], coord[0], hsurf->GetData()[i]);
	}
	vtkSmartPointer<vtkPolyData> hfield_polydata = vtkSmartPointer<vtkPolyData>::New();
	hfield_polydata->SetPoints(hfield_points);
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(hfield_polydata);
	vtkSmartPointer<vtkPolyData> landscape = delaunay->GetOutput();
	delaunay->Update();// TODO this takes a while
	delete hsurf;
	
	//Read UVW and convert into RegVectorfield
	vtkNew<vtkXMLImageDataReader> imageReader;
	std::string filename = "../cmd/UVW.vti";
	cout << "Reading file " << filename << endl;
	imageReader->SetFileName(filename.c_str());
	imageReader->Update();
	vtkSmartPointer<vtkImageData> imageData = imageReader->GetOutput();
	cout << "Got imageData\n";
	int* dims = imageData->GetDimensions();
	cout << "Dimensions: " << dims[0] << " x " << dims[1] << " x " << dims[2] << endl;
	int num_components = imageData->GetNumberOfScalarComponents();
	cout << "Number of scalar components: " << num_components << endl;
	double* bounds = imageData->GetBounds();
	cout << "Bounds: " << bounds[0] << " - " << bounds[1] << "\t" << bounds[2] << " x " << bounds[3] << "\t" << bounds[4] << " x " << bounds[5] << endl;
	Vec3f midpoint(0.5*(bounds[0] + bounds[1]), 0.5*(bounds[2] + bounds[3]), 0.5*(bounds[4] + bounds[5]));

	vector<vector<Vec3f>> paths;
	NetCDF::ReadPaths("../cmd/imp_paths_first.nc", paths);
	int nPaths = paths.size();
	std::vector<Vec3f> pathColors(nPaths);
	for (int i = 0; i < nPaths; ++i) {
		pathColors[i] = Vec3f(1 - ((float)i / nPaths), 1, 1);
	}

	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkSphereSource> sphere;
	//sphere->SetRadius(bb_size);
	//sphere->SetRadius(bounds[1] - bounds[0]);
	//sphere->SetCenter(midpoint[0], midpoint[1], midpoint[2]);
	sphere->Update();
	vtkNew<vtkOutlineFilter> outline;
	outline->SetInputData(sphere->GetOutput());
	vtkNew<vtkPolyDataMapper> outlineMapper;
	outlineMapper->SetInputConnection(outline->GetOutputPort());
	vtkNew<vtkActor> outlineActor;
	outlineActor->SetMapper(outlineMapper);
	double scale[] = { bounds[1] - bounds[0], bounds[3] - bounds[2], (bounds[5] - bounds[4]) };
	outlineActor->SetScale(scale);
	outlineActor->SetPosition(midpoint[0], midpoint[1], midpoint[2]);

	vtkNew<vtkPolyDataMapper> landscapeMapper;
	landscapeMapper->SetInputData(landscape);
	vtkNew<vtkActor> landscapeActor;
	landscapeActor->SetMapper(landscapeMapper);

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(landscapeActor);
	renderer->AddActor(outlineActor);
	for (int i = 0; i < paths.size(); ++i) {
		renderer->AddActor(createLineActor(paths[i], pathColors[i]));
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
