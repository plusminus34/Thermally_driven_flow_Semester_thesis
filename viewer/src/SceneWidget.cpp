#include "SceneWidget.h"
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>

#include "core//Math.hpp"
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vector>
#include "core/ParticleTracer.hpp"


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

void SceneWidget::CreateTestScene()
{
	const Vec3i res(10, 10, 10);
	const Vec3d p1(-1, -1, -1);
	const Vec3d p2(1, 1, 1);
	const BoundingBox3d bb(p1, p2);
	RegVectorField3f field(res, bb);

	const int numCells = res[0] * res[1] * res[2];
	for (size_t i = 0; i < numCells; ++i)
	{
		Vec3i gridCoord = field.GetGridCoord(i);
		Vec3d coord = field.GetCoordAt(gridCoord);
		Vec3f field_val(coord[1], -coord[0], -0.5*coord[2]);
		field.SetVertexDataAt(gridCoord, field_val);
	}

	ParticleTracer<Vec3f, 3> tracer;

	int nSteps = 128;
	std::vector<Vec3f> dots(nSteps+1);
	Vec3f start(1, 0, 1);
	double dt = 0.1;
	dots[0] = (start);
	for (int i = 0; i < nSteps; ++i) {
		Vec3f res = tracer.traceParticle(field, dots[i], dt);
		dots[i+1] = res;
	}


	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkSphereSource> sphereSource;

	vtkNew<vtkPolyDataMapper> sphereMapper;
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

	vtkNew<vtkActor> sphereActor;
	sphereActor->SetMapper(sphereMapper);
	sphereActor->GetProperty()->SetColor(colors->GetColor4d("Tomato").GetData());

	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(sphereActor);
	renderer->AddActor(createLineActor(dots, Vec3f(1, 0, 0)));
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
