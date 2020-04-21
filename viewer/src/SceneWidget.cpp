#include "SceneWidget.h"
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>

#include "core/Math.hpp"
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vector>
#include "core/ParticleTracer.hpp"
#include "core/ISampleField.hpp"

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

	const int sampled_field_resolution = 2;
	const double t0 = 0.0;
	const double t1 = 2*PI;
	const double dt = 0.2;

	// 

	ParticleTracer<Vec2f, 2> tracer2;
	ParticleTracer<Vec3f, 3> tracer3;
	ParticleTracer<Vec4f, 4> tracer4;

	ABCFlow abcFlow;
	CenterField centerField_analytic;

	const Vec2i res2(sampled_field_resolution, sampled_field_resolution);
	const BoundingBox2d bb2d(Vec2d(-2, -2), Vec2d(2, 2));
	RegVectorField2f centerField_sampled(res2, bb2d);
	for (size_t i = 0; i < res2[0] * res2[1]; ++i) {
		Vec2i gridCoord = centerField_sampled.GetGridCoord(i);
		Vec2d coord = centerField_sampled.GetCoordAt(gridCoord);
		centerField_sampled.SetVertexDataAt(gridCoord, centerField_analytic.Sample(coord));
	}

	const int nSteps = ceil((t1 - t0) / dt);

	int nPaths = 6;
	std::vector<std::vector<Vec3f>> paths(nPaths);
	std::vector<Vec3f> pathColors(nPaths);
	for (int i = 0; i < nPaths; ++i) paths[i].resize(nSteps + 1);

	{
		int path = 0;
		pathColors[path] = Vec3f(1, 0, 0);
		Vec2f position(0.5f, 0);
		paths[path][0] = Vec3f(position[0], position[1], 0);
		for (int i = 0; i < nSteps; ++i) {
			position = tracer2.traceParticle(centerField_analytic, position, dt);
			paths[path][i+1] = Vec3f(position[0], position[1], 0);
		}
	}
	{
		int path = 3;
		pathColors[path] = Vec3f(1, 0, 0.4f);
		Vec2f position(0.505f, 0);
		paths[path][0] = Vec3f(position[0], position[1], 0);
		for (int i = 0; i < nSteps; ++i) {
			position = tracer2.traceParticle(centerField_sampled, position, dt);
			paths[path][i + 1] = Vec3f(position[0], position[1], 0);
		}
	}

	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(createLineActor(paths[0], pathColors[0]));
	renderer->AddActor(createLineActor(paths[3], pathColors[3]));

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
