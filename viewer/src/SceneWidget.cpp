#include "SceneWidget.h"
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkOutlineFilter.h>

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

	const int sampled_field_resolution = 30;
	const double t0 = 0.0;
	const double t1 = 2*PI;
	const double dt = 0.25;

	const double bb_size = 7.0;

	// 

	ParticleTracer<Vec2f, 2> tracer2;
	ParticleTracer<Vec3f, 3> tracer3;
	ParticleTracer<Vec4f, 4> tracer4;

	ABCFlow abcFlow_analytic;
	CenterField centerField_analytic;

	const Vec2i res2(sampled_field_resolution, sampled_field_resolution);
	const BoundingBox2d bb2d(Vec2d(-2, -2), Vec2d(2, 2));
	RegVectorField2f centerField_sampled(res2, bb2d);
	for (size_t i = 0; i < res2[0] * res2[1]; ++i) {
		Vec2i gridCoord = centerField_sampled.GetGridCoord(i);
		Vec2d coord = centerField_sampled.GetCoordAt(gridCoord);
		centerField_sampled.SetVertexDataAt(gridCoord, centerField_analytic.Sample(coord));
	}

	const Vec4i res4(sampled_field_resolution, sampled_field_resolution, sampled_field_resolution, sampled_field_resolution);
	const BoundingBox<Vec4d> bb4d(Vec4d(-bb_size, -bb_size, -bb_size, -bb_size), Vec4d(bb_size, bb_size, bb_size, bb_size));
	RegularGrid<Vec4f,4> abcFlow_sampled(res4, bb4d);
	for (size_t i = 0; i < res4[0] * res4[1] * res4[2] * res4[3]; ++i) {
		Vec4i gridCoord = abcFlow_sampled.GetGridCoord(i);
		Vec4d coord = abcFlow_sampled.GetCoordAt(gridCoord);
		abcFlow_sampled.SetVertexDataAt(gridCoord, abcFlow_analytic.Sample(coord));
	}

	const int nSteps = ceil((t1 - t0) / dt);

	int nPaths = 12;
	std::vector<std::vector<Vec3f>> paths(nPaths);
	std::vector<Vec3f> pathColors(nPaths);
	for (int i = 0; i < nPaths; ++i) paths[i].resize(nSteps + 1);

	for (int p = 0; p < 6;++p) {
		int path = p;
		pathColors[path] = Vec3f(1, 0, 0);
		if (p % 2 == 0) pathColors[path][2] = 0.7f;
		Vec2f position(0.5f + 0.1f*p, 0);
		paths[path][0] = Vec3f(position[0], position[1], 0);
		for (int i = 0; i < nSteps; ++i) {
			Vec2d positiond(position[0], position[1]);
			if (p % 2 == 0) {
				position = tracer2.traceParticle(centerField_sampled, positiond, dt);
			}
			else {
				position = tracer2.traceParticle(centerField_analytic, positiond, dt);
			}
			paths[path][i+1] = Vec3f(position[0], position[1], 0);
		}
	}
	for (int p = 6; p < 12; ++p) {
		int path = p;
		pathColors[path] = Vec3f(0.2f, 0.9f, 0.1f);
		if (p % 2 == 0)pathColors[path][0] = 0.8f;
		Vec4f position(0.5f, -0.2f + 0.05f*p, 0.1f, 0);
		paths[path][0] = Vec3f(position[0], position[1], position[2]);
		for (int i = 0; i < nSteps; ++i) {
			Vec4d positiond(position[0], position[1], position[2], position[3]);
			if (p % 2 == 0) {
				position = tracer4.traceParticle(abcFlow_analytic, positiond, dt);
			}
			else {
				position = tracer4.traceParticle(abcFlow_sampled, positiond, dt);
			}
			paths[path][i + 1] = Vec3f(position[0], position[1], position[2]);
		}
	}

	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkSphereSource> sphere;
	sphere->SetRadius(bb_size);
	sphere->Update();
	vtkNew<vtkOutlineFilter> outline;
	outline->SetInputData(sphere->GetOutput());
	vtkNew<vtkPolyDataMapper> outlineMapper;
	outlineMapper->SetInputConnection(outline->GetOutputPort());
	vtkNew<vtkActor> outlineActor;
	outlineActor->SetMapper(outlineMapper);

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(outlineActor);
	for (int i = 0; i < paths.size(); ++i) {
		renderer->AddActor(createLineActor(paths[i], pathColors[i]));
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
