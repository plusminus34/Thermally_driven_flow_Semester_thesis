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

const double M_PI = 3.14159265358979;

class ABCFlow : public ISampleField<Vec4f,4> {
public:
	virtual Vec4f Sample(const Vec<double, 4>& coord) const override {
		return Vec4f(
			((1 - exp(-0.1*coord[3]))*sin(2.0*M_PI*coord[3]) + sqrt(3))*sin(coord[2]) + cos(coord[1]),
			((1 - exp(-0.1*coord[3]))*sin(2.0*M_PI*coord[3]) + sqrt(3))*cos(coord[2]) + sqrt(2)*sin(coord[0]),
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
	const int sampled_field_resolution = 20;

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

	CenterField centerField_analytic;
	const Vec2i res2(sampled_field_resolution, sampled_field_resolution);
	const BoundingBox2d bb2d(Vec2d(-2, -2), Vec2d(2, 2));
	RegVectorField2f centerField_sampled(res2, bb2d);
	for (size_t i = 0; i < res2[0] * res2[1]; ++i) {
		Vec2i gridCoord = centerField_sampled.GetGridCoord(i);
		Vec2d coord = centerField_sampled.GetCoordAt(gridCoord);
		centerField_sampled.SetVertexDataAt(gridCoord, centerField_analytic.Sample(coord));
	}

	ParticleTracer<Vec3f, 3> tracer;

	ParticleTracer<Vec4f, 4Ui64> tracer4;
	ABCFlow abcFlow;

	int nSteps = 128;
	std::vector<Vec3f> dots(nSteps+1);

	std::vector<Vec4f> dots2(nSteps + 1);
	std::vector<Vec3f> dots2_vec3(dots2.size());
	dots2[0] = Vec4f(1, 1, 1, 0);

	std::vector<Vec2f> dots3(nSteps + 1);
	std::vector<Vec3f> dots3_vec3(dots3.size());
	dots3[0] = Vec2f(1, 0);
	ParticleTracer<Vec2f, 2> tracer2;

	Vec3f start(1, 0, 1);
	double dt = 0.1;
	dots[0] = (start);
	for (int i = 0; i < nSteps; ++i) {
		Vec3f res = tracer.traceParticle(field, dots[i], dt);
		dots[i+1] = res;
		dots2[i + 1] = tracer4.traceParticle(abcFlow, dots2[i], dt);
		dots3[i + 1] = tracer2.traceParticle(centerField_analytic, dots3[i], dt);
	}
	for (int i = 0; i < dots2.size();++i) {
		dots2_vec3[i] = Vec3f(dots2[i][0], dots2[i][1], dots2[i][2]);
	}
	for (int i = 0; i < dots3.size(); ++i) {
		dots3_vec3[i] = Vec3f(dots3[i][0], dots3[i][1], 0.f);
	}


	vtkNew<vtkNamedColors> colors;

//	vtkNew<vtkSphereSource> sphereSource;

//	vtkNew<vtkPolyDataMapper> sphereMapper;
//	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

//	vtkNew<vtkActor> sphereActor;
//	sphereActor->SetMapper(sphereMapper);
//	sphereActor->GetProperty()->SetColor(colors->GetColor4d("Tomato").GetData());

	vtkNew<vtkRenderer> renderer;
//	renderer->AddActor(sphereActor);
	renderer->AddActor(createLineActor(dots, Vec3f(1, 0, 0)));
	renderer->AddActor(createLineActor(dots2_vec3, Vec3f(1, 0, 1)));
	renderer->AddActor(createLineActor(dots3_vec3, Vec3f(1, 1, 1)));
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
