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

	int nPaths = 100;//12;
	std::vector<std::vector<Vec3f>> paths(nPaths);
	std::vector<Vec3f> pathColors(nPaths);
	for (int i = 0; i < nPaths; ++i) paths[i].resize(nSteps + 1);

	/*
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
	}*/
	
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

		Vec3i res(dims[0], dims[1], dims[2]);
		Vec3d bb_min(bounds[0], bounds[2], bounds[3]);
		Vec3d bb_max(bounds[1], bounds[3], bounds[5]);
		RegVectorField3f field(res, BoundingBox3d(bb_min, bb_max));

		vtkSmartPointer<vtkPointData> pointData = imageData->GetPointData();
		//vtkSmartPointer<vtkDataArray> dataArray = pointData->GetArray("U");// apparently it's still called W or U ... TODO change that
		vtkSmartPointer<vtkDataArray> dataArray = pointData->GetArray(0);

		int nPts = imageData->GetNumberOfPoints();
		assert(dataArray->GetNumberOfTuples() == nPts);
		cout << "Filling field ...\n";
		for (int i = 0; i < nPts; ++i) {
			double* pt = dataArray->GetTuple3(i);
			//cout << "Pt " << i << ": " << pt[0] << "\t" << pt[1] << "\t" << pt[2] << endl;
			Vec3i coord = field.GetGridCoord(i);
			field.SetVertexDataAt(coord, Vec3f(pt[0], pt[1], pt[2]));
			if (i % (nPts / 10) == 0) {
				cout << "  " << i << " of " << nPts << endl;
			}
		}

		Vec3f midpoint(0.5*(bounds[0] + bounds[1]), 0.5*(bounds[2] + bounds[3]), 0.5*(bounds[4] + bounds[5]));
		Vec3f eee((bounds[1] - bounds[0]) / dims[0], (bounds[3] - bounds[2]) / dims[1], (bounds[5] - bounds[4]) / dims[2]);
		for (int i = 0; i < 10; ++i) {
			for (int j = 0; j < 10; ++j) {
				int path = i * 10 + j;
				pathColors[path] = Vec3f(0.7f + i*0.01f, 0.9f + j*0.01f, 0.3f);
				Vec3f position = midpoint + Vec3f(i, j, 0)*eee;
				float z0 = position[2];
				paths[path][0] = Vec3f(position[0], position[1], position[2]);
				for (int k = 0; k < nSteps; ++k) {
					Vec3d positiond(position[0], position[1], position[2]);
					position = tracer3.traceParticle(field, positiond, dt);
					position[2] = z0;
					paths[path][k + 1] = Vec3f(position[0], position[1], position[2]);
					if (position[0] < bounds[0] || position[0] > bounds[1] ||
						position[1] < bounds[2] || position[1] > bounds[3] ||
						position[2] < bounds[4] || position[2] > bounds[5]) {
						paths[path].resize(k + 2);
						break;
					}
				}
			}
			Vec3f eyo = field.Sample(Vec3d(0.5*(bounds[0] + bounds[1]), 0.5*(bounds[3] + bounds[2]), 0.5*(bounds[4] + bounds[5])));
			cout << "sampled " << eyo[0] << " " << eyo[1] << " " << eyo[2] << endl;
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
	double scale[] = { bounds[1] - bounds[0], bounds[3] - bounds[2],bounds[5] - bounds[4] };
	outlineActor->SetScale(scale);
	outlineActor->SetPosition(midpoint[0], midpoint[1], midpoint[2]);

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(outlineActor);
	for (int i = 0; i < paths.size(); ++i) {
		renderer->AddActor(createLineActor(paths[i], pathColors[i]));
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
