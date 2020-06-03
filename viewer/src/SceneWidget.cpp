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
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>

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

vtkSmartPointer<vtkActor> createTrajectoryActor(const TrajectoryData& td, int trajectory_id)
{
	int lon_id = td.get_var_id("lon");
	int lat_id = td.get_var_id("lat");
	int z_id = td.get_var_id("z");
	float val = 0; int jojo = 1;
	// create polyline vertices
	vtkNew<vtkPoints> points;
	for (int i = 0; i < td.points_per_trajectory; ++i) {
		points->InsertNextPoint(td.get_value(lon_id, trajectory_id, i), td.get_value(lat_id, trajectory_id, i), td.get_value(z_id, trajectory_id, i)*ZSCALE);
		val += td.get_value(jojo, trajectory_id, i);
	}
	val /= td.points_per_trajectory;

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
	actor->GetProperty()->SetColor(1, 0.5*trajectory_id / td.num_trajectories, val/td.min_values[jojo]);
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
	vtkNew<vtkPolyDataMapper> landscapeMapper;
	bool build_landscape = false;// TODO check if file exists automatically
	if (build_landscape) {
		std::string constantsfile = "../../../lfff00000000c.nc";
		vtkSmartPointer<vtkPoints> hfield_points = vtkSmartPointer<vtkPoints>::New();
		RegScalarField3f* hsurf = NetCDF::ImportScalarField3f(constantsfile, "HSURF", "rlon", "rlat", "time");
		for (int i = 0; i < hsurf->GetData().size(); ++i) {
			Vec3i gridCoord = hsurf->GetGridCoord(i);
			/*
			for some reason delaunay fails (stack overflow) when there are too many points
			therefore: skip a few
			*/
			if (gridCoord[0] % 100 == 1 || gridCoord[1] % 100 == 1) continue;

			Vec3d coord = hsurf->GetCoordAt(gridCoord);
			double lon, lat;
			CoordinateTransform::RlatRlonToLatLon(coord[1], coord[0], lat, lon);
			hfield_points->InsertNextPoint(lon, lat, hsurf->GetVertexDataAt(gridCoord) * ZSCALE);
		}
		delete hsurf;
		vtkSmartPointer<vtkPolyData> hfield_polydata = vtkSmartPointer<vtkPolyData>::New();
		hfield_polydata->SetPoints(hfield_points);
		vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
		delaunay->SetInputData(hfield_polydata);
		vtkSmartPointer<vtkPolyData> landscape = delaunay->GetOutput();
		delaunay->Update();
		landscapeMapper->SetInputData(landscape);
		// write it down
		vtkNew<vtkXMLPolyDataWriter> landscapeWriter;
		landscapeWriter->SetFileName("landscape_reduced.vtp");
		landscapeWriter->SetInputData(landscape);
		landscapeWriter->Write();
	}
	else {
		vtkNew<vtkXMLPolyDataReader> landscapeReader;
		landscapeReader->SetFileName("landscape_reduced.vtp");
		landscapeReader->Update();
		landscapeMapper->SetInputConnection(landscapeReader->GetOutputPort());
	}
	vtkNew<vtkActor> landscapeActor;
	landscapeActor->SetMapper(landscapeMapper);
	
	string file_1 = "../../../outputs/lag_trajectory_compare_2.nc";
	string file_2 = "../../../outputs/imp_trajectory_compare_2_T.nc";
	//string file_1 = "../../../outputs/lag_trajectory_spread_dense.nc";
	//string file_2 = "../../../outputs/imp_trajectory_spread_dense.nc";
	
	TrajectoryData td;
	NetCDF::ReadTrajectoryData(file_1, td);
	assert(td.num_trajectories > 0);
	assert(td.points_per_trajectory > 0);
	int x_id = td.get_var_id("lon");
	int y_id = td.get_var_id("lat");
	int z_id = td.get_var_id("z");
	if (x_id < 0 || y_id < 0 || z_id < 0) return;
	vector<vector<Vec3f>> paths(td.num_trajectories);
	vector<Vec3f> pathColors(paths.size());
	for (int i = 0; i < paths.size(); ++i) {
		paths[i].resize(td.points_per_trajectory);
		for (int j = 0; j < paths[i].size(); ++j) {
			paths[i][j] = Vec3f(td.get_value(x_id, i, j), td.get_value(y_id, i, j), td.get_value(z_id, i, j)*ZSCALE);
		}
		pathColors[i] = Vec3f(1 - ((float)i / paths.size()), 1, 1);
	}

	TrajectoryData td2;
	NetCDF::ReadTrajectoryData(file_2, td2);
	assert(td2.num_trajectories > 0);
	assert(td2.points_per_trajectory > 0);
	x_id = td2.get_var_id("lon");
	y_id = td2.get_var_id("lat");
	z_id = td2.get_var_id("z");
	if (x_id < 0 || y_id < 0 || z_id < 0) return;
	vector<vector<Vec3f>> paths2(td2.num_trajectories);
	vector<Vec3f> pathColors2(paths2.size());
	for (int i = 0; i < paths2.size(); ++i) {
		paths2[i].resize(td2.points_per_trajectory);
		for (int j = 0; j < paths2[i].size(); ++j) {
			paths2[i][j] = Vec3f(td2.get_value(x_id, i, j), td2.get_value(y_id, i, j), td2.get_value(z_id, i, j)*ZSCALE);
		}
		pathColors2[i] = Vec3f(1, 1 - ((float)i / paths2.size()), 0);
	}

	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	bool debugSpheres = false;
	if(debugSpheres){
		vtkNew<vtkSphereSource> sphere;
		sphere->SetRadius(0.5);
		sphere->SetCenter(0, 0, 0);
		sphere->Update();
		vtkNew<vtkPolyDataMapper> sphereMapper;
		sphereMapper->SetInputConnection(sphere->GetOutputPort());

		/*
		vtkNew<vtkActor> sphereActor;
		sphereActor->SetMapper(sphereMapper);
		sphereActor->GetProperty()->SetColor(0.5, 0.5, 0.5);
		renderer->AddActor(sphereActor);

		vtkNew<vtkActor> sphereActorX;
		sphereActorX->SetMapper(sphereMapper);
		sphereActorX->GetProperty()->SetColor(1, 0, 0);
		sphereActorX->SetScale(5, 0.1, 0.1);
		renderer->AddActor(sphereActorX);
		vtkNew<vtkActor> sphereActorY;
		sphereActorY->SetMapper(sphereMapper);
		sphereActorY->GetProperty()->SetColor(0, 1, 0);
		sphereActorY->SetScale(0.1, 5, 0.1);
		renderer->AddActor(sphereActorY);
		vtkNew<vtkActor> sphereActorZ;
		sphereActorZ->SetMapper(sphereMapper);
		sphereActorZ->GetProperty()->SetColor(0, 0, 1);
		sphereActorZ->SetScale(0.1, 0.1, 5);
		renderer->AddActor(sphereActorZ);
		*/

		double rlon, rlat, lon, lat;

		rlon = 0; rlat = 0;
		CoordinateTransform::RlatRlonToLatLon(rlat, rlon, lat, lon);
		vtkNew<vtkActor> sphereActor_loco;
		sphereActor_loco->SetMapper(sphereMapper);
		sphereActor_loco->GetProperty()->SetColor(1, 0, 1);
		sphereActor_loco->SetScale(0.3, 0.3, 1);
		sphereActor_loco->SetPosition(lon, lat, 0);
		renderer->AddActor(sphereActor_loco);

		rlon = 1; rlat = 0;
		CoordinateTransform::RlatRlonToLatLon(rlat, rlon, lat, lon);
		vtkNew<vtkActor> sphereActor_lon;
		sphereActor_lon->SetMapper(sphereMapper);
		sphereActor_lon->GetProperty()->SetColor(1, 0.5, 0);
		sphereActor_lon->SetScale(0.5, 0.2, 0.8);
		sphereActor_lon->SetPosition(lon, lat, 0);
		renderer->AddActor(sphereActor_lon);

		rlon = 0; rlat = 1;
		CoordinateTransform::RlatRlonToLatLon(rlat, rlon, lat, lon);
		vtkNew<vtkActor> sphereActor_lat;
		sphereActor_lat->SetMapper(sphereMapper);
		sphereActor_lat->GetProperty()->SetColor(0, 0.5, 1);
		sphereActor_lat->SetScale(0.2, 0.5, 0.8);
		sphereActor_lat->SetPosition(lon, lat, 0);
		renderer->AddActor(sphereActor_lat);


		for (int i = 0; i <= 10; ++i) {
			for (int j = 0; j <= 5; ++j) {
				double rlon = i * 0.1;
				double rlat = j * 0.1;
				double lon, lat;
				CoordinateTransform::RlatRlonToLatLon(rlat, rlon, lat, lon);
				vtkNew<vtkActor> sphereActor_ij;
				sphereActor_ij->SetMapper(sphereMapper);
				sphereActor_ij->GetProperty()->SetColor(rlon, rlat, 0.5);
				sphereActor_ij->SetScale(0.03, 0.03, 0.1);
				sphereActor_ij->SetPosition(lon, lat, 1);
				renderer->AddActor(sphereActor_ij);
			}
		}
	}

	renderer->AddActor(landscapeActor);
	for (int i = 0; i < paths.size(); ++i) {
		renderer->AddActor(createLineActor(paths[i], pathColors[i]));
	}
	for (int i = 0; i < paths2.size(); ++i) {
		//renderer->AddActor(createLineActor(paths2[i], pathColors2[i]));
		renderer->AddActor(createTrajectoryActor(td2, i));
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
