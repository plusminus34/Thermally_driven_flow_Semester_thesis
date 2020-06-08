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

SceneWidget::SceneWidget()
{
	mRenderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
	SetRenderWindow(mRenderWindow);

	CreateTestScene();
}

void SceneWidget::CreateTestScene()
{
	// Settings

	bool build_landscape = false;// TODO check if file exists automatically
	bool use_rotated = true;//use rlonrlat (if false, lonlat are used)

	// display surface
	vtkNew<vtkPolyDataMapper> landscapeMapper;
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
			if (use_rotated) {
				hfield_points->InsertNextPoint(coord[0], coord[1], hsurf->GetVertexDataAt(gridCoord) * ZSCALE);
			}
			else {
				double lon, lat;
				CoordinateTransform::RlatRlonToLatLon(coord[1], coord[0], lat, lon);
				hfield_points->InsertNextPoint(lon, lat, hsurf->GetVertexDataAt(gridCoord) * ZSCALE);
			}
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
		if (use_rotated) landscapeWriter->SetFileName("landscape_reduced_r.vtp");
		else landscapeWriter->SetFileName("landscape_reduced.vtp");
		
		landscapeWriter->SetInputData(landscape);
		landscapeWriter->Write();
	}
	else {
		vtkNew<vtkXMLPolyDataReader> landscapeReader;
		if(use_rotated) landscapeReader->SetFileName("landscape_reduced_r.vtp");
		else landscapeReader->SetFileName("landscape_reduced.vtp");
		landscapeReader->Update();
		landscapeMapper->SetInputConnection(landscapeReader->GetOutputPort());
	}
	vtkNew<vtkActor> landscapeActor;
	landscapeActor->SetMapper(landscapeMapper);
	
	//string file_1 = "../../../outputs/lagranto_demo_605.nc";
	string file = "../../../outputs/trajectory_demo_605.nc";
	//string file_1 = "../../../outputs/testoutput.4";
	//string file_2 = "../../../outputs/trajectory_TEST.nc";
	
	TrajectoryData td;

	NetCDF::ReadTrajectoryData(file, td);
	assert(td.num_trajectories > 0);
	assert(td.points_per_trajectory > 0);

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

	// Create trajectory actor
	if(true){
		assert(td.num_trajectories > 0);
		assert(td.points_per_trajectory > 0);
		vtkNew<vtkPoints> points;

		int x_id, y_id, z_id;
		if (use_rotated) {
			x_id = td.get_var_id("rlon");
			y_id = td.get_var_id("rlat");
		}
		else {
			x_id = td.get_var_id("lon");
			y_id = td.get_var_id("lat");
		}
		z_id = td.get_var_id("z");
		assert(x_id > -1);
		assert(y_id > -1);
		assert(z_id > -1);

		int T_id = td.get_var_id("T"); assert(T_id > -1);
		int P_id = td.get_var_id("P");
		int hum_id = td.get_var_id("RELHUM");

		vtkNew<vtkUnsignedCharArray> colors;
		colors->SetName("T_color");
		colors->SetNumberOfComponents(3);
		colors->SetNumberOfTuples(td.num_trajectories*td.points_per_trajectory);

		for (int i = 0; i < td.num_trajectories; ++i) {
			for (int j = 0; j < td.points_per_trajectory; ++j) {
				//points->InsertNextPoint(i, j, i*j);
				points->InsertNextPoint(td.val(x_id, i, j), td.val(y_id, i, j), td.val(z_id, i, j) * ZSCALE);

				float sat = (td.val(T_id, i, j) - td.min_values[T_id]) / (td.max_values[T_id] - td.min_values[T_id]);
				//float sat = (float)j / td.points_per_trajectory;
				//float sat = (td.val(T_id, i, j) - 200) / (td.max_values[T_id] - 200);
				//float sat = (td.val(hum_id, i, j) - td.min_values[hum_id]) / (td.max_values[hum_id] - td.min_values[hum_id]);
				int ij = i * td.points_per_trajectory + j;
				//colors->InsertTuple3(ij, (int)(100 + 155 * sat), (int)(150 - 140 * sat), (int)(255 - 255 * sat));
				colors->InsertTuple3(ij, 1, 255*j/td.points_per_trajectory, 255*i/td.num_trajectories);
				//colors->InsertTuple3(ij, td.val(T_id, i, j)-100, 20 * j / td.points_per_trajectory, 25 * i / td.num_trajectories);
			}
		}
		//size_t nPts = points->GetNumberOfPoints();

		vtkNew<vtkCellArray> cells;
		for (int i = 0; i < td.num_trajectories; ++i) {
			vtkNew<vtkPolyLine> polyLine;
			polyLine->GetPointIds()->SetNumberOfIds(td.points_per_trajectory);
			for (unsigned int j = 0; j < td.points_per_trajectory; j++)
				polyLine->GetPointIds()->SetId(j, j + i * td.points_per_trajectory);

			cells->InsertNextCell(polyLine);
		}

		vtkNew<vtkPolyData> polyData;
		polyData->SetPoints(points);	// Add the points to the dataset
		polyData->SetLines(cells);			// Add the lines to the dataset

		polyData->GetPointData()->AddArray(colors);

		vtkNew<vtkPolyDataMapper> mapper;
		mapper->SetInputData(polyData);
		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointFieldData();
		mapper->SelectColorArray("T_color");

		vtkNew<vtkActor> actor;
		actor->SetMapper(mapper);
		actor->GetProperty()->SetLineWidth(3);

		renderer->AddActor(actor);
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
