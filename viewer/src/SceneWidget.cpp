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
	bool use_rotated = false;//use rlonrlat (if false, lonlat are used)

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

	vector<string> files(0);
	files.push_back("../../../outputs/trajectory_TEST_correctish.nc");
	files.push_back("../../../outputs/lag_TEST.4");
	//files.push_back("../../../outputs/trajectory_demo_605.nc");
	//files.push_back("../../../outputs/testoutput.4");
	
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(landscapeActor);

	// Create trajectory actor
	for (int file_i = 0; file_i < files.size(); ++file_i) {
		TrajectoryData td;
		NetCDF::ReadTrajectoryData(files[file_i], td);
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

		//int T_id = td.get_var_id("T"); assert(T_id > -1);
		//int P_id = td.get_var_id("P");
		//int hum_id = td.get_var_id("RELHUM");

		vtkNew<vtkUnsignedCharArray> colors;
		colors->SetName("T_color");
		colors->SetNumberOfComponents(3);
		colors->SetNumberOfTuples(td.num_trajectories*td.points_per_trajectory);

		for (int i = 0; i < td.num_trajectories; ++i) {
			for (int j = 0; j < td.points_per_trajectory; ++j) {
				//points->InsertNextPoint(i, j, i*j);
				points->InsertNextPoint(td.val(x_id, i, j), td.val(y_id, i, j), td.val(z_id, i, j) * ZSCALE);
				//points->InsertNextPoint(i,j*0.01,file_i);

				//float sat = (td.val(T_id, i, j) - td.min_values[T_id]) / (td.max_values[T_id] - td.min_values[T_id]);
				//float sat = (float)j / td.points_per_trajectory;
				//float sat = (td.val(T_id, i, j) - 200) / (td.max_values[T_id] - 200);
				//float sat = (td.val(hum_id, i, j) - td.min_values[hum_id]) / (td.max_values[hum_id] - td.min_values[hum_id]);
				int ij = i * td.points_per_trajectory + j;
				//colors->InsertTuple3(ij, (int)(100 + 155 * sat), (int)(150 - 140 * sat), (int)(255 - 255 * sat));
				colors->InsertTuple3(ij, 100*file_i, 255*j/td.points_per_trajectory, 255*i/td.num_trajectories);
				//colors->InsertTuple3(ij, 100 * file_i, 240-100*file_i, 20);
				//colors->InsertTuple3(ij, td.val(T_id, i, j)-100, 20 * j / td.points_per_trajectory, 25 * i / td.num_trajectories);
				//colors->InsertTuple3(ij, 200*file_i, 255 - 250*file_i, 255);
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
