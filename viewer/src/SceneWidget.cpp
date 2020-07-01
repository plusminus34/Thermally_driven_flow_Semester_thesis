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
	int trajectory_set = 9;
	int colordisplay = 3;//0: RGB, 1: Temperature, 2: Pressure, 3: relative humidity 

	// potentially modified by trajectory set
	int line_limit = 100000;
	int line_increment = 1;// only draw every x-th line
	int maxtime = 999999;

	// constants for mapping variables to colors?
	float T_min = 263.15f; int r_T_min = 50; int g_T_min = 50; int b_T_min = 255;
	float T_0 = 273.15f;   int r_T_0 = 222;  int g_T_0 = 200; int b_T_0 = 222;
	float T_max = 293.15f; int r_T_max = 255; int g_T_max = 100; int b_T_max = 50;
	float P_min =  65000; int r_P_min = 180; int g_P_min = 255; int b_P_min = 210;
	float P_max = 105000; int r_P_max = 180; int g_P_max = 20; int b_P_max = 20;
	float RELHUM_min = 0; int r_RELHUM_min = 250; int g_RELHUM_min = 200; int b_RELHUM_min = 40;
	float RELHUM_max = 100; int r_RELHUM_max = 50; int g_RELHUM_max = 0; int b_RELHUM_max = 255;

	bool T_in_kelvin = true;
	if (!T_in_kelvin) { T_min -= T_0; T_max -= T_0; T_0 -= T_0; }


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

	int r[] = { 255,0,0, 255,255,0 };
	int g[] = { 0,255,0, 255,0,255 };
	int b[] = { 0,0,255, 0,255,255 };
	float thickness[] = { 1,2,3,4,5,6 };

	if (trajectory_set == 0) {
		//jepp, die sind gleich
		files.push_back("../../../outputs/trajectory_TEST_correctish.nc");
		files.push_back("../../../outputs/lag_TEST.4");
		thickness[0] = thickness[1] = 2;
	}else if(trajectory_set == 1) {
		//konstantes Feld, RK vs Euler
		// sind verschieden, siehe Westen: Meiner geht eher um Berge als durch
		files.push_back("../../../outputs/trajectory_Lcomparison_const.nc");
		//files.push_back("../../../outputs/trajectory_comparison_const.nc");
		files.push_back("../../../outputs/trajectory_Lcomparison_const.4");
		thickness[0] = thickness[1] = 2;
		line_limit = 13;
	}
	else if (trajectory_set == 2) {
		//Feld nicht konstant: sehr veschieden
		// diesmal ist rot mein altes Feld, blau das Lagrantoartige
		// Ein Problem von Lagranto: Es benutzt nur 2 Felder statt 6, das sie im Stundentakt statt Minutentakt liest
		line_limit = 10;
		files.push_back("../../../outputs/trajectory_comparison.nc");
		files.push_back("../../../outputs/trajectory_Lcomparison.4");
		files.push_back("../../../outputs/trajectory_Lcomparison.nc");
		thickness[0] = thickness[2] = 1.5f;
	}
	else if (trajectory_set == 3) {
		//Rueckwaerts tracen
		// nicht sicher, was mit dem roten los ist: Blau hat die gleichen Einstellungen und funktioniert
		files.push_back("../../../outputs/trajectory_backwardtest_TI.nc");
		files.push_back("../../../outputs/trajectory_backwardtest.4");
		files.push_back("../../../outputs/trajectory_backwardtest_TI1.nc");
		thickness[0] = 1; thickness[1] = thickness[2] = 2;
	}
	else if (trajectory_set == 4) {
		//Verschiedene Zeitschritte
		string base = "../../../outputs/trajectory_compare_";
		string end = ".nc";
		string timestep[] = { "dt60","dt120","dt300" };
		string integrator[] = {"rk_", "ie_"};//RungeKutta,IterativeEuler
		line_limit = 10;
		for (int i = 0; i < 6; ++i) {
			int aa = i / 3;
			int bb = i % 3;
			files.push_back(base + integrator[aa] + timestep[bb] + end);
			r[i] = 25 * aa * (2 + 4 * bb);
			g[i] = 25 * (1-aa) * (2 + 4 * bb);
			b[i] = 50*bb;
			thickness[i] = 1 + 0.5*bb;
		}
	}
	else if (trajectory_set == 5 || trajectory_set == 7) {
		// Ein paar Orte, an denen ich mal war und die vielleicht gute Beispiele sind
		if (trajectory_set == 5) {
			files.push_back("../../../outputs/trajectory_luzern_wohin.nc");
			files.push_back("../../../outputs/trajectory_luzern_woher.nc");
			files.push_back("../../../outputs/trajectory_thusis_wohin.nc");
			files.push_back("../../../outputs/trajectory_thusis_woher.nc");
			files.push_back("../../../outputs/trajectory_zermatt_wohin.nc");
			files.push_back("../../../outputs/trajectory_zermatt_woher.nc");
			line_limit = 1000;
			line_increment = 10;
		}
		else {
			files.push_back("../../../outputs/trajectory_wohin_a_623_2.nc");
			files.push_back("../../../outputs/trajectory_woher_a_623_2.nc");
			files.push_back("../../../outputs/trajectory_wohin_b_623_2.nc");
			files.push_back("../../../outputs/trajectory_woher_b_623_2.nc");
			files.push_back("../../../outputs/trajectory_wohin_c_623_2.nc");
			files.push_back("../../../outputs/trajectory_woher_c_623_2.nc");
			
			line_limit = 5000;
			line_increment = 97;
		}
		for (int i = 0; i < files.size(); ++i) {
			thickness[i] = thickness[0];
			r[i] = i > 3 ? 200 : 70;
			g[i] = i < 2 ? 200 : 70;
			b[i] = i == 2 || i==3 ? 200 : 70;
			if (i % 2) { r[i] += 55; b[i] += 55; g[i] += 55; }
			else { r[i] -= 70; b[i] -= 70; g[i] -= 70; }
		}
	}
	else if (trajectory_set == 6) {
		//Test
		files.push_back("../../../outputs/trajectory_test__lagrantolike.nc");
		files.push_back("../../../outputs/trajectory_test__lagranto_betterW.nc");
		files.push_back("../../../outputs/trajectory_test__mine.nc");
		files.push_back("../../../outputs/trajectory_test__mine.nc");

		line_limit = 25;
		line_increment = 5;
		r[0] = 0;   g[0] = 255; b[0] = 0;   thickness[0] = 2;
		r[1] = 0;   g[1] = 255; b[1] = 255; thickness[1] = 1;
		r[2] = 0;   g[2] = 0;   b[2] = 255; thickness[2] = 1.5;
		r[3] = 255; g[3] = 0;   b[3] = 0;   thickness[3] = 1;
	}
	else if (trajectory_set == 8) {
		//"numbers"
		files.push_back("../../../outputs/trajectory_numbers3_meins_3h_dt1.nc");
		files.push_back("../../../outputs/trajectory_numbers3_lagranto_3h_dt1.4");
		files.push_back("../../../outputs/trajectory_numbers3_lagrantolike_3h_dt1.nc");
		files.push_back("../../../outputs/trajectory_numbers3_rungekutta_3h_dt1.nc");
		files.push_back("../../../outputs/trajectory_numbers3_betterW_3h_dt1.nc");
		files.push_back("../../../outputs/trajectory_numbers3_4pillars_3h_dt1.nc");
		//files.push_back("../cmd/trajectory_testr.nc");

		for (int i = 0; i < 6; ++i)thickness[i] = 1;
		line_increment = 11;
		line_limit = 7000;
		maxtime = 36000;
	}
	else if (trajectory_set == 9) {
		//files.push_back("../../../outputs/trajectory_626_breit3_unten_wohin.nc");
		//files.push_back("../../../outputs/trajectory_626_breit3_unten_woher.nc");
		//files.push_back("../../../outputs/trajectory_626_breit3_oben_wohin.nc");
		//files.push_back("../../../outputs/trajectory_626_breit3_oben_woher.nc");
		r[0] = 150; g[0] = 0; b[0] = 100;
		r[1] = 220; g[1] = 100; b[1] = 0;
		r[2] = 100; g[2] = 155; b[2] = 0;
		r[3] = 0; g[3] = 155; b[3] = 100;
		//files.push_back("../../../outputs/trajectory_626_nord_woher.nc");
		//files.push_back("../../../outputs/trajectory_626_sued_wohin.nc");
		//files.push_back("../../../outputs/trajectory_626_taltest.4");
		files.push_back("../../../outputs/trajectory_626_taltest_wohin.nc");
		//files.push_back("../../../outputs/trajectory_626_taltest_wohin_lagrantolike.nc");
		files.push_back("../../../outputs/trajectory_626_taltest_woher.nc");

		for (int i = 0; i < 6; ++i)thickness[i] = 1;//2 - (i/2);
		line_increment = 1;
	}

	
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());

	renderer->AddActor(landscapeActor);

	std::vector<TrajectoryData> td(files.size());
	for (int file_i = 0; file_i < files.size(); ++file_i) {
		NetCDF::ReadTrajectoryData(files[file_i], td[file_i]);
		assert(td.num_trajectories > 0);
		assert(td.points_per_trajectory > 0);

		int T_id = td[file_i].get_var_id("T");
		assert(T_id > -1);
		if (td[file_i].min_values[T_id] < T_min) T_min = td[file_i].min_values[T_id];
		else if (td[file_i].max_values[T_id] > T_max) T_max = td[file_i].max_values[T_id];

		int P_id = td[file_i].get_var_id("P");
		assert(p_id > -1);
		if (td[file_i].min_values[P_id] < P_min) P_min = td[file_i].min_values[P_id];
		else if (td[file_i].max_values[P_id] > P_max) P_max = td[file_i].max_values[P_id];

		int hum_id = td[file_i].get_var_id("RELHUM");
		assert(hum_id > -1);
		if (td[file_i].min_values[hum_id] < RELHUM_min) RELHUM_min = td[file_i].min_values[hum_id];
		else if (td[file_i].max_values[hum_id] > RELHUM_max) RELHUM_max = td[file_i].max_values[hum_id];
	}

	// Create trajectory actor
	for (int file_i = 0; file_i < files.size(); ++file_i) {
		vtkNew<vtkPoints> points;

		int x_id, y_id, z_id;
		if (use_rotated) {
			x_id = td[file_i].get_var_id("rlon");
			y_id = td[file_i].get_var_id("rlat");
		}
		else {
			x_id = td[file_i].get_var_id("lon");
			y_id = td[file_i].get_var_id("lat");
		}
		z_id = td[file_i].get_var_id("z");
		assert(x_id > -1);
		assert(y_id > -1);
		assert(z_id > -1);

		int T_id = td[file_i].get_var_id("T"); assert(T_id > -1);
		int P_id = td[file_i].get_var_id("P");
		int hum_id = td[file_i].get_var_id("RELHUM");

		vtkNew<vtkUnsignedCharArray> colors;
		colors->SetName("color");
		colors->SetNumberOfComponents(3);
		colors->SetNumberOfTuples(td[file_i].num_trajectories * td[file_i].points_per_trajectory);

		for (int i = 0; i < td[file_i].num_trajectories; ++i) {
			for (int j = 0; j < td[file_i].points_per_trajectory; ++j) {
				points->InsertNextPoint(td[file_i].val(x_id, i, j), td[file_i].val(y_id, i, j), td[file_i].val(z_id, i, j) * ZSCALE);

				const int ij = i * td[file_i].points_per_trajectory + j;
				if (colordisplay == 0) {
					//color
					if (td[file_i].times[j] < maxtime)
					colors->InsertTuple3(ij, r[file_i], g[file_i], b[file_i]);
					else colors->InsertTuple3(ij, 111, 123, 135);
				}
				else if (colordisplay == 1) {
					// temperature
					const float T = td[file_i].val(T_id, i, j);
					if (T < T_min) {
						colors->InsertTuple3(ij, r_T_min, g_T_min, b_T_min);
					}
					else if (T < T_0) {
						const float alpha = (T - T_min) / (T_0 - T_min);
						const int r = r_T_min + alpha * (r_T_0 - r_T_min);
						const int g = g_T_min + alpha * (g_T_0 - g_T_min);
						const int b = b_T_min + alpha * (b_T_0 - b_T_min);
						colors->InsertTuple3(ij, r, g, b);
					}
					else if (T < T_max) {
						const float alpha = (T - T_0) / (T_max - T_0);
						const int r = r_T_0 + alpha * (r_T_max - r_T_0);
						const int g = g_T_0 + alpha * (g_T_max - g_T_0);
						const int b = b_T_0 + alpha * (b_T_max - b_T_0);
						colors->InsertTuple3(ij, r, g, b);
					}
					else {
						colors->InsertTuple3(ij, r_T_max, g_T_max, b_T_max);
					}
				}
				else if (colordisplay == 2) {
					// pressure
					float P = td[file_i].val(P_id, i, j);
					if (files[file_i][files[file_i].size() - 1] == '4') P *= 100;
					if (P < P_min) {
						colors->InsertTuple3(ij, r_P_min, g_P_min, b_P_min);
					}
					else if (P < P_max) {
						const float alpha = (P - P_min) / (P_max - P_min);
						const int r = r_P_min + alpha * (r_P_max - r_P_min);
						const int g = g_P_min + alpha * (g_P_max - g_P_min);
						const int b = b_P_min + alpha * (b_P_max - b_P_min);
						colors->InsertTuple3(ij, r, g, b);
					}
					else {
						colors->InsertTuple3(ij, r_P_max, g_P_max, b_P_max);
					}
				}
				else if (colordisplay == 3) {
					// relative humidity
					const float hum = td[file_i].val(hum_id, i, j);
					if (hum < RELHUM_min) {
						colors->InsertTuple3(ij, r_RELHUM_min, g_RELHUM_min, b_RELHUM_min);
					}
					else if (hum < RELHUM_max) {
						const float alpha = (hum - RELHUM_min) / (RELHUM_max - RELHUM_min);
						const int r = r_RELHUM_min + alpha * (r_RELHUM_max - r_RELHUM_min);
						const int g = g_RELHUM_min + alpha * (g_RELHUM_max - g_RELHUM_min);
						const int b = b_RELHUM_min + alpha * (b_RELHUM_max - b_RELHUM_min);
						colors->InsertTuple3(ij, r, g, b);
					}
					else {
						colors->InsertTuple3(ij, r_RELHUM_max, g_RELHUM_max, b_RELHUM_max);
					}
				}
			}
		}
		//size_t nPts = points->GetNumberOfPoints();

		vtkNew<vtkCellArray> cells;
		for (int i = 0; i < line_limit*line_increment && i < td[file_i].num_trajectories; i+=line_increment) {
			vtkNew<vtkPolyLine> polyLine;
			polyLine->GetPointIds()->SetNumberOfIds(td[file_i].points_per_trajectory);
			for (unsigned int j = 0; j < td[file_i].points_per_trajectory; j++)
				polyLine->GetPointIds()->SetId(j, j + i * td[file_i].points_per_trajectory);

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
		mapper->SelectColorArray("color");

		vtkNew<vtkActor> actor;
		actor->SetMapper(mapper);
		actor->GetProperty()->SetLineWidth(thickness[file_i]);

		renderer->AddActor(actor);
	}

	GetRenderWindow()->AddRenderer(renderer);
	GetRenderWindow()->SetWindowName("RenderWindowNoUIFile");
}
