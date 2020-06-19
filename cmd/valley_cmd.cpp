#include "core/NetCDF.hpp"

#include <iostream>
#include <fstream>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include "ImportantPart.hpp"
#include "core/LagrantoUVW.hpp"

// ---------------------------------------
// Entry point
// ---------------------------------------
int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "No file to process!" << std::endl;
		return 0;
	}
	int input = 0;
	bool custom = false;
	std::cout << "Choose option:\n";
	std::cout << "0:\tRead and interpolate UVW field\n";
	std::cout << "1:\tRead custom field\n";
	std::cout << "2:\tRead UVW and trace particles\n";
	std::cout << "3:\tCompare trajectories\n";
	std::cout << "4:\tDebug output\n> ";
	cin >> input;
	if (input == 4) {
		return 0;
	}
	else if (input == 3) {
		string file_a = "trajectory_test_longtime.nc";
		string file_b = "trajectory_test_longtime2.nc";

		TrajectoryData td_a, td_b;
		NetCDF::ReadTrajectoryData(file_a, td_a);
		NetCDF::ReadTrajectoryData(file_b, td_b);

		size_t n_tra = td_a.num_trajectories;
		size_t n_tim = td_b.points_per_trajectory;

		if (n_tra != td_b.num_trajectories) {
			cout << "Error: Different number of trajectories\n";
			return 1;
		}
		if (n_tim != td_b.points_per_trajectory) {
			cout << "Error: Different number of timesteps\n";
			return 1;
		}

		cout << "timessize " << td_a.times.size()<<endl;

		int lon_id = td_a.get_var_id("lon");
		int lat_id = td_a.get_var_id("lat");
		int rlon_id = td_a.get_var_id("rlon");
		int rlat_id = td_a.get_var_id("rlat");
		int z_id = td_a.get_var_id("z");
		//TODO check that a and b have the same ids

		for (int i = 0; i < n_tra; ++i) {
			cout << "Comparing trajectory " << i << endl;
			double sum_horiz = 0;
			double sum_vert = 0;
			double int_horiz = 0;
			double int_vert = 0;

			for (int j = 0; j < n_tim; ++j) {
				float dx = td_a.val(rlon_id, i, j) - td_b.val(rlon_id, i, j);
				float dy = td_a.val(rlat_id, i, j) - td_b.val(rlat_id, i, j);
				float dz = td_a.val(z_id, i, j) - td_b.val(z_id, i, j);
				sum_horiz += sqrt(dx*dx + dy * dy);
				sum_vert += abs(dz);
				if (j > 0) {
					double dt = td_a.times[j] - td_a.times[j-1];
					int_horiz += sum_horiz*dt;
					int_vert += sum_vert*dt;
				}
				if (j < n_tim - 1) {
					double dt = td_a.times[j + 1] - td_a.times[j];
					int_horiz += sum_horiz * dt;
					int_vert += sum_vert * dt;
				}
			}
			int_horiz *= 0.5;
			int_vert *= 0.5;

			cout << "  Sum of horizontal distance differences: " << sum_horiz << endl;
			cout << "  Sum of vertical distance differences: " << sum_vert << endl;
			cout << "  Integrated horizontal distance differences: " << int_horiz << endl;
			cout << "  Integrated vertical distance differences: " << int_vert << endl;
		}

		return 0;
	}
	else if (input == 2) {

		string path(argv[1]);
		ImportantPart imp;
		TrajectoryData td;

		// Get user input
		NetCDF::Info info;
		if (!NetCDF::ReadInfo(path, info)) { return -1; }

		double tracing_bounds[6];
		double spacing[3];
		int paths_dim[3];
		int nPaths;
		double t0, t1, dt;
		int nSteps;
		vector<int> extra_variables(0);
		bool confirmed = false;
		while (!confirmed) {
			cout << "Input bounds for tracing in (lon,lat,h) (Format: min_x max_x min_y max_y min_z max_z)\n> ";
			for (int i = 0; i < 6; ++i) {
				cin >> tracing_bounds[i];
				/*
				if (i % 2 == 0) tracing_bounds[i] = std::max(tracing_bounds[i], bounds[i]);
				else {
					tracing_bounds[i] = std::min(tracing_bounds[i], bounds[i]);
					if (tracing_bounds[i - 1] > tracing_bounds[i]) swap(tracing_bounds[i - 1], tracing_bounds[i]);
				}
				*/
			}
			cout << "Input spacing between individual starting points  (Format: dx dy dz)> ";
			for (int i = 0; i < 3; ++i) {
				cin >> spacing[i];
				paths_dim[i] = floor((tracing_bounds[2 * i + 1] - tracing_bounds[2 * i]) / spacing[i]) + 1;
				assert(paths_dim[i] > 0);
			}
			nPaths = paths_dim[0] * paths_dim[1] * paths_dim[2];
			cout << "Input start and end time (Format: t0 t1)> "; cin >> t0; cin >> t1;
			cout << "Input timestep> "; cin >> dt;
			if (t1 < t0 && dt < 0) {
				cout << "Sorry, backtracing has not been implemented yet.\n";
				//continue;
			}
			else if (t1 > t0 && dt > 0) {
				//standard case
			}
			else {
				cout << "This doesn't work, start over please.\n";
				continue;
			}
			nSteps = ceil((t1 - t0) / dt);

			int more = 0;
			cout << "How many variables beyond position should be tracked?> ";
			cin >> more;
			if (more > 0) {
				extra_variables.resize(more);
				cout << "Possible variables include:\n";
				for (int i = 0; i < info.NumVariables; ++i) {
					try
					{
						NetCDF::Info::Dimension dim0 = info.Variables[i].GetDimensionByName("rlon");
						NetCDF::Info::Dimension dim1 = info.Variables[i].GetDimensionByName("rlat");
						NetCDF::Info::Dimension dim2 = info.Variables[i].GetDimensionByName("level");

						cout << "  " << i << ": " << info.Variables[i].GetName() << endl;
					}
					catch (const string& str) {}
				}
				cout << "Choose which " << more << " to use> ";
				for (int i = 0; i < more; ++i) {
					cin >> extra_variables[i];
				}
			}
			more = 0;
			cout << "Change integrator settings? (0/1)> ";
			cin >> more;
			if (more > 0) {
				cout << "Jump back into domain if leaving? (0/1)> ";
				cin >> input;
				imp.setJumpIntoDomain((bool)input);
				cout << "Use Lagranto-style UVW? (0/1)> ";
				cin >> input;
				if (input == 1) {
					imp.setUseLagrantoUVW(true);
					cout << "Use level interpolation on four columns? (0/1)> ";
					cin >> input;
					imp.setZInterpolationFlag((bool)input);
					if (input == 0) {
						cout << "Use staggered HHL for unstaggered W? (0/1)> ";
						cin >> input;
						imp.setUseWrongW((bool)input);
					}
				}
				cout << "Use which integrator?\n 0) Runge-Kutta 4\n 1) Iterative Euler\n> ";
				cin >> input;
				if (input == 1) imp.setIntegratorToIterativeEuler();
				cout << "Use time-invariant UVW? (0/1)> ";
				cin >> input;
				if (input == 1) imp.setTimeInvariant(true);
			}

			cout << "This will trace a total of " << nPaths << " trajectories (" << paths_dim[0] << " x " << paths_dim[1] << " x " << paths_dim[2] << ")\n";
			cout << "There will be up to " << nSteps << " steps per trajectory, storing "<<extra_variables.size()<<" extra variables\n";

			cout << "Confirm (0/1)> "; cin >> confirmed;
		}
		string name = "";
		cout << "Enter a name for points/trajectory> ";
		cin >> name;

		// Prepare tracing
		td.num_trajectories = nPaths;
		td.points_per_trajectory = nSteps + 1;
		td.time_begin = t0;
		td.time_end = t1;
		td.varnames = { "rlon", "rlat", "z", "lon", "lat" };
		for (int i = 0; i < extra_variables.size(); ++i) {
			td.varnames.push_back(info.Variables[extra_variables[i]].GetName());
		}

		string basefile = path;
		basefile = basefile.substr(0, basefile.size() - 11);
		imp.setBaseFileName(basefile);

		imp.setTimeBoundaries(t0, t1);
		imp.setTimestep(dt);

		imp.trajectories.resize(nPaths);
		cout << "Printing initial points\n";

		// Write initial points into file
		ofstream points_file;
		points_file.open("points_" + name + ".txt");
		//points_file.open(name + "_points" + ".txt");

		string hhmm = "";
		int h = td.time_begin / 3600;
		int m = (td.time_begin - 3600 * h) / 60;
		if (h < 10) hhmm.append("0");
		hhmm.append(to_string(h) + ".");
		if (m < 10) hhmm.append("0");
		hhmm.append(to_string(m));
		for (int i = 0; i < paths_dim[0]; ++i) {
			for (int j = 0; j < paths_dim[1]; ++j) {
				for (int k = 0; k < paths_dim[2]; ++k) {
					int path = i * paths_dim[1] * paths_dim[2] + j * paths_dim[2] + k;
					Vec3f position = Vec3f(tracing_bounds[0] + i * spacing[0], tracing_bounds[2] + j * spacing[1], tracing_bounds[4] + k * spacing[2]);
					points_file << hhmm << " " << position[0] << " " << position[1] << " " << position[2] << endl;
					double lat = position[1];
					double lon = position[0];
					double rlat, rlon;
					CoordinateTransform::LatLonToRlanRlon(lat, lon, rlat, rlon);
					position[0] = rlon; position[1] = rlat;
					imp.trajectories[path].resize(1);
					imp.trajectories[path][0] = position;
				}
			}
		}
		points_file.close();

		// Do the important part
		imp.computeTrajectoryData(td);

		string trajectory_file = "trajectory_" + name + ".nc";
		cout << "Writing trajectories to " << trajectory_file << "...\n";
		NetCDF::WriteTrajectoryData(trajectory_file, td);
		cout << "Finished" << endl;
	}
	else if (input < 2) {
		std::string varname;
		std::string dimname_0, dimname_1, dimname_2;
		custom = input;
		if (!custom) {
			varname = "U";
			dimname_0 = "srlon";
			dimname_1 = "rlat";
			dimname_2 = "level";
		}
		std::string path(argv[1]);

		std::cout << "Reading info object..." << std::endl;
		// read info object
		NetCDF::Info info;
		if (!NetCDF::ReadInfo(path, info)) { return -1; }

		if (custom) {
			int confirm = 0;
			while (confirm == 0) {
				for (int i = 0; i < info.NumDimensions; ++i) {
					std::cout << "dimension " << i << ": " << info.Dimensions[i].GetName() << std::endl;
				}
				for (int i = 0; i < info.NumVariables; ++i) {
					std::cout << "Variable " << i << ": " << info.Variables[i].GetName() << std::endl;
				}
				int input = 0;
				std::cout << "Selected variable is "; std::cin >> input;
				varname = info.Variables[input].GetName();
				const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);
				for (int i = 0; i < variable.Dimensions.size(); ++i) {
					std::cout << varname << " dimension " << i << ": " << variable.Dimensions[i].GetName() << std::endl;
				}
				for (int i = 0; i < variable.Attributes.size(); ++i) {
					std::cout << varname << " attribute " << i << ": " << variable.Attributes[i].GetName() << " : " << variable.Attributes[i].GetValue() << std::endl;
				}
				std::cout << "X dimension is "; std::cin >> input;
				dimname_0 = variable.Dimensions[input].GetName();
				std::cout << "Y dimension is "; std::cin >> input;
				dimname_1 = variable.Dimensions[input].GetName();
				std::cout << "Z dimension is "; std::cin >> input;
				dimname_2 = variable.Dimensions[input].GetName();

				cout << "Confirm? (0/1)> "; cin >> confirm;
			}
		}
		std::cout << "Importing field..." << std::endl;
		// import field
		RegScalarField3f* field = NetCDF::ImportScalarField3f(path, varname, dimname_0, dimname_1, dimname_2);
		std::cout << "Checking values afterwards...\n";
		int64_t field_size = (int64_t)field->GetResolution()[0] * field->GetResolution()[1] * field->GetResolution()[2];
		int64_t checkpoint = field_size / 50;

		if (!custom) {
			RegVectorField3f* vecField = UVWFromNCFile(path);
			int64_t vecField_size = (int64_t)vecField->GetResolution()[0] * vecField->GetResolution()[1] * vecField->GetResolution()[2];

			int input = 0;
			std::cout << "write into vti? (0/1)> "; cin >> input;
			if (input) {
				delete field;


				vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
				imageData->SetDimensions(vecField->GetResolution().ptr());
				double spacing[] = {
					(vecField->GetDomain().GetMax()[0] - vecField->GetDomain().GetMin()[0]) / (vecField->GetResolution()[0] - 1.),
					(vecField->GetDomain().GetMax()[1] - vecField->GetDomain().GetMin()[1]) / (vecField->GetResolution()[1] - 1.),
					(vecField->GetDomain().GetMax()[2] - vecField->GetDomain().GetMin()[2]) / (vecField->GetResolution()[2] - 1.),
				};
				imageData->SetSpacing(spacing);
				imageData->SetNumberOfScalarComponents(3, imageData->GetInformation());

				std::cout << "Filling image data..." << std::endl;
				// fill image data
				vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
				floatArray->SetNumberOfComponents(3);
				int numTuples = vecField_size;
				int tuple_checkpoint = numTuples / 50;
				floatArray->SetNumberOfTuples(numTuples);
				floatArray->SetName(varname.c_str());
				for (int linearIndex = 0; linearIndex < numTuples; ++linearIndex)
				{
					Vec3i gridCoord = vecField->GetGridCoord(linearIndex);
					Vec3d dataVec = vecField->GetVertexDataAt(gridCoord);
					floatArray->SetTuple3(linearIndex, dataVec[0], dataVec[1], dataVec[2]);
					if (linearIndex % tuple_checkpoint == 0) { std::cout << linearIndex << " of " << numTuples << std::endl; }
				}
				imageData->GetPointData()->AddArray(floatArray);

				double origin[3] = { vecField->GetDomain().GetMin()[0], vecField->GetDomain().GetMin()[1], vecField->GetDomain().GetMin()[2] };
				imageData->SetOrigin(origin);

				delete vecField;
				std::cout << "Writing to file..." << std::endl;
				// write to file
				vtkSmartPointer<vtkXMLImageDataWriter> imageWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
				std::string filename = "UVW.vti";
				imageWriter->SetFileName(filename.c_str());
				imageWriter->SetInputData(imageData);
				imageWriter->Update();

				std::cout << "Done." << std::endl;
				// delete resources and return
				return 0;
			}
			else {
				delete vecField;
			}
		}
		std::cout << "Allocating image data object..." << std::endl;
		// allocate image data object
		vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
		imageData->SetDimensions(field->GetResolution().ptr());
		double spacing[] = {
			(field->GetDomain().GetMax()[0] - field->GetDomain().GetMin()[0]) / (field->GetResolution()[0] - 1.),
			(field->GetDomain().GetMax()[1] - field->GetDomain().GetMin()[1]) / (field->GetResolution()[1] - 1.),
			(field->GetDomain().GetMax()[2] - field->GetDomain().GetMin()[2]) / (field->GetResolution()[2] - 1.),
		};
		imageData->SetSpacing(spacing);

		std::cout << "Filling image data..." << std::endl;
		// fill image data
		vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
		floatArray->SetNumberOfComponents(1);
		int numTuples = field->GetResolution()[0] * field->GetResolution()[1] * field->GetResolution()[2];
		int tuple_checkpoint = numTuples / 50;
		floatArray->SetNumberOfTuples(numTuples);
		floatArray->SetName(varname.c_str());
		for (int linearIndex = 0; linearIndex < numTuples; ++linearIndex)
		{
			Vec3i gridCoord = field->GetGridCoord(linearIndex);
			float value = field->GetVertexDataAt(gridCoord);
			floatArray->SetValue(linearIndex, value);
			if (linearIndex % tuple_checkpoint == 0) { std::cout << linearIndex << " of " << numTuples << std::endl; }
		}
		imageData->GetPointData()->AddArray(floatArray);
		double origin[3] = { field->GetDomain().GetMin()[0], field->GetDomain().GetMin()[1], field->GetDomain().GetMin()[2] };
		imageData->SetOrigin(origin);

		std::cout << "Writing to file..." << std::endl;
		// write to file
		vtkSmartPointer<vtkXMLImageDataWriter> imageWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
		std::string filename = varname + "_" + dimname_0 + "_" + dimname_1 + "_" + dimname_2 + ".vti";
		imageWriter->SetFileName(filename.c_str());
		imageWriter->SetInputData(imageData);
		imageWriter->Update();

		std::cout << "Done." << std::endl;
		// delete resources and return
		delete field;
	}
	return 0;
}