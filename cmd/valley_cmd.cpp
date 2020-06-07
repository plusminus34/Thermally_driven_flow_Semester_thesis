#include "core/NetCDF.hpp"

#include <iostream>
#include <fstream>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include "ImportantPart.hpp"

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
	std::cout << "3:\tDebug output\n> ";
	cin >> input;
	if (input == 3) {
		std::string path(argv[1]);
		RegScalarField3f* U = nullptr;
		RegScalarField3f* V = nullptr;
		RegScalarField3f* W = nullptr;
		SeparateUVWFromNCFile(path, U, V, W);
		assert(U != nullptr);
		assert(U);
		Vec3i gc = U->GetGridCoord(67);
		cout << "ggg " << gc[0] << " " << gc[1] << " " << gc[2] << endl;
		Vec3i dim = U->GetResolution();
		cout << "U dim: " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
		dim = V->GetResolution();
		cout << "V dim: " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
		dim = W->GetResolution();
		cout << "W dim: " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
		Vec3i wo(111, 111, 33);
		cout << "at ijk " << wo[0] << " " << wo[1] << " " << wo[2] << endl;
		cout << "  U " << U->GetVertexDataAt(wo)<<endl;
		cout << "  V " << V->GetVertexDataAt(wo) << endl;
		cout << "  W " << W->GetVertexDataAt(wo) << endl;
		Vec3d wa(-2.0837, -0.960101426, 4000);
		Vec3f trg(15.15658, 13.92212799, -2.81133294);
		//wa = U->GetCoordAt(wo);
		for (double i = -2; i < 2; i+=10.001) {
			for (double j = -2; j < 2; j+=10.01) {
				for (double k = -2; k < 2; k+=100.01) {
					Vec3d ww = wa + Vec3d(i, j, k);
					Vec3f uvw(U->Sample(ww), V->Sample(ww), W->Sample(ww));
					if ((uvw - trg).length() < 0.1) {
						cout << "pos " << ww[0] << " " << ww[1] << " " << ww[2] << endl;
						cout << "  uvw   " << uvw[0] << " " << uvw[1] << " " << uvw[2] << endl;
					}
				}
			}
		}
		cout << "at pos " << wa[0] << " " << wa[1] << " " << wa[2] << endl;
		cout << "  U " << U->Sample(wa) << endl;
		cout << "  V " << V->Sample(wa) << endl;
		cout << "  W " << W->Sample(wa) << endl;
		string constantsfile = path;
		for (int i = 0; i < constantsfile.size(); ++i) {
			if (constantsfile[i] == '0' && constantsfile[i + 1] == '.')
			{
				constantsfile.append("c");
				constantsfile[i + 1] = 'c';
				constantsfile[i + 2] = '.';
				constantsfile[i + 3] = 'n';
				break;
			}
		}
		cout << "Reading hhl from " << constantsfile << endl;
		RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");


		//TODO  test-integrator similar to input 2


		Vec3d coord(-2.08377838, -0.960101426, 4000);
		//Vec3d coord(-2.04549885, 0.0391660146, 4000);
		//Vec3d coord(0, 0, 4000);

		ImportantPart imp;
		int nPaths = 1;

		TrajectoryData td;
		td.num_trajectories = nPaths;
		td.points_per_trajectory = 61;
		td.time_begin = 0;
		td.time_end = 3600;
		td.varnames = { "rlon", "rlat", "z", "lon", "lat" };

		string basefile = path;
		basefile = basefile.substr(0, basefile.size() - 11);
		imp.setBaseFileName(basefile);

		imp.setTimeBoundaries(td.time_begin, td.time_end);
		imp.setTimestep(60);

		imp.trajectories.resize(nPaths);
		for (int path = 0; path < nPaths; ++path) {
			imp.trajectories[path].resize(1);
			imp.trajectories[path][0] = coord;
		}

		// Do the important part
		imp.computeTrajectoryDataTEST(td, U, V, W);



		return 0;// The end for the working part

		cout << "Domains\n";
		cout << "U min x " << U->GetDomain().GetMin()[0] << endl;
		cout << "V min x " << V->GetDomain().GetMin()[0] << endl;
		cout << "W min x " << W->GetDomain().GetMin()[0] << endl;
		cout << "hhlmn x " << hhl->GetDomain().GetMin()[0] << endl;

		cout << "U max x " << U->GetDomain().GetMax()[0] << endl;
		cout << "V max x " << V->GetDomain().GetMax()[0] << endl;
		cout << "W max x " << W->GetDomain().GetMax()[0] << endl;
		cout << "h max x " << hhl->GetDomain().GetMax()[0] << endl;

		cout << "U min y " << U->GetDomain().GetMin()[1] << endl;
		cout << "V min y " << V->GetDomain().GetMin()[1] << endl;
		cout << "W min y " << W->GetDomain().GetMin()[1] << endl;
		cout << "h min y " << hhl->GetDomain().GetMin()[1] << endl;

		cout << "U max y " << U->GetDomain().GetMax()[1] << endl;
		cout << "V max y " << V->GetDomain().GetMax()[1] << endl;
		cout << "W max y " << W->GetDomain().GetMax()[1] << endl;
		cout << "h max y " << hhl->GetDomain().GetMax()[1] << endl;

		//cout << "U min z " << U->GetDomain().GetMin()[2] << endl;
		cout << "U max z " << U->GetDomain().GetMax()[2] << endl;
		//cout << "V min z " << V->GetDomain().GetMin()[2] << endl;
		cout << "V max z " << V->GetDomain().GetMax()[2] << endl;
		//cout << "W min z " << W->GetDomain().GetMin()[2] << endl;
		cout << "W max z " << W->GetDomain().GetMax()[2] << endl;
		//cout << "h min z " << hhl->GetDomain().GetMin()[2] << endl;
		cout << "h max z " << hhl->GetDomain().GetMax()[2] << endl;

		cout << "U res: " << U->GetResolution()[0] << " " << U->GetResolution()[1] << " " << U->GetResolution()[2] << endl;
		cout << "v res: " << V->GetResolution()[0] << " " << V->GetResolution()[1] << " " << V->GetResolution()[2] << endl;
		cout << "w res: " << W->GetResolution()[0] << " " << W->GetResolution()[1] << " " << W->GetResolution()[2] << endl;
		cout << "hhl res: " << hhl->GetResolution()[0] << " " << hhl->GetResolution()[1] << " " << hhl->GetResolution()[2] << endl;

		cout << "making uvw\n";
		RegVectorField3f uvw(W->GetResolution(), W->GetDomain());
		for (int i = 0; i < uvw.GetData().size(); ++i) {
			Vec3i gc = uvw.GetGridCoord(i);
			Vec3d val(U->GetVertexDataAt(gc), V->GetVertexDataAt(gc), W->GetVertexDataAt(gc));
			uvw.SetVertexDataAt(gc, val);
		}
		cout << "uvw done\n";

		for (double haa = 49.2; haa < 49.3; haa += 0.001) {
			//-2.08378 - 0.960101 49.241
			Vec3d p = coord; p[2] = haa;
			Vec3f was = uvw.Sample(p);
			if (abs(was[2] - trg[2]) < 0.001) {
				cout << "hm? p   " << p[0] << " " << p[1] << " " << p[2] << endl;
				cout << "   uvw  " << was[0] << " " << was[1] << " " << was[2] << endl;
				for (double hb = 48; hb < 52; hb += 0.001) {
					//no success :(			 what can i do?
					Vec3d p2 = coord; p[2] = hb;
					Vec3f wat = uvw.Sample(p2);
					if (abs(wat[0] - trg[0]) < 0.001 && abs(wat[1] - trg[1]) < 0.001) {
						cout << "      p2  " << p2[0] << " " << p2[1] << " " << p2[2] << endl;
						cout << "      uvw " << wat[0] << " " << wat[1] << " " << wat[2] << endl;
					}
				}
			}
		}


		for (int i = 110; i < 15; ++i) {
			for (int j = 10; j < 15; ++j) {
				int k = 20;
				Vec3i ij(i, j, k);
				cout << "HHL " << i << " " << j << " " << k << ": " << hhl->GetVertexDataAt(ij)<<endl;
			}
		}

		return 0;
		/*
		Vec3i res(8, 8, 10);
		Vec3d bmin(-1, -2, 0);
		Vec3d bmax(1, 1, 1);
		BoundingBox3d bb(bmin, bmax);
		RegScalarField3f* jojo = new RegScalarField3f(res,bb);
		RegScalarField3f hhl(res, bb);
		for (int i = 0; i < res[0]; ++i) {
			for (int j = 0; j < res[1]; ++j) {
				for (int k = 0; k < res[2]; ++k) {
					jojo->SetVertexDataAt(Vec3i(i, j, k), k);
					hhl.SetVertexDataAt(Vec3i(i, j, k), (res[2] - 1)*(res[2] - 1) - k * k);
				}
			}
		}
		for (int i = 0; i < res[2]; ++i) {
			cout << "lvl " << i << ": h = " << hhl.GetVertexDataAt(Vec3i(0, 0, i));
			cout << "\tand jojo = " << jojo->GetVertexDataAt(Vec3i(0, 0, i)) << endl;
		}
		RlonRlatHField field(jojo, &hhl);
		for (int i = 0; i < 50; ++i) {
			Vec3d pos(0, 0, i*1.2);
			cout << " h" << i << " = " << pos[2];
			float lvl = field.SampleXYiHd(0, 0, pos[2]);
			cout << "   has level " << lvl << endl;
			cout << "     and sample " << field.Sample(pos) << endl;
		}

		return 0;
		*/
	}
	else if (input == 2) {

		string path(argv[1]);
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
			cout << "Input bounds for tracing in (lat,lon,h) (Format: min_x max_x min_y max_y min_z max_z)\n> ";
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
			assert(t1 > t0);
			cout << "Input timestep> "; cin >> dt;
			assert(dt > 0);
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

			cout << "This will trace a total of " << nPaths << " trajectories (" << paths_dim[0] << " x " << paths_dim[1] << " x " << paths_dim[2] << ")\n";
			cout << "There will be up to " << nSteps << " steps per trajectory, storing "<<extra_variables.size()<<" extra variables\n";

			cout << "Confirm (0/1)> "; cin >> confirmed;
		}
		string name = "";
		cout << "Enter a filename> ";
		cin >> name;

		ImportantPart imp;

		TrajectoryData td;
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

		// 2 0 0.3 2 2.3 -10 8000 0.1 0.1 300 11 1111 44 1

		imp.trajectories.resize(nPaths);
		cout << "Printing initial points\n";

		ofstream points_file;
		points_file.open("points_" + name + ".txt");

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

		NetCDF::WriteTrajectoryData("trajectory_" + name + ".nc", td);
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
		/*
	#ifdef NDEBUG
	#pragma omp parallel for schedule(dynamic,16)
	#endif
		for (int64_t linearIndex = 0; linearIndex < field_size ; ++linearIndex)
		{
			Vec3i gridCoord = field->GetGridCoord(linearIndex);
			float value = field->GetVertexDataAt(gridCoord);
			// do something with value
			//if (linearIndex % checkpoint == 0) { std::cout << linearIndex << " of " << field_size<<std::endl; }
		}
		*/

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