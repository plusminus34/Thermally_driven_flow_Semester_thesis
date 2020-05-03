#include "core/NetCDF.hpp"

#include <iostream>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkFloatArray.h>

#include "core/CoordinateTransform.hpp"
#include "core/ParticleTracer.hpp"
#include "core/UVWReader.hpp"

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
	std::cout << "2:\tRead UVW and trace particles\n> ";
	cin >> input;
	if (input < 2){
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
				std::cout << varname << " attribute " << i << ": " << variable.Attributes[i].GetName() << std::endl;
				std::cout << varname << " attribute " << i << ": " << variable.Attributes[i].GetValue() << std::endl;
			}
			std::cout << "X dimension is "; std::cin >> input;
			dimname_0 = variable.Dimensions[input].GetName();
			std::cout << "Y dimension is "; std::cin >> input;
			dimname_1 = variable.Dimensions[input].GetName();
			std::cout << "Z dimension is "; std::cin >> input;
			dimname_2 = variable.Dimensions[input].GetName();
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
			// init vector field
			const Vec3i uX(1, 0, 0);
			const Vec3i uY(0, 1, 0);
			const Vec3i uZ(0, 0, 1);
			RegVectorField3f* vecField = new RegVectorField3f(field->GetResolution() - uX - uY, field->GetDomain());
			int64_t vecField_size = (int64_t)vecField->GetResolution()[0] * vecField->GetResolution()[1] * vecField->GetResolution()[2];

			// U
			std::cout << varname << ": " << dimname_0 << " " << field->GetResolution()[0] << " x "
				<< dimname_1 << " " << field->GetResolution()[1] << " x "
				<< dimname_2 << " " << field->GetResolution()[2] << std::endl;
			std::cout << "  domain: " << field->GetDomain().GetMin()[0] << " to " << field->GetDomain().GetMax()[0] << "\t"
				<< field->GetDomain().GetMin()[1] << " to " << field->GetDomain().GetMax()[1] << "\t"
				<< field->GetDomain().GetMin()[2] << " to " << field->GetDomain().GetMax()[2] << std::endl;

			Vec3f min(field->GetDomain().GetMin());
			Vec3f max(field->GetDomain().GetMax());
			std::cout << "Interpolating U...\n";
			for (int64_t i = 0; i < vecField_size; ++i)
			{
				Vec3i gridCoord = vecField->GetGridCoord(i);
				float value = (field->GetVertexDataAt(gridCoord + uY) + field->GetVertexDataAt(gridCoord + uX + uY)) * 0.5f;
				Vec3f v(value, 0, 0);
				vecField->SetVertexDataAt(gridCoord, v);
			}


			// V
			varname = "V"; dimname_0 = "rlon"; dimname_1 = "srlat";
			field = NetCDF::ImportScalarField3f(path, varname, dimname_0, dimname_1, dimname_2);
			std::cout << varname << ": " << dimname_0 << " " << field->GetResolution()[0] << " x "
				<< dimname_1 << " " << field->GetResolution()[1] << " x "
				<< dimname_2 << " " << field->GetResolution()[2] << std::endl;
			std::cout << "  domain: " << field->GetDomain().GetMin()[0] << " to " << field->GetDomain().GetMax()[0] << "\t"
				<< field->GetDomain().GetMin()[1] << " to " << field->GetDomain().GetMax()[1] << "\t"
				<< field->GetDomain().GetMin()[2] << " to " << field->GetDomain().GetMax()[2] << std::endl;

			Vec3f minV(field->GetDomain().GetMin());
			Vec3f maxV(field->GetDomain().GetMax());
			for (int i = 0; i < 3; ++i) {
				min[i] = (min[i] > minV[i] ? min[i] : minV[i]);
				max[i] = (max[i] < maxV[i] ? max[i] : maxV[i]);
			}
			std::cout << "Interpolating V...\n";
			for (int64_t i = 0; i < vecField_size; ++i)
			{
				Vec3i gridCoord = vecField->GetGridCoord(i);
				float value = (field->GetVertexDataAt(gridCoord + uX) + field->GetVertexDataAt(gridCoord + uX + uY)) * 0.5f;
				Vec3f v = vecField->GetVertexDataAt(gridCoord); v[1] = value;
				vecField->SetVertexDataAt(gridCoord, v);
			}

			//rescale U and V to rlon/rlat per second
			cout << "Rescaling U and V...\n";
			for (int64_t i = 0; i < vecField_size; ++i)
			{
				Vec3i gridCoord = vecField->GetGridCoord(i);
				Vec3d coord = vecField->GetCoordAt(gridCoord);
				Vec3f v = vecField->GetVertexDataAt(gridCoord);
				double lon, lat;
				RlatRlonToLatLon(coord[1], coord[0], lat, lon);
				double dlon, dlat;
				degreeLengths(lat, dlat, dlon);
				v[1] /= dlon;
				v[1] /= dlat;
				vecField->SetVertexDataAt(gridCoord, v);
			}

			// W
			delete field;
			varname = "W"; dimname_1 = "rlat"; dimname_2 = "level1";
			field = NetCDF::ImportScalarField3f(path, varname, dimname_0, dimname_1, dimname_2);
			std::cout << varname << ": " << dimname_0 << " " << field->GetResolution()[0] << " x "
				<< dimname_1 << " " << field->GetResolution()[1]
				<< " x " << dimname_2 << " " << field->GetResolution()[2] << std::endl;
			std::cout << "  domain: " << field->GetDomain().GetMin()[0] << " to " << field->GetDomain().GetMax()[0] << "\t"
				<< field->GetDomain().GetMin()[1] << " to " << field->GetDomain().GetMax()[1] << "\t"
				<< field->GetDomain().GetMin()[2] << " to " << field->GetDomain().GetMax()[2] << std::endl;

			Vec3f minW(field->GetDomain().GetMin());
			Vec3f maxW(field->GetDomain().GetMax());
			for (int i = 0; i < 3; ++i) {
				min[i] = (min[i] > minW[i] ? min[i] : minW[i]);
				max[i] = (max[i] < maxW[i] ? max[i] : maxW[i]);
			}
			BoundingBox3d domain(min, max);
			vecField->UpdateDomain(domain);
			std::cout << "Interpolating W...\n";
			for (int64_t i = 0; i < vecField_size; ++i)
			{
				Vec3i gridCoord = vecField->GetGridCoord(i);
				float value = (field->GetVertexDataAt(gridCoord + uX + uY) + field->GetVertexDataAt(gridCoord + uX + uY + uZ)) * 0.5f;
				Vec3f v = vecField->GetVertexDataAt(gridCoord); v[2] = value;
				vecField->SetVertexDataAt(gridCoord, v);
			}

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
		int input = 0;
		std::cout << "continue? "; cin >> input; if (!input) return 0;
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
	else {
		RegVectorField3f* field = UVWFromVTIFile("UVW.vti");
		field->GetDomain();
		double bounds[6];
		for (int i = 0; i < 3; ++i) {
			bounds[2 * i] = field->GetDomain().GetMin()[i];
			bounds[2 * i + 1] = field->GetDomain().GetMax()[i];
		}
		cout << "Bounds: x " << bounds[0] << " - " << bounds[1] << "\ty " << bounds[2] << " - " << bounds[3] << "\t z " << bounds[4] << " - " << bounds[5] << endl;
		double tracing_bounds[6];
		double spacing;
		int paths_dim[3];
		int nPaths;
		double t0, t1, dt;
		int nSteps;
		bool confirmed = false;
		while (!confirmed) {
			cout << "Input bounds for tracing (Format: min_x max_x min_y max_y min_z max_z)\n> ";
			for (int i = 0; i < 6; ++i) {
				cin >> tracing_bounds[i];
				if (i % 2 == 0) tracing_bounds[i] = std::max(tracing_bounds[i], bounds[i]);
				else tracing_bounds[i] = std::min(tracing_bounds[i], bounds[i]);
			}
			cout << "Input spacing between individual starting points> ";
			cin >> spacing;
			for (int i = 0; i < 3; ++i) {
				paths_dim[i] = floor((tracing_bounds[2 * i + 1] - tracing_bounds[2 * i]) / spacing);
				assert(paths_dim[i] > 0);
			}
			nPaths = paths_dim[0] * paths_dim[1] * paths_dim[2];
			cout << "Input start and end time (Format: t0 t1)> "; cin >> t0; cin >> t1;
			assert(t1 > t0);
			cout << "Input timestep> "; cin >> dt;
			assert(dt > 0);
			nSteps = ceil((t1 - t0) / dt);
			cout << "This will trace a total of " << nPaths << " trajectories (" << paths_dim[0] << " x " << paths_dim[1] << " x " << paths_dim[2] << ")\n";
			cout << "There will be up to "<<nSteps<<" steps per trajectory\n";
			cout << "Confirm (0/1)> "; cin >> confirmed;
		}

		std::vector<std::vector<Vec3f>> paths(nPaths);

		ParticleTracer<Vec3f, 3> tracer;
		for (int i = 0; i < paths_dim[0];++i) {
			for (int j = 0; j < paths_dim[1]; ++j) {
				for (int k = 0; k < paths_dim[2]; ++k) {
					int path = i * paths_dim[1] * paths_dim[2] + j * paths_dim[2] + k;
					Vec3f position = Vec3f(tracing_bounds[0] + i*spacing, tracing_bounds[2] + j * spacing, tracing_bounds[3] + k * spacing);
					paths[path].push_back(position);
					for (int l = 0; l < nSteps; ++l) {
						Vec3d pos_d(position[0], position[1], position[2]);
						position = tracer.traceParticle(*field, pos_d, dt);
						paths[path].push_back(position);
					}
				}
			}
		}

		//TODO
		cout << "THE END\n";
		cout << "Trajectories are discarded\n";
		delete field;
	}
	return 0;
}