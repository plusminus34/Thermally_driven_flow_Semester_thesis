#include "core/NetCDF.hpp"

#include <iostream>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkXMLImageDataWriter.h>

// ---------------------------------------
// Entry point
// ---------------------------------------
int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "No file to process!" << std::endl;
		return 0;
	}
	bool custom = false;
	std::string varname;
	std::string dimname_0, dimname_1, dimname_2;
	if (argc >= 3) {
		varname = argv[2];
		custom = true;
	}
	else {
		varname = "U";
		dimname_0 = "srlon";
		dimname_1 = "rlat";
		dimname_2 = "level";
	}
	std::string path(argv[1]);

	std::cout << "Variable is: " << varname << std::endl;

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
		std::cout << "X dimension is "; std::cin >> input;
		dimname_0 = variable.Dimensions[input].GetName();
		std::cout << "Y dimension is "; std::cin >> input;
		dimname_1 = variable.Dimensions[input].GetName();
		std::cout << "Z dimension is "; std::cin >> input;
		dimname_2 = variable.Dimensions[input].GetName();
	}

	std::cout << "Importing field..." << std::endl;
	// import field
	RegScalarField3f* field = NetCDF::ImportScalarField3f(path, varname, dimname_0, dimname_1, dimname_2 );
	std::cout << "Checking values afterwards...\n";
	int64_t cap = (int64_t)field->GetResolution()[0] * field->GetResolution()[1] * field->GetResolution()[2];
	int64_t checkpoint = cap / 50;
#ifdef NDEBUG
#pragma omp parallel for schedule(dynamic,16)
#endif
	for (int64_t linearIndex = 0; linearIndex < cap; ++linearIndex)
	{
		Vec3i gridCoord = field->GetGridCoord(linearIndex);
		float value = field->GetVertexDataAt(gridCoord);
		// do something with value
		//if (linearIndex % checkpoint == 0) { std::cout << linearIndex << " of " << cap<<std::endl; }
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
	return 0;
}