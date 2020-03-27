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
	std::string path(argv[1]);

	// read info object
 	NetCDF::Info info;
	if (!NetCDF::ReadInfo(path, info)) { return -1; }

	// import field
	RegScalarField3f* U = NetCDF::ImportScalarField3f(path, "U", "srlon", "rlat", "level" );

#ifdef NDEBUG
#pragma omp parallel for schedule(dynamic,16)
#endif
	for (int64_t linearIndex = 0; linearIndex < (int64_t)U->GetResolution()[0] * U->GetResolution()[1] * U->GetResolution()[2]; ++linearIndex)
	{
		Vec3i gridCoord = U->GetGridCoord(linearIndex);
		float value = U->GetVertexDataAt(gridCoord);
		// do something with value
	}

	// allocate image data object
	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->SetDimensions(U->GetResolution().ptr());
	double spacing[] = {
		(U->GetDomain().GetMax()[0] - U->GetDomain().GetMin()[0]) / (U->GetResolution()[0] - 1.),
		(U->GetDomain().GetMax()[1] - U->GetDomain().GetMin()[1]) / (U->GetResolution()[1] - 1.),
		(U->GetDomain().GetMax()[2] - U->GetDomain().GetMin()[2]) / (U->GetResolution()[2] - 1.),
	};
	imageData->SetSpacing(spacing);

	// fill image data
	vtkSmartPointer<vtkFloatArray> floatArray = vtkSmartPointer<vtkFloatArray>::New();
	floatArray->SetNumberOfComponents(1);
	int numTuples = U->GetResolution()[0] * U->GetResolution()[1] * U->GetResolution()[2];
	floatArray->SetNumberOfTuples(numTuples);
	floatArray->SetName("U");
	for (int linearIndex = 0; linearIndex < numTuples; ++linearIndex)
	{
		Vec3i gridCoord = U->GetGridCoord(linearIndex);
		float value = U->GetVertexDataAt(gridCoord);
		floatArray->SetValue(linearIndex, value);
	}
	imageData->GetPointData()->AddArray(floatArray);

	// write to file
	vtkSmartPointer<vtkXMLImageDataWriter> imageWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	imageWriter->SetFileName("test.vti");
	imageWriter->SetInputData(imageData);
	imageWriter->Update();

	// delete resources and return
	delete U;
	return 0;
}