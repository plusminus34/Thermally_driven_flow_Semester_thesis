#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

using namespace std;

RegVectorField3f* UVWFromVTIFile(string filename = "UVW.vti") {
	//Read UVW and convert into RegVectorfield
	vtkNew<vtkXMLImageDataReader> imageReader;
	cout << "Reading file " << filename << endl;
	imageReader->SetFileName(filename.c_str());
	imageReader->Update();
	vtkSmartPointer<vtkImageData> imageData = imageReader->GetOutput();

	cout << "Got imageData\n";
	int* dims = imageData->GetDimensions();
	double* bounds = imageData->GetBounds();

	Vec3i res(dims[0], dims[1], dims[2]);
	Vec3d bb_min(bounds[0], bounds[2], bounds[3]);
	Vec3d bb_max(bounds[1], bounds[3], bounds[5]);
	RegVectorField3f* uvw = new RegVectorField3f(res, BoundingBox3d(bb_min, bb_max));

	vtkSmartPointer<vtkPointData> pointData = imageData->GetPointData();
	vtkSmartPointer<vtkDataArray> dataArray = pointData->GetArray("W");// apparently it's still called W ... TODO change that

	int nPts = imageData->GetNumberOfPoints();
	assert(dataArray->GetNumberOfTuples() == nPts);
	cout << "Filling field ...\n";
	for (int i = 0; i < nPts; ++i) {
		double* pt = dataArray->GetTuple3(i);
		Vec3i coord = uvw->GetGridCoord(i);
		uvw->SetVertexDataAt(coord, Vec3f(pt[0], pt[1], pt[2]));
	}
	return uvw;
}