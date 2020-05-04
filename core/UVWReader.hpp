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

RegVectorField3f* UVWFromNCFile(string filename) {
	std::string varname[3];
	std::string dimname[3][3];
	varname[0] = "U"; dimname[0][0] = "srlon"; dimname[0][1] = "rlat"; dimname[0][2] = "level";
	varname[1] = "V"; dimname[1][0] = "rlon"; dimname[1][1] = "srlat"; dimname[1][2] = "level";
	varname[2] = "W"; dimname[2][0] = "rlon"; dimname[2][1] = "rlat"; dimname[2][2] = "level1";

	NetCDF::Info info;
	if (!NetCDF::ReadInfo(filename, info)) { return nullptr; }

	// import field
	RegScalarField3f* field = NetCDF::ImportScalarField3f(filename, varname[0], dimname[0][0], dimname[0][1], dimname[0][2]);
	//int64_t field_size = (int64_t)field->GetResolution()[0] * field->GetResolution()[1] * field->GetResolution()[2];

	// init vector field
	const Vec3i uX(1, 0, 0);
	const Vec3i uY(0, 1, 0);
	const Vec3i uZ(0, 0, 1);
	RegVectorField3f* vecField = new RegVectorField3f(field->GetResolution() - uX - uY, field->GetDomain());
	int64_t vecField_size = (int64_t)vecField->GetResolution()[0] * vecField->GetResolution()[1] * vecField->GetResolution()[2];

	// U
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
	field = NetCDF::ImportScalarField3f(filename, varname[1], dimname[1][0], dimname[1][1], dimname[1][2]);
	Vec3f minV(field->GetDomain().GetMin());
	Vec3f maxV(field->GetDomain().GetMax());
	for (int i = 0; i < 3; ++i) {
		min[i] = (min[i] > minV[i] ? min[i] : minV[i]);
		max[i] = (max[i] < maxV[i] ? max[i] : maxV[i]);
	}
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
	field = NetCDF::ImportScalarField3f(filename, varname[2], dimname[2][0], dimname[2][1], dimname[2][2]);

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

	delete field;
	return vecField;
}