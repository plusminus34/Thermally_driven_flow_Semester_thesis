#pragma once
#include "UVWReader.hpp"
#include <vector>
#include "NetCDF.hpp"
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

RegVectorField3f* UVWFromVTIFile(string filename) {
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
	Vec3d bb_min(bounds[0], bounds[2], bounds[4]);
	Vec3d bb_max(bounds[1], bounds[3], bounds[5]);
	RegVectorField3f* uvw = new RegVectorField3f(res, BoundingBox3d(bb_min, bb_max));

	vtkSmartPointer<vtkPointData> pointData = imageData->GetPointData();
	vtkSmartPointer<vtkDataArray> dataArray = pointData->GetArray("U");// apparently it's still called W ... TODO change that

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
	cout << "Importing UVW from " << filename << endl;
	// Constants: variable names and unit vectors
	std::string varname[3];
	std::string dimname[3][3];
	varname[0] = "U"; dimname[0][0] = "srlon"; dimname[0][1] = "rlat"; dimname[0][2] = "level";
	varname[1] = "V"; dimname[1][0] = "rlon"; dimname[1][1] = "srlat"; dimname[1][2] = "level";
	varname[2] = "W"; dimname[2][0] = "rlon"; dimname[2][1] = "rlat"; dimname[2][2] = "level1";
	const Vec3i uX(1, 0, 0); const Vec3i uY(0, 1, 0); const Vec3i uZ(0, 0, 1);

	// import U for the start
	RegScalarField3f* field = NetCDF::ImportScalarField3f(filename, varname[0], dimname[0][0], dimname[0][1], dimname[0][2]);
	if (field == NULL) return nullptr;

	// init vector field based on dimensions of U
	RegVectorField3f* uvw = new RegVectorField3f(field->GetResolution() - uX - uY, field->GetDomain());
	int64_t vecField_size = (int64_t)uvw->GetResolution()[0] * uvw->GetResolution()[1] * uvw->GetResolution()[2];

	// U
	Vec3f min(field->GetDomain().GetMin());
	Vec3f max(field->GetDomain().GetMax());
	for (int64_t i = 0; i < vecField_size; ++i)
	{
		Vec3i gridCoord = uvw->GetGridCoord(i);
		float value = (field->GetVertexDataAt(gridCoord + uY) + field->GetVertexDataAt(gridCoord + uX + uY)) * 0.5f;
		Vec3f v(value, 0, 0);
		uvw->SetVertexDataAt(gridCoord, v);
	}
	delete field;

	// V
	field = NetCDF::ImportScalarField3f(filename, varname[1], dimname[1][0], dimname[1][1], dimname[1][2]);
	if (field == NULL){
		delete uvw;
		return nullptr;
	}
	Vec3f minV(field->GetDomain().GetMin());
	Vec3f maxV(field->GetDomain().GetMax());
	for (int i = 0; i < 3; ++i) {
		min[i] = (min[i] > minV[i] ? min[i] : minV[i]);
		max[i] = (max[i] < maxV[i] ? max[i] : maxV[i]);
	}
	for (int64_t i = 0; i < vecField_size; ++i)
	{
		Vec3i gridCoord = uvw->GetGridCoord(i);
		float value = (field->GetVertexDataAt(gridCoord + uX) + field->GetVertexDataAt(gridCoord + uX + uY)) * 0.5f;
		Vec3f v = uvw->GetVertexDataAt(gridCoord); v[1] = value;
		uvw->SetVertexDataAt(gridCoord, v);
	}
	delete field;

	//rescale U and V to rlon/rlat per second
	for (int64_t i = 0; i < vecField_size; ++i)
	{
		Vec3i gridCoord = uvw->GetGridCoord(i);
		Vec3d coord = uvw->GetCoordAt(gridCoord);
		Vec3f v = uvw->GetVertexDataAt(gridCoord);
		double lon, lat;
		CoordinateTransform::RlatRlonToLatLon(coord[1], coord[0], lat, lon);
		double dlon, dlat;
		CoordinateTransform::degreeLengthsSimple(lat, dlat, dlon);
		if (dlon != 0) v[0] /= dlon;
		else v[0] = 0;
		v[1] /= dlat;
		uvw->SetVertexDataAt(gridCoord, v);
	}
	// W
	field = NetCDF::ImportScalarField3f(filename, varname[2], dimname[2][0], dimname[2][1], dimname[2][2]);
	if (field == NULL) {
		delete uvw;
		return nullptr;
	}
	// also here: finish defining the domain (TODO do it in a nicer way)
	Vec3f minW(field->GetDomain().GetMin());
	Vec3f maxW(field->GetDomain().GetMax());
	for (int i = 0; i < 3; ++i) {
		min[i] = (min[i] > minW[i] ? min[i] : minW[i]);
		max[i] = (max[i] < maxW[i] ? max[i] : maxW[i]);
	}
	BoundingBox3d domain(min, max);
	uvw->UpdateDomain(domain);
	// interpolate W values
	for (int64_t i = 0; i < vecField_size; ++i)
	{
		Vec3i gridCoord = uvw->GetGridCoord(i);
		float value = (field->GetVertexDataAt(gridCoord + uX + uY) + field->GetVertexDataAt(gridCoord + uX + uY + uZ)) * 0.5f;
		Vec3f v = uvw->GetVertexDataAt(gridCoord); v[2] = value;
		uvw->SetVertexDataAt(gridCoord, v);
	}
	delete field;

	return uvw;
}

void SeparateUVWFromNCFile(string filename, RegScalarField3f* &U, RegScalarField3f* &V, RegScalarField3f * &W)
{
	cout << "Importing UVW from " << filename << " (Lagranto style)" << endl;
	// Constants: variable names and unit vectors
	const Vec3i uX(1, 0, 0); const Vec3i uY(0, 1, 0); const Vec3i uZ(0, 0, 1);
	cout << "reading U\n";
	RegScalarField3f* U0 = NetCDF::ImportScalarField3f(filename, "U", "srlon", "rlat", "level");
	U = new RegScalarField3f(U0->GetResolution(), U0->GetDomain());
	for (int i = 0; i < U->GetData().size(); ++i) {
		Vec3i gc = U->GetGridCoord(i);
		if (gc[0] > 0) {
			float val = 0.5*(U0->GetVertexDataAt(gc) + U0->GetVertexDataAt(gc - uX));
			U->SetVertexDataAt(gc, val);
		}
		else {
			// inaccurate slice at the end: just like in lagranto
			U->SetVertexDataAt(gc, U0->GetVertexDataAt(gc));
		}
	}
	delete U0;
	cout << "reading V\n";
	RegScalarField3f* V0 = NetCDF::ImportScalarField3f(filename, "V", "rlon", "srlat", "level");
	V = new RegScalarField3f(V0->GetResolution(), V0->GetDomain());
	for (int i = 0; i < V->GetData().size(); ++i) {
		Vec3i gc = V->GetGridCoord(i);
		if (gc[1] > 0) {
			float val = 0.5*(V0->GetVertexDataAt(gc) + V0->GetVertexDataAt(gc - uY));
			V->SetVertexDataAt(gc, val);
		}
		else {
			V->SetVertexDataAt(gc, V0->GetVertexDataAt(gc));
		}
	}
	delete V0;
	cout << "reading W\n";
	RegScalarField3f* W0 = NetCDF::ImportScalarField3f(filename, "W", "rlon", "rlat", "level1");
	Vec3d wmin = W0->GetDomain().GetMin();
	Vec3d wmax = W0->GetDomain().GetMax();
	wmax[2] -= 1;
	W = new RegScalarField3f(W0->GetResolution() - uZ, BoundingBox3d(wmin, wmax));
	for (size_t i = 0; i < W->GetData().size(); ++i) {
		Vec3i gridCoord = W->GetGridCoord(i);
		float val = 0.5*(W0->GetVertexDataAt(gridCoord) + W0->GetVertexDataAt(gridCoord + uZ));
		W->SetVertexDataAt(gridCoord, val);
	}
	delete W0;
}
