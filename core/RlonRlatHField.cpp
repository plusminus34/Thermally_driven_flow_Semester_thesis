#pragma once

#include "RlonRlatHField.hpp"

RlonRlatHField::RlonRlatHField(RegVectorField3f * velocityField, RegScalarField3f * levelToHeight)
{
	uvw = velocityField;
	hhl = levelToHeight;
	lvl_top = hhl->GetResolution()[2] - 1;
	dx = uvw->GetVoxelSize()[0];
	dy = uvw->GetVoxelSize()[1];
}

Vec3f RlonRlatHField::Sample(const Vec3d & coord) const
{
	int rlon_i = floor(coord[0]);
	int rlat_i = floor(coord[1]);
	float w0 = (coord[0] - rlon_i * dx) / dx;
	float w1 = (coord[1] - rlat_i * dy) / dy;
	Vec3f tmp0 = SampleXYiHd(rlon_i, rlat_i, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i, coord[2]) * w0;
	Vec3f tmp1 = SampleXYiHd(rlon_i, rlat_i + 1, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i + 1, coord[2]) * w0;
	return tmp0 * (1 - w1) + tmp1 * w1;
}

Vec3f RlonRlatHField::SampleXYiHd(int rlon_i, int rlat_i, double h) const
{
	//binary search at rlon_i,rlat_i to find level for h
	int lower = 0;
	int upper = lvl_top;
	if (h < hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lower))) return Vec3f(0,0,0);
	if (h > hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, upper))) return Vec3f(0,0,0);
	while (upper > lower + 1) {
		int half = (upper + lower) / 2;
		if (hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, half)) > h)
			upper = half;
		else lower = half;
	}
	const float h1 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, upper));
	const float h0 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lower));
	double lvl = (h - h0) / (h1 - h0);
	// and then sample uvw
	Vec3d samplePoint = uvw->GetCoordAt(Vec3i(rlon_i, rlat_i, 0));
	samplePoint[2] = lvl;
	return uvw->Sample(samplePoint);
}