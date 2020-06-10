#pragma once
#include "LagrantoUVW.hpp"

#include "UVWReader.hpp"

LagrantoUVW::LagrantoUVW(string filename, RegScalarField3f * levelToHeight)
{
	SeparateUVWFromNCFile(filename, U, V, W);
	hhl = levelToHeight;

	lvl_top = U->GetResolution()[2] - 1;
	res_x = hhl->GetResolution()[0];
	res_y = hhl->GetResolution()[1];
	dx = hhl->GetVoxelSize()[0];
	dy = hhl->GetVoxelSize()[1];
	x0 = hhl->GetDomain().GetMin()[0];
	y0 = hhl->GetDomain().GetMin()[1];
	x1 = hhl->GetDomain().GetMax()[0];
	y1 = hhl->GetDomain().GetMax()[1];
}

LagrantoUVW::~LagrantoUVW()
{
	delete U;
	delete V;
	delete W;
}

Vec3f LagrantoUVW::Sample(const Vec3d & coord) const
{
	Vec3f uvw(0, 0, 0);

	double rlon_d = (coord[0] - x0) / (x1 - x0) * (res_x - 1);
	double rlat_d = (coord[1] - y0) / (y1 - y0) * (res_y - 1);

	double h = coord[2];
	double lvl_d[2];//[0] rlon,rlat,level	[1] rlon,rlat,level1
	for (int stag = 0; stag < 2; ++stag) {
		int lvl_0 = 0;
		int lvl_1 = lvl_top + stag;

		float h0, h1, hh;

		h0 = sampleHHL(rlon_d, rlat_d, lvl_0, !stag);
		h1 = sampleHHL(rlon_d, rlat_d, lvl_1, !stag);

		if (h > h0 || h < h1) return Vec3f(0, 0, 0);

		while (lvl_1 > lvl_0 + 1) {
			int half = (lvl_1 + lvl_0) / 2;
			//cout << "half " << half << " between " << lvl_1 << " and " << lvl_0 << endl;
			hh = sampleHHL(rlon_d, rlat_d, half, !stag);
			if (hh < h) {
				lvl_1 = half;
				h1 = hh;
			}
			else {
				lvl_0 = half;
				h0 = hh;
			}
		}

		float alpha = (h - h0) / (h1 - h0);
		//cout << "localized: h " << h << " belongs between h0 and h1 " << h0 << " " << h1 << endl;
		//cout << "   which means level1s " << lvl_0 << " " << lvl_1 << endl;
		lvl_d[stag] = lvl_0 + alpha;
	}
	// The fact that these make the result more Lagranto-like is worrying
	lvl_d[0] -= 0.0;//1.0;
	lvl_d[1] -= 1.0;//2.0;

	const int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 1);
	const int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 1);
	const float wx = rlon_d - rlon_i;
	const float wy = rlat_d - rlat_i;

	for (int stag = 0; stag < 2; ++stag) {
		int lvl_i = lvl_d[stag];
		double wz = lvl_d[stag] - lvl_i;

		Vec3i gcs[8];
		double ws[8];
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				for (int k = 0; k < 2; ++k) {
					ws[4 * i + 2 * j + k] = 0;
				}
			}
		}
		ws[0] = (1 - wx) * (1 - wy) * (1 - wz); gcs[0] = Vec3i(rlon_i, rlat_i, lvl_i);
		ws[1] = wx * (1 - wy) * (1 - wz); gcs[1] = Vec3i(rlon_i + 1, rlat_i, lvl_i);
		ws[2] = (1 - wx) * wy * (1 - wz); gcs[2] = Vec3i(rlon_i, rlat_i + 1, lvl_i);
		ws[3] = wx * wy * (1 - wz); gcs[3] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i);
		ws[4] = (1 - wx) * (1 - wy) * wz; gcs[4] = Vec3i(rlon_i, rlat_i, lvl_i + 1);
		ws[5] = wx * (1 - wy) * wz; gcs[5] = Vec3i(rlon_i + 1, rlat_i, lvl_i + 1);
		ws[6] = (1 - wx) * wy * wz; gcs[6] = Vec3i(rlon_i, rlat_i + 1, lvl_i + 1);
		ws[7] = wx * wy * wz; gcs[7] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i + 1);

		for (int i = 0; i < 8; ++i) {
			if (stag == 0) {
				uvw[0] += U->GetVertexDataAt(gcs[i])*ws[i];
				uvw[1] += V->GetVertexDataAt(gcs[i])*ws[i];
			}
			else if (stag == 1) {
				uvw[2] += W->GetVertexDataAt(gcs[i])*ws[i];
			}
		}
	}

	/*
	// rescale to rlon/s rlat/s m/s
	double lon, lat;
	CoordinateTransform::RlatRlonToLatLon(coord[1], coord[0], lat, lon);
	double dlon, dlat;
	CoordinateTransform::degreeLengthsSimple(lat, dlat, dlon);
	if (dlon != 0) uvw[0] /= dlon;
	else uvw[0] = 0;
	uvw[1] /= dlat;
	*/
	//cout << "uvw unmodified " << uvw[0] << " " << uvw[1] << " " << uvw[2] << endl;
	float deltay = 111200;
	uvw[1] /= deltay;
	uvw[0] /= deltay * cos(coord[1] * 3.1415926535 / 180.0);

	return uvw;
}

float LagrantoUVW::sampleHHL(double rlon_d, double rlat_d, int level, bool destagger) const
{
	float res = 0;
	const int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 1);
	const int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 1);
	const float w0 = rlon_d - rlon_i;
	const float w0o = 1 - w0;
	const float w1 = rlat_d - rlat_i;
	const float w1o = 1 - w1;
	for (int ix = 0; ix < 2; ++ix) {
		float wx = ix == 0 ? w0o : w0;
		for (int iy = 0; iy < 2; ++iy) {
			float wy = iy == 0 ? w1o : w1;
			res += wx * wy * hhl->GetVertexDataAt(Vec3i(rlon_i + ix, rlat_i + iy, level));
			if(destagger) res += wx * wy * hhl->GetVertexDataAt(Vec3i(rlon_i + ix, rlat_i + iy, level+1));
		}
	}
	if (destagger) res *= 0.5f;
	return res;
}
