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

	// assorted constants
	const double rlon_d = (coord[0] - x0) / (x1 - x0) * (res_x - 1);
	const double rlat_d = (coord[1] - y0) / (y1 - y0) * (res_y - 1);
	const int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 2);
	const int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 2);
	const float w0 = rlon_d - rlon_i;
	const float w0o = 1 - w0;
	const float w1 = rlat_d - rlat_i;
	const float w1o = 1 - w1;

	if (z_sampling_alternative) {
		// Sample each column separately, then combine via bilinear interpolation
		double h = coord[2];
		for (int i = rlon_i; i < rlon_i + 2; ++i) {
			const float wx = i > rlon_i ? rlon_d - rlon_i : rlon_i + 1 - rlon_d;
			for (int j = rlat_i; j < rlat_i + 2; ++j) {
				const float wy = j > rlat_i ? rlat_d - rlat_i : rlat_i + 1 - rlat_d;
				int lvl_0 = 0;
				int lvl_1 = hhl->GetResolution()[2] - 2;

				float h0, h1, hh;

				h0 = (hhl->GetVertexDataAt(Vec3i(i, j, lvl_0)) + hhl->GetVertexDataAt(Vec3i(i, j, lvl_0 + 1)))*0.5f;
				h1 = (hhl->GetVertexDataAt(Vec3i(i, j, lvl_1)) + hhl->GetVertexDataAt(Vec3i(i, j, lvl_1 + 1)))*0.5f;
				if (h > h0 || h < h1) return Vec3f(0, 0, 0);


				while (lvl_1 > lvl_0 + 1) {
					int half = (lvl_1 + lvl_0) / 2;
					hh = (hhl->GetVertexDataAt(Vec3i(i, j, half)) + hhl->GetVertexDataAt(Vec3i(i, j, half + 1)))*0.5f;
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
				uvw[0] += U->GetVertexDataAt(Vec3i(i, j, lvl_0)) * wx*wy*(1 - alpha);
				uvw[0] += U->GetVertexDataAt(Vec3i(i, j, lvl_1)) * wx*wy*alpha;
				uvw[1] += V->GetVertexDataAt(Vec3i(i, j, lvl_0)) * wx*wy*(1 - alpha);
				uvw[1] += V->GetVertexDataAt(Vec3i(i, j, lvl_1)) * wx*wy*alpha;
				uvw[2] += W->GetVertexDataAt(Vec3i(i, j, lvl_0)) * wx*wy*(1 - alpha);
				uvw[2] += W->GetVertexDataAt(Vec3i(i, j, lvl_1)) * wx*wy*alpha;
			}
		}
	}
	else {
		// Construct a local HHL-column by bilinearly interpolating heights at each level
		// Not suggested because equal levels at different rlon_i,rlat_i are not actually
		//   guaranteed to be at similar heights
		double h = coord[2];
		double lvl_d[2];//[0] is for rlon,rlat,level	[1] rlon,rlat,level1
		for (int stag = 0; stag < 2; ++stag) {
			if (use_correct_zw && stag > 0) continue;
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
		// The fact that this makes the result more Lagranto-like is worrying
		if (!use_correct_zw) lvl_d[1] -= 1.0;

		for (int stag = 0; stag < 2; ++stag) {
			if (use_correct_zw && stag > 0) continue;
			int lvl_i = lvl_d[stag];
			if (lvl_i<0 || lvl_i> U->GetResolution()[2] - 2) continue;
			const double w2 = lvl_d[stag] - lvl_i;
			const double w2o = 1.0 - w2;

			Vec3i gcs[8];
			double ws[8];

			ws[0] = w0o * w1o * w2o; gcs[0] = Vec3i(rlon_i, rlat_i, lvl_i);
			ws[1] = w0 * w1o * w2o; gcs[1] = Vec3i(rlon_i + 1, rlat_i, lvl_i);
			ws[2] = w0o * w1 * w2o; gcs[2] = Vec3i(rlon_i, rlat_i + 1, lvl_i);
			ws[3] = w0 * w1 * w2o; gcs[3] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i);
			ws[4] = w0o * w1o * w2; gcs[4] = Vec3i(rlon_i, rlat_i, lvl_i + 1);
			ws[5] = w0 * w1o * w2; gcs[5] = Vec3i(rlon_i + 1, rlat_i, lvl_i + 1);
			ws[6] = w0o * w1 * w2; gcs[6] = Vec3i(rlon_i, rlat_i + 1, lvl_i + 1);
			ws[7] = w0 * w1 * w2; gcs[7] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i + 1);

			for (int i = 0; i < 8; ++i) {
				if (stag == 0) {
					uvw[0] += U->GetVertexDataAt(gcs[i])*ws[i];
					uvw[1] += V->GetVertexDataAt(gcs[i])*ws[i];
					if (use_correct_zw) uvw[2] += W->GetVertexDataAt(gcs[i])*ws[i];
				}
				else if (stag == 1) {
					uvw[2] += W->GetVertexDataAt(gcs[i])*ws[i];
				}
			}
		}
	}

	// rescale to rlon/s rlat/s m/s
	//cout << "uvw unmodified " << uvw[0] << " " << uvw[1] << " " << uvw[2] << endl;
	const float deltay = 111200;// length of one degree latitude in meters
	uvw[1] /= deltay;
	uvw[0] /= deltay * cos(coord[1] * 3.1415926535 / 180.0);
	//cout << "uvw: " << uvw[0] << " " << uvw[1] << " " << uvw[2]<<endl;
	return uvw;
}

float LagrantoUVW::sampleHHL(double rlon_d, double rlat_d, int level, bool destagger) const
{
	float res = 0;
	const int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 2);
	const int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 2);
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