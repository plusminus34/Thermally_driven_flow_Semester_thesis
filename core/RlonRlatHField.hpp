#pragma once

#include "RegularGrid.hpp"

#include <iostream>
#include <cmath>

using namespace std;

template<typename TValueType>
class RlonRlatHField : public ISampleField<TValueType, 3>{
public:
	RlonRlatHField(RegularGrid<TValueType, 3>* velocityField, RegScalarField3f* levelToHeight) {
		uvw = velocityField;
		hhl = levelToHeight;
		hhl_offset = Vec3i(1, 1, 0);
		lvl_top = hhl->GetResolution()[2] - 1;
		dx = uvw->GetVoxelSize()[0];
		dy = uvw->GetVoxelSize()[1];
		x0 = uvw->GetDomain().GetMin()[0];
		y0 = uvw->GetDomain().GetMin()[1];
		x1 = uvw->GetDomain().GetMax()[0];
		y1 = uvw->GetDomain().GetMax()[1];
		res_x = uvw->GetResolution()[0];
		res_y = uvw->GetResolution()[1];
	}
	~RlonRlatHField() { delete uvw; }

	virtual TValueType Sample(const Vec3d& coord) const override {
		double rlon_d = (coord[0] - x0) / (x1 - x0) * (res_x-1);
		double rlat_d = (coord[1] - y0) / (y1 - y0) * (res_y - 1);
		int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 1);
		int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 1);
		double w0 = rlon_d - rlon_i;
		double w1 = rlat_d - rlat_i;
		if (w0 < 0) { cout << "ERROR: w0 = " << w0 << " < 0" << endl; w0 = 0; }
		else if (w0 > 1) { cout << "ERROR: w0 = " << w0 << " > 1" << endl; w0 = 1; }
		if (w1 < 0) { cout << "ERROR: w1 = " << w1 << " < 0" << endl; w1 = 0; }
		else if (w1 > 1) { cout << "ERROR: w1 = " << w1 << " > 1" << endl; w1 = 1; }
		TValueType tmp0 = SampleXYiHd(rlon_i, rlat_i, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i, coord[2]) * w0;
		TValueType tmp1 = SampleXYiHd(rlon_i, rlat_i + 1, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i + 1, coord[2]) * w0;
		return tmp0 * (1 - w1) + tmp1 * w1;
	}
	TValueType SampleXYiHd(int rlon_i, int rlat_i, double h) const {
		/*
		seems	level	0		80
				height	22000	800-something
		*/
		assert(rlon_i > -1 && rlon_i < res_x);
		assert(rlat_i > -1 && rlat_i < res_y);

		//binary search at rlon_i,rlat_i to find level for h
		int lower = 0;
		int upper = lvl_top;
		if (h > hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lower) + hhl_offset)) return TValueType();
		if (h < hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, upper) + hhl_offset)) return TValueType();

		while (upper > lower + 1) {
			int half = (upper + lower) / 2;
			//cout << "half " << half << " between " << upper << " and " << lower << endl;
			if (hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, half) + hhl_offset) < h)
				upper = half;
			else lower = half;
		}
		const float h1 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, upper) + hhl_offset);
		const float h0 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lower) + hhl_offset);
		// and then sample uvw
		float alpha = (h - h0) / (h1 - h0);
		Vec3i smp0(rlon_i, rlat_i, lower);
		Vec3i smp1(rlon_i, rlat_i, upper);
		TValueType res = TValueType();
		if (lower > 0) res += uvw->GetVertexDataAt(smp0)*(1 - alpha);
		if (upper < uvw->GetResolution()[2]) res += uvw->GetVertexDataAt(smp1)*alpha;
		return res;
	}
private:
	RegularGrid<TValueType, 3>* uvw;
	RegScalarField3f* hhl;

	Vec3i hhl_offset;

	int lvl_top, res_x, res_y;
	double x0, y0, x1, y1, dx, dy;
};

typedef RlonRlatHField<Vec3f> RlonRlatHField_Vec3f;
typedef RlonRlatHField<float> RlonRlatHField_f;