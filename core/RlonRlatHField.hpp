#pragma once

#include "RegularGrid.hpp"

#include <iostream>

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
	}
	~RlonRlatHField() { delete uvw; }

	virtual TValueType Sample(const Vec3d& coord) const override {
		double x0 = uvw->GetDomain().GetMin()[0];
		double y0 = uvw->GetDomain().GetMin()[1];
		int rlon_i = floor((coord[0] - x0) / dx);
		int rlat_i = floor((coord[1] - y0) / dy);
		float w0 = (coord[0] - x0 - rlon_i * dx) / dx;
		float w1 = (coord[1] - y0 - rlat_i * dy) / dy;
		//cout << "Sample at "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<endl;
		//if(w0<0||w0>1||w1<0||w1>1)
		//cout << "  w0,w1 = " << w0 << "," << w1 << endl;
		if (w0 < 0) { cout << "ERROR: w0 = " << w0 << " < 0" << endl; w0 = 0; }
		else if (w0 > 1) { cout << "ERROR: w0 = " << w0 << " > 1" << endl; w0 = 1; }
		if (w1 < 0) { cout << "ERROR: w1 = " << w1 << " < 0" << endl; w1 = 0; }
		else if (w1 > 1) { cout << "ERROR: w1 = " << w1 << " > 1" << endl; w1 = 1; }
		//cout << "RlonRlat sample 2 for tmp0\n";
		TValueType tmp0 = SampleXYiHd(rlon_i, rlat_i, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i, coord[2]) * w0;
		//cout << "RlonRlat sample 2 for tmp1\n";
		TValueType tmp1 = SampleXYiHd(rlon_i, rlat_i + 1, coord[2]) * (1 - w0) + SampleXYiHd(rlon_i + 1, rlat_i + 1, coord[2]) * w0;
		//cout << "  tmp0 " << tmp0[0] << " " << tmp0[2] << "   tmp1 " << tmp1[0] << " " << tmp1[2] << endl;
		return tmp0 * (1 - w1) + tmp1 * w1;
	}
	TValueType SampleXYiHd(int rlon_i, int rlat_i, double h) const {
		//cout << "samplexyihd " << rlon_i << " " << rlat_i << " " << h << endl;
		/*
		seems	level	0		80
				height	22000	800-something
		*/
		//TODO figure out why this is necessary
		if (rlon_i < 0) {
			//cout << "rloni = " << rlon_i << endl;
			rlon_i = 0;
		}
		else if (rlon_i >= uvw->GetResolution()[0]) {
			//cout << "rloni = " << rlon_i << endl;
			rlon_i = uvw->GetResolution()[0] - 1;
		}
		if (rlat_i < 0) {
			//cout << "rlati = " << rlat_i << endl;
			rlat_i = 0;
		}
		else if (rlat_i >= uvw->GetResolution()[1]) {
			//cout << "rlati = " << rlat_i << endl;
			rlat_i = uvw->GetResolution()[1] - 1;
		}

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

	int lvl_top;
	float dx, dy;
};

typedef RlonRlatHField<Vec3f> RlonRlatHField_Vec3f;
typedef RlonRlatHField<float> RlonRlatHField_f;