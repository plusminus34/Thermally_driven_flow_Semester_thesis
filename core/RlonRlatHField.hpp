#pragma once

#include "RegularGrid.hpp"

using namespace std;

class RlonRlatHField : public ISampleField<Vec3f, 3>{
public:
	RlonRlatHField(RegVectorField3f* velocityField, RegScalarField3f* levelToHeight);
	RlonRlatHField() { delete uvw; }

	virtual Vec3f Sample(const Vec3d& coord) const override;
	Vec3f SampleXYiHd(int rlon_i, int rlat_i, double h) const;
private:
	RegVectorField3f* uvw;
	RegScalarField3f* hhl;

	Vec3i hhl_offset;

	int lvl_top;
	float dx, dy;
};