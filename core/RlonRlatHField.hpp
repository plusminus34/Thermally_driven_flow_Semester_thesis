#pragma once

#include "RegularGrid.hpp"

using namespace std;

class RlonRlatHField : public ISampleField<Vec3f, 3>{
public:
	RlonRlatHField(RegVectorField3f* velocityField, RegScalarField3f* levelToHeight);
	RlonRlatHField() { delete uvw; }

	virtual Vec3f Sample(const Vec3d& coord) const override;
	Vec3f SampleXYiHd(int rlon_i, int rlat_i, double h) const;
	RegScalarField3f* hhl;
private:
	RegVectorField3f* uvw;
	int lvl_top;
	float dx, dy;
};