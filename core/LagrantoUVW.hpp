#pragma once

#include "RegularGrid.hpp"

#include <iostream>
#include <cmath>

using namespace std;

class LagrantoUVW : public ISampleField<Vec3f, 3>{
public:
	LagrantoUVW(string filename, RegScalarField3f* levelToHeight);
	~LagrantoUVW();

	virtual Vec3f Sample(const Vec3d& coord) const override;

private:
	RegScalarField3f* hhl;
	RegScalarField3f* U;
	RegScalarField3f* V;
	RegScalarField3f* W;

	int lvl_top, res_x, res_y;
	double x0, y0, x1, y1, dx, dy;

	float sampleHHL(double rlon_d, double rlat_d, int level, bool destagger) const;
};
