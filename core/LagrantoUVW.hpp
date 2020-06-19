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

	// Settings that probably make the results worse, but match what LAGRANTO does
	bool getZSamplingAlternative() const { return z_sampling_alternative; }
	void setZSamplingAlternative(bool value) { z_sampling_alternative = value; }
	bool isCorrectZW() const { return use_correct_zw; }
	void setCorrectZW(bool value) { use_correct_zw = value; }

private:
	RegScalarField3f* hhl;
	RegScalarField3f* U;
	RegScalarField3f* V;
	RegScalarField3f* W;

	int lvl_top, res_x, res_y;
	double x0, y0, x1, y1, dx, dy;

	float sampleHHL(double rlon_d, double rlat_d, int level, bool destagger) const;

	bool z_sampling_alternative = false;
	bool use_correct_zw = true;
};
