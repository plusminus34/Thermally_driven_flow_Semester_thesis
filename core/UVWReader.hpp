#pragma once

#include <string>
#include <vector>
#include "RlonRlatHField.hpp"
#include "CoordinateTransform.hpp"

using namespace std;

RegVectorField3f* UVWFromVTIFile(string filename = "UVW.vti");

RegVectorField3f* UVWFromNCFile(string filename);

//read UVW Lagranto-style
void SeparateUVWFromNCFile(string filename, RegScalarField3f* &U, RegScalarField3f* &V, RegScalarField3f* &W);