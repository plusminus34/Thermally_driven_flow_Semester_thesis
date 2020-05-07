#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"
#include "CoordinateTransform.hpp"

using namespace std;

RegVectorField3f* UVWFromVTIFile(string filename = "UVW.vti");

RegVectorField3f* UVWFromNCFile(string filename);