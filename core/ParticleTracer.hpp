#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"

//float, 3d, time invariant
Vec3f traceParticle(RegVectorField3f field, Vec3f start, float dt);
//float, 3d, field0 at time t0, field1 at time to+dt, linear interpolation in between
Vec3f traceParticle(RegVectorField3f field0, RegVectorField3f field1, Vec3f start, float dt);