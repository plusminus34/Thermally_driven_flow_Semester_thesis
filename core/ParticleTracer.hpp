#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"

//float, 3d, time invariant
Vec3f traceParticle(const RegVectorField3f& field, Vec3f start, float dt);
//time varying, field0 at time t0, field1 at time t1 = t0+dt, linear interpolation in between
Vec3f traceParticle(const RegVectorField3f& field0, const RegVectorField3f& field1, Vec3f start, float dt);
//as before, but t0 <= start_t < start_t+dt <= t1
Vec3f traceParticle(const RegVectorField3f& field0, float t0, const RegVectorField3f& field1, float t1, Vec3f start, float start_t, float dt);