#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"

template<typename TValueType, size_t TDimensions>
class ParticleTracer {
public:
	using field_t = ISampleField<TValueType, TDimensions>;

	//time invariant
	TValueType traceParticle(const field_t& field, const TValueType& start, double dt) {
		TValueType k1 = field.Sample(start);
		TValueType k2 = field.Sample(start + k1 * dt*0.5);
		TValueType k3 = field.Sample(start + k2 * dt*0.5);
		TValueType k4 = field.Sample(start + k3 * dt);
		return start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}
	//between field0 at t0 and field1 at t1
	//TValueType traceParticle(const field_t& field0, double t0, const field_t& field1, double t1, const TValueType& start, double start_t, float dt);
};

//time varying, field0 at time t0, field1 at time t1 = t0+dt, linear interpolation in between
Vec3f traceParticle(const RegVectorField3f& field0, const RegVectorField3f& field1, Vec3f start, float dt);
//as before, but t0 <= start_t < start_t+dt <= t1
Vec3f traceParticle(const RegVectorField3f& field0, float t0, const RegVectorField3f& field1, float t1, Vec3f start, float start_t, float dt);