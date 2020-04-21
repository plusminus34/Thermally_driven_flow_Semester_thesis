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
	TValueType traceParticle(const field_t& field0, double t0, const field_t& field1, double t1, const TValueType& start, double start_t, float dt) {
		const double iT = 1.0/ (t1 - t0);
		const double t_a = (start_t - t0) * iT;
		const double t_b = (start_t + 0.5*dt - t0) * iT;
		const double t_c = (start_t + dt - t0) * iT;
		TValueType k1 = field0.Sample(start) * (1 - t_a) + field1.Sample(start) * t_a;
		TValueType k2 = field0.Sample(start + k1 * dt*0.5) * (1 - t_b) + field1.Sample(start + k1 * dt*0.5) * t_b;
		TValueType k3 = field0.Sample(start + k2 * dt*0.5) * (1 - t_b) + field1.Sample(start + k2 * dt*0.5) * t_b;
		TValueType k4 = field0.Sample(start) * (1 - t_c) + field1.Sample(start) * t_c;
		return start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}
};