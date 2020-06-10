#pragma once

#include <string>
#include <vector>
#include "UVWReader.hpp"
#include "core/TimeRelatedFields.hpp"

template<typename TValueType, size_t TDimensions>
class ParticleTracer {
public:
	using field_t = ISampleField<TValueType, TDimensions>;
	using time_field_t = TimeRelatedFields<TValueType, TDimensions>;

	//--- time invariant case ---
	TValueType traceParticle(const field_t& field, const TValueType& start, double dt) {
		TValueType k1 = field.Sample(start);
		TValueType k2 = field.Sample(start + k1 * dt*0.5);
		TValueType k3 = field.Sample(start + k2 * dt*0.5);
		TValueType k4 = field.Sample(start + k3 * dt);
		return start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}

	//--- Variants with more than one field ---

	//between field0 at t0 and field1 at t1
	TValueType traceParticle(const field_t& field0, double t0, const field_t& field1, double t1, const TValueType& start, double start_t, float dt) {
		const double iT = 1.0/ (t1 - t0);
		const double t_a = (start_t - t0) * iT;
		const double t_b = (start_t + 0.5*dt - t0) * iT;
		const double t_c = (start_t + dt - t0) * iT;
		TValueType k1 = field0.Sample(start) * (1 - t_a) + field1.Sample(start) * t_a;
		TValueType k2 = field0.Sample(start + k1 * dt*0.5) * (1 - t_b) + field1.Sample(start + k1 * dt*0.5) * t_b;
		TValueType k3 = field0.Sample(start + k2 * dt*0.5) * (1 - t_b) + field1.Sample(start + k2 * dt*0.5) * t_b;
		TValueType k4 = field0.Sample(start + k3 * dt) * (1 - t_c) + field1.Sample(start + k3 * dt) * t_c;
		return start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}

	// field0 at field_t0, field2 at field_t2, field1 at their midpoint; start between the first two fields, end possibly between the second and third field
	TValueType traceParticle(const field_t& field0, const field_t& field1, const field_t& field2, double field_t0, double field_t2, const TValueType& start, double start_t, float dt) {
		assert(start_t >= field_t0 && start_t <= 0.5*(field_t0 + field_t2));
		const double field_t1 = 0.5*(field_t0 + field_t2);
		const double iT = 2.0 / (field_t2 - field_t0);
		double alpha = (start_t - field_t0) * iT;
		TValueType k1 = field0.Sample(start) * (1 - alpha) + field1.Sample(start) * alpha;
		TValueType k2, k3, k4;
		const double th = start_t + 0.5*dt;
		if (th < field_t1) {
			alpha = (th - field_t0) * iT;
			k2 = field0.Sample(start + k1 * dt*0.5) * (1 - alpha) + field1.Sample(start + k1 * dt*0.5) * alpha;
			k3 = field0.Sample(start + k2 * dt*0.5) * (1 - alpha) + field1.Sample(start + k2 * dt*0.5) * alpha;
		}
		else {
			alpha = (th - field_t1) * iT;
			k2 = field1.Sample(start + k1 * dt*0.5) * (1 - alpha) + field2.Sample(start + k1 * dt*0.5) * alpha;
			k3 = field1.Sample(start + k2 * dt*0.5) * (1 - alpha) + field2.Sample(start + k2 * dt*0.5) * alpha;
		}
		const double t1 = start_t + dt;
		if (t1 < field_t1) {
			alpha = (t1 - field_t0) * iT;
			k4 = field0.Sample(start + k3 * dt) * (1 - alpha) + field1.Sample(start + k3 * dt) * alpha;
		}
		else {
			alpha = (t1 - field_t1) * iT;
			k4 = field1.Sample(start + k3 * dt) * (1 - alpha) + field2.Sample(start + k3 * dt) * alpha;
		}
		return start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}

	// between field0 at t0 and field1 at t1 using Explicit Euler instead of Runge-Kutta
	TValueType traceParticleEE(const field_t& field0, double t0, const field_t& field1, double t1, const TValueType& start, double start_t, float dt) {
		const double iT = 1.0 / (t1 - t0);
		const double t_a = (start_t - t0) * iT;
		TValueType v = field0.Sample(start) * (1 - t_a) + field1.Sample(start) * t_a;
		return start + v * dt;
	}

	// Lagranto actually uses an iterative Euler scheme (nIter=2 for explicit trapezoidal)
	TValueType traceParticleIterativeEuler(const field_t& field0, double t0, const field_t& field1, double t1, const TValueType& start, double start_t, float dt, int nIter) {
		const double iT = 1.0 / (t1 - t0);
		const double t_a = (start_t - t0) * iT;
		const double t_b = (start_t + dt - t0) * iT;
		TValueType v0 = field0.Sample(start) * (1 - t_a) + field1.Sample(start) * t_a;
		TValueType res = start;
		for (int i = 0; i < nIter; ++i) {
			TValueType v1 = field0.Sample(res) * (1 - t_b) + field1.Sample(res) * t_b;
			res = start + (v0 + v1)*dt*0.5;
		}
		return res;
	}

	//Variants using TimeRelatedFields

	TValueType traceParticleRungeKutta(const time_field_t& field, const TValueType& start_pos, double start_t, float dt) {
		TValueType k1 = field.Sample(start_pos, start_t);
		TValueType k2 = field.Sample(start_pos + k1 * dt*0.5, start_t + dt * 0.5);
		TValueType k3 = field.Sample(start_pos + k2 * dt*0.5, start_t + dt * 0.5);
		TValueType k4 = field.Sample(start_pos + k3 * dt, start_t + dt);
		return start_pos + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
	}

	TValueType traceParticleIterativeEuler(const time_field_t& field, const TValueType& start_pos, double start_t, float dt, int nIter) {
		TValueType v0 = field.Sample(start_pos, start_t);
		TValueType res = start_pos;
		for (int i = 0; i < nIter; ++i) {
			TValueType v1 = field.Sample(res, start_t + dt);
			res = start_pos + (v0 + v1)*dt*0.5;
		}
		return res;
	}
};