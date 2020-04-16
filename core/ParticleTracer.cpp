#include "ParticleTracer.hpp"



Vec3f traceParticle(RegVectorField3f field, Vec3f start, float dt) {
	Vec3f k1 = field.Sample(start);
	Vec3f k2 = field.Sample(start + k1 * dt*0.5f);
	Vec3f k3 = field.Sample(start + k2 * dt*0.5f);
	Vec3f k4 = field.Sample(start + k3 * dt);
	return Vec3f(start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.f));
}

Vec3f traceParticle(RegVectorField3f field0, RegVectorField3f field1, Vec3f start, float dt){
	Vec3f k1 = field0.Sample(start);
	Vec3f k2 = (field0.Sample(start + k1 * dt*0.5f) + field1.Sample(start + k1 * dt*0.5f)) *0.5f;
	Vec3f k3 = (field0.Sample(start + k2 * dt*0.5f) + field1.Sample(start + k2 * dt*0.5f)) *0.5f;
	Vec3f k4 = field1.Sample(start + k3 * dt);
	return Vec3f(start + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.f));
}