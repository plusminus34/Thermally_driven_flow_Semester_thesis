#pragma once

#include <string>
#include <vector>
//#include "core/Math.hpp"
#include "core/ParticleTracer.hpp"
//#include "core/UVWReader.hpp"

class ImportantPart
{
public:
	std::vector<std::vector<Vec3f>> trajectories;// our "output"

	ImportantPart();

	// convert the DDHHMMSS string to int (number of seconds)
	int DDHHMMSSToDouble(std::string input) const;

	void doStuff();

private:
	std::vector<RegVectorField3f*> ringbuffer;
	std::string basefilename;

};