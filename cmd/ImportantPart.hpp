#pragma once

#include <string>
#include <vector>
#include "core/ParticleTracer.hpp"

class ImportantPart
{
public:
	std::vector<std::vector<Vec3f>> trajectories;// our "output"

	ImportantPart();

	// Where trajectories are computed
	void doStuff();

	// Set parameters
	void setBaseFileName(string name) { basefilename = name; }
	void setFileOriginTime(int t0) { file_t0 = t0; }
	void setFileTimeStep(int dt) { file_dt = dt; }

	// convert the DDHHMMSS string to int (number of seconds)
	int DDHHMMSSToInt(std::string input) const;
	// and the other way around
	string IntToDDHHMMSS(int seconds) const;


private:
	vector<RegVectorField3f*> ringbuffer;
	string basefilename;
	int file_t0 = 0;
	int file_dt = 600;
};