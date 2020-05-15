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
	// Start and end time of tracing
	void setTimeBoundaries(double t0, double t1) { start_t = t0; end_t = t1; }
	// Define step size
	void setTimestep(double timestep);
	void setNumberOfTimesteps(int num_steps);
	// Information about which files to use
	void setBaseFileName(string name) { basefilename = name; }
	void setFileOriginTime(int t0) { file_t0 = t0; }
	void setFileTimestep(int dt) { file_dt = dt; }

	// for resizing hhl
	void helpWithStuff(RegVectorField3f* uvw);

	// convert the DDHHMMSS string to int (number of seconds)
	int DDHHMMSSToInt(std::string input) const;
	// and the other way around
	string IntToDDHHMMSS(int seconds) const;


private:
	string basefilename;
	string fileending = ".nc";
	int file_t0 = 0;
	int file_dt = 600;
	double start_t = 0.0;
	double end_t = 1200.0;
	double dt = 1; int nSteps = 1200;

	//kinda stupid way to find out how big hhl should be
	Vec3d boundsmin, boundsmax;
	Vec3i resolution;
};