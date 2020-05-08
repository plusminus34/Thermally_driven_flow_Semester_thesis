#pragma once

#include <iostream>
#include "ImportantPart.hpp"

ImportantPart::ImportantPart(): ringbuffer(3), trajectories(1) {
	basefilename = "lfff";//followed by DDHHMMSS, possibly 'c', and .nc
}

int ImportantPart::DDHHMMSSToInt(std::string input) const
{
	assert(input.size() >= 8);
	int days = stoi(input.substr(0, 2));
	int hours = stoi(input.substr(2, 2));
	int minutes = stoi(input.substr(4, 2));
	int seconds = stoi(input.substr(6, 2));
	return days * 86400 + hours * 3600 + minutes * 60 + seconds;
}

string ImportantPart::IntToDDHHMMSS(int seconds) const {
	int days = seconds / 86400;
	seconds -= 86400 * days;
	int hours = seconds / 3600;
	seconds -= 3600 * hours;
	int minutes = seconds / 60;
	seconds -= 60 * minutes;
	if (days > 99) {
		days = 99;
		cout << "Warning: DDHHMMSS with more than 99 days\n";
	}
	string res = "";
	if (days < 10) res += "0"; res += to_string(days);
	if (hours < 10) res += "0"; res += to_string(hours);
	if (minutes < 10) res += "0"; res += to_string(minutes);
	if (seconds < 10) res += "0"; res += to_string(seconds);
	return res;

}

void ImportantPart::setTimestep(double timestep) {
	dt = timestep;
	nSteps = ceil((end_t - start_t) / dt);
}

void ImportantPart::SetNumberOfTimesteps(int num_steps) {
	nSteps = num_steps;
	dt = (end_t - start_t) / nSteps;
}

void ImportantPart::doStuff() {
	std::cout << "I swear by God that the ImportantPart is doing stuff!\n";
	//--------------------- Initialization
	// name of the file containing constants (mainly HHL)
	string constantsfile = basefilename + IntToDDHHMMSS(file_t0) + 'c' + fileending;

	// Create list of relevant files
	vector<string> files;// the list
	vector<int> file_t;// the time of each file
	int file_i_0 = floor((start_t - file_t0) / file_dt);
	int file_i_1 = ceil((end_t - file_t0) / file_dt);
	for (int i = file_i_0; i < file_i_1; ++i) {
		int t = file_t0 + i * file_dt;
		files.push_back(basefilename + IntToDDHHMMSS(t) + fileending);
		file_t.push_back(t);
	}

	// Yep, that's a particle tracer
	ParticleTracer<Vec3f, 3> tracer;
	
	// setup ringbuffer
	if (files.size() < 3) {
		cout << "Sorry, at least 3 files are currently required.\n";
		return;
	}
	int file_i = 0;
	ringbuffer[0] = UVWFromNCFile(files[0]);
	ringbuffer[1] = UVWFromNCFile(files[1]);
	ringbuffer[2] = UVWFromNCFile(files[2]);

	// store a list of trajectories that haven't left the bounding box
	std::vector<std::vector<Vec3f>*> active_paths(trajectories.size());
	for (int i = 0; i < trajectories.size(); ++i) active_paths[i] = &trajectories[i];

	//--------------------- Where particles are traced and paths filled
	std::cout << "THE IMPORTANT PART\n";
	// compute trajectories
	for (double t = start_t; t < end_t; t += dt) {
		//if ( t reaches ring_next) {
		//	delete ringbuffer[ring_i];
		//	ring_next = ring_third;
		//	ring_i = ring_next;
		//	ring_third = (ring_next + 1) % 3;
		//}
		//loop over active paths
		for (int i = 0; i < active_paths.size(); ++i) {
			//compute next position
			std::vector<Vec3f>* path = active_paths[i];
			//Vec3f* pos_fp = &(*path)[path->size() - 1];//maybe not use so many pointers?
			//Vec3d pos_d((*pos_fp)[0], (*pos_fp)[1], (*pos_fp)[2]);
			//TODO path->push_back(trace(???)); 
			//make inactive if it leaves the bounds of the field
		}
	}
	//--------------------- the end
	for (int i = 0; i < 3; ++i) delete ringbuffer[i];
	std::cout << "Paths are finished\n";
}