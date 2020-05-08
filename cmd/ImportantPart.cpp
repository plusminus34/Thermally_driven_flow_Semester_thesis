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

void ImportantPart::doStuff() {
	std::cout << "I swear by God that the ImportantPart is doing stuff!\n";
	//--------------------- beginning
	// get/set important variables
	double t0=0, t1=1;//beginning and end time
	double dt=0.1; int nSteps;//timestep size and number

	int ring_i = 0;//index of the current earliest element in the ringbuffer
	int ring_next = 1;
	int ring_third = 2;

	std::string filepath;//path to where files lie
	std::vector<std::string> files;//list of relevant files
	std::string constantsfile;//name of the file containing constants (mainly HHL)

	//ParticleTracer<Vec3f, 3> tracer;
	
	std::cout << "Ringbuffer ftw\n";
	// setup ringbuffer
	std::string filename = basefilename + "00000000.nc";//placeholder

	//ringbuffer[0] = UVWFromNCFile(filename);//from the first file timestep
	//ringbuffer[1] = UVWFromNCFile(filename);//actually from the next file timestep
	//ringbuffer[2] = UVWFromNCFile(filename);//actually from the file timestep after that

	// store a list of trajectories that haven't left the bounding box
	std::vector<std::vector<Vec3f>*> active_paths(trajectories.size());
	for (int i = 0; i < trajectories.size(); ++i) active_paths[i] = &trajectories[i];

	//--------------------- the important part
	std::cout << "THE IMPORTANT PART\n";
	// compute trajectories
	for (double t = t0; t < t1; t += dt) {
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
	std::cout << "fin\n";
	//for (int i = 0; i < 3; ++i) delete ringbuffer[i];
	//return paths somehow
}