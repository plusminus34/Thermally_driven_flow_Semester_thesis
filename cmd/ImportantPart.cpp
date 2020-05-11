#pragma once

#include <iostream>
#include "ImportantPart.hpp"
#include "core/NetCDF.hpp"

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

void ImportantPart::setNumberOfTimesteps(int num_steps) {
	nSteps = num_steps;
	dt = (end_t - start_t) / nSteps;
}

void ImportantPart::doStuff() {
	//--------------------- Initialization
	// name of the file containing constants (mainly HHL)
	string constantsfile = basefilename + IntToDDHHMMSS(file_t0) + 'c' + fileending;

	// Create list of relevant files
	vector<string> files;// the list
	vector<int> file_t;// the time of each file
	int file_i_0 = floor((start_t - file_t0) / file_dt);
	int file_i_1 = ceil((end_t - file_t0) / file_dt);
	for (int i = file_i_0; i <= file_i_1; ++i) {
		int t = file_t0 + i * file_dt;
		files.push_back(basefilename + IntToDDHHMMSS(t) + fileending);
		file_t.push_back(t);
	}

	// Yep, that's a particle tracer
	ParticleTracer<Vec3f, 3> tracer;

	// Helping: store a 1d-array mapping level1 to height (TODO is this correct?)
	cout << "Will import HHL from " << constantsfile << endl;
	RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");
	vector<float> lv_to_h(hhl->GetResolution()[2]);
	for (int i = 0; i < lv_to_h.size(); ++i) {
		Vec3i gridCoord(0, 0, i);
		lv_to_h[i] = hhl->GetVertexDataAt(gridCoord);
	}
	delete hhl;
	
	// setup ringbuffer
	int file_i = 0;
	ringbuffer[0] = UVWFromNCFile(files[0], lv_to_h);
	ringbuffer[1] = UVWFromNCFile(files[1], lv_to_h);
	ringbuffer[2] = UVWFromNCFile(files[2], lv_to_h);
	RegVectorField3f& field0 = *ringbuffer[0];
	RegVectorField3f& field1 = *ringbuffer[1];
	RegVectorField3f& field2 = *ringbuffer[2];

	// store a list of trajectories that haven't left the bounding box
	std::vector<std::vector<Vec3f>*> active_paths(trajectories.size());
	for (int i = 0; i < trajectories.size(); ++i) active_paths[i] = &trajectories[i];

	// domain is relevant for checking
	const BoundingBox3d bb = field0.GetDomain();

	//--------------------- Where particles are traced and paths filled
	// compute trajectories
	for (double t = start_t; t < end_t; t += dt) {
		cout << "time " << t << endl;
		if (t >= file_t[file_i + 1]) {
			cout << "delete\n";
			delete ringbuffer[file_i % 3];
			cout << "switch 0=1\n";
			field0 = field1;
			cout << "switch 1=2\n";
			field1 = field2;
			cout << "read 2 = ringbuffernew [" << (file_i % 3) << "] file " << files[file_i + 3] << endl;
			ringbuffer[file_i % 3] = UVWFromNCFile(files[file_i + 3], lv_to_h);
			cout << "set 2\n";
			field2 = *ringbuffer[file_i % 3];
			cout << "++\n";
			++file_i;
			cout << "switched\n";
		}
		for (int i = 0; i < trajectories.size(); ++i) {
			cout << "  trajectories: " << trajectories.size() << " am at " << i << endl;
			//cout << "  traj[i]: " << trajectories[i].size() << " am at " << trajectories[i].size() - 1 << endl;
			Vec3f pos_f = trajectories[i][trajectories[i].size() - 1];
			//cout << "  pos_f " << pos_f[0] << " " << pos_f[1] << " " << pos_f[2];
			pos_f = tracer.traceParticle(field0, field1, field2, file_t[file_i], file_t[file_i + 2], pos_f, t, dt);
			trajectories[i].push_back(pos_f);
		}
		// loop over active paths
		/*
		for (int i = 0; i < active_paths.size(); ++i) {
			//compute next position
			vector<Vec3f>*  path = active_paths[i];
			Vec3f pos_f = (*path)[path->size() - 1];
			pos_f = tracer.traceParticle(field0, field1, field2, file_t[file_i], file_t[file_i + 2], pos_f, t, dt);
			path->push_back(pos_f);

			//TODO changing active_paths while iterating over it is not nice
			// remove from active paths pos_f is outside the domain
			Vec3d pos_d(pos_f[0], pos_f[1], pos_f[2]);
			if (!bb.Contains(pos_d)) {
				swap(active_paths[i], active_paths[active_paths.size() - 1]);
				--i;
				active_paths.pop_back();
			}
		}
		*/
	}
	//--------------------- the end
	for (int i = 0; i < 3; ++i) if (ringbuffer[i] != nullptr)delete ringbuffer[i];
	std::cout << "Paths are finished\n";
}