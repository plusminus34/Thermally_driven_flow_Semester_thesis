#pragma once

#include <iostream>
#include "ImportantPart.hpp"
#include "core/NetCDF.hpp"

ImportantPart::ImportantPart(): trajectories(1) {
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
	std::cout << "doStuff is deprecated, use computeTrajectoryData instead\n";
}

// Computes the data of td
// Assumes that the other members of td have been initialized
void ImportantPart::computeTrajectoryData(TrajectoryData& td)
{
	//--------------------- Initialization
	// see that data has correct sizes
	td.data.resize(td.varnames.size());
	for (int i = 0; i < td.varnames.size(); ++i) {
		td.data[i].resize(td.num_trajectories*td.points_per_trajectory);
	}

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

	// Helping: Map of level to acutal height in meters
	RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");

	// setup ringbuffer
	int file_i = 0;
	vector<RlonRlatHField<Vec3f>*> ringbuffer(3);
	ringbuffer[0] = new RlonRlatHField_Vec3f(UVWFromNCFile(files[0]), hhl);
	ringbuffer[1] = new RlonRlatHField_Vec3f(UVWFromNCFile(files[1]), hhl);
	ringbuffer[2] = new RlonRlatHField_Vec3f(UVWFromNCFile(files[2]), hhl);
	int ri0 = 0, ri1 = 1, ri2 = 2;

	// store current position of each trajectory
	vector<Vec3f> position(td.num_trajectories);
	for (int i = 0; i < position.size(); ++i) {
		position[i] = trajectories[i][0];
	}

	// domain is relevant for checking
	const BoundingBox3d bb = hhl->GetDomain();//TODO that's probably incorrect
	BoundingBox<Vec3f> bb_f(Vec3f(bb.GetMin()[0], bb.GetMin()[1], bb.GetMin()[1]), Vec3f(bb.GetMax()[0], bb.GetMax()[1], bb.GetMax()[1]));
	// and an extra variable to mark the final part where only 2 fields are used
	bool finalPart = false;

	// Figure out which variables are stored where
	int rlon_id, rlat_id, z_id, lon_id, lat_id;
	vector<int> other_vars;//TODO trace other variables
	for (int i = 0; i < td.varnames.size(); ++i) {
		if (td.varnames[i] == "rlon") rlon_id = i;
		else if (td.varnames[i] == "rlat") rlat_id = i;
		else if (td.varnames[i] == "z") z_id = i;
		else if (td.varnames[i] == "lon") lon_id = i;
		else if (td.varnames[i] == "lat") lat_id = i;
		else other_vars.push_back(i);
	}

	// need another two fields per other variable
	vector<vector<RlonRlatHField_f*>> other_fields(other_vars.size());
	for (int i = 0; i < other_fields.size(); ++i) {
		other_fields[i].resize(2);
		for (int j = 0; j < other_fields[i].size(); ++j) {
			RegScalarField3f* field = NetCDF::ImportScalarField3f(files[j], td.varnames[other_vars[i]], "rlon", "rlat", "level");
			other_fields[i][j] = new RlonRlatHField_f(field, hhl);
		}
	}
	int other_ri = 0;

	//oh, and write the initial points
	double lon, lat;
	for(int i=0;i<td.num_trajectories;++i){
		td.val(rlon_id, i, 0) = position[i][0];
		td.val(rlat_id, i, 0) = position[i][1];
		td.val(z_id, i, 0) = position[i][2];
		CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
		td.val(lon_id, i, 0) = lon;
		td.val(lat_id, i, 0) = lat;
		for (int j = 0; j < other_fields.size(); ++j) {
			Vec3f coord(position[i][0], position[i][1], position[i][2]);
			td.val(other_vars[j], i, 0) = other_fields[j][0]->Sample(coord);//TODO sample correctly
		}
	}

	//--------------------- Where particles are traced and paths filled
	// compute trajectories
	int step_i = 1;
	for (double t = td.time_begin; t < td.time_end; t += dt) {
		cout << "begin step " << step_i<<"   at time "<<t << endl;
		if (t >= file_t[file_i + 1]) {
			delete ringbuffer[file_i % 3];
			ri0 = ri1;
			ri1 = ri2;
			if (files.size() > file_i + 3) {
				ringbuffer[file_i % 3] = new RlonRlatHField_Vec3f(UVWFromNCFile(files[file_i + 3]), hhl);
				for (int i = 0; i < other_fields.size(); ++i) {
					delete other_fields[i][other_ri];
					cout << "Loading " << td.varnames[other_vars[i]] << "-field from " << files[file_i + 2];
					RegScalarField3f* field = NetCDF::ImportScalarField3f(files[file_i + 2], td.varnames[other_vars[i]], "rlon", "rlat", "level");
					other_fields[i][other_ri] = new RlonRlatHField_f(field, hhl);
				}
			}
			else {
				ringbuffer[file_i % 3] = nullptr;
				finalPart = true;
			}
			ri2 = file_i % 3;
			other_ri = file_i % 2;
			++file_i;
		}
		for (int i = 0; i < td.num_trajectories; ++i) {
			if (true) {//if (bb_f.Contains(position[i])) { TODO domain matters
				if (!finalPart) {
					position[i] = tracer.traceParticle(*ringbuffer[ri0], *ringbuffer[ri1], *ringbuffer[ri2], file_t[file_i], file_t[file_i + 2], position[i], t, dt);
				}
				else {
					position[i] = tracer.traceParticle(*ringbuffer[ri0], file_t[file_i], *ringbuffer[ri1], file_t[file_i + 1], position[i], t, dt);
				}
			}
			//TODO write td
			td.val(rlon_id, i, step_i) = position[i][0];
			td.val(rlat_id, i, step_i) = position[i][1];
			td.val(z_id, i, step_i) = position[i][2];
			CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
			td.val(lon_id, i, step_i) = lon;
			td.val(lat_id, i, step_i) = lat;
			for (int j = 0; j < other_fields.size(); ++j) {
				Vec3f coord(position[i][0], position[i][1], position[i][2]);
				td.val(other_vars[j], i, step_i) = other_fields[j][0]->Sample(coord);//TODO sample correctly
			}
		}
		++step_i;
	}
	//--------------------- the end
	for (int i = 0; i < 3; ++i) if (ringbuffer[i] != nullptr) delete ringbuffer[i];
	std::cout << "Trajectories have been computed\n";
}

void ImportantPart::helpWithStuff(RegVectorField3f* uvw)
{
	boundsmin = uvw->GetDomain().GetMin();
	boundsmax = uvw->GetDomain().GetMax();
	resolution = uvw->GetResolution();
}