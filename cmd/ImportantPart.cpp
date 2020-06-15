#pragma once

#include <iostream>
#include "ImportantPart.hpp"
#include "core/NetCDF.hpp"
#include "core/LagrantoUVW.hpp"

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
	if (dt > 0)
		nSteps = ceil((end_t - start_t) / dt);
	else if (dt < 0)
		nSteps = floor((end_t - start_t) / dt);
}

void ImportantPart::setNumberOfTimesteps(int num_steps) {
	nSteps = num_steps;
	if(nSteps > 0) dt = (end_t - start_t) / nSteps;
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
	td.times.resize(td.points_per_trajectory);

	// name of the file containing constants (mainly HHL)
	string constantsfile = basefilename + IntToDDHHMMSS(file_t0) + 'c' + fileending;

	// check if tracing forward or backward
	bool backward = false;
	if (start_t > end_t) backward = true;
	if (backward && dt > 0) dt *= -1;
	else if (!backward && dt < 0) dt *= -1;

	// Create list of relevant files
	vector<string> files;// the list
	vector<int> file_t;// the time of each file
	int file_i_min = floor((std::min(start_t, end_t) - file_t0) / file_dt);
	int file_i_max = ceil((std::max(start_t, end_t) - file_t0) / file_dt);
	int nFiles = file_i_max - file_i_min + 1;
	files.resize(nFiles);
	file_t.resize(nFiles);
	for (int i = file_i_min; i <= file_i_max; ++i) {
		int t = file_t0 + i * file_dt;
		int j = backward ? file_i_max - i : i - file_i_min;
		files[j] = basefilename + IntToDDHHMMSS(t) + fileending;
		if (time_invariant) files[j] = basefilename + IntToDDHHMMSS(file_t0) + fileending;//TODO this is a quick and dirty way to do this
		file_t[j] = t;
	}

	// Yep, that's a particle tracer
	ParticleTracer<Vec3f, 3> tracer;

	// Helping: Map of level to actual height in meters
	RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");

	// setup ringbuffer
	int file_i = 0;
	TimeRelatedFields<Vec3f, 3> UVW(3);
	for (int i = 0; i < 3; ++i) {
		if (!use_lagranto_uvw)
			UVW.InsertNextField(new RlonRlatHField_Vec3f(UVWFromNCFile(files[i]), hhl), file_t[i]);
		else
			UVW.InsertNextField(new LagrantoUVW(files[i], hhl), file_t[i]);
	}
	if (backward) UVW.setBackward(true);

	// store current position of each trajectory
	vector<Vec3f> position(td.num_trajectories);
	for (int i = 0; i < position.size(); ++i) {
		position[i] = trajectories[i][0];//TODO trajectories is no longer used otherwise
	}

	// domain is relevant for checking active trajectories
	const BoundingBox3d bb = hhl->GetDomain();
	const float rlon_min = bb.GetMin()[0];
	const float rlon_max = bb.GetMax()[0];
	const float rlat_min = bb.GetMin()[1];
	const float rlat_max = bb.GetMax()[1];

	// and an extra variable to mark the final part where only 2 fields are used
	bool finalPart = false;

	// Figure out which variables are stored where
	int rlon_id, rlat_id, z_id, lon_id, lat_id;
	vector<int> other_vars;
	for (int i = 0; i < td.varnames.size(); ++i) {
		if (td.varnames[i] == "rlon") rlon_id = i;
		else if (td.varnames[i] == "rlat") rlat_id = i;
		else if (td.varnames[i] == "z") z_id = i;
		else if (td.varnames[i] == "lon") lon_id = i;
		else if (td.varnames[i] == "lat") lat_id = i;
		else other_vars.push_back(i);
	}

	//oh, and write the initial points
	for (int i = 0; i < td.num_trajectories; ++i) {
		td.val(rlon_id, i, 0) = position[i][0];
		td.val(rlat_id, i, 0) = position[i][1];
		td.val(z_id, i, 0) = position[i][2];
		double lon, lat;
		CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
		td.val(lon_id, i, 0) = lon;
		td.val(lat_id, i, 0) = lat;
		td.times[0] = td.time_begin;
	}

	//--------------------- Where particles are traced and paths filled
	// compute trajectories
	int nSteps = ceil((td.time_end - td.time_begin) / dt);
	double t = td.time_begin;
	for (size_t step_i = 1; step_i < td.points_per_trajectory; ++step_i) {
		cout << "Begin step " << step_i<<"   at time "<< t << endl;
		double dt_real = step_i == td.points_per_trajectory - 1 ? td.time_end - t : dt;
		// insert new files if necessary
		bool need_new_file = backward ? (t <= file_t[file_i + 1]) : t >= file_t[file_i + 1];
		if (need_new_file) {
			if (files.size() > file_i + 3) {
				if (!use_lagranto_uvw)
					UVW.InsertNextField(new RlonRlatHField_Vec3f(UVWFromNCFile(files[file_i + 3]), hhl), file_t[file_i + 3]);
				else
					UVW.InsertNextField(new LagrantoUVW(files[file_i + 3], hhl), file_t[file_i + 3]);
			}
			++file_i;
		}
		// propagate all active trajectories by one timestep
		#pragma omp parallel for
		for (int i = 0; i < td.num_trajectories; ++i) {
			if (position[i][0]>=rlon_min && position[i][0] <= rlon_max
				&& position[i][1] >= rlat_min && position[i][1] <= rlat_max) {
				if (integrator==0)
					position[i] = tracer.traceParticleRungeKutta(UVW, position[i], t, dt_real);
				else
					position[i] = tracer.traceParticleIterativeEuler(UVW, position[i], t, dt_real, 3);
			}
			const size_t data_i = step_i * td.num_trajectories + i;
			td.data[rlon_id][data_i] = position[i][0];
			td.data[rlat_id][data_i] = position[i][1];
			td.data[z_id][data_i] = position[i][2];
			double lon, lat;
			CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
			td.data[lon_id][data_i] = lon;
			td.data[lat_id][data_i] = lat;
		}
		t += dt_real;
		td.times[step_i] = t;
	}
	//--------------------- get scalar extra variables
	for (int var_i = 0; var_i < other_vars.size(); ++var_i) {
		cout << "now sampling " << td.varnames[other_vars[var_i]] << endl;
		// initialize current fields and variables
		TimeRelatedFields<float, 3 > fields(2);
		for (int j = 0; j < 2; ++j) {
			RegScalarField3f* field = NetCDF::ImportScalarField3f(files[j], td.varnames[other_vars[var_i]], "rlon", "rlat", "level");
			fields.InsertNextField(field, file_t[var_i]);
		}
		if (backward) fields.setBackward(true);
		int file_i = 0;
		int id = td.get_var_id(td.varnames[other_vars[var_i]]);

		// go over all points
		for (size_t step_i = 0; step_i < td.points_per_trajectory; ++step_i) {
			double t = td.times[step_i];
			bool need_new_file = backward ? (t < file_t[file_i + 1]) : t > file_t[file_i + 1];
			if (need_new_file) {
				RegScalarField3f* field = NetCDF::ImportScalarField3f(files[file_i + 2], td.varnames[other_vars[var_i]], "rlon", "rlat", "level");
				fields.InsertNextField(field, file_t[file_i + 2]);
				++file_i;
			}

			for (size_t traj_i=0; traj_i < td.num_trajectories; ++traj_i) {
				Vec3d coord(td.val(rlon_id,traj_i,step_i), td.val(rlat_id, traj_i, step_i), td.val(z_id, traj_i, step_i));
				td.val(id, traj_i, step_i) = fields.Sample(coord, t);
				//if(traj_i == 0)
				//cout << ": tdval "<< td.val(id, traj_i, step_i) <<endl;
			}
		}

	}
	//--------------------- the end
	std::cout << "Trajectories have been computed\n";
}

void ImportantPart::computeTrajectoryDataTEST(TrajectoryData & td)
{
	//--------------------- Initialization
	// see that data has correct sizes
	td.data.resize(td.varnames.size());
	for (int i = 0; i < td.varnames.size(); ++i) {
		td.data[i].resize(td.num_trajectories*td.points_per_trajectory);
	}

	// name of the file containing constants (read: HHL)
	string constantsfile = basefilename + IntToDDHHMMSS(file_t0) + 'c' + fileending;

	ParticleTracer<Vec3f, 3> tracer;

	// Helping: Map of level to acutal height in meters
	RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");

	//This is the only field used for now
	LagrantoUVW luv(basefilename + IntToDDHHMMSS(file_t0) + fileending, hhl);

	// store current position of each trajectory
	vector<Vec3f> position(td.num_trajectories);
	for (int i = 0; i < position.size(); ++i) {
		position[i] = trajectories[i][0];
	}

	// Figure out which variables are stored where
	int rlon_id, rlat_id, z_id, lon_id, lat_id;
	for (int i = 0; i < td.varnames.size(); ++i) {
		if (td.varnames[i] == "rlon") rlon_id = i;
		else if (td.varnames[i] == "rlat") rlat_id = i;
		else if (td.varnames[i] == "z") z_id = i;
		else if (td.varnames[i] == "lon") lon_id = i;
		else if (td.varnames[i] == "lat") lat_id = i;
	}

	//oh, and write the initial points
	double lon, lat;
	for (int i = 0; i < td.num_trajectories; ++i) {
		td.val(rlon_id, i, 0) = position[i][0];
		td.val(rlat_id, i, 0) = position[i][1];
		td.val(z_id, i, 0) = position[i][2];
		CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
		td.val(lon_id, i, 0) = lon;
		td.val(lat_id, i, 0) = lat;
	}

	for (int traj = 0; traj < position.size();++traj) {
		//cout << "coord "<<traj<<" begin: " << position[traj][0] << " " << position[traj][1] << " " << position[traj][2] << endl;
		int step_i = 1;
		for (double t = td.time_begin; t < td.time_end; t += dt) {
			position[traj] = tracer.traceParticleIterativeEuler(luv, start_t, luv, end_t, position[traj], t, dt, 3);

			td.val(rlon_id, traj, step_i) = position[traj][0];
			td.val(rlat_id, traj, step_i) = position[traj][1];
			td.val(z_id, traj, step_i) = position[traj][2];
			CoordinateTransform::RlatRlonToLatLon(position[traj][1], position[traj][0], lat, lon);
			td.val(lon_id, traj, step_i) = lon;
			td.val(lat_id, traj, step_i) = lat;
			//cout << "coord t " << t + dt << ": " << position[traj][0] << " " << position[traj][1] << " " << position[traj][2] << endl;
			//cout << "coord t " << t + dt << " lonlat: " << lon << " " << lat << " " << position[traj][2] << endl;
			++step_i;
		}
	}
}
