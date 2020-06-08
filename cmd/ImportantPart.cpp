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

float hhl_gdata(RegScalarField3f* hhl, Vec3i gc) {
	return 0.5f*(hhl->GetVertexDataAt(gc) + hhl->GetVertexDataAt(gc + Vec3i(0, 0, 1)));
}

Vec3f ImportantPart::sampleUVWTEST(Vec3d coord, RegScalarField3f * U, RegScalarField3f * V, RegScalarField3f * W, RegScalarField3f* hhl)
{
	Vec3f uvw(0, 0, 0);


// constants
const double x0 = hhl->GetDomain().GetMin()[0];
const double x1 = hhl->GetDomain().GetMax()[0];
const double y0 = hhl->GetDomain().GetMin()[1];
const double y1 = hhl->GetDomain().GetMax()[1];
const int res_x = hhl->GetResolution()[0];
const int res_y = hhl->GetResolution()[1];
const double dx = hhl->GetVoxelSize()[0];
const double dy = hhl->GetVoxelSize()[1];

	double rlon_d = (coord[0] - x0) / (x1 - x0) * (res_x - 1);
	double rlat_d = (coord[1] - y0) / (y1 - y0) * (res_y - 1);
	int rlon_i = std::min(std::max(0, (int)std::floor(rlon_d)), res_x - 1);
	int rlat_i = std::min(std::max(0, (int)std::floor(rlat_d)), res_y - 1);
	double w0 = rlon_d - rlon_i;
	double w1 = rlat_d - rlat_i;

	//cout << "point " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
	//cout << " has rlon_i, rlat_i " << rlon_i << ", " << rlat_i << endl;
	//cout << " and rlon_d, rlat_d " << rlon_d << ", " << rlat_d << endl;

	double h = coord[2];
	double lvl_d = 0;
	{
		int lvl_0 = 0;
		int lvl_1 = 79;
		if (h > hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lvl_0))) lvl_d = -1;
		if (h < hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lvl_1))) lvl_d = -1;

		while (lvl_1 > lvl_0 + 1) {
			int half = (lvl_1 + lvl_0) / 2;
			//cout << "half " << half << " between " << lvl_1 << " and " << lvl_0 << endl;
			//if (hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, half)) < h)
			if (hhl_gdata(hhl,Vec3i(rlon_i, rlat_i, half)) < h)
				lvl_1 = half;
			else lvl_0 = half;
		}
		//const float h1 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lvl_1));
		//const float h0 = hhl->GetVertexDataAt(Vec3i(rlon_i, rlat_i, lvl_0));
		const float h1 = hhl_gdata(hhl,Vec3i(rlon_i, rlat_i, lvl_1));
		const float h0 = hhl_gdata(hhl,Vec3i(rlon_i, rlat_i, lvl_0));
		
		// and then sample uvw
		float alpha = (h - h0) / (h1 - h0);
		Vec3i smp0(rlon_i, rlat_i, lvl_0);
		Vec3i smp1(rlon_i, rlat_i, lvl_1);
		//cout << "localized: h " << h << " belongs between h0 and h1 " << h0 << " " << h1 << endl;
		//cout << "   which means level1s " << lvl_0 << " " << lvl_1 << endl;
		double res = 0;
		if (lvl_0 > 0) res += hhl->GetVertexDataAt(smp0)*(1 - alpha);
		if (lvl_1 < hhl->GetResolution()[2]) res += hhl->GetVertexDataAt(smp1)*alpha;
		lvl_d = lvl_0 + alpha;
	}
	//cout << "lvl_d is " << lvl_d << endl;

	double wx = w0;
	double wy = w1;

	double lvl_d_0 = lvl_d;
	for (int count = 0; count <= 3; ++count) {
		lvl_d = lvl_d_0 - count * 0.5;
		int lvl_i = lvl_d;
		double wz = lvl_d - lvl_i;

		Vec3i gcs[8];
		double ws[8];
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				for (int k = 0; k < 2; ++k) {
					ws[4 * i + 2 * j + k] = 0;
				}
			}
		}
		ws[0] = (1 - wx) * (1 - wy) * (1 - wz); gcs[0] = Vec3i(rlon_i, rlat_i, lvl_i);
		ws[1] = wx * (1 - wy) * (1 - wz); gcs[1] = Vec3i(rlon_i + 1, rlat_i, lvl_i);
		ws[2] = (1 - wx) * wy * (1 - wz); gcs[2] = Vec3i(rlon_i, rlat_i + 1, lvl_i);
		ws[3] = wx * wy * (1 - wz); gcs[3] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i);
		ws[4] = (1 - wx) * (1 - wy) * wz; gcs[4] = Vec3i(rlon_i, rlat_i, lvl_i + 1);
		ws[5] = wx * (1 - wy) * wz; gcs[5] = Vec3i(rlon_i + 1, rlat_i, lvl_i + 1);
		ws[6] = (1 - wx) * wy * wz; gcs[6] = Vec3i(rlon_i, rlat_i + 1, lvl_i + 1);
		ws[7] = wx * wy * wz; gcs[7] = Vec3i(rlon_i + 1, rlat_i + 1, lvl_i + 1);

		float uu = 0, vv = 0, ww = 0;
		for (int i = 0; i < 8; ++i) {
			if (count == 1) {
				uvw[0] += U->GetVertexDataAt(gcs[i])*ws[i];
				uvw[1] += V->GetVertexDataAt(gcs[i])*ws[i];
			}
			else if (count == 2) {
				uvw[2] += W->GetVertexDataAt(gcs[i])*ws[i];
			}
			uu += U->GetVertexDataAt(gcs[i])*ws[i];
			vv += V->GetVertexDataAt(gcs[i])*ws[i];
			ww += W->GetVertexDataAt(gcs[i])*ws[i];
		}

		//cout << "uvw found at rlon_d " << rlon_d << " rlat_d" << rlat_d << " lvl_d " << lvl_d << endl;
		//cout << "   is " << uu << " " << vv << " " << ww << endl;
		//if (count == 1)cout << "UV from lvl_d " << lvl_d << ":   " << uu << " " << vv << endl;
		//if (count == 2) cout << " W from lvl_d " << lvl_d << ":   " << ww << endl;
	}

	/*
	// rescale to rlon/s rlat/s m/s
	double lon, lat;
	CoordinateTransform::RlatRlonToLatLon(coord[1], coord[0], lat, lon);
	double dlon, dlat;
	CoordinateTransform::degreeLengthsSimple(lat, dlat, dlon);
	if (dlon != 0) uvw[0] /= dlon;
	else uvw[0] = 0;
	uvw[1] /= dlat;
	*/
	//cout << "uvw unmodified " << uvw[0] << " " << uvw[1] << " " << uvw[2] << endl;
	float deltay = 111200;
	uvw[1] /= deltay;
	uvw[0] /= deltay * cos(coord[1] * 3.1415926535 / 180.0);

	return uvw;
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
			float val_0 = other_fields[j][other_ri]->Sample(coord);
			float val_1 = other_fields[j][(other_ri + 1) % 2]->Sample(coord);
			float alpha = (td.time_begin - file_t[file_i]) / (file_t[file_i + 1] - file_t[file_i]);
			td.val(other_vars[j], i, 0) = (1 - alpha)*val_0 + alpha * val_1;
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
			}
			else {
				ringbuffer[file_i % 3] = nullptr;
				finalPart = true;
			}
			ri2 = file_i % 3;
			++file_i;
		}
		if (t + dt > file_t[file_i + 1]) {
			for (int i = 0; i < other_fields.size(); ++i) {
				delete other_fields[i][other_ri];
				cout << "Loading " << td.varnames[other_vars[i]] << "-field from " << files[file_i + 2] << endl;
				RegScalarField3f* field = NetCDF::ImportScalarField3f(files[file_i + 2], td.varnames[other_vars[i]], "rlon", "rlat", "level");
				other_fields[i][other_ri] = new RlonRlatHField_f(field, hhl);
			}
			other_ri = file_i % 2;
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
			td.val(rlon_id, i, step_i) = position[i][0];
			td.val(rlat_id, i, step_i) = position[i][1];
			td.val(z_id, i, step_i) = position[i][2];
			CoordinateTransform::RlatRlonToLatLon(position[i][1], position[i][0], lat, lon);
			td.val(lon_id, i, step_i) = lon;
			td.val(lat_id, i, step_i) = lat;
			for (int j = 0; j < other_fields.size(); ++j) {
				// hope this works correctly now
				Vec3f coord(position[i][0], position[i][1], position[i][2]);
				float val_0 = other_fields[j][other_ri]->Sample(coord);
				float val_1 = other_fields[j][(other_ri + 1) % 2]->Sample(coord);
				float alpha = (t - file_t[file_i]) / (file_t[file_i + 1] - file_t[file_i]);
				td.val(other_vars[j], i, step_i) = (1 - alpha)*val_0 + alpha * val_1;
			}
		}
		++step_i;
	}
	//--------------------- the end
	for (int i = 0; i < 3; ++i) if (ringbuffer[i] != nullptr) delete ringbuffer[i];
	std::cout << "Trajectories have been computed\n";
}

void ImportantPart::computeTrajectoryDataTEST(TrajectoryData & td, RegScalarField3f * U, RegScalarField3f * V, RegScalarField3f * W)
{
	//--------------------- Initialization
	// see that data has correct sizes
	td.data.resize(td.varnames.size());
	for (int i = 0; i < td.varnames.size(); ++i) {
		td.data[i].resize(td.num_trajectories*td.points_per_trajectory);
	}

	// name of the file containing constants (mainly HHL)
	string constantsfile = basefilename + IntToDDHHMMSS(file_t0) + 'c' + fileending;

	// Not sure if tracer is used ... I'll probably do it "by hand"
	ParticleTracer<Vec3f, 3> tracer;

	// Helping: Map of level to acutal height in meters
	RegScalarField3f* hhl = NetCDF::ImportScalarField3f(constantsfile, "HHL", "rlon", "rlat", "level1");

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

	//cout << "hhl domain: x " << x0 << " to " << x1 << "\ty " << y0 << " to " << y1 << endl;
	//cout << "    voxelsize " << dx << " " << dy << "\tres " << res_x << " " << res_y << endl;

	for (int traj = 0; traj < position.size();++traj) {
		cout << "coord "<<traj<<" begin: " << position[traj][0] << " " << position[traj][1] << " " << position[traj][2] << endl;
		int step_i = 1;
		for (double t = td.time_begin; t < td.time_end; t += dt) {

			Vec3f uvw0 = sampleUVWTEST(position[traj], U, V, W, hhl);
			Vec3f res = position[traj];
			for (int i = 0; i < 3; ++i) {
				Vec3f uvw1 = sampleUVWTEST(res, U, V, W, hhl);
				Vec3f uvw = (uvw0 + uvw1)*0.5;
				res = position[traj] + uvw*dt;
				//cout << "p0 " << position[traj][0] << " " << position[traj][1] << " " << position[traj][2] << endl;
				//cout << "uvw " << uvw[0] << " " << uvw[1] << " " << uvw[2] << endl;
				cout << "p1 " << res[0] << " " << res[1] << " " << res[2] << endl;
			}
			position[traj] = res;
			td.val(rlon_id, traj, step_i) = position[traj][0];
			td.val(rlat_id, traj, step_i) = position[traj][1];
			td.val(z_id, traj, step_i) = position[traj][2];
			CoordinateTransform::RlatRlonToLatLon(position[traj][1], position[traj][0], lat, lon);
			td.val(lon_id, traj, step_i) = lon;
			td.val(lat_id, traj, step_i) = lat;
			//cout << "coord t " << t + dt << ": " << position[traj][0] << " " << position[traj][1] << " " << position[traj][2] << endl;
			cout << "coord t " << t + dt << " lonlat: " << lon << " " << lat << " " << position[traj][2] << endl;
			++step_i;
		}
	}
	cout << "jo\n";
}
