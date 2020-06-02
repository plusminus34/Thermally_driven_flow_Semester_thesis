#pragma once

#include <vector>
#include <string>
#include "CoordinateTransform.hpp"

// Struct for the kind of data we want as output
struct TrajectoryData {
	// info about trajectories
	int num_trajectories;
	int points_per_trajectory;
	float time_begin, time_end;// in seconds

	// info that is constant across my input data
	int ref_year = 2016;
	int ref_month = 11;
	int ref_day = 21;
	int ref_hour = 0;
	int ref_min = 0;
	float pole_lat = POLE_PHI;
	float pole_lon = POLE_LAMBDA;

	// name of variables (generally, the first two will be rlon and rlat)
	std::vector<std::string> varnames;

	// the data itself, stored as one array per variable
	/* Format:
	data[i] for variable [i]:
		(t0,p0), (t0,p1), ..., (t0,pN),
		(t1,p0), (t1,p1), ..., (t1,pN),
		...
		(tT,p0), (tT,p1), ..., (tT,pN)
	*/
	std::vector<std::vector<float>> data;

	// min and max values of each variable in data
	std::vector<float> min_values;
	std::vector<float> max_values;
	void compute_value_bounds() {
		min_values.resize(data.size());
		max_values.resize(data.size());
		for (int i = 0; i < data.size(); ++i) {
			min_values[i] = data[i][0];
			max_values[i] = data[i][0];
			for (int j = 1; j < data[i].size(); ++j) {
				if (data[i][j] < min_values[i]) min_values[i] = data[i][j];
				else if (data[i][j] > max_values[i]) max_values[i] = data[i][j];
			}
		}
	}

	int get_var_id(std::string varname) const {
		for (int i = 0; i < varnames.size(); ++i)
			if (varnames[i] == varname) return i;
		return -1;
	}

	float get_value(int var_id, int trajectory_id, int point_i) const {
		return data[var_id][trajectory_id + point_i * num_trajectories];
	}

	float& val(int var_id, int trajectory_id, int point_i) {
		return data[var_id][trajectory_id + point_i * num_trajectories];
	}
};