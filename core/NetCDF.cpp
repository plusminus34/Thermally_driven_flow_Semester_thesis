#include "NetCDF.hpp"
//#include <netcdf.h>
#include "C:\Program Files\netCDF 4.7.4\include\netcdf.h"
//#include "C:/Users/Linus/Documents/Semester_thesis/netCDF4.7.4/include/netcdf.h"
#include "RegularGrid.hpp"

const NetCDF::Info::Attribute& NetCDF::Info::Variable::GetAttributeByName(const std::string& name) const
{
	for (size_t i = 0; i < Attributes.size(); ++i)
		if (Attributes[i].GetName() == name)
			return Attributes[i];
	throw "Attribute name '" + name + "' not found.";
}

const NetCDF::Info::Dimension& NetCDF::Info::Variable::GetDimensionByName(const std::string& name) const
{
	for (size_t i = 0; i < Dimensions.size(); ++i)
		if (Dimensions[i].GetName() == name)
			return Dimensions[i];
	throw "Dimension name '" + name + "' not found.";
}

const NetCDF::Info::Dimension& NetCDF::Info::GetDimensionByName(const std::string& name) const
{
	for (size_t i = 0; i < Dimensions.size(); ++i)
		if (Dimensions[i].GetName() == name)
			return Dimensions[i];
	throw "Dimension name '" + name + "' not found.";
}

const NetCDF::Info::Variable& NetCDF::Info::GetVariableByName(const std::string& name) const
{
	for (size_t i = 0; i < Variables.size(); ++i)
		if (Variables[i].GetName() == name)
			return Variables[i];
	throw "Variable name '" + name + "' not found.";
}


bool NetCDF::Info::HasDimension(const std::string& name) const {
	for (size_t i = 0; i < Dimensions.size(); ++i)
		if (Dimensions[i].GetName() == name)
			return true;
	return false;
}

bool NetCDF::Info::HasVariable(const std::string& name) const {
	for (size_t i = 0; i < Variables.size(); ++i)
		if (Variables[i].GetName() == name)
			return true;
	return false;
}

bool NetCDF::ReadInfo(const std::string& path, Info& info)
{
	int status, ncid, unlimdimid;

	// open the file
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return false; } //handle_error(status);

	// read basic counters
	status = nc_inq(ncid, &info.NumDimensions, &info.NumVariables, &info.NumAttributes, &unlimdimid);
	if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);

	// read all dimensions
	for (int dimid = 0; dimid < info.NumDimensions; ++dimid)
	{
		char dim_name[NC_MAX_NAME + 1] = { 0 };
		size_t lenp;
		nc_inq_dim(ncid, dimid, dim_name, &lenp);

		NetCDF::Info::Dimension dim(dim_name, dimid, lenp);
		info.Dimensions.push_back(dim);
	}

	// read all variables
	for (int varid = 0; varid < info.NumVariables; ++varid)
	{
		char var_name[NC_MAX_NAME + 1] = { 0 };
		nc_type var_type;
		int var_ndims;
		int var_dimids[NC_MAX_VAR_DIMS] = { 0 };
		int var_numatts;

		status = nc_inq_var(ncid, varid, var_name, &var_type, &var_ndims, var_dimids, &var_numatts);
		if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);

		Info::Variable var(std::string(var_name), varid, (Info::EType)var_type);
		for (int i = 0; i < var_ndims; ++i)
			var.Dimensions.push_back(info.Dimensions[var_dimids[i]]);
			
		// read the attributes of that variable
		for (int attid = 0; attid < var_numatts; ++attid)
		{
			char att_name[NC_MAX_NAME + 1] = { 0 };
			nc_type att_type;
			size_t att_lenp;

			nc_inq_attname(ncid, varid, attid, att_name);
			if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);
			nc_inq_atttype(ncid, varid, att_name, &att_type);
			if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);
			nc_inq_attlen(ncid, varid, att_name, &att_lenp);
			if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);
				
			switch (att_type)
			{
			case (nc_type)Info::EType::None:	break;
			case (nc_type)Info::EType::Byte: break;
			case (nc_type)Info::EType::Char: att_lenp *= sizeof(char); break;
			case (nc_type)Info::EType::SHORT: att_lenp *= sizeof(short); break;
			case (nc_type)Info::EType::INT: att_lenp *= sizeof(int); break;
			case (nc_type)Info::EType::FLOAT: att_lenp *= sizeof(float); break;
			case (nc_type)Info::EType::DOUBLE: att_lenp *= sizeof(double); break;
			case (nc_type)Info::EType::UBYTE: break;
			case (nc_type)Info::EType::USHORT: att_lenp *= sizeof(unsigned short); break;
			case (nc_type)Info::EType::UINT: att_lenp *= sizeof(unsigned int); break;
			case (nc_type)Info::EType::INT64: att_lenp *= sizeof(int64_t); break;
			case (nc_type)Info::EType::UINT64: att_lenp *= sizeof(uint64_t); break;
			case (nc_type)Info::EType::STRING: att_lenp = (att_lenp + 1) * sizeof(char); break;
			}
				
			Info::Attribute attr(att_name, attid, (Info::EType)att_type, att_lenp);
			nc_get_att(ncid, varid, att_name, attr.GetValue());
			if (status != NC_NOERR) { nc_close(ncid); return false; } //handle_error(status);
			var.Attributes.push_back(attr);
		}
			
		info.Variables.push_back(var);
	}
	nc_close(ncid);
	return true;
}

RegScalarField2f* NetCDF::ImportScalarField2f(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return NULL;

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);
	size_t resX = variable.GetDimensionByName(dimXname).GetLength();
	size_t resY = variable.GetDimensionByName(dimYname).GetLength();

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return NULL;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return NULL; }

	// read the bounds if not provided
	BoundingBox2d domain;
	std::vector<float> dimX, dimY;
	ImportFloatArray(path, dimXname, dimX);
	ImportFloatArray(path, dimYname, dimY);
	domain.ExpandByPoint(Vec2d({ dimX.front(), dimY.front() }));
	domain.ExpandByPoint(Vec2d({ dimX.back(), dimY.back() }));

	// allocate the scalar field
	RegScalarField2f* field = new RegScalarField2f(Vec2i({ (int)resX, (int)resY }), domain);

	if (vartype == Info::EType::FLOAT)
	{
		float* rawdata = field->GetData().data();
		status = nc_get_var_float(ncid, varid, rawdata);
		if (status != NC_NOERR) { delete field; nc_close(ncid); return NULL; }
	}
	else {
		printf("Incompatible format.\n");
		delete field;
		nc_close(ncid); 
		return NULL;
	}
	nc_close(ncid);
	return field;
}

RegScalarField2d* NetCDF::ImportScalarField2d(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return NULL;

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);
	size_t resX = variable.GetDimensionByName(dimXname).GetLength();
	size_t resY = variable.GetDimensionByName(dimYname).GetLength();

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return NULL;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return NULL; }

	// read the bounds if not provided
	BoundingBox2d domain;
	std::vector<float> dimX, dimY;
	ImportFloatArray(path, dimXname, dimX);
	ImportFloatArray(path, dimYname, dimY);
	domain.ExpandByPoint(Vec2d({ dimX.front(), dimY.front() }));
	domain.ExpandByPoint(Vec2d({ dimX.back(), dimY.back() }));

	// allocate the scalar field
	RegScalarField2d* field = new RegScalarField2d(Vec2i({ (int)resX, (int)resY }), domain);

	if (vartype == Info::EType::DOUBLE)
	{
		double* rawdata = field->GetData().data();
		status = nc_get_var_double(ncid, varid, rawdata);
		if (status != NC_NOERR) { delete field; nc_close(ncid); return NULL; }
	}
	else {
		printf("Incompatible format.\n");
		delete field;
		nc_close(ncid);
		return NULL;
	}
	nc_close(ncid);
	return field;
}

RegScalarField3f* NetCDF::ImportScalarField3f(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname, const std::string& dimZname)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return NULL;

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);
	size_t resX = variable.GetDimensionByName(dimXname).GetLength();
	size_t resY = variable.GetDimensionByName(dimYname).GetLength();
	size_t resZ = variable.GetDimensionByName(dimZname).GetLength();

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return NULL;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return NULL; }

	// read the bounds if not provided
	BoundingBox3d domain;
	std::vector<float> dimX, dimY, dimZ;
	if (!ImportFloatArray(path, dimXname, dimX)) { dimX.push_back(0); dimX.push_back(resX - 1); }
	if (!ImportFloatArray(path, dimYname, dimY)) { dimY.push_back(0); dimY.push_back(resY - 1); }
	if (!ImportFloatArray(path, dimZname, dimZ)) { dimZ.push_back(0); dimZ.push_back(resZ - 1); }
	domain.ExpandByPoint(Vec3d({ dimX.front(), dimY.front(), dimZ.front() }));
	domain.ExpandByPoint(Vec3d({ dimX.back(), dimY.back(), dimZ.back() }));

	// allocate the scalar field
	RegScalarField3f* field = new RegScalarField3f(Vec3i({ (int)resX, (int)resY, (int)resZ }), domain);
	
	if (vartype == Info::EType::FLOAT)
	{
		float* rawdata = field->GetData().data();
		status = nc_get_var_float(ncid, varid, rawdata);
		if (status != NC_NOERR) { delete field; nc_close(ncid); return NULL; }
	}
	else {
		printf("Incompatible format.\n");
		delete field;
		nc_close(ncid);
		return NULL;
	}
	nc_close(ncid);
	return field;
}
	
RegScalarField3d* NetCDF::ImportScalarField3d(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname, const std::string& dimZname)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return NULL;

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);
	size_t resX = variable.GetDimensionByName(dimXname).GetLength();
	size_t resY = variable.GetDimensionByName(dimYname).GetLength();
	size_t resZ = variable.GetDimensionByName(dimZname).GetLength();

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return NULL;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return NULL; }

	// read the bounds if not provided
	BoundingBox3d domain;
	std::vector<float> dimX, dimY, dimZ;
	ImportFloatArray(path, dimXname, dimX);
	ImportFloatArray(path, dimYname, dimY);
	ImportFloatArray(path, dimZname, dimZ);
	domain.ExpandByPoint(Vec3d({ dimX.front(), dimY.front(), dimZ.front() }));
	domain.ExpandByPoint(Vec3d({ dimX.back(), dimY.back(), dimZ.back() }));

	// allocate the scalar field
	RegScalarField3d* field = new RegScalarField3d(Vec3i({ (int)resX, (int)resY, (int)resZ }), domain);

	if (vartype == Info::EType::DOUBLE)
	{
		double* rawdata = field->GetData().data();
		status = nc_get_var_double(ncid, varid, rawdata);
		if (status != NC_NOERR) { delete field; nc_close(ncid); return NULL; }
	}
	else {
		printf("Incompatible format.\n");
		delete field;
		nc_close(ncid);
		return NULL;
	}
	nc_close(ncid);
	return field;
}

bool NetCDF::ImportFloat(const std::string& path, const std::string& varname, float& output)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return false;

	if (!info.HasVariable(varname)) { return false; }

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return false;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return false; }
	
	// allocate the scalar field
	if (vartype == Info::EType::FLOAT) {
		status = nc_get_var_float(ncid, varid, &output);
		if (status != NC_NOERR) { nc_close(ncid); return false; }
	}
	else if (vartype == Info::EType::DOUBLE) {
		double output_double;
		status = nc_get_var_double(ncid, varid, &output_double);
		if (status != NC_NOERR) { nc_close(ncid); return false; }
		output = static_cast<float>(output_double);
	}

	nc_close(ncid);
	return true;
}


bool NetCDF::ImportFloatArray(const std::string& path, const std::string& varname, std::vector<float>& floatArray)
{
	// get the info object
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return false;

	if (!info.HasVariable(varname)) { return false; }

	// read the resolution from the info object
	const NetCDF::Info::Variable& variable = info.GetVariableByName(varname);

	// get meta information on the variable
	int varid = variable.GetID();
	Info::EType vartype = variable.GetType();
	if (vartype != Info::EType::FLOAT && vartype != Info::EType::DOUBLE) {
		printf("Unsupported format!");
		return false;
	}

	// open the file
	int status, ncid;
	status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) { return false; } 

	// allocate the scalar field
	if (vartype == Info::EType::FLOAT) {
		floatArray.resize(variable.Dimensions[0].GetLength());
		float* rawdata = floatArray.data();
		status = nc_get_var_float(ncid, varid, rawdata);
		if (status != NC_NOERR) { nc_close(ncid); return false; } 
	}
	else if (vartype == Info::EType::DOUBLE) {
		floatArray.resize(variable.Dimensions[0].GetLength());
		float* rawdata = floatArray.data();
		double* rawdbl = new double[variable.Dimensions[0].GetLength()];
		status = nc_get_var_double(ncid, varid, rawdbl);
		for (size_t i = 0; i < variable.Dimensions[0].GetLength(); ++i)
			rawdata[i] = static_cast<float>(rawdbl[i]);
		if (status != NC_NOERR) { nc_close(ncid); return false; }
	}
	nc_close(ncid);
	return true;
}

bool NetCDF::ReadPaths(const std::string& path, std::vector<std::vector<Vec3f>>& paths) {
	NetCDF::Info info;
	if (!ReadInfo(path, info)) return false;

	int nPaths = info.GetVariableByName("X").GetDimensionByName("pathID").GetLength();
	int pathLength = info.GetVariableByName("X").GetDimensionByName("pathT").GetLength();

	paths.resize(nPaths);
	for (int i = 0; i < nPaths; ++i) paths[i].resize(pathLength);
	int x_id = info.GetVariableByName("X").GetID();
	int y_id = info.GetVariableByName("Y").GetID();
	int z_id = info.GetVariableByName("Z").GetID();

	int ncid;
	if (nc_open(path.c_str(), NC_NOWRITE, &ncid)) return false;
	std::vector<float> data_x(nPaths*pathLength);
	std::vector<float> data_y(nPaths*pathLength);
	std::vector<float> data_z(nPaths*pathLength);
	if (nc_get_var_float(ncid, x_id, &data_x[0])) return false;
	if (nc_get_var_float(ncid, y_id, &data_y[0])) return false;
	if (nc_get_var_float(ncid, z_id, &data_z[0])) return false;
	for (int i = 0; i < nPaths; ++i) {
		for (int j = 0; j < pathLength; ++j) {
			paths[i][j][0] = data_x[i*pathLength + j];
			paths[i][j][1] = data_y[i*pathLength + j];
			paths[i][j][2] = data_z[i*pathLength + j];
		}
	}

	nc_close(ncid);

	return true;
}
bool NetCDF::WritePaths(const std::string& path, const std::vector<std::vector<Vec3f>>& paths) {
	size_t nPaths = paths.size();
	size_t pathLength = paths[0].size();
	for (int i = 1; i < nPaths; ++i) if (paths[i].size() != pathLength) return false;
	const size_t nPoints = nPaths * pathLength;
	std::vector<float> pathsData(nPoints * 3);
	for (int i = 0; i < nPaths; ++i) {
		for (int j = 0; j < pathLength; ++j) {
			pathsData[i*pathLength + j] = paths[i][j][0];
			pathsData[i*pathLength + j + nPoints] = paths[i][j][1];
			pathsData[i*pathLength + j + 2 * nPoints] = paths[i][j][2];
		}
	}

	int status, ncid;
	if (status = nc_create(path.c_str(), NC_CLOBBER, &ncid)) { return false; }

	int pathID_id, pathT_id;
	char pathIDName[] = "Path ID"; char pathTName[] = "Path T";
	char long_name[] = "long_name", standard_name[] = "standard_name", units[] = "units";
	char axis[] = "axis";
	char xname[] = "X", yname[] = "Y", zname[] = "Z";

	//important
	//define variables and dimensions and attributes

	if (status = nc_def_dim(ncid, "pathID", nPaths, &pathID_id)) return false;
	if (status = nc_def_dim(ncid, "pathT", pathLength, &pathT_id)) return false;

	int fieldIDs[] = { pathID_id,pathT_id };
	int xID, yID, zID;
	int pathID_var, pathT_var;
	if (nc_def_var(ncid, pathIDName, NC_FLOAT, 1, &pathID_id, &pathID_var)) return false;
	if (nc_def_var(ncid, pathTName, NC_FLOAT, 1, &pathT_id, &pathT_var)) return false;
	if (nc_def_var(ncid, xname, NC_FLOAT, 2, fieldIDs, &xID)) return false;
	if (nc_def_var(ncid, yname, NC_FLOAT, 2, fieldIDs, &yID)) return false;
	if (nc_def_var(ncid, zname, NC_FLOAT, 2, fieldIDs, &zID)) return false;

	if (status = nc_put_att_text(ncid, xID, standard_name, strlen(xname), xname)) return false;
	if (status = nc_put_att_text(ncid, yID, standard_name, strlen(yname), yname)) return false;
	if (status = nc_put_att_text(ncid, zID, standard_name, strlen(zname), zname)) return false;

	if (nc_enddef(ncid))return false;

	//fill in data
	size_t startp[] = { 0, 0 };
	size_t countp[] = { nPaths, pathLength };
	for (int i = 0; i < nPaths; ++i)
		for (int j = 0; j < pathLength; ++j) {
			if (nc_put_vara_float(ncid, xID, startp, countp, &pathsData[0])) return false;
			if (nc_put_vara_float(ncid, yID, startp, countp, &pathsData[nPoints])) return false;
			if (nc_put_vara_float(ncid, zID, startp, countp, &pathsData[2 * nPoints])) return false;
		}

	if (status = nc_close(ncid)) return false;

	return true;
}

bool NetCDF::WriteTrajectoryData(const std::string & path, const TrajectoryData & td)
{
	int status, ncid;
	if (status = nc_create(path.c_str(), NC_CLOBBER, &ncid)) return false;

	printf("Writing dimensions\n");
	// dimensions
	int ntra_id, ntim_id;
	if (status = nc_def_dim(ncid, "ntra", td.num_trajectories, &ntra_id)) return false;
	if (status = nc_def_dim(ncid, "ntim", td.points_per_trajectory, &ntim_id)) return false;
	int fieldIDs[] = { ntra_id, ntim_id };

	printf("Writing variables\n");
	// variables
	int time_id;
	std::vector<int> var_id(td.varnames.size());
	if (nc_def_var(ncid, "time", NC_DOUBLE, 2, fieldIDs, &time_id)) return false;
	printf("  time there\n");
	for (int i = 0; i < td.varnames.size(); ++i) {
		if (nc_def_var(ncid, td.varnames[i].c_str(), NC_FLOAT, 2, fieldIDs, &var_id[i])) return false;
	}

	printf("Writing attributes\n");
	// attributes
	int global_id = -1;//TODO is this how you write global attributes?
	if (status = nc_put_att_int(ncid, global_id, "ref_year", NC_INT, 1, &td.ref_year)) return false;
	if (status = nc_put_att_int(ncid, global_id, "ref_month", NC_INT, 1, &td.ref_month)) return false;
	if (status = nc_put_att_int(ncid, global_id, "ref_day", NC_INT, 1, &td.ref_day)) return false;
	if (status = nc_put_att_int(ncid, global_id, "ref_hour", NC_INT, 1, &td.ref_hour)) return false;
	if (status = nc_put_att_int(ncid, global_id, "ref_min", NC_INT, 1, &td.ref_min)) return false;
	int duration = floor(td.time_end - td.time_begin) / 60;//TODO duration in minutes?
	if (status = nc_put_att_int(ncid, global_id, "duration", NC_INT, 1, &duration)) return false;
	if (status = nc_put_att_float(ncid, global_id, "pollon", NC_FLOAT, 1, &td.pole_lon)) return false;
	if (status = nc_put_att_float(ncid, global_id, "pollat", NC_FLOAT, 1, &td.pole_lat)) return false;
	
	/*
TODO
data:
	varname:
		(tra0, t0), (tra1, t0), ..., (tran, t0),
		(tra0, t1), (tra1, t1), ..., (tran, t1),
		...
		(tra0, tT), (tra1, tT), ..., (tran, tT),
	(note: time has format hh.mm)
*/
	printf("Writing data\n");
	size_t startp[] = { 0, 0 };
	size_t countp[] = { td.num_trajectories, td.points_per_trajectory };
	for (int i = 0; i < td.points_per_trajectory; ++i) {
		printf(("i " + std::to_string(i) + "\n").c_str());
		double time = 0;//TODO format hh.mm
		for (int j = 0; j < td.num_trajectories; ++j) {
			printf(("  j " + std::to_string(j) + "\n").c_str());
			if (nc_put_vara_double(ncid, time_id, startp, countp, &time)) return false;
		}
	}
	for (int i = 0; i < td.varnames.size(); ++i) {
		//TODO write data
	}
	/*up there
	for (int i = 0; i < nPaths; ++i)
		for (int j = 0; j < pathLength; ++j) {
			if (nc_put_vara_float(ncid, xID, startp, countp, &pathsData[0])) return false;
			if (nc_put_vara_float(ncid, yID, startp, countp, &pathsData[nPoints])) return false;
			if (nc_put_vara_float(ncid, zID, startp, countp, &pathsData[2 * nPoints])) return false;
		}
	*/

	printf("actually, forget it\n");
	if (status = nc_close(ncid)) return false;

	return true;
}

bool NetCDF::ReadTrajectoryData(const std::string & path, TrajectoryData & td)
{
	return false;
}