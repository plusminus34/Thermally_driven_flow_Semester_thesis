#include "NetCDF.hpp"
//#include <netcdf.h>
//#include "C:/Program Files/netCDF 4.6.2/include/netcdf.h"
#include "C:/Users/Linus/Documents/Semester_thesis/netCDF4.7.4/include/netcdf.h"
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
