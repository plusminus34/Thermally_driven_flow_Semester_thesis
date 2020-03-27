#pragma once

#include <string>
#include <vector>
#include "RegularGrid.hpp"

class NetCDF
{
public:
	class Info
	{
	public:
		// an enumeration of all supported types
		enum class EType
		{
			None,
			Byte,
			Char,
			SHORT,
			INT,
			FLOAT,
			DOUBLE,
			UBYTE,
			USHORT,
			UINT,
			INT64,
			UINT64,
			STRING
		};

		// class that stores a dimension
		class Dimension
		{
		public:
			Dimension(const std::string& name, const int& ID, const size_t& length) : mName(name), mID(ID), mLength(length) {}
			const std::string& GetName() const { return mName; }
			const int& GetID() const { return mID; }
			const size_t& GetLength() const { return mLength; }
		private:
			std::string mName;
			int mID;
			size_t mLength;
		};

		// class that stores an attribute
		class Attribute
		{
		public:
			Attribute(const std::string& name, const int& id, const EType& type, size_t& length) : mName(name), mID(id), mType(type), mLength(length) { mValue.resize(length, 0); }
			~Attribute() {}
			const std::string& GetName() const { return mName; }
			const int& GetID() const { return mID; }
			const EType& GetType() const { return mType; }

			char* GetValue() { return mValue.data(); }
			const char* GetValue() const { return mValue.data(); }
			const int* GetValueAsInt() const { return (int*)mValue.data(); }
			const float* GetValueAsFloat() const { return (float*)mValue.data(); }
			const double* GetValueAsDouble() const { return (double*)mValue.data(); }
			const char* GetValueAsChar() const { return (char*)mValue.data(); }
		private:
			std::string mName;
			int mID;
			EType mType;
			std::vector<char> mValue;
			size_t mLength;
		};

		// class that stores a variable
		class Variable
		{
		public:
			Variable(const std::string& name, const int& ID, const EType& type) : mName(name), mID(ID), mType(type) {}
				
			const std::string& GetName() const { return mName; }
			const int& GetID() const { return mID; }
			const EType& GetType() const { return mType; }

			const Attribute& GetAttributeByName(const std::string& name) const;
			const Dimension& GetDimensionByName(const std::string& name) const;

			std::vector<Dimension> Dimensions;
			std::vector<Attribute> Attributes;

		private:
			std::string mName;
			int mID;
			EType mType;
		};

		Info() : NumDimensions(0), NumAttributes(0), NumVariables(0) {}

		int NumDimensions;
		int NumAttributes;
		int NumVariables;

		std::vector<Dimension> Dimensions;
		std::vector<Variable> Variables;

		const Dimension& GetDimensionByName(const std::string& name) const;
		const Variable& GetVariableByName(const std::string& name) const;

		bool HasDimension(const std::string& name) const;
		bool HasVariable(const std::string& name) const;
	};
		
	// reads the info object, desccribing the nc file
	static bool ReadInfo(const std::string& path, Info& info);

	// imports a steady 2d scalar field from an nc file
	static RegScalarField2f* ImportScalarField2f(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname);
	// imports a steady 3d scalar field from an nc file. Providing the bounding box is optional. If it is not provided, this functions reads the dimensions to get the bounds itself.
	static RegScalarField3f* ImportScalarField3f(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname, const std::string& dimZname);		

	// imports a steady 2d scalar field from an nc file
	static RegScalarField2d* ImportScalarField2d(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname);
	// imports a steady 3d scalar field from an nc file. Providing the bounding box is optional. If it is not provided, this functions reads the dimensions to get the bounds itself.
	static RegScalarField3d* ImportScalarField3d(const std::string& path, const std::string& varname, const std::string& dimXname, const std::string& dimYname, const std::string& dimZname);
		
	// imports a float value
	static bool ImportFloat(const std::string& path, const std::string& varname, float& output);
	// imports a float array
	static bool ImportFloatArray(const std::string& path, const std::string& varname, std::vector<float>& output);
};
