#pragma once

#include "Math.hpp"
#include <vector>
#include "ISampleField.hpp"

// Base class for a field on a regular grid
template<typename TValueType, size_t TDimensions>
class RegularGrid : public ISampleField<TValueType, TDimensions>
{
private:
	using TGridCoord = Vec<int, TDimensions>;
	static constexpr size_t Dimensions = TDimensions;
	using TValue = TValueType;
	using TBoundingBox = BoundingBox<Vec<double, TDimensions>>;
	using TDomainCoord = Vec<double, TDimensions>;

public:
	
	// Sampling methods.
	enum class ESampler {
		Nearest,
		Linear
	};

	// Constructor.
	RegularGrid(const TGridCoord& res, const TBoundingBox& domain) :
		mResolution(res),
		mDomain(domain),
		mSampler(ESampler::Linear)
	{
		int numElements = 1;
		for (int d = 0; d < this->Dimensions; ++d)
			numElements *= res[d];
		mData.resize(numElements);
	}

	// Disable copy-constructor.
	RegularGrid(const RegularGrid&) = delete;
	// Destructor.
	~RegularGrid() {}

	// Samples the field.
	virtual TValue Sample(const TDomainCoord& coord) const override
	{
		TDomainCoord position = this->mDomain.ClampToDomain(coord);
		TDomainCoord vfTex = (position - this->mDomain.GetMin()) / (this->mDomain.GetMax() - this->mDomain.GetMin());
		TDomainCoord vfSample = vfTex * (static_cast<TDomainCoord>(mResolution) - TDomainCoord::ones());
		switch (mSampler) {
		case ESampler::Nearest:
		{
			TGridCoord gridCoord;
			for (int i = 0; i < TDomainCoord::Dimensions; ++i) {
				gridCoord[i] = std::min(std::max(0, (int)(vfSample[i] + typename TDomainCoord::TScalar(0.5))), mResolution[i] - 1);
			}
			return GetVertexDataAt(gridCoord);
		}
		case ESampler::Linear:
		{
			TGridCoord viSampleBase0, viSampleBase1;
			for (int i = 0; i < TDomainCoord::Dimensions; ++i) {
				viSampleBase0[i] = std::min(std::max(0, (int)vfSample[i]), mResolution[i] - 1);
				viSampleBase1[i] = std::min(std::max(0, viSampleBase0[i] + 1), mResolution[i] - 1);
			}
			TDomainCoord vfSampleInterpol = vfSample - static_cast<TDomainCoord>(viSampleBase0);

			size_t numCorners = (size_t)std::pow(2, TDomainCoord::Dimensions);
			TValue result({ 0 });
			for (size_t i = 0; i < numCorners; ++i) {
				typename TDomainCoord::TScalar weight(1);
				TGridCoord gridCoord;
				for (size_t d = 0; d < TDomainCoord::Dimensions; ++d) {
					if (i & (size_t(1) << (TDomainCoord::Dimensions - size_t(1) - d))) {
						gridCoord[d] = viSampleBase1[d];
						weight *= vfSampleInterpol[d];
					}
					else {
						gridCoord[d] = viSampleBase0[d];
						weight *= 1 - vfSampleInterpol[d];
					}
				}
				result += static_cast<TValueType>(GetVertexDataAt(gridCoord) * weight);
			}
			return result;
		}
		}
		return TValue();
	}

	// Gets the linear array index based on a grid coordinate index.
	int GetLinearIndex(const TGridCoord& gridCoord) const
	{
		int stride = 1;
		int linearIndex = gridCoord[0];
		for (int d = 1; d < this->Dimensions; ++d) {
			stride *= mResolution[d - 1];
			linearIndex += gridCoord[d] * stride;
		}
		return linearIndex;
	}

	// Gets the spatial location of a grid vertex.
	TDomainCoord GetCoordAt(const TGridCoord& gridCoord) const
	{
		TDomainCoord s;
		for (int i = 0; i < TDomainCoord::Dimensions; ++i) {
			s[i] = mResolution[i] < 2 ? 0.5 : gridCoord[i] / (mResolution[i] - typename TDomainCoord::TScalar(1.));
		}
		return TDomainCoord::lerp(this->mDomain.GetMin(), this->mDomain.GetMax(), s);
	}

	// Gets the grid coordinate based on the linear array index.
	TGridCoord GetGridCoord(const size_t& linearIndex) const
	{
		TGridCoord result;
		int stride = 1;
		for (int d = 0; d < TGridCoord::Dimensions - 1; ++d)
			stride *= mResolution[d];

		int t = (int)linearIndex;
		for (int d = TGridCoord::Dimensions - 1; d >= 0; --d) {
			result[d] = t / stride;
			t = t % stride;
			if (d > 0)
				stride /= mResolution[d - 1];
		}
		return result;
	}

	// Gets the vertex data stored at a certain grid coordinate.
	TValue GetVertexDataAt(const TGridCoord& gridCoord) const
	{
		for (int dim = 0; dim < this->Dimensions; ++dim)
			assert(gridCoord[dim] >= 0 && gridCoord[dim] < mResolution[dim]);

		int addr = GetLinearIndex(gridCoord);
		return mData[addr];
	}

	// Sets the vertex data at a certain grid coordinate.
	void SetVertexDataAt(const TGridCoord& gridCoord, const TValue& value)
	{
		for (int dim = 0; dim < this->Dimensions; ++dim)
			assert(gridCoord[dim] >= 0 && gridCoord[dim] < mResolution[dim]);

		int addr = GetLinearIndex(gridCoord);
		mData[addr] = value;
	}

	// Gets the grid resolution.
	const TGridCoord& GetResolution() const { return mResolution; }
	// Gets the domain.
	const TBoundingBox& GetDomain() const { return mDomain; }

	// Gets the sampling method.
	const ESampler& GetSamplingMethod() const { return mSampler; }
	// Sets the sampling method.
	void SetSamplingMethod(const ESampler& sampler) { mSampler = sampler; }

	std::vector<TValue>& GetData() { return mData; }
	const std::vector<TValue>& GetData() const { return mData; }

	TDomainCoord GetVoxelSize() const {
		return (this->mDomain.GetMax() - this->mDomain.GetMin()) / (static_cast<TDomainCoord>(mResolution) - TDomainCoord::ones());
	}

private:
	std::vector<TValue> mData;
	TGridCoord mResolution;
	TBoundingBox mDomain;
	ESampler mSampler;
};

typedef RegularGrid<float, 2> RegScalarField2f;
typedef RegularGrid<float, 3> RegScalarField3f;
typedef RegularGrid<double, 2> RegScalarField2d;
typedef RegularGrid<double, 3> RegScalarField3d;

typedef RegularGrid<Vec2f, 2> RegVectorField2f;
typedef RegularGrid<Vec3f, 3> RegVectorField3f;
typedef RegularGrid<Vec2d, 2> RegVectorField2d;
typedef RegularGrid<Vec3d, 3> RegVectorField3d;