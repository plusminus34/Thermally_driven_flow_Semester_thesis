#pragma once

#include <assert.h>
#include <algorithm>

template<typename TScalarType, size_t TDimensions>
class Vec
{
public:
	using TScalar = TScalarType;
	static constexpr size_t Dimensions = TDimensions;

	// Default constructor. All components are zero.
	Vec() { std::memset(mData, 0, TDimensions * sizeof(TScalarType)); }

	// Initializes the vector from another vector of equal size
	template<class TOther>
	Vec(const TOther& other) { for (size_t i = 0; i < TDimensions; ++i) mData[i] = static_cast<TScalarType>(other[i]); }

	// Initializes a vector with initializer list.
	explicit Vec(std::initializer_list<TScalarType> list) {
		int i = 0;
		for (auto elem : list) {
			if (i >= TDimensions) return;
			mData[i++] = elem;
		}
		for (; i < TDimensions; ++i) mData[i] = 0;
	}

	Vec(TScalarType x, TScalarType y, TScalarType z) {
		assert(TDimensions >= 3);
		mData[0] = x;
		mData[1] = y;
		mData[2] = z;
		for (int = 3; i < TDimensions; ++i) mData[i] = 0;
	}

	const TScalarType& operator[](size_t index) const {
		assert(index >= 0 && index < Dimensions);
		return mData[index];
	}
	TScalarType& operator[](size_t index) {
		assert(index >= 0 && index < Dimensions);
		return mData[index];
	}

	TScalarType* ptr() { return mData; }
	const TScalarType* ptr() const { return mData; }

	// Adds a matrix element-wise.
	Vec operator+(const Vec& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] + other.ptr()[i];
		return result;
	}

	// Subtracts a matrix element-wise.
	Vec operator-(const Vec& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] - other.ptr()[i];
		return result;
	}

	// Adds a scalar element-wise.
	Vec operator+(const TScalar& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] + other;
		return result;
	}

	// Subtracts a scalar element-wise.
	Vec operator-(const TScalar& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] - other;
		return result;
	}

	// Multiplies by a scalar component-wise.
	Vec operator*(const TScalar& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] * other;
		return result;
	}

	// Divides by a scalar component-wise.
	Vec operator/(const TScalar& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = self.ptr()[i] / other;
		return result;
	}

	// Negates a vector component-wise.
	Vec operator-() const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = -self.ptr()[i];
		return result;
	}

	// Adds a matrix element-wise.
	Vec& operator+=(const Vec& other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] += other.ptr()[i];
		return self;
	}

	// Subtracts a matrix element-wise.
	Vec& operator-=(const Vec& other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] -= other.ptr()[i];
		return self;
	}

	// Adds a scalar component-wise.
	Vec& operator+=(TScalar other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] += other;
		return self;
	}

	// Subtracts a scalar component-wise.
	Vec& operator-=(TScalar other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] -= other;
		return self;
	}

	// Multiplies by a scalar component-wise.
	Vec& operator*=(TScalar other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] *= other;
		return self;
	}

	// Divides by a scalar component-wise.
	Vec& operator/=(TScalar other)
	{
		Vec& self = *static_cast<Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			self.ptr()[i] /= other;
		return self;
	}

	// Multiplies a vector element-wise.
	Vec operator*(const Vec& other) const
	{
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = mData[i] * other.ptr()[i];
		return result;
	}

	// Divides a vector element-wise.
	Vec operator/(const Vec& other) const
	{
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = mData[i] / other.ptr()[i];
		return result;
	}

	// Equality operator
	bool operator==(const Vec& other) const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		return memcmp(self.ptr(), other.ptr(), Dimensions * sizeof(TScalar)) == 0;
	}

	// Inequality operator
	bool operator!=(const Vec& other) const
	{
		return !operator==(other);
	}

	// Sets the matrix to a zero matrix.
	void setZero() {
		Vec& self = *static_cast<Vec*>(this);
		std::memset(self.ptr(), 0, Dimensions * sizeof(TScalar));
	}

	// Checks if the matrix is exactly the zero-matrix.
	bool isZero() const
	{
		const Vec& self = *static_cast<const Vec*>(this);
		for (int i = 0; i < Dimensions; ++i)
			if (self.ptr()[i])
				return false;
		return true;
	}

	// Computes the component-wise minimum of two vectors
	static Vec min(const Vec& a, const Vec& b) {
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = std::min(a.ptr()[i], b.ptr()[i]);
		return result;
	}

	// Computes the component-wise maximum of two vectors
	static Vec max(const Vec& a, const Vec& b) {
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = std::max(a.ptr()[i], b.ptr()[i]);
		return result;
	}

	// Linealy interpolates between two vectors.
	static Vec lerp(const Vec& x, const Vec& y, TScalar t) {
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = x.ptr()[i] + (y.ptr()[i] - x.ptr()[i]) * t;
		return result;
	}

	// Linealy interpolates between two vectors.
	static Vec lerp(const Vec& x, const Vec& y, const Vec& t) {
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = x.ptr()[i] + (y.ptr()[i] - x.ptr()[i]) * t.ptr()[i];
		return result;
	}

	// Returns a zero matrix.
	static Vec zeros() { return Vec(); }

	// Returns a one matrix.
	static Vec ones() {
		Vec result;
		for (int i = 0; i < Dimensions; ++i)
			result.ptr()[i] = TScalar(1);
		return result;
	}

	// Computes the length of the vector (L1 norm)
	TScalarType length() const
	{
		TScalarType result = (TScalarType)0;
		for (int i = 0; i < Dimensions; ++i)
			result += mData[i] * mData[i];
		return (TScalarType)std::sqrt(result);
	}

	// Computes the squared length of the vector (L2 norm)
	TScalarType lengthSquared() const
	{
		TScalarType result = (TScalarType)0;
		for (int i = 0; i < Dimensions; ++i)
			result += mData[i] * mData[i];
		return result;
	}

	// Normalizes the vector and optionally returns the length of the vector.
	const Vec& normalize(TScalarType* len = nullptr)
	{
		TScalarType l = length();
		if (len)
			*len = l;
		if (l)
			*this *= (TScalarType)(1.0 / l);
		return *(Vec*)this;
	}

	// Sorts the vector components in ascending order.
	Vec& sortAscend() {
		std::sort(mData, mData + Dimensions);
		Vec& self = *static_cast<Vec*>(this);
		return self;
	}

	// Sorts the vector components in descending order.
	Vec& sortDescend() {
		std::sort(mData, mData + Dimensions, std::greater<TScalarType>());
		Vec& self = *static_cast<Vec*>(this);
		return self;
	}

	// Computes a dot product between this vector and another one.
	TScalarType dot(const Vec& other) const {
		TScalarType result(0);
		for (int i = 0; i < Dimensions; ++i)
			result += mData[i] * other.ptr()[i];
		return result;
	}

	// Computes the dot product between two vectors.
	static TScalarType dot(const Vec& a, const Vec& b) { return a.dot(b); }

private:
	TScalarType mData[Dimensions];
};

// Base class for bounding boxes of any dimension. The template argument is a vector valued type.
template<typename VecType>
class BoundingBox
{
public:
	typedef VecType Vec;

	// Default constructor. The bounding box is set to infinite with clamp boundary handling.
	BoundingBox() {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			mMin[i] = std::numeric_limits<typename Vec::TScalar>::max();
			mMax[i] = -std::numeric_limits<typename Vec::TScalar>::max();
		}
	}

	// Constructor. Creates a bounding box from two points.
	BoundingBox(const Vec& Pnt1, const Vec& Pnt2) {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			mMin[i] = std::min(Pnt1[i], Pnt2[i]);
			mMax[i] = std::max(Pnt1[i], Pnt2[i]);
		}
	}

	// Copy-constructor.
	BoundingBox(const BoundingBox& other) {
		this->mMin = other.mMin;
		this->mMax = other.mMax;
	}

	// Gets the min corner.
	const Vec& GetMin() const { return mMin; }
	// Gets the max corner.
	const Vec& GetMax() const { return mMax; }
	
	// Expands the bounding box.
	void ExpandByBox(const BoundingBox& other) {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			this->mMin[i] = std::min(this->mMin[i], other.mMin[i]);
			this->mMax[i] = std::max(this->mMax[i], other.mMax[i]);
		}
	}

	// Expands the box by a point.
	void ExpandByPoint(const Vec& point) {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			this->mMin[i] = std::min(this->mMin[i], point[i]);
			this->mMax[i] = std::max(this->mMax[i], point[i]);
		}
	}

	// Checks if the box is infinite.
	bool IsInfinite() const {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			if (mMin[i] == std::numeric_limits<typename Vec::TScalar>::max() ||
				mMax[i] == -std::numeric_limits<typename Vec::TScalar>::max())
				return true;
		}
		return false;
	}

	// Makes the box inverse infinite.
	void Reset() {
		for (int i = 0; i < Vec::Dimensions; ++i) {
			mMin[i] = std::numeric_limits<typename Vec::TScalar>::max();
			mMax[i] = -std::numeric_limits<typename Vec::TScalar>::max();
		}
	}

	// Checks if a point is inside the box.
	bool Contains(const Vec& x) const {
		for (size_t i = 0; i < Vec::Dimensions; ++i) {
			if (x[i] < mMin[i] || mMax[i] < x[i])
				return false;
		}
		return true;
	}

	// Clamps a point by the boundary of the domain.
	template<typename TVecOther>
	TVecOther ClampToDomain(const TVecOther& pnt) const {
		TVecOther result = pnt;
		for (int i = 0; i < std::min(Vec::Dimensions, TVecOther::Dimensions); ++i) {
			result[i] = std::min(std::max(mMin[i], pnt[i]), mMax[i]);
		}
		return result;
	}

protected:
	// min corner
	Vec mMin;
	// max corner
	Vec mMax;
};

typedef Vec<float, 2> Vec2f; 
typedef Vec<float, 3> Vec3f;
typedef Vec<double, 2> Vec2d;
typedef Vec<double, 3> Vec3d;
typedef Vec<int, 2> Vec2i;
typedef Vec<int, 3> Vec3i;

// Computes the cross product of two vectors
template<typename TScalarType>
Vec<TScalarType, 3> cross(const Vec<TScalarType, 3>& v1, const Vec<TScalarType, 3>& v2) {
	return Vec<TScalarType, 3>({
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0] });
}


typedef BoundingBox<Vec2d> BoundingBox2d;
typedef BoundingBox<Vec3d> BoundingBox3d;