#pragma once

#include "Math.hpp"
#include <vector>
#include <cmath>

template<typename TValueType, size_t TDimensions>
class ISampleField
{
	virtual TValueType Sample(const Vec<double, TDimensions>& coord) const = 0;
};