#pragma onceä

#include "ISampleField.hpp"

template<typename TValueType, size_t TDimensions>
class TimeRelatedFields: public ISampleField<TValueType, TDimensions>
{
public:
	TimeRelatedFields(int size) {
		assert(size > 0);
		current_time = 0.0;
		fields.resize(size);
		times.resize(size);
		for (int i = 0; i < size; ++i) {
			fields[i] = nullptr;
			times[i] = 0;
		}
		current_index = 0;
		times_decreasing = false;
	}
	~TimeRelatedFields() {
		for (int i = 0; i < fields.size(); ++i)
			if (fields[i] != nullptr) delete fields[i];
	}

	// sample at any point in time
	virtual TValueType Sample(const Vec<double, TDimensions>& coord, double t) const {
		int ia = -1, ib;
		double wa, wb;
		for (int b = 1; b < fields.size(); ++b) {
			if ( (!times_decreasing && t >= times[ri(b - 1)] && t <= times[ri(b)])
				|| (times_decreasing && t <= times[ri(b - 1)] && t >= times[ri(b)])) {
				ib = ri(b);
				ia = ri(b - 1);
				wb = (t - times[ia]) / (times[ib] - times[ia]);
				wa = 1 - wb;
				break;
			}
		}
		if (ia == -1) return TValueType();
		return fields[ia]->Sample(coord)*wa + fields[ib]->Sample(coord)*wb;
	}
	// sample at the currently stored time
	virtual TValueType Sample(const Vec<double, TDimensions>& coord) const override {
		return Sample(coord, current_time);
	}

	virtual void InsertNextField(ISampleField<TValueType, TDimensions>* field, double t) {
		if (fields[current_index] != nullptr) delete fields[current_index];
		fields[current_index] = field;
		times[current_index] = t;
		current_index = ri(1);
	}

	void setCurrentTime(double t) {
		current_time = t;
		//curernt field
	}
	double getCurrentTime() const { return current_time; }

	bool isBackward() const { return times_decreasing; }
	void setBackward(bool val) { times_decreasing = val; }

private:
	std::vector<ISampleField<TValueType, TDimensions>*> fields;
	std::vector<double> times;

	double current_time;
	int current_index;
	bool times_decreasing;

	inline int ri(int adjust) const { return (current_index + fields.size() + adjust) % fields.size(); }
};