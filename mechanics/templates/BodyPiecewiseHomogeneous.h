#ifndef BODYPWH_H
#define BODYPWH_H
#include <vector>
#include "Inclusion.h"
#include "Region.h"
#include <ostream>
#include <istream>
#include <limits>

template<typename T>
class BodyPiecewiseHomogeneous
{
public:
	std::vector<Inclusion<T> > inclusions;
	Region<T> areaHi;
	Region<T> areaLow;

	BodyPiecewiseHomogeneous() = default;
	~BodyPiecewiseHomogeneous() = default;
	BodyPiecewiseHomogeneous(const BodyPiecewiseHomogeneous&) = default;
	BodyPiecewiseHomogeneous& operator= (const BodyPiecewiseHomogeneous&) = default;
	BodyPiecewiseHomogeneous(BodyPiecewiseHomogeneous&&) = default;
	BodyPiecewiseHomogeneous& operator= (BodyPiecewiseHomogeneous&&) = default;

	template<typename U>
	friend std::ostream& operator<< (std::ostream& os_, const BodyPiecewiseHomogeneous<U>& b_);
	template<typename U>
	friend std::istream& operator>> (std::istream& is_, BodyPiecewiseHomogeneous<U>& b_);
};

template<typename T>
std::ostream& operator<< (std::ostream& os_, const BodyPiecewiseHomogeneous<T>& b_)
{
	os_ << b_.areaHi;
	os_ << b_.areaLow;
	os_ << "inclusion_count = " <<  b_.inclusions.size() << "\n";
	for (size_t i = 0; i < b_.inclusions.size(); ++i) 
	{
		os_ << b_.inclusions[i];
	}
	return os_;
}

template<typename T>
std::istream& operator>> (std::istream& is_, BodyPiecewiseHomogeneous<T>& b_)
{

	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};
	
	if (is_.fail()) 
		return is_;

	Inclusion<T> incl_;
	is_ >> b_.areaHi;
	is_ >> b_.areaLow;

	int incl_count = 0;
	readValueFromSymbol(incl_count);

	for(;incl_count != 0; --incl_count)
	{
		is_ >> incl_;		
		b_.inclusions.push_back(incl_);
	}
	return is_;
}

#endif