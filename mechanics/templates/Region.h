#ifndef REGION_H
#define REGION_H

#include "BasicAnisotropy.h"
#include "Stress.h"
#include <ostream>
#include <istream>

template<typename T>
class Region
{
public:
	BasicAnisotropy<T> material;
	Stress<T> syz_isxz;

	Region() = default;
	~Region() = default;
	Region(const Region&) = default;
	Region& operator= (const Region&) = default;
	Region(Region&&) = default;
	Region& operator= (Region&&) = default;

	template<typename U>
	friend std::ostream& operator<<(std::ostream& os_, const Region<U>& reg_);
	template<typename U>
	friend std::istream& operator>> (std::istream& is_, Region<U>& reg_);

};

template<typename T>
std::ostream& operator<< (std::ostream& os_, const Region<T>& reg_)
{
	os_ << reg_.material;
	os_ << reg_.syz_isxz;
	return os_; 
}

template<typename T>
std::istream& operator>> (std::istream& is_, Region<T>& reg_)
{
	if (is_.fail()) 
		return is_;
	is_ >> reg_.material;
	is_ >> reg_.syz_isxz;
	return is_;
}

#endif