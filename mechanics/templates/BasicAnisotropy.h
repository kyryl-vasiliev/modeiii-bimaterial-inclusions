#ifndef BASICANISOTROPY_H
#define BASICANISOTROPY_H
#include <ostream>
#include <istream>
#include <limits>

template<typename T>
struct BasicAnisotropy 
{
	T a44;
	T a45;
	T a55;
	BasicAnisotropy() = default;
	~BasicAnisotropy() = default;
	BasicAnisotropy(const BasicAnisotropy&) = default;
	BasicAnisotropy& operator= (const BasicAnisotropy&) = default;
	BasicAnisotropy(BasicAnisotropy&&) = default;
	BasicAnisotropy& operator= (BasicAnisotropy&&) = default;

	template<typename U>
    friend std::ostream& operator<<(std::ostream& os_, const BasicAnisotropy<U>& ba_);
	template<typename U>
	friend std::istream& operator>>(std::istream& is_, BasicAnisotropy<U>& ba_);

};

template<typename U>
std::ostream& operator<< (std::ostream& os_, const BasicAnisotropy<U>& ba_)
{
	os_ << "a44 = " << ba_.a44 << "\n";
	os_ << "a45 = " << ba_.a45 << "\n";
	os_ << "a55 = " << ba_.a55 << "\n";
	return os_;
}

template<typename U>
std::istream& operator>> (std::istream& is_, BasicAnisotropy<U>& ba_)
{
	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};
	if (is_.fail()) 
		return is_;
	readValueFromSymbol(ba_.a44);		
	readValueFromSymbol(ba_.a45);		
	readValueFromSymbol(ba_.a55);		
	return is_;
}

#endif