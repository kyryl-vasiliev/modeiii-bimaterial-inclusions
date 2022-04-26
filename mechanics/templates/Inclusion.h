#ifndef INCLUSION_H
#define INCLUSION_H
#include "AnisotropyAngle.h"
#include "../h/CalculationParams.h"
#include <ostream>
#include <istream>
#include <limits>

template<typename T>
class Inclusion
{
public:
	struct XY 
	{
		T x;
		T y;
		XY(T x_, T y_)
		{
			x = x_;
			y = y_;
		}
		XY() :x(T()), y(T()) 
		{	
		};
		XY(const XY&) = default;
		XY& operator=(const XY&) = default;
		~XY() = default;
	};
public:
	XY center;
	T phi;
	T a;
	T h;
	T tau_plus;
	T tau_minus;
	T w_plus;
	T w_minus;
	
	BasicAnisotropy<T> material;
	CalculationParams calc_params;

	Inclusion() = default;
	~Inclusion() = default;
	Inclusion(const Inclusion &in) = default;
	Inclusion& operator=(const Inclusion &in) = default;
	Inclusion(Inclusion &&in) = default;
	Inclusion& operator=(Inclusion &&in) = default;

	void outInclusionGeometryForce(std::ostream& os_) const;
	void  inInclusionGeometryForce(std::istream& is_);
	template<typename U>
	friend std::ostream& operator<< (std::ostream& os_, const Inclusion<U>& incl_);
	template<typename U>
	friend std::istream& operator>> (std::istream& is_, Inclusion<U>& incl_);
};

template<typename T>
void  Inclusion<T>::outInclusionGeometryForce(std::ostream& os_) const
{
	os_ << "x0 = " << center.x << "\n";
	os_ << "y0 = " << center.y << "\n";
	os_ << "phi = " << phi << "\n";
	os_ << "a = " << a << "\n";
	os_ << "h = " << h << "\n";
	os_ << "tau+ = " << tau_plus << "\n";
	os_ << "tau- = " << tau_minus << "\n";
	os_ << "w+ = " << w_plus << "\n";
	os_ << "w- = " << w_minus << "\n";
}

template<typename T>
void  Inclusion<T>::inInclusionGeometryForce(std::istream& is_)
{
	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};
	
	readValueFromSymbol(center.x);
	readValueFromSymbol(center.y);
	readValueFromSymbol(phi);
	readValueFromSymbol(a);
	readValueFromSymbol(h);
	readValueFromSymbol(tau_plus);
	readValueFromSymbol(tau_minus);
	readValueFromSymbol(w_plus);
	readValueFromSymbol(w_minus);
}

template<typename T>
std::ostream& operator<< (std::ostream& os_, const Inclusion<T>& incl_)
{
	incl_.outInclusionGeometryForce(os_);
	os_ << incl_.material;
	os_ << incl_.calc_params;
	return os_; 
}

template<typename T>
std::istream& operator>> (std::istream& is_, Inclusion<T>& incl_)
{
	if (is_.fail()) 
		return is_;
	incl_.inInclusionGeometryForce(is_);
	is_ >> incl_.material;
	is_ >> incl_.calc_params;
	return is_;
}

#endif