#ifndef STRESS_H
#define STRESS_H

#include <complex>
#include <ostream>
#include <istream>
#include <limits>

using namespace std::complex_literals;

template<typename T>
class Stress 
{
public:
	Stress() = default;
	~Stress() = default;
	Stress(	const Stress&) = default;
	Stress& operator= (const Stress&) = default;
	Stress(Stress&&) = default;
	Stress& operator= (Stress&&) = default;

	Stress(T syz_, T sxz_) : Syz0_plus_iSxz0 { syz_, sxz_ }
	{		
	}
	void init(T syz_, T sxz_) 
	{ 
		Syz0_plus_iSxz0 = { syz_, sxz_ }; 
	}
	std::complex<T> getSnz_plus_iSsz_(T alpha_) const 
	{
		const std::complex<T> I = std::complex<T>(1.0i);
		return Syz0_plus_iSxz0 * exp(I*alpha_);
	}

	T getSnz(T alpha_) const 
	{
		const std::complex<T> I = std::complex<T>(1.0i);
		return (Syz0_plus_iSxz0 * exp(I*alpha_)).real();
	}

	T getSsz(T alpha_) const 
	{
		const std::complex<T> I = std::complex<T>(1.0i);
		return (Syz0_plus_iSxz0 * exp(I*alpha_)).imag();
	}

	template<typename U>
	friend std::ostream& operator<<(std::ostream& os_, const Stress<U>& str_);
	template<typename U>
	friend std::istream& operator>> (std::istream& is_, Stress<U>& str_);

private:
	std::complex<T> Syz0_plus_iSxz0;
};

template<typename T>
std::ostream& operator<< (std::ostream& os_, const Stress<T>& str_)
{
	os_ << "syz_odn = " << str_.Syz0_plus_iSxz0.real() << "\n";
	os_ << "sxz_odn = " << str_.Syz0_plus_iSxz0.imag() << "\n";
	return os_; 
}

template<typename T>
std::istream& operator>> (std::istream& is_, Stress<T>& str_)
{
	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};

	if (is_.fail()) 
		return is_;
	T syz = 0, sxz = 0;
	readValueFromSymbol(syz);
	readValueFromSymbol(sxz);
	str_.init(syz, sxz);
	return is_;
}

#endif