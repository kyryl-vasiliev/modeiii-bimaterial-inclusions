#ifndef ANISOTROPYANGLE_H
#define ANISOTROPYANGLE_H

#include <complex>
#include <math.h>

#include "BasicAnisotropy.h"

template<typename T>
class AnisotropyAngle 
{
public:
	T get_a44j() { return a44j; }
	T get_a45j() { return a45j; }
	T get_a55j() { return a55j; }
	T get_alphaj() { return alphaj; }
	T get_betaj() { return betaj; }

	std::complex<T> get_gpj() { return gpj; }
	std::complex<T> get_gmj() { return gmj; }

	void init(T a44, T a45, T a55, T phij) 
	{
		T cos_phij = cosl(phij);
		T sin_phij = sinl(phij);

		a44j = a44 * cos_phij*cos_phij - 2.0*a45*sin_phij*cos_phij + a55 * sin_phij*sin_phij;
		a45j = (a44 - a55)*sin_phij*cos_phij + a45 * (cos_phij*cos_phij - sin_phij * sin_phij);
		a55j = a44 * sin_phij*sin_phij + 2.0*a45*sin_phij*cos_phij + a55 * cos_phij*cos_phij;

		betaj = a45j / a55j;
		alphaj = sqrtl(a44j*a55j - a45j * a45j) / a55j;
		gpj = std::complex<T>(betaj, 1.0 + alphaj);
		gmj = std::complex<T>(betaj, 1.0 - alphaj);
	}
	
	void init(const BasicAnisotropy<T> &an, T phij) 
	{
		init(an.a44, an.a45, an.a55, phij);
	}
	
	AnisotropyAngle() = default;
	~AnisotropyAngle() = default;
	AnisotropyAngle(const AnisotropyAngle& ) = default;
	AnisotropyAngle& operator= (const AnisotropyAngle& ) = default;
	AnisotropyAngle(AnisotropyAngle&& ) = default;
	AnisotropyAngle& operator= (AnisotropyAngle&& ) = default;

private:
	T a44j;
	T a45j;
	T a55j;
	T alphaj;
	T betaj;
	std::complex<T> gpj;
	std::complex<T> gmj;
};

#endif