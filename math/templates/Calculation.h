#ifndef CALCULATION_H
#define CALCULATION_H
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>
#include <vector>
#include <complex>
#include <algorithm>

template<typename T>
T Tn(int n, T x) 
{
	return cosl(n*acosl(x));
}

template<typename T>
T Un(int n, T x) 
{
	return sinl( (n + 1.0)*acosl(x) ) / sinl(acosl(x));
}

template<typename T, typename X, typename ... Y>
T gaussQuadrature(int N, T (*f)(X t, const Y& ... y1), const Y& ... y1) 
{
	std::vector<X> tk(N, X());
	for (int k = 0; k < N; ++k) 
	{
		tk[k] = cosl((2.0*(k + 1.0) - 1.0)*M_PI / 2 / N);
	}
	return T(M_PI / N) * std::accumulate(std::begin(tk), std::end(tk), T(), 
							[&](T x1, X x2) { return x1 + f(x2, y1...); });
}

template<typename T> 
long double distance(T t1, T t2) 
{
	return fabsl(t1 - t2);
}

template<typename T>
long double distance(std::complex<T> t1, std::complex<T> t2) 
{
	return sqrtl((t1.real() - t2.real())*(t1.real() - t2.real()) +
		(t1.imag() - t2.imag())*(t1.imag() - t2.imag()));
}


template<typename T, typename X, typename ... Y>
T gaussQuadratureEpsilon(double eps, size_t N_max, T(*f)(X t, const Y& ... y1), const Y& ... y1) 
{
	std::vector<long double> v_epsilon;
	std::vector<T> v_result;
	
	size_t N = 4;
	long double dist = 0;
	T r1, r2;
	r1 = gaussQuadrature<T>(N, f, y1 ...);
	do 
	{
		r2 = gaussQuadrature<T>(2*N, f, y1 ...);
		dist = distance(r1, r2);
		if (  dist < eps) 
		{
			return r2;
		}
		else 
		{
			v_epsilon.push_back(dist);
			v_result.push_back(r2);
			N = 2 * N;
			r1 = r2;
		}
	} while (N <= N_max*2);
	
	int min_element_index = std::min_element(v_epsilon.begin(),v_epsilon.end()) - v_epsilon.begin();
	return v_result[min_element_index];
}

template<typename T, typename ... Y>
T integralTrapetion(int N, T a, T b, T(*f)(T t, const Y& ... y), const Y& ... y) 
{
	T integral = 0.0;
	T h = (b - a) / N;

	for (int i = 1; i < N; ++i)
	{
		integral = integral + f(a + h * i, y...);
	}
	integral = h * (integral + (f(a, y...) + f(b, y...)) / 2);
	return integral;
}

template<typename T, typename ... Y>
T integralSimpson(int N, T a, T b, T(*f)(T t, const Y& ... y), const Y& ... y) 
{
	T integral = 0;
	T h = (b - a) / N / 2;

	for (int i = 1; i <= N; ++i)
	{
		integral = integral + 4.0 * f(a + h * (2.0 * i - 1), y...);
	}

	for (int i = 2; i <= N; ++i)
	{
		integral = integral + 2.0 * f(a + h * (2.0 * i - 2.0), y...);
	}

	integral = h * (integral + f(a, y...) + f(b, y...)) / 3;
	return integral;
}

template<typename T, typename ... Y>
T integral38(int N, T a, T b, T(*f)(T t, const Y& ... y), const Y& ... y) 
{
	T integral = 0.0;
	T h = (b - a) / N / 3;

	for (int i = 1; i <= N; ++i)
	{
		integral = integral + 3.0 * f(a + h * (3.0 * i - 2.0), y...);
	}

	for (int i = 1; i <= N; ++i)
	{
		integral = integral + 3.0 * f(a + h * (3.0 * i - 1.0), y...);
	}

	for (int i = 2; i <= N; ++i)
	{
		integral = integral + 2.0 * f(a + h * (3.0 * i - 3.0), y...);
	}
	integral = 3.0 * h*(integral + f(a, y...) + f(b, y...)) / 8;
	return integral;
}

#endif