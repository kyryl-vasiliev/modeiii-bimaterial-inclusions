#ifndef FORMATE_MATRIX_H
#define FORMATE_MATRIX_H

#include <numeric>
#include "../../math/templates/Slae.h"
#include "../../mechanics/templates/BodyPiecewiseHomogeneous.h"
#include "InclusionIntegrals.h"
template <typename T>
void formateMatrix(Slae<T> & s_, const BodyPiecewiseHomogeneous<T> &b_) 
{
	int vkl_count = b_.inclusions.size();
	int dim = std::accumulate(std::begin(b_.inclusions), std::end(b_.inclusions), 0,
		[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
	int dim2 = dim - 2 * vkl_count;
	Matrix<T> m(dim, dim);
	Row<T> m_b(dim);
	std::vector<T> s_f5, s_f6;
	for (int h_vkl = 0; h_vkl < vkl_count; h_vkl++) 
	{
		int N1h_f5 = b_.inclusions[h_vkl].calc_params.N1_f5;
		int N1h_f6 = b_.inclusions[h_vkl].calc_params.N1_f6;
		s_f5.resize(N1h_f5);
		for (int i = 1; i <= N1h_f5; ++i) 
		{
			s_f5[i - 1] = cosl(M_PI*i / (N1h_f5 + 1.0));
		}

		s_f6.resize(N1h_f6);
		for (int i = 1; i <= N1h_f6; ++i) {
			s_f6[i - 1] = cosl(M_PI*i / (N1h_f6 + 1.0));
		}

		int currHvkl_f5 = std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + h_vkl, 0,
			[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + incl.calc_params.N1_f6; });
		int currHvkl_f6 = currHvkl_f5 + b_.inclusions[h_vkl].calc_params.N1_f5;

		for (int w_vkl = 0; w_vkl < vkl_count; ++w_vkl) 
		{
			int N1w_f5 = b_.inclusions[w_vkl].calc_params.N1_f5 + 1;
			int N1w_f6 = b_.inclusions[w_vkl].calc_params.N1_f6 + 1;
			int currWvkl_f5 = std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + w_vkl, 0,
				[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
			int currWvkl_f6 = currWvkl_f5 + b_.inclusions[w_vkl].calc_params.N1_f5 + 1;
			for (int h_N1_f5 = 0; h_N1_f5 < N1h_f5; ++h_N1_f5) 
			{
				for (int w_N1 = 0; w_N1 < N1w_f5; ++w_N1) 
				{
					if (w_vkl == h_vkl) 
					{
						m[currHvkl_f5 + h_N1_f5][currWvkl_f5 + w_N1] =
							mainDiagonalSyz_f5<T>(w_N1, s_f5[h_N1_f5], h_vkl, b_);
					}
					else 
					{
						m[currHvkl_f5 + h_N1_f5][currWvkl_f5 + w_N1] =
							nonMainDiagonalSyz_f5<T>(w_N1, s_f5[h_N1_f5], h_vkl, w_vkl, b_); 
					}
				}
				for (int w_N1 = 0; w_N1 < N1w_f6; ++w_N1) 
				{	
					if (w_vkl == h_vkl) 
					{
						m[currHvkl_f5 + h_N1_f5][currWvkl_f6 + w_N1] =
							mainDiagonalSyz_f6<T>(w_N1, s_f5[h_N1_f5], h_vkl, b_);
					}
					else 
					{
						m[currHvkl_f5 + h_N1_f5][currWvkl_f6 + w_N1] =
							nonMainDiagonalSyz_f6<T>(w_N1, s_f5[h_N1_f5], h_vkl, w_vkl, b_);
					}
				}
				
				const auto & a44_incl = b_.inclusions[h_vkl].material.a44;
				const auto & a45_incl = b_.inclusions[h_vkl].material.a45;
				const auto & a55_incl = b_.inclusions[h_vkl].material.a55;
				const auto & phi_incl = b_.inclusions[h_vkl].phi;
				const auto & h_incl = b_.inclusions[h_vkl].h;

				const auto & tau_plus_incl = b_.inclusions[h_vkl].tau_plus;
				const auto & tau_minus_incl = b_.inclusions[h_vkl].tau_minus;

				T ssz0, snz0, snz_area_phi;
				AnisotropyAngle<T> Consts;
				if (b_.inclusions[h_vkl].center.y > 0)
				{
					ssz0 = b_.areaHi.syz_isxz.getSsz(phi_incl);
					snz0 = b_.areaHi.syz_isxz.getSnz(phi_incl);
					Consts.init(b_.areaHi.material, phi_incl);
					snz_area_phi = b_.areaHi.syz_isxz.getSnz(phi_incl);
				}
				else
				{
					ssz0 = b_.areaLow.syz_isxz.getSsz(phi_incl);
					snz0 = b_.areaLow.syz_isxz.getSnz(phi_incl);
					Consts.init(b_.areaLow.material, phi_incl);
					snz_area_phi = b_.areaLow.syz_isxz.getSnz(phi_incl);
				}

				const auto  a44_mat = Consts.get_a44j();
				const auto  a45_mat = Consts.get_a45j();
				const auto  a55_mat = Consts.get_a55j();

				auto  ssz0_mid = ssz0* a44_mat / fmaxl(a44_mat, a44_incl);
				auto w_ = 2 * h_incl*( (a44_mat*(snz0 + (tau_plus_incl + tau_minus_incl)/2.0 ) + a45_mat * ssz0_mid) )* fminl(a44_mat, a44_incl) / a44_mat;
			
				m_b[currHvkl_f5 + h_N1_f5] = - 2.0 * a45_incl / a44_incl * ssz0_mid + w_/h_incl/a44_incl
					- 2 * snz_area_phi - (tau_plus_incl + tau_minus_incl)
					- nonSymmetricalLoadSyz_f5(s_f5[h_N1_f5], h_vkl, b_);
			}
			for (int h_N1_f6 = 0; h_N1_f6 < N1h_f6; ++h_N1_f6) 
			{
				for (int w_N1 = 0; w_N1 < N1w_f5; ++w_N1) 
				{
					if (w_vkl == h_vkl) 
					{
						m[currHvkl_f6 + h_N1_f6][currWvkl_f5 + w_N1] =
							mainDiagonalW_f5<T>(w_N1, s_f6[h_N1_f6], h_vkl, b_);
					}
					else 
					{
						m[currHvkl_f6 + h_N1_f6][currWvkl_f5 + w_N1] = 
							nonMainDiagonalW_f5<T>(w_N1, s_f6[h_N1_f6], h_vkl, w_vkl, b_);
					}
				}
				for (int w_N1 = 0; w_N1 < N1w_f6; ++w_N1) 
				{
					if (w_vkl == h_vkl) 
					{
						m[currHvkl_f6 + h_N1_f6][currWvkl_f6 + w_N1] =
							mainDiagonalW_f6<T>(w_N1, s_f6[h_N1_f6], h_vkl, b_);
					}
					else 
					{
						m[currHvkl_f6 + h_N1_f6][currWvkl_f6 + w_N1] =
							nonMainDiagonalW_f6<T>(w_N1, s_f6[h_N1_f6], h_vkl, w_vkl, b_);
					}
				}
				const auto & a44_incl = b_.inclusions[h_vkl].material.a44;
				const auto & a45_incl = b_.inclusions[h_vkl].material.a45;
				const auto & a55_incl = b_.inclusions[h_vkl].material.a55;
				const auto & phi_incl = b_.inclusions[h_vkl].phi;
				const auto & w_plus_incl = b_.inclusions[h_vkl].w_plus;
				const auto & w_minus_incl = b_.inclusions[h_vkl].w_minus;
				const auto & h_incl = b_.inclusions[h_vkl].h;
		
				T ssz0, snz0, snz_area_phi;
				AnisotropyAngle<T> Consts;
				if (b_.inclusions[h_vkl].center.y > 0)
				{
					ssz0 = b_.areaHi.syz_isxz.getSsz(phi_incl);
					snz0 = b_.areaHi.syz_isxz.getSnz(phi_incl);
					Consts.init(b_.areaHi.material, phi_incl);
				}
				else
				{
					ssz0 = b_.areaLow.syz_isxz.getSsz(phi_incl);
					snz0 = b_.areaLow.syz_isxz.getSnz(phi_incl);
					Consts.init(b_.areaLow.material, phi_incl);
				}

				const auto a44_mat = Consts.get_a44j();
				const auto a45_mat = Consts.get_a45j();
				const auto a55_mat = Consts.get_a55j();

				T r0_2 = fabsl(a44_incl * a55_incl - a45_incl * a45_incl);

				const auto & tau_plus_incl = b_.inclusions[h_vkl].tau_plus;
				const auto & tau_minus_incl = b_.inclusions[h_vkl].tau_minus;

				auto  ssz0_mid = ssz0* a44_mat / fmaxl(a44_mat, a44_incl);
				auto w_ = 2 * h_incl*( (a44_mat*(snz0 + (tau_plus_incl + tau_minus_incl)/2.0 ) + a45_mat * ssz0_mid) )*fminl(a44_mat, a44_incl) / a44_mat;

				m_b[currHvkl_f6 + h_N1_f6] = 2.0 * r0_2 / a44_incl * ssz0_mid + a45_incl / a44_incl * w_ / h_incl 
					- 2 * (snz0 * a45_mat + ssz0 * a55_mat)
					- nonSymmetricalLoadW_f5(s_f6[h_N1_f6], h_vkl, b_);
			}
		}
		int currWvkl_f5 = std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + h_vkl, 0,
			[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
		int currWvkl_f6 = currWvkl_f5 + b_.inclusions[h_vkl].calc_params.N1_f5 + 1;

		for (int i = 0; i < N1h_f5 + 1; ++i) {
			m[dim2 + h_vkl * 2][currWvkl_f5 + i] = int_ax(i, b_.inclusions[h_vkl].a, b_.inclusions[h_vkl].a);
		}
		m_b[dim2 + h_vkl * 2] = - 2 * b_.inclusions[h_vkl].a *( b_.inclusions[h_vkl].tau_minus - b_.inclusions[h_vkl].tau_plus );

		for (int i = 0; i < N1h_f6 + 1; ++i) {
			m[dim2 + h_vkl * 2 + 1][currWvkl_f6 + i] = int_ax(i, b_.inclusions[h_vkl].a, b_.inclusions[h_vkl].a);
		}
		m_b[dim2 + h_vkl * 2 + 1] = b_.inclusions[h_vkl].h * 0;
	}
	s_.init(std::move(m), std::move(m_b) );
}

template<typename T>
std::vector< std::vector<T> > calculateSIF(const BodyPiecewiseHomogeneous<T>& b_, const Row<T> &res_) {
	std::vector< std::vector<T> > ret(b_.inclusions.size(), std::vector<T>(4, T()));
	std::vector<T> tmp = std::vector<T>(res_);
	for (size_t i = 0; i < b_.inclusions.size(); ++i) 
	{
		int currWvkl_f5 = std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + i, 0,
			[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
		int currWvkl_f6 = currWvkl_f5 + b_.inclusions[i].calc_params.N1_f5 + 1;
		T K32_B = std::accumulate(std::begin(tmp) + currWvkl_f5,
			std::begin(tmp) + currWvkl_f5 + b_.inclusions[i].calc_params.N1_f5 + 1, T(),
			[](T sum, T x1) { return sum + x1; });
		T K31_B = std::accumulate(std::begin(tmp) + currWvkl_f6,
			std::begin(tmp) + currWvkl_f6 + b_.inclusions[i].calc_params.N1_f6 + 1, T(),
			[](T sum, T x1) { return sum + x1; });

		int v = -1;
		T K32_A = std::accumulate(std::begin(tmp) + currWvkl_f5,
			std::begin(tmp) + currWvkl_f5 + b_.inclusions[i].calc_params.N1_f5 + 1, T(),
			[&v](T sum, T x1) { v = -v; return sum + v * x1; });
		v = -1;
		T K31_A = std::accumulate(std::begin(tmp) + currWvkl_f6,
			std::begin(tmp) + currWvkl_f6 + b_.inclusions[i].calc_params.N1_f6 + 1, T(),
			[&v](T sum, T x1) { v = -v;  return sum + v * x1; });

		AnisotropyAngle<T> ConstsHi, ConstsLow ;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[i].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[i].phi);

		ret[i][2] = K32_A / 2.0;
		ret[i][3] = - K32_B / 2.0;
		if (b_.inclusions[i].center.y > 0)
		{
			ret[i][0] = -K31_A / ConstsHi.get_a55j() / ConstsHi.get_alphaj() / 2.0;
			ret[i][1] = K31_B / ConstsHi.get_a55j() / ConstsHi.get_alphaj() / 2.0;
		}
		else
		{
			ret[i][0] = -K31_A / ConstsLow.get_a55j() / ConstsLow.get_alphaj() / 2.0;
			ret[i][1] = K31_B / ConstsLow.get_a55j() / ConstsLow.get_alphaj() / 2.0;
		}
	}
	return ret;
}

template<typename T>
bool calculateStressInArea(const BodyPiecewiseHomogeneous<T>& b_, const Row<T>& row_, std::istream& is_, std::ostream& os_)
{
	if (is_.fail())
		return false;
	std::vector<T> Af5f6 = row_;
	std::vector<T> x;
	std::vector<T> y;
	auto processToVector = [&is_](std::vector<T>& vctr_) {
		T val;
		std::string line;
		getline(is_, line);
		std::istringstream iss(line);
		while (iss >> val)
			vctr_.push_back(val);
	};
	processToVector(x);
	processToVector(y);
	for (const auto& xx : x)
	{
		for (const auto& yy : y)
		{
			auto res = Int_Sum_Snzj_plus_ISszj(xx, yy, b_, Af5f6);
			if (res)
			{
				os_ << xx << " " << yy << " " << (*res).real() << " " << (*res).imag() << '\n';
			}
		}
	}
	return true;
}

#endif