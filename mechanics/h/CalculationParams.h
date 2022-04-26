#ifndef CALCULATIONPARAMS_H
#define CALCULATIONPARAMS_H
#include <ostream>
#include <istream>
#include <limits>

struct CalculationParams 
{
	CalculationParams();
	CalculationParams(unsigned int n1f5,unsigned int n1f6,
					unsigned int qn1f5,unsigned intqn1f6,double eps_,int ue_);
	~CalculationParams() = default;
	CalculationParams(const CalculationParams& ) = default;
	CalculationParams& operator= (const CalculationParams& ) = default; 
	CalculationParams(CalculationParams&& ) = default;
	CalculationParams& operator= (CalculationParams&& ) = default;
	
	unsigned int N1_f5;
	unsigned int N1_f6;
	unsigned int N2_f5;
	unsigned int N2_f6;
	double eps;
	int use_epsilon;
	friend std::ostream& operator<<(std::ostream& os_, const CalculationParams& cp_);
	friend std::istream& operator>>(std::istream& is_, CalculationParams& cp_);

}; 

CalculationParams::CalculationParams():N1_f5{40},N1_f6{40},
					N2_f5{512},N2_f6{512},
					eps{0.0001},use_epsilon(1)
{
	
}

CalculationParams::CalculationParams(unsigned int n1f5,unsigned int n1f6,
				unsigned int qn2f5,unsigned int qn2f6,double eps_,int ue_):N1_f5{n1f5},N1_f6{n1f6},
					N2_f5{qn2f5},N2_f6{qn2f6},
					eps{eps_},use_epsilon(ue_)
{
	
}

std::ostream& operator<< (std::ostream& os_, const CalculationParams& cp_)
{
	os_ << "N1f5 = " << cp_.N1_f5 << "\n";
	os_ << "N1f6 = " << cp_.N1_f6 << "\n";
	os_ << "N2f5_QUADRATURE = " << cp_.N2_f5 << "\n";
	os_ << "N2f6_QUADRATURE = " << cp_.N2_f6 << "\n";
	os_ << "epsilon_QUADRATURE = " << cp_.eps << "\n";
	os_ << "use_epsilon = " << cp_.use_epsilon << "\n";
	return os_;
}

std::istream& operator>> (std::istream& is_, CalculationParams& cp_)
{
	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};

	if (is_.fail()) 
		return is_;

	readValueFromSymbol(cp_.N1_f5);
	readValueFromSymbol(cp_.N1_f6);
	readValueFromSymbol(cp_.N2_f5);
	readValueFromSymbol(cp_.N2_f6);
	readValueFromSymbol(cp_.eps);
	readValueFromSymbol(cp_.use_epsilon);
	return is_;
}

#endif