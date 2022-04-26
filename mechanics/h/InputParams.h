#ifndef INPUTPARAMS_H
#define INPUTPARAMS_H
#include <ostream>
#include <istream>
#include <limits>
#include <vector>

class InputParams 
{
	public:
		std::string changed_param;	
		std::vector<double> points_;
		int inclusion_index;		
	private:
		friend std::ostream& operator<<(std::ostream& os_, const InputParams& ip_);
		friend std::istream& operator>>(std::istream& is_, InputParams& cp_);

}; 

std::ostream& operator<< (std::ostream& os_, const InputParams& ip_)
{
	for (const auto& i:ip_.points_)
		os_<< i << '\n';
	return os_;
}


std::istream& operator>> (std::istream& is_, InputParams& ip_)
{
	auto readValueFromSymbol = [&is_](auto &x, char delim = '=')
	{ 
		is_.ignore(std::numeric_limits<std::streamsize>::max(), '=');
		is_ >> x;		
	};
	
	if (is_.fail()) 
		return is_;
	
	double pnt;
	std::string s_begin_end;
	readValueFromSymbol(ip_.inclusion_index);
	readValueFromSymbol(ip_.changed_param);
	is_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	is_ >> s_begin_end;
	if (s_begin_end != "begin")
	{
		is_.setstate(is_.badbit);
		return is_;
	};

	is_>> pnt;
	while(is_.good())
	{
		ip_.points_.push_back(pnt);		
		is_>> pnt;
	}
	is_.clear(is_.goodbit);
	is_ >> s_begin_end;
	if (s_begin_end != "end")
	{	
		is_.setstate(is_.badbit);
		return is_;

	};
	
	return is_;
}

#endif