#include <iostream>
#include <fstream>
#include <numeric>
#include "../../math/templates/Slae.h"
#include "../../mechanics/templates/BodyPiecewiseHomogeneous.h"
#include "../../mechanics/templates/BasicAnisotropy.h"
#include "../../mechanics/templates/Inclusion.h"
#include "../../mechanics/h/InputParams.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "../templates/FormateSlaeInclusion.h"

using std::vector;

bool loadBodyFromFile(const char*filename_, BodyPiecewiseHomogeneous< long double> &b)
{
	BodyPiecewiseHomogeneous< long double> tmp;
	std::ifstream f;
	f.open(filename_);
	if (!f) 
	{
		return false;
	}
	f >> tmp;
	if (f.good())
	{
		b = std::move(tmp);
		f.close();
		return true;
	}
	else
	{
		f.close();
		return false;		
	}
}
void saveBodyToFile(const char*filename_, BodyPiecewiseHomogeneous< long double>& b)
{
	std::ofstream f;
	f.open(filename_);
	f << "material hi a44 a55 a45" << '\n';
	f << b.areaHi.material.a44 <<" "<<b.areaHi.material.a55 <<" "<<b.areaHi.material.a45<<'\n';
	f << "material low a44 a55 a45" << '\n';
	f << b.areaLow.material.a44 << " " << b.areaLow.material.a55 << " " << b.areaLow.material.a45 << '\n';

	for (size_t i = 0; i < b.inclusions.size(); ++i) {
		f << "inclusion number _ " << i << "_" << '\n';
		f << "phi " << b.inclusions[i].phi << " x0 " << b.inclusions[i].center.x << " y0 " << b.inclusions[i].center.y << '\n';
		f << "a " << b.inclusions[i].a << " h " << b.inclusions[i].h << '\n';
		f << "tau+ " << b.inclusions[i].tau_plus << " tau- " << b.inclusions[i].tau_minus<< '\n';
		f << "Nf5 " << b.inclusions[i].calc_params.N1_f5 << " Nf6 " << b.inclusions[i].calc_params.N1_f6 << '\n';
		f << "inclusion a44 a55 a45" << '\n';
		f << b.inclusions[i].material.a44 << " " << b.inclusions[i].material.a55 << " " << b.inclusions[i].material.a45 << '\n' << '\n';
	}
	f.close();
	
}

bool loadInputFromFile(const char*filename_, InputParams &ip)
{	
	InputParams tmp;
	std::ifstream f;
	f.open(filename_);
	if (!f) 
	{
		return false;
	}
	f >> tmp;
	if (f.good())
	{
		ip = std::move(tmp);
		f.close();
		return true;
	}
	else
	{
		f.close();
		return false;		
	}

}

void madeCalculation(BodyPiecewiseHomogeneous< long double> &b, std::vector< std::vector<long double> >& sif)
{

	Slae< long double> slae_;
	formateMatrix(slae_, b);
	slae_.sortTriangular(); 
	slae_.solveTriangular();

	if (slae_.singleSolution())
	sif = calculateSIF(b, slae_.getResult());
	std::cout <<"solution:\n";
	for (size_t i = 0; i < b.inclusions.size(); ++i) 
	{
		std::cout << sif[i][0] << " " << sif[i][1] << " " 
				  << sif[i][2] << " "<< sif[i][3] << '\n';
	}
	return;
}

void saveSif(std::vector< std::vector<long double> > sif, std::string basicFilename)
{
	std::ofstream file;
	std::string fname = "";


	for (size_t i = 0; i < sif.size(); ++i)
	{
		fname = basicFilename + std::to_string(int(i)) + std::string(".txt");
		file.open(fname, std::ios::app);

		file << sif[i][0] << " " << sif[i][1] << " "
			<< sif[i][2] << " " << sif[i][3] << '\n';
		file.close();
	}
}

int main(int argc, char** argv)
{
	BodyPiecewiseHomogeneous< long double> b;
	InputParams ip;
	BasicAnisotropy< long double> space_material{ 1, 0, 1 };

	const char*constDataFilename_ = "one.txt";
	const char*changeDataFilename_ = "input.txt";

	if(loadBodyFromFile(constDataFilename_, b))
	{	
		std::cout << b;
		long double * changed_param1 = nullptr;
		long double * changed_param2 = nullptr;

		if (loadInputFromFile(changeDataFilename_, ip) )
		{
			int indx = ip.inclusion_index;
			if (ip.changed_param == "a")
			{
				changed_param1 = &(b.inclusions[indx].a);
			}
			else
			if (ip.changed_param == "h")
			{
				changed_param1 = &(b.inclusions[indx].h);
			}
			else
			if (ip.changed_param == "x0")
			{
				changed_param1 = &(b.inclusions[indx].center.x);
			}
			else
			if (ip.changed_param == "y0")
			{
				changed_param1 = &(b.inclusions[indx].center.y);
			}
			else
			if (ip.changed_param == "phi")
			{
				changed_param1 = &(b.inclusions[indx].phi);
			}
			else
			if (ip.changed_param == "a44")
			{
				changed_param1 = &b.inclusions[indx].material.a44;
			}
			else
			if (ip.changed_param == "a45")
			{
				changed_param1 = &b.inclusions[indx].material.a45;
			}
			else
			if (ip.changed_param == "a55")
			{
				changed_param1 = &b.inclusions[indx].material.a55;
			}
			else
			if (ip.changed_param == "a44a55")
			{
				changed_param1 = &b.inclusions[indx].material.a44;
				changed_param2 = &b.inclusions[indx].material.a55;
			}
			else
			if (ip.changed_param == "ma44hi")
			{
				changed_param1 = &b.areaHi.material.a44;
			}
			else
			if (ip.changed_param == "ma45hi")
			{
				changed_param1 = &b.areaHi.material.a45;
			}
			else
			if (ip.changed_param == "ma55hi")
			{
				changed_param1 = &b.areaHi.material.a55;
			}
			else
			if (ip.changed_param == "ma44a55hi")
			{
				changed_param1 = &b.areaHi.material.a44;
				changed_param2 = &b.areaHi.material.a55;
			}
			else
			if (ip.changed_param == "ma44low")
			{
				changed_param1 = &b.areaLow.material.a44;
			}
			else
			if (ip.changed_param == "ma45low")
			{
				changed_param1 = &b.areaLow.material.a45;
			}
			else
			if (ip.changed_param == "ma55low")
			{
				changed_param1 = &b.areaLow.material.a55;
			}
			else
			if (ip.changed_param == "ma44a55low")
			{
				changed_param1 = &b.areaLow.material.a44;
				changed_param2 = &b.areaLow.material.a55;
			}
			else
			if (ip.changed_param == "a1a2")
			{
				changed_param1 = &(b.inclusions[1].a);
				changed_param2 = &(b.inclusions[2].a);
			}
			else
			if (ip.changed_param == "x12")
			{

			}
			else
			if (ip.changed_param == "wedgea12")
			{

			}
			else
			if (ip.changed_param == "wedgeeps12")
			{
			}
			else
			{
				std::cout << "error loop parameter name";
				exit(1);
			}
		}
			
		for(const auto &i: ip.points_)
		{
			if (changed_param1!= nullptr)
				*changed_param1 = i;
			if (changed_param2 != nullptr)
			{
				*changed_param2 = i;
			}
			if (ip.changed_param == "x12")
			{
				b.inclusions[1].center.x = i;
				b.inclusions[2].center.x = -i;
			}
			
			if (ip.changed_param == "wedgea12") 
			{
				double eps = 0.001;
				b.inclusions[1].a = i;
				b.inclusions[2].a = i;

				b.inclusions[1].center.x = (b.inclusions[1].a + eps)*cosl(b.inclusions[1].phi);
				b.inclusions[1].center.y = (b.inclusions[1].a + eps) * sinl(b.inclusions[1].phi);
				
				b.inclusions[2].center.x = b.inclusions[1].center.x;
				b.inclusions[2].center.y = - b.inclusions[1].center.y;
				b.inclusions[2].phi = -b.inclusions[1].phi;
			}

			if (ip.changed_param == "wedgeeps12") 
			{
				double eps = i;
				b.inclusions[1].center.x = (b.inclusions[1].a + eps) * cosl(b.inclusions[1].phi);
				b.inclusions[1].center.y = (b.inclusions[1].a + eps) * sinl(b.inclusions[1].phi);
				b.inclusions[2].center.x = b.inclusions[1].center.x;
				b.inclusions[2].center.y = -b.inclusions[1].center.y;
				b.inclusions[2].phi = -b.inclusions[1].phi;
			}


			std::vector< std::vector<long double> >  sif;
			madeCalculation(b, sif);
			saveSif(sif, "inclusion_"+ ip.changed_param);
		}
	}
	return 0;
}