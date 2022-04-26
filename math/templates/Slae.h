#ifndef SLAE_H
#define SLAE_H
#include "Matrix.h"

	template<typename T>
	class Slae 
	{
	public:
		Slae();
		Slae(size_t range_);
		void init(const Matrix<T> &a_, const Row<T> &b_);
		void initMatrix(const Matrix<T> &a_);
		void initFreeCol(const Row<T> &b_);

		void init(Matrix<T> &&a_, Row<T> &&b_);
		void initMatrix(Matrix<T> &&a_);
		void initFreeCol(Row<T> &&b_);

		void clear();
		void sortTriangular();
		void solveTriangular();
		void setEps(long double eps_);
		long double getEps();
		bool singleSolution();
		bool isSorted();
		
		Row<T> getResult();
	private:
		Row<T> result_; 
		bool singleSolution_;
		void slaeInit();
		std::vector<std::pair<size_t, size_t> > getSwapColHistory();
		void swapRows(size_t row_one, size_t row_two);
		void swapCols(size_t col_one, size_t col_two);
		bool checkMatrix();
		bool is_sorted;
		Matrix<T> a;
		Row<T> b;
		std::vector<std::pair<size_t, size_t> > swap_col_history;
		long double eps;
		template<typename T2>
		friend std::ostream &operator<<(std::ostream &out, const Slae<T2> &slae_);
		void zeroNearZero();
	};

	template<typename T>
	bool Slae<T>::singleSolution()
	{
		return singleSolution_;
	}

	template<typename T>
	bool Slae<T>::isSorted()
	{
		return is_sorted;
	}

	template<typename T>
	Row<T> Slae<T>::getResult()
	{
			return result_;
	}

	template<typename T>
	void Slae<T>::zeroNearZero()
	{
		for (size_t i = 0; i < a.size(); ++i ) 
		{
			for (size_t j = 0; j < a[i].size(); ++j ) 
			{
				if ( (a[i][j]>T()?a[i][j]:-a[i][j]) < eps)
					a[i][j] = T();
			}
		}

		for (size_t i = 0; i < b.size(); ++i ) 
		{
			if ( (b[i]>T()?b[i]:-b[i]) < eps)
				b[i] = T();		
		}

	};
	
	template<typename T>
	void Slae<T>::setEps(long double eps_)
	{
		eps = eps_;
	}
	
	template<typename T>
	long double Slae<T>::getEps()
	{
		return eps;
	}

	template<typename T>
	void Slae<T>::slaeInit()
	{
		eps = 1.e-15;	
		is_sorted = false;
		singleSolution_ = false;
	}

	template<typename T>
	Slae<T>::Slae()
	{
		slaeInit();
	}

	template<typename T>
	Slae<T>::Slae(size_t range_) : a(range_, range_), 
									b(range_) 
	{ 
		slaeInit();
	};

	template<typename T>
	void Slae<T>::init(const Matrix<T> &a_, const Row<T> &b_) 
	{
		a = a_;
		b = b_;
		is_sorted = false;
	};

	template<typename T>
	void Slae<T>::initMatrix(const Matrix<T> &a_) 
	{
		a = a_; 
		is_sorted = false; 
	};

	template<typename T>
	void Slae<T>::initFreeCol(const Row<T> &b_) 
	{
		b = b_; 
		is_sorted = false; 
	};
	
	template<typename T>
	void Slae<T>::init(Matrix<T> &&a_, Row<T> &&b_)
	{
		a = std::move(a_);
		b = std::move(b_);
		is_sorted = false;		
	}

	template<typename T>	
	void Slae<T>::initMatrix(Matrix<T> &&a_)
	{
		a = std::move(a_); 
		is_sorted = false; 		
	}

	template<typename T>
	void Slae<T>::initFreeCol(Row<T> &&b_)
	{
		b = std::move(b_); 
		is_sorted = false; 		
	}

	template<typename T>
	bool Slae<T>::checkMatrix() 
	{
		if ((a.size() == 0)||(b.size() == 0 )) 
		{
			return false;
		}
		if ( a.size() != b.size()) 
		{
			return false;
		}
		for (size_t i = 0; i < a.size(); ++i) 
		{
			if (a[i].size() != a.size()) 
			{
				return false;
			}
		}
		return true;
	};

	template<typename T>
	void Slae<T>::swapRows(size_t row_one, size_t row_two) 
	{
		
		std::swap(a[row_one], a[row_two]);
		std::swap(b[row_one], b[row_two]);
		is_sorted = false;
	}

	template<typename T>
	void Slae<T>::swapCols(size_t col_one, size_t col_two) 
	{
		for (size_t i = 0; i < a.size(); ++i ) 
		{
			std::swap(a[i][col_one], a[i][col_two]);
		}
		swap_col_history.push_back({ col_one, col_two });
		is_sorted = false;
	}

	template<typename T>
	std::vector<std::pair<size_t, size_t> >  Slae<T>::getSwapColHistory() 
	{
		return swap_col_history;
	}

	template<typename T>
	void Slae<T>::sortTriangular() 
	{
		if (checkMatrix()) 
		{
			const size_t &max_y = a.size();
			const size_t &max_x = max_y;

			for (size_t i = 0; i < max_y; ++i) 
			{			
				auto pos = a.findMaxElementByModule({ i, 0 }, { max_y, max_x });
				
				if ( (a[pos.first][pos.second]>T()?a[pos.first][pos.second]:-a[pos.first][pos.second] ) < eps)
				{
					break;
				}

				if (pos.first != i) 
				{
					swapRows(i, pos.first);
				}
				if (pos.second != i) 
				{
					swapCols(i, pos.second);
				}
			
				auto a_ii = a[i] / (-a[i][i]);
				T b_i = b[i] / (-a[i][i]);
				for (size_t j = i + 1; j < max_y; ++j) 
				{
					T a_ji = a[j][i];
					b[j] = b_i * a_ji + b[j];
					a[j] = a_ii * a_ji + a[j];
				}
			}
			is_sorted = true;
		}
		else 
		{
			is_sorted = false;
		}	
	}

	template<typename T>
	void Slae<T>::solveTriangular() 
	{
		if (is_sorted) 
		{
			size_t range = a.size();
 
			result_ = std::move(Row<T>(range));

			size_t i = range - 1;
			if ( (a[i][i]>T()?a[i][i]:-a[i][i]) < eps) 
			{
				singleSolution_ = false;
				if ( (b[i]>T()?b[i]:-b[i])< eps)
				{					
					for (size_t j = 0; j< result_.size(); ++j)
						result_[j] = 1;
				}
				return;
			}
			
			result_[i] = b[i] / a[i][i];
			do 
			{
				--i;
				result_[i] = b[i];
				for (size_t j = i + 1; j < range; ++j) 
				{
					result_[i] -= result_[j] * a[i][j];
				}
				
				if ( (a[i][i]>T()?a[i][i]:-a[i][i]) < eps ) 
				{
					return;
				}
				result_[i] /= a[i][i];
				
			} while ( i != 0 );
			result_ = swapByPatternBackward(result_, swap_col_history);
			singleSolution_ = true;
			return;
		}
		else 
		{
			return;
		}
	}

	template<typename T2>
	std::ostream &operator<<(std::ostream &out, const Slae<T2> &slae_) 
	{
		for (size_t i = 0; i < slae_.a.size(); ++i) 
		{
			out << slae_.a[i] << ' ';
			out << slae_.b[i] << '\n';
		}
		return out;
	}
 
#endif