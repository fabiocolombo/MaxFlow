#ifndef FLOATUTILS_HPP
#define FLOATUTILS_HPP
#include <cmath>
#include <string>
#include <sstream>

#ifndef EPS
#define EPS 1e-6
#endif
namespace util{
	namespace floatUtils{
		template <class RealType> bool le(RealType a, RealType b, RealType eps=EPS){
			return a+eps < b;
		}

		template <class RealType> bool gr(RealType a, RealType b, RealType eps=EPS){
			return a-eps > b;
		}

		template <class RealType> bool leq(RealType a, RealType b, RealType eps=EPS){
			return !gr(a,b,eps);
		}

		template <class RealType> bool geq(RealType a, RealType b, RealType eps=EPS){
			return !le(a,b,eps);
		}

		template <class RealType> bool eq(RealType a, RealType b, RealType eps=EPS){
			return b-a < eps && a-b<eps;
		}

		template <class RealType> int roundInt(RealType a){
			return (int) floor(a + 0.5);
		}

		template <class NumType> NumType sq(NumType num){
			return num*num;
		}

		template <class RealType> double l2(RealType x1, RealType y1, RealType x2, RealType y2){
			return sqrt(static_cast<double>(sq(x1-x2) + sq(y1-y2)));
		}

		template <class RealType> double trunc(RealType x, int n){
			double mul=pow(10.0,n);
			return static_cast<int>((x*mul))/mul;
		}

		template <class RealType> std::string truncString(const RealType& x, int n){
			std::stringstream sstr;
			sstr.setf(std::ios::fixed,std::ios::floatfield);
			sstr.precision(2);
			sstr<<x;
			return sstr.str();
		}
	}
}




#endif //end FLOATUTILS_HPP
