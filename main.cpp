#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random.hpp>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
	array3D r0 = {{1.12, 2.32, 0.123}};
	array3D v0 = {{1/sqrt(2),1/sqrt(2),0.23}};
	//array3D v0 = {{0, 1, 0}};
	const array3D u = {{0.1,0.2,0.3}};

	double m0=1000;
	std::cout << r0 << v0 << m0 << std::endl;
	double t = M_PI_2;
	propagate_taylor_jorba(r0,v0,m0,u,t,2.1,1.45,-30,-30);
	std::cout << r0 << v0 << m0 << std::endl;
	propagate_taylor_jorba(r0,v0,m0,u,-t,2.1,1.45,-30,-30);
	std::cout << r0 << v0 << m0 <<  std::endl<<std:: endl<<std::endl;

	propagate_taylor(r0,v0,m0,u,t,2.1,1.45,-30,-30);
	std::cout << r0 << v0 << m0 << std::endl;
	propagate_taylor(r0,v0,m0,u,-t,2.1,1.45,-30,-30);
	std::cout << r0 << v0 << m0 <<  std::endl<<std:: endl<<std::endl;

	return 0;
}
