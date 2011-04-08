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

//	std::string line("99942   19.2   0.15 K107N 202.49545  126.41859  204.43202    3.33173  0.1911104  1.11267324   0.9223398  1 MPO164109  1397   2 2004-2008 0.40 M-v 3Eh MPCAPO     C802  (99942) Apophis            20080109");
	planet_mpcorb();
	return 0;
}
