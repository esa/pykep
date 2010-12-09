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
    array3D r0 = {{1, 0, 0}};
    array3D v0 = {{0, 1, 0}};
    array3D u = {{1,0,0}};
    double m0=1000;
	std::cout << r0 << v0 << m0 << std::endl;
    double t = M_PI_2;
    propagate_taylor(r0,v0,m0,u,t,1,1,-12,-12);
	std::cout << r0 << v0 << m0 << std::endl;
    propagate_taylor(r0,v0,m0,u,-t,1.1,1,-12,-12);
    std::cout << r0 << v0 << m0 <<  std::endl;
    return 0;
}
