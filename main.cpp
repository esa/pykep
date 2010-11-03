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
    array3D r0 = {{-134510944015.16229, 1544909621.3945236, 1673333377.1862454}};
    array3D v0 = {{-277730.23684193398, -451656.73287408578, 170099.07655321722}};
    double t = -7.01253e+06;
    propagate_lagrangian(r0,v0,t,ASTRO_MU_SUN);
    propagate_lagrangian(r0,v0,-t,ASTRO_MU_SUN);
    std::cout << r0 << v0 << std::endl;
    return 0;
}
