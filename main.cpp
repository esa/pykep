#include <iostream>
#include <iomanip>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    planet_ss pl1("mars");
    planets::jpl_low_precision pl2("earth");
    array3D r,v,r1,v1;
 	pl1.get_eph(kep_toolbox::epoch(1.23),r ,v);
    std::cout << r << v << std::endl;
 	pl2.eph(kep_toolbox::epoch(1.23), r1, v1);
    std::cout << r1 << v1 << std::endl;

    return 0;
}
