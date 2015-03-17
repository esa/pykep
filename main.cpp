#include <iostream>
#include <iomanip>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    asteroid_gtoc2 pl1(33);
    planets::gtoc2 pl2(33);
    array3D r,v,r1,v1;
 	pl1.get_eph(kep_toolbox::epoch(1.23),r ,v);
    std::cout << r << v << std::endl;
 	pl2.eph(kep_toolbox::epoch(1.23), r1, v1);
    std::cout << r1 << v1 << std::endl;

    std::cout << pl1 << std::endl;
    std::cout << pl2 << std::endl;

    std::cout << pl1.get_elements(kep_toolbox::epoch(0)) << std::endl;
    std::cout << pl2.compute_elements() << std::endl;

    std::cout << pl1.compute_period() << std::endl;
    std::cout << pl2.compute_period() << std::endl;
    return 0;
}
