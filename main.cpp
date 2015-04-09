#include <iostream>
#include <iomanip>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
	planet::gtoc5 pl2(33);
	array3D r1, v1;

	pl2.eph(kep_toolbox::epoch(1.23), r1, v1);
	std::cout << r1 << v1 << std::endl;
	
	
	std::cout << pl2 << std::endl;
	
	
	std::cout << pl2.compute_elements() << std::endl;
	
	
	std::cout << pl2.compute_period() << std::endl;
	return 0;
}
