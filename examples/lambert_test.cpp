#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random.hpp>

#include "../src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    // Preamble
    array3D r1,r2;
    double tof;
    boost::mt19937 rng;
    boost::uniform_int<> dist(0, 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_bit(rng, dist);
    boost::uniform_real<> dist1(-2,2);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > drng(rng, dist1);
    double acc=0,err_max=0;
    int count=0;

    // Experiment Settings
    unsigned int Ntrials = 43800;

    // Start Experiment
    for (unsigned int i = 0; i<Ntrials; ++i){
        //1 - generate a random problem geometry
        r1[0] = drng() * 2; r1[1] = drng() * 2; r1[2] = drng() * 2;
        r2[0] = drng() * 2; r2[1] = drng() * 2; r2[2] = drng() * 2;
        tof = (drng () + 2) * 50 + 0.1;

        //2 - Solve the lambert problem
        lambert_problem lp(r1,r2,tof);

        //3 - Check its precision using propagate_lagrangian
        for (unsigned int i = 0; i < lp.get_v1().size(); ++i){
            array3D r1_p(r1),v1_p,err;
            v1_p = lp.get_v1()[i];
            propagate_lagrangian(r1_p,v1_p,tof,1.0);
            diff(err,r2,r1_p);
            err_max = std::max(err_max,norm(err));
            acc += norm(err);

        }
        count += (lp.get_Nmax()*2+1);
    }
    std::cout << "Max error: " << err_max <<std::endl;
    std::cout << "Average Error: " << acc / count <<std::endl;
    std::cout << "Number of Problems Solved: " << count << std::endl;
    return 0;
}
