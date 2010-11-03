#ifndef CODING_H
#define CODING_H

#include <numeric>
#include <boost/iterator/counting_iterator.hpp>
#include "../astro_constants.h"
#include "../epoch.h"
#include "../serialization.h"


namespace kep_toolbox {

/**
 * indicates that the variables are stored in a cartesian format
 */
class cartesian_coding_tag {};

/**
* @brief An instance of this class represents a format which
 * can be used to read the chromosome. In order to read some
 * data from the chromosome, instead of having to somehow know
 * in what position the data is and read the appropriate
 * variable, you can just ask this class for the data,
 * whithout having to know its location or the format used to
 * express it. For example, to read the start epoch of leg
 * number 2, you could use, given the vector x
 *
 * @verbatim
coding.leg_start_epoch(2, x.begin());
 * @endverbatim
 *
 * This allows you to build a trajectory without having to
 * know which position, the various elements have, which
 * format and which physical units are used to express
 * them. More importantly, it allows to change all of that
 * without having to modify the code which handles the
 * trajectories.
 *
 * The base_format class represents the following chromosome format
 *
 * starting epoch expressed in MJD2000
 * then for each leg:
 * start velocity expressed in cartesian coodinates in km/s relative to the planet
 * all the impulses expressed in cartesian coordinates as a fraction of the maximum allowed impulse
 * end velocity expressed in cartesian coodinates in km/s relative to the planet
 * leg end epoch expressed in MJD2000
 *
 * Mass is not in the chromosome, and is considered fixed for the whole trajectory.
 */
class base_format : public cartesian_coding_tag {

public:
	/// Constructor.
	/**
	 * @param n_legs The number of legs in the trajectory
	 * @param n_thrust The number of thrusts in each leg (which is assumed to be the same for all the legs)
	 * @param mass The mass of the spacecraft in kg, which is assumed to be constant for the whole trajectory.
	 */
	template<typename it_type>
	base_format(it_type begin, it_type end, double mass)
	    : n_thrusts_m(begin, end), mass_m(mass) {}

	base_format(int n_legs, int n_segments, double mass)
	    : n_thrusts_m(n_legs, n_segments), mass_m(mass) {}

	std::vector<int> leg_start_velocity_i(int leg_n) const {
	    int index = leg_start(leg_n);
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> leg_end_velocity_i(int leg_n) const {
		int index = leg_start(leg_n) + 3 * (n_thrusts_m[leg_n] + 1);
		return std::vector<int>(boost::counting_iterator<int>(index),
					boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> segment_thrust_i(int leg_n, int thrust_n) const {
		assert(thrust_n < n_thrusts_m[leg_n]);
		int index = leg_start(leg_n) + 3 * (thrust_n + 1);
		return std::vector<int>(boost::counting_iterator<int>(index),
					boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> leg_start_epoch_i(int leg_n) const {
		return std::vector<int>(1, leg_start(leg_n) - 1);
	}

	std::vector<int> leg_end_epoch_i(int leg_n) const {
		return leg_start_epoch_i(leg_n + 1);
	}

	std::vector<int> segment_start_epoch_i(int leg_n, int segment_n) const {
		(void) segment_n;
		std::vector<int> indexes(2);
		indexes[0] = leg_start(leg_n) - 1;
		indexes[1] = leg_start(leg_n + 1) - 1;
		return indexes;
	}

	std::vector<int> segment_end_epoch_i(int leg_n, int segment_n) const {
		return segment_start_epoch_i(leg_n, segment_n);
	}

	std::vector<int> leg_start_mass_i(int leg_n) const {
		(void)leg_n;
		return std::vector<int>(0);
	}

	std::vector<int> leg_end_mass_i(int leg_n) const {
		(void)leg_n;
		return std::vector<int>(0);
	}

	std::vector<int> leg_i(int leg_n) const {
	    int index = leg_start(leg_n) - 1;
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + leg_size(leg_n) + 1));

	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the leg start velocity of leg leg_n in cartesian coortinates in m/s relative to the starting object
	 */
	template<typename it_type>
	array3D leg_start_velocity(int leg_n, it_type begin) const {
		array3D x;
		std::vector<int> indexes(leg_start_velocity_i(leg_n));
		for(size_t i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]] * 1000.;
		return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the leg end velocity of leg leg_n in cartesian coortinates in m/s relative to the starting object
	 */
	template<typename it_type>
	array3D leg_end_velocity(int leg_n, it_type begin) const {
	    array3D x;
	    std::vector<int> indexes(leg_end_velocity_i(leg_n));
	    for(size_t i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]] * 1000.;

	    return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the starting epoch of leg leg_n
	 */
	template<typename it_type>
	epoch leg_start_epoch(int leg_n, it_type begin) const {
	    return epoch(begin[leg_start(leg_n) - 1], epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the end epoch of leg leg_n
	 */
	template<typename it_type>
	epoch leg_end_epoch(int leg_n, it_type begin) const {
	    return leg_start_epoch(leg_n + 1, begin);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the flight time of leg leg_n in seconds
	 */
	template<typename it_type>
	double leg_flight_time(int leg_n, it_type begin) const {
	    return (leg_end_epoch(leg_n, begin).mjd2000() -
		    leg_start_epoch(leg_n, begin).mjd2000()) * 60. * 60. * 24.;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the thrust magnitude with respect to the maximum thrust allowed
	 * during that segment in cartesian coortinates
	 */
	template<typename it_type>
	array3D segment_thrust(int leg_n, int segment_n, it_type begin) const {
	    assert(segment_n < n_thrusts_m[leg_n]);
	    array3D x;
	    std::vector<int> indexes(segment_thrust_i(leg_n, segment_n));
	    for(size_t i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]];

	    return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the start epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_start_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() +
			 (leg_end_epoch(leg_n, begin).mjd2000() - leg_start_epoch(leg_n, begin).mjd2000())
			 / double(n_segments(leg_n)) * double(segment_n), epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the end epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_end_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() +
			 (leg_end_epoch(leg_n, begin).mjd2000() - leg_start_epoch(leg_n, begin).mjd2000())
			 / double(n_segments(leg_n)) * double(segment_n + 1), epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the spacecraft mass in kg at the beginning of leg leg_n
	 */
	template<typename it_type>
	double leg_start_mass(int, it_type) const {
	    return mass_m;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the spacecraft mass in kg at the end of leg leg_n
	 */
	template<typename it_type>
	double leg_end_mass(int leg_n, it_type begin) const {
		(void)leg_n;
		(void) begin;
		return mass_m;
	}

	/**
	 * @param leg_n the index of the leg
	 * @returns the number of segments in leg leg_n
	 */
	int n_segments(int leg_n) const {
	    return n_thrusts_m[leg_n];
	}

	const std::vector<int>& n_segments() const {
	    return n_thrusts_m;
	}

	/**
	 * @returns the number of legs in the trajectory
	 */
	int n_legs() const {
	    return n_thrusts_m.size();
	}

	/**
	 * @returns the number of variables of the chromosome
	 * which represents the trajectory
	 */
	int size() const {
	    return leg_start(n_thrusts_m.size());
	}



    private:
	/**
	 * @param leg_n the index of the leg
	 * @returns the number of variables used to represent the leg
	 */
	int leg_size(int leg_n) const {
	    return (n_thrusts_m[leg_n] + 2) * 3 + 1;
	}

	int leg_start(int leg_n) const {
	    return std::accumulate(n_thrusts_m.begin(), n_thrusts_m.begin() + leg_n, 0) * 3 + 7 * leg_n + 1;
	}
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & n_thrusts_m;
		ar & mass_m;
	}

	std::vector<int> n_thrusts_m;
	double mass_m;
    };

    /**
     * The base_format class represents the following chromosome format
     *
     * starting epoch expressed in MJD2000
     * starting mass in kg
     * then for each leg:
     * start velocity expressed in cartesian coodinates in km/s relative to the planet
     * all the impulses expressed in cartesian coordinates as a fraction of the maximum allowed impulse
     * end velocity expressed in cartesian coodinates in km/s relative to the planet
     * leg end epoch expressed in MJD2000
     * leg end mass expressed in kg
     */
    class mass_format : public cartesian_coding_tag {
    public:
	template<typename it_type>
	mass_format(it_type begin, it_type end)
	    : n_thrusts_m(begin, end) {}


	std::vector<int> leg_start_velocity_i(int leg_n) const {
	    int index = leg_start(leg_n);
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> leg_end_velocity_i(int leg_n) const {
	    int index = leg_start(leg_n) + 3 * (n_thrusts_m[leg_n] + 1);
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> segment_thrust_i(int leg_n, int thrust_n) const {
	    assert(thrust_n < n_thrusts_m[leg_n]);
	    int index = leg_start(leg_n) + 3 * (thrust_n + 1);
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + 3));
	}

	std::vector<int> leg_start_epoch_i(int leg_n) const {
	    return std::vector<int>(1, leg_start(leg_n) - 2);
	}

	std::vector<int> leg_end_epoch_i(int leg_n) const {
		return leg_start_epoch_i(leg_n + 1);
	}

	std::vector<int> segment_start_epoch_i(int leg_n, int segment_n) const {
		(void)segment_n;
		std::vector<int> indexes(2);
		indexes[0] = leg_start(leg_n) - 2;
		indexes[1] = leg_start(leg_n + 1) - 2;
		return indexes;
	}

	std::vector<int> segment_end_epoch_i(int leg_n, int segment_n) const {
		return segment_start_epoch_i(leg_n, segment_n);
	}

	std::vector<int> leg_start_mass_i(int leg_n) const {
	    return std::vector<int>(1, leg_start(leg_n) - 1);
	}

	std::vector<int> leg_end_mass_i(int leg_n) const {
	    return leg_start_mass_i(leg_n + 1);
	}

	std::vector<int> leg_i(int leg_n) const {
	    int index = leg_start(leg_n) - 2;
	    return std::vector<int>(boost::counting_iterator<int>(index),
				    boost::counting_iterator<int>(index + leg_size(leg_n) + 2));

	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the leg start velocity of leg leg_n in cartesian coortinates in m/s relative to the starting object
	 */
	template<typename it_type>
	array3D leg_start_velocity(int leg_n, it_type begin) const {
	    array3D x;
	    std::vector<int> indexes(leg_start_velocity_i(leg_n));
	    for(int i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]] * 1000.;

	    return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the leg end velocity of leg leg_n in cartesian coortinates in m/s relative to the starting object
	 */
	template<typename it_type>
	array3D leg_end_velocity(int leg_n, it_type begin) const {
	    array3D x;
	    std::vector<int> indexes(leg_end_velocity_i(leg_n));
	    for(int i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]] * 1000.;

	    return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the starting epoch of leg leg_n
	 */
	template<typename it_type>
	epoch leg_start_epoch(int leg_n, it_type begin) const {
	    return epoch(begin[leg_start(leg_n) - 2], epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the end epoch of leg leg_n
	 */
	template<typename it_type>
	epoch leg_end_epoch(int leg_n, it_type begin) const {
	    return leg_start_epoch(leg_n + 1, begin);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the flight time of leg leg_n in seconds
	 */
	template<typename it_type>
	double leg_flight_time(int leg_n, it_type begin) const {
	    return (leg_end_epoch(leg_n, begin).mjd2000() -
		    leg_start_epoch(leg_n, begin).mjd2000()) * 60. * 60. * 24.;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the thrust magnitude with respect to the maximum thrust allowed
	 * during that segment in cartesian coortinates
	 */
	template<typename it_type>
	array3D segment_thrust(int leg_n, int segment_n, it_type begin) const {
	    assert(segment_n < n_thrusts_m[leg_n]);
	    array3D x;
	    std::vector<int> indexes(segment_thrust_i(leg_n, segment_n));
	    for(int i = 0; i < indexes.size(); ++i)
		x[i] = begin[indexes[i]];

	    return x;
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the start epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_start_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() +
			 (leg_end_epoch(leg_n, begin).mjd2000() - leg_start_epoch(leg_n, begin).mjd2000())
			 / double(n_segments(leg_n)) * double(segment_n), epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the end epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_end_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() +
			 (leg_end_epoch(leg_n, begin).mjd2000() - leg_start_epoch(leg_n, begin).mjd2000())
			 / double(n_segments(leg_n)) * double(segment_n + 1), epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the spacecraft mass in kg at the beginning of leg leg_n
	 */
	template<typename it_type>
	double leg_start_mass(int leg_n, it_type begin) const {
	    return begin[leg_start(leg_n) - 1];
	}

	/**
	 * @param leg_n the index of the leg
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the spacecraft mass in kg at the end of leg leg_n
	 */
	template<typename it_type>
	double leg_end_mass(int leg_n, it_type begin) const {
	    return leg_start_mass((leg_n + 1), begin);
	}

	/**
	 * @param leg_n the index of the leg
	 * @returns the number of segments in leg leg_n
	 */
	int n_segments(int leg_n) const {
	    return n_thrusts_m[leg_n];
	}

	const std::vector<int>& n_segments() const {
	    return n_thrusts_m;
	}

	/**
	 * @returns the number of legs in the trajectory
	 */
	int n_legs() const {
		return n_thrusts_m.size();
	}

	/**
	 * @returns the number of variables of the chromosome
	 * which represents the trajectory
	 */
	int size() const {
		return leg_start(n_thrusts_m.size());
	}


	protected:
	/**
	 * @param leg_n the index of the leg
	 * @returns the number of variables used to represent the leg
	 */
	int leg_size(int leg_n) const {
		return (n_thrusts_m[leg_n] + 2) * 3 + 2;
	}

	int leg_start(int leg_n) const {
		return std::accumulate(n_thrusts_m.begin(), n_thrusts_m.begin() + leg_n, 0) * 3 + 8 * leg_n + 2;
	}

	std::vector<int> n_thrusts_m;
	};






	/**
	* Same as the mass format, but forces the engines to be off for a
	* fixed amount of time before the end of each leg
	*/
	template<typename base_format_type>
	class coast_arcs : public base_format_type {
	public:

	/**
	 * @param days the number of days where the engines must be
	 * switched off
	 *
	 * @note arc times are in days
	 */
	template<typename dep_it_type, typename arr_it_type>
	coast_arcs(const base_format_type& base,
		   dep_it_type departure_arc_b, dep_it_type departure_arc_e,
		   arr_it_type arrival_arc_b, arr_it_type arrival_arc_e)
	    : base_format_type(base),
	      departure_arc_times_m(departure_arc_b, departure_arc_e),
	      arrival_arc_times_m(arrival_arc_b, arrival_arc_e)
	{
	    assert(base_format_type::n_legs() == departure_arc_times_m.size());
	    assert(base_format_type::n_legs() == arrival_arc_times_m.size());
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the start epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_start_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() + departure_arc_times_m[leg_n] +
			 (std::max(0., leg_end_epoch(leg_n, begin).mjd2000() - departure_arc_times_m[leg_n]
				   - arrival_arc_times_m[leg_n] - leg_start_epoch(leg_n, begin).mjd2000()))
			 / double(base_format_type::n_segments(leg_n)) * double(segment_n), epoch::MJD2000);
	}

	/**
	 * @param leg_n the index of the leg
	 * @param segment_n the index of the segment
	 * @param begin an iterator pointing at the beginning of the chromosome
	 * @returns the end epoch of segment segment_n of leg leg_n
	 */
	template<typename it_type>
	epoch segment_end_epoch(int leg_n, int segment_n, it_type begin) const {
	    return epoch(leg_start_epoch(leg_n, begin).mjd2000() + departure_arc_times_m[leg_n] +
			 (std::max(0., leg_end_epoch(leg_n, begin).mjd2000() - departure_arc_times_m[leg_n]
				   - arrival_arc_times_m[leg_n] - leg_start_epoch(leg_n, begin).mjd2000()))
			 / double(base_format_type::n_segments(leg_n)) * double(segment_n + 1), epoch::MJD2000);
	}

    private:
	std::vector<double> departure_arc_times_m;
	std::vector<double> arrival_arc_times_m;
    };
}

#endif
