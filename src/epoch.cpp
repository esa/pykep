/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <boost/date_time/gregorian/gregorian.hpp>
#include <iomanip>
#include <cmath>
#include <iostream>

#include <keplerian_toolbox/astro_constants.hpp>
#include <keplerian_toolbox/core_functions/convert_dates.hpp>
#include <keplerian_toolbox/epoch.hpp>

namespace kep_toolbox
{
using namespace boost::gregorian;
using namespace boost::posix_time;

/// Constructor.
/**
* Constructs an epoch from a non-gregorian date.
* \param[in] epoch_in A double indicating the non-gregorian date
* \param[in] epoch_type One of [epoch::MJD2000, epoch::MJD, epoch::JD]
*/
epoch::epoch(const double &epoch_in, type epoch_type) : mjd2000_m(epoch_in)
{
    switch (epoch_type) {
        case MJD2000:
            break;
        case MJD:
            mjd2000_m = mjd2mjd2000(epoch_in);
            break;
        case JD:
            mjd2000_m = jd2mjd2000(epoch_in);
            break;
    }
}

/// Constructor.
/**
* Constructs an epoch from a gregorian date. The time of the day is assumed to be midnight.
* \param[in] year The gregorian year
* \param[in] month The month of the year
* \param[in] day The day of the month
*/
epoch::epoch(const greg_year &year, const greg_month &month, const greg_day &day)
{
    set_posix_time(ptime(date(year, month, day)));
}

/// Constructor.
/**
* Constructs an epoch from a boost ptime object (posix time)
* \param[in] posix_time The posix_time
*/
epoch::epoch(const boost::posix_time::ptime &posix_time)
{
    time_duration dt = posix_time - ptime(date(2000, 1, 1));
    bool flag = false;
    if (dt.is_negative()) {
        flag = true;
        dt = dt.invert_sign();
    }
    double fr_secs = static_cast<double>(dt.fractional_seconds()) * BOOST_DATE_PRECISION;
    mjd2000_m = static_cast<double>(dt.hours()) / 24.0 + static_cast<double>(dt.minutes()) / 1440.0 + (static_cast<double>(dt.seconds()) + fr_secs) / 86400.0;
    if (flag) mjd2000_m = -mjd2000_m;
}

/// jd getter.
/**
* Returns the julian date
*
* @return double containing the julian date
*
*/
double epoch::jd() const
{
    return mjd20002jd(mjd2000_m);
}

/// mjd getter.
/**
* Returns the modified julian date
*
* @return double containing the modified julian date
*
*/
double epoch::mjd() const
{
    return mjd20002mjd(mjd2000_m);
}

/// mjd2000 getter.
/**
* Gets the modified julian date 2000
* @return const reference to mjd2000
*/
double epoch::mjd2000() const
{
    return mjd2000_m;
}

/// Extracts the posix time
/**
* Returns the posix_time representation of the epoch. The method evaluates
* from the mjd2000 the number of days, months, seconds and
* micro/nano seconds passed since the 1st of January 2000 and uses this information
* to build the posix time
*
* @return ptime containing the posix time
*
*/
ptime epoch::get_posix_time() const
{
    long hrs, min, sec, fsec;
    bool flag = false;
    double copy = mjd2000_m;
    if (copy < 0) {
        copy = -copy;
        flag = true;
    }
    hrs = static_cast<long>(copy * 24);
    min = static_cast<long>((copy * 24 - static_cast<double>(hrs)) * 60);
    sec = static_cast<long>((((copy * 24 - static_cast<double>(hrs)) * 60) - static_cast<double>(min)) * 60);
    double dblfsec = ((((copy * 24 - static_cast<double>(hrs)) * 60) - static_cast<double>(min)) * 60) - static_cast<double>(sec);
    std::ostringstream fsecstr;
    fsecstr << std::setiosflags(std::ios::fixed) << std::setprecision(-std::log10(BOOST_DATE_PRECISION)) << dblfsec;
    fsec = boost::lexical_cast<long>(fsecstr.str().substr(2, -std::log10(BOOST_DATE_PRECISION) + 1));
    ptime retval;
    if (flag)
        retval = ptime(date(2000, 1, 1), time_duration(-hrs, -min, -sec, -fsec));
    else
        retval = ptime(date(2000, 1, 1), time_duration(hrs, min, sec, fsec));
    return retval;
}

/// Sets the epoch from a posix time
/**
* Sets the epoch to a particular posix_time.
*
* \param[in] posix_time containing the posix time
*
*/
void epoch::set_posix_time(const boost::posix_time::ptime &posix_time)
{

    mjd2000_m = epoch(posix_time).mjd2000();
}

/// Returns an epoch constructed from a delimited string containing a date
/**
 *  Builds an epoch from a delimited string. Excess digits in fractional seconds will be dropped. Ex:
 * "1:02:03.123456999" => "1:02:03.123456".
 *  This behavior depends on the precision defined in astro_constant.h used to compile
 *
 * Example:
 * 	std::string ts("2002-01-20 23:59:54.003");
 * 	epoch e(epoch_from_string(ts))
 *
 */
epoch epoch_from_string(const std::string date)
{
    return epoch(boost::posix_time::ptime(boost::posix_time::time_from_string(date)));
}

/// Returns an epoch constructed from a non delimited iso string containing a date
/**
 *  Builds an epoch from a non delimited iso string containing a date.
 *
 * Example:
 * 	std::string ts("20020131T235959");
 * 	epoch e(epoch_from_iso_string(ts))
 *
 */
epoch epoch_from_iso_string(const std::string date)
{
    return epoch(boost::posix_time::ptime(boost::posix_time::from_iso_string(date)));
}

} // end of namespace kep_toolbox

/// Overload the stream operator for kep_toolbox::epoch
/**
 * Streams out a date in the format 2000-Jan-01 00:12:30.123457
 *
 * \param[in] s stream to which the epoch will be sent
 * \param[in] now epoch to be sent to stream
 *
 * \return reference to s
 *
 */
std::ostream &kep_toolbox::operator<<(std::ostream &s, const kep_toolbox::epoch &now)
{
    s << now.get_posix_time();
    return s;
}
