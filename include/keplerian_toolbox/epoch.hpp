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

#ifndef KEP_TOOLBOX_EPOCH_H
#define KEP_TOOLBOX_EPOCH_H

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

#include <keplerian_toolbox/detail/visibility.hpp>
#include <keplerian_toolbox/serialization.hpp>

/// Keplerian Toolbox
/**
 * This namespace contains astrodynamics and space flight mechanics routines that are related to
 * keplerian motions or models.
 */
namespace kep_toolbox
{

/// epoch class.
/**
 * This class defines and contains a non-gregorian date (i.e. a date expressed in julian form). It also provides the
 * user with an
 * interface to boost gregorian dates (see boost documentation at
 * http://www.boost.org/doc/libs/1_42_0/doc/html/date_time.html)
 * using the posix time.
 * To achieve higher performance the date is defined in MJD2000 (double) as a private member
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
class KEP_TOOLBOX_DLL_PUBLIC epoch
{
public:
    /** Types of non gregorian dates supported. Julian Date (JD) is the number of days passed since
     * January 1, 4713 BC at noon. Modified Julian Date (MJD) is the number of days passed since
     * November 17, 1858 at 00:00 am. The Modified Julian Date 2000 (MJD2000) is the number of days passed since
     * Juanuary 1, 2000 at 00:00am.
     */
    enum type { MJD2000, MJD, JD };

    /** @name Constructors */
    //@{
    epoch(const double &epoch_in = 0, type epoch_type = MJD2000);
    epoch(const boost::gregorian::greg_year &year, const boost::gregorian::greg_month &month,
          const boost::gregorian::greg_day &day);
    epoch(const boost::posix_time::ptime &posix_time);
    //@}

    /** @name Computing non-gregorian dates */
    //@{
    double mjd2000() const;
    double jd() const;
    double mjd() const;
    //@}

    /** @name Interface to boost::posix_time::ptime */
    //@{
    boost::posix_time::ptime get_posix_time() const;
    void set_posix_time(const boost::posix_time::ptime &);
    //@}

private:
    // Serialization code
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &mjd2000_m;
    }
    // Serialization code (END)

    /// the modified julian date 2000 stored in a double
    double mjd2000_m;
};

std::ostream &operator<<(std::ostream &s, const epoch &epoch_in);

KEP_TOOLBOX_DLL_PUBLIC epoch epoch_from_string(const std::string date);

KEP_TOOLBOX_DLL_PUBLIC epoch epoch_from_iso_string(const std::string date);

} // end of namespace kep_toolbox

#endif // KEP_TOOLBOX_EPOCH_H
