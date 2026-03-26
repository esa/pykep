// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_CONVERT_JULIAN_DATES_HPP
#define kep3_CONVERT_JULIAN_DATES_HPP

namespace kep3
{
inline double jd2mjd(double in)
{
    return (in - 2400000.5);
}
inline double jd2mjd2000(double in)
{
    return (in - 2451544.5);
}
inline double mjd2jd(double in)
{
    return (in + 2400000.5);
}
inline double mjd2mjd2000(double in)
{
    return (in - 51544);
}
inline double mjd20002jd(double in)
{
    return (in + 2451544.5);
}
inline double mjd20002mjd(double in)
{
    return (in + 51544);
}

} // namespace kep3

#endif // kep3_CONVERT_JULIAN_DATES_HPP
