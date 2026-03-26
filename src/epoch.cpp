// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <fmt/chrono.h>
#include <fmt/core.h>

#include <kep3/epoch.hpp>

namespace kep3
{

/**
 * @brief Constructs a default epoch .
 */
epoch::epoch() = default;

/**
 * @brief Constructs an epoch from a non-Gregorian date.
 *
 * @param[in] julian_epoch_in A double indicating the number of days
                              since day 0 in the specified calendar.
 * @param[in] epoch_type epoch::julian_type
 */
epoch::epoch(double julian_epoch_in, julian_type epoch_type)
{
    switch (epoch_type) {
        case julian_type::MJD2000:
            m_tp = epoch::tp_from_days(julian_epoch_in);
            break;
        case julian_type::MJD:
            m_tp = epoch::tp_from_days(julian_epoch_in) - mjd_offset;
            break;
        case julian_type::JD:
            m_tp = epoch::tp_from_days(julian_epoch_in) - jd_offset;
            break;
        default:
            throw std::invalid_argument(
                fmt::format("An unsupported julian_type enumerator with value {} was used in the epoch constructor",
                            static_cast<std::underlying_type_t<julian_type>>(epoch_type)));
    }
}

/**
 * @brief Constructs an epoch from offsets relative to 0 MJD2000.
 *
 * @param[in] y The number of years.
 * @param[in] d The number of days.
 * @param[in] h The number of hours.
 * @param[in] min The number of minutes.
 * @param[in] s The number of seconds.
 * @param[in] ms The number of milliseconds.
 * @param[in] us The number of microseconds.
 */
epoch::epoch(int y, unsigned mon, unsigned d, const std::int32_t h, const std::int32_t min, const std::int32_t s,
             const std::int32_t ms, const std::int32_t us)
    : m_tp{make_tp(y, mon, d, h, min, s, ms, us)}
{
}

// Epoch constructor from string
epoch::epoch(const std::string &in, string_format sf)
{
    // NOTE: only ISO format supported so far.
    if (sf != string_format::ISO) {
        throw std::invalid_argument(
            fmt::format("An unsupported string_format enumerator with value {} was used in the epoch constructor",
                        static_cast<std::underlying_type_t<string_format>>(sf)));
    }

    // We assume: 1980-10-17T11:36:21.121841 and allow crops such as 1980-10.
    constexpr std::array<decltype(in.size()), 11> allowed_lenghts{7, 10, 13, 16, 19, 21, 22, 23, 24, 25, 26};
    const auto len = in.size();
    auto foo = std::find(std::begin(allowed_lenghts), std::end(allowed_lenghts), len); // NOLINT
    if (foo == std::end(allowed_lenghts)) {
        throw std::logic_error(
            "Malformed input string when constructing an epoch. Must be 'YYYY-MM-DD HH:MM:SS:XXXXXX'. "
            "D,H,M,S and X can be missing incrementally.");
    }
    unsigned d = 1;
    std::int32_t h = 0, min = 0, s = 0, us = 0;
    int const y = std::stoi(in.substr(0, 4));
    auto mon = boost::numeric_cast<unsigned>(std::stoi(in.substr(5, 2)));
    if (len >= 10) {
        d = boost::numeric_cast<unsigned>(std::stoi(in.substr(8, 2)));
        if (len >= 13) {
            h = boost::numeric_cast<std::int32_t>(std::stoi(in.substr(11, 2)));
            if (len >= 16) {
                min = boost::numeric_cast<std::int32_t>(std::stoi(in.substr(14, 2)));
                if (len >= 19) {
                    s = boost::numeric_cast<std::int32_t>(std::stoi(in.substr(17, 2)));
                    if (len >= 21) {
                        std::string rest = in.substr(20);
                        us = boost::numeric_cast<std::int32_t>(std::stoi(rest));
                        for (decltype(rest.size()) i = 0; i < 6u - rest.size(); ++i) {
                            us *= 10;
                        }
                    }
                }
            }
        }
    }
    m_tp = make_tp(y, mon, d, h, min, s, 0, us);
}

/**
 * @brief Constructs an epoch from a time point.
 *
 * @param[in] time_point Self-explanatory.
 */
epoch::epoch(time_point tp) : m_tp{tp} {}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
time_point epoch::make_tp(int y, unsigned mon, unsigned d, const std::int32_t h, const std::int32_t min,
                          const std::int32_t s, const std::int32_t ms, const std::int32_t us)
{
    // Construct and validate year, month and day.
    const std::chrono::year chr_y(y);
    if (!chr_y.ok()) {
        throw std::invalid_argument(
            fmt::format("An invalid year value of {} was specified in the construction of a time point", y));
    }

    const std::chrono::month chr_mon(mon);
    if (!chr_mon.ok()) {
        throw std::invalid_argument(
            fmt::format("An invalid month value of {} was specified in the construction of a time point", mon));
    }

    const std::chrono::day chr_d(d);
    if (!chr_d.ok()) {
        throw std::invalid_argument(
            fmt::format("An invalid day value of {} was specified in the construction of a time point", d));
    }

    // Build the year_month_day object.
    auto ymd = std::chrono::year_month_day(chr_y, chr_mon, chr_d);

    // Normalise year and month (see https://en.cppreference.com/w/cpp/chrono/year_month_day/operator_days).
    ymd += std::chrono::months{0};

    // Convert ymd to a std::chrono::system_clock::time_point. This also normalises the day.
    const auto chr_tp = static_cast<std::chrono::sys_days>(ymd);

    // Convert chr_tp to our time point.
    // NOTE: this is guaranteed to work because the conversion factor is 86400 * 10**6,
    // which is always representable by std::int_max_t. See the docs for the constructor
    // of std::time_point and the constructor for std::duration:
    // https://en.cppreference.com/w/cpp/chrono/time_point/time_point
    // https://en.cppreference.com/w/cpp/chrono/duration/duration
    time_point tp = chr_tp;

    // Add the rest.
    tp += ensure_64bit<std::chrono::hours>(h);
    tp += ensure_64bit<std::chrono::minutes>(min);
    tp += ensure_64bit<std::chrono::seconds>(s);
    tp += ensure_64bit<std::chrono::milliseconds>(ms);
    tp += ensure_64bit<std::chrono::microseconds>(us);

    return tp;
}

bool operator>(const epoch &c1, const epoch &c2)
{
    return c1.m_tp > c2.m_tp;
}

bool operator<(const epoch &c1, const epoch &c2)
{
    return c1.m_tp < c2.m_tp;
}

bool operator>=(const epoch &c1, const epoch &c2)
{
    return c1.m_tp >= c2.m_tp;
}

bool operator<=(const epoch &c1, const epoch &c2)
{
    return c1.m_tp <= c2.m_tp;
}

bool operator==(const epoch &c1, const epoch &c2)
{
    return c1.m_tp == c2.m_tp;
}

bool operator!=(const epoch &c1, const epoch &c2)
{
    return c1.m_tp != c2.m_tp;
}

/**
 * @brief Creates time point from the number of days since 0 MJD2000.
 *
 * @return A time point
 */
time_point epoch::tp_from_days(double days)
{
    return y2k + std::chrono::duration_cast<microseconds>(std::chrono::duration<double, seconds_day_ratio>(days));
}

/**
 * @brief Returns a time point formatted as a date/time string
 * in the in the format 2000-12-31T12:34:56.123456.
 *
 * @param tp The time point.
 * @return A formatted date/time string.
 */
std::string epoch::as_utc_string(const time_point &tp)
{
    // Format it with fmt.
    // NOTE: fmt is seemingly able to format all system_clock time points,
    // even ours which will typically have a different duration.
    return fmt::format("{:%FT%T}", tp);
}

std::string epoch::as_utc_string() const
{
    return epoch::as_utc_string(m_tp);
}

/**
 * @brief Streams out an epoch as a UTC string.
 *
 * @param[in] s Stream to which the epoch will be sent.
 * @param[in] ep Epoch to be sent to the stream.
 *
 * @return Reference to s.
 */
std::ostream &operator<<(std::ostream &s, const epoch &ep)
{
    s << ep.as_utc_string();
    return s;
}

microseconds operator-(const epoch &lhs, const epoch &rhs)
{
    return lhs.m_tp - rhs.m_tp;
}

time_point epoch::get_tp() const
{
    return m_tp;
}

epoch epoch::now()
{
    // NOTE: need explicit conversion here as system_clock might have a different
    // duration than microseconds.
    auto tp = std::chrono::time_point_cast<microseconds>(std::chrono::system_clock::now());

    return epoch(tp);
}

} // namespace kep3
