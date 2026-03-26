// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_EPOCH_HPP
#define kep3_EPOCH_HPP

#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <ratio>
#include <type_traits>

#include <fmt/ostream.h>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>

/// Keplerian Toolbox
/**
 * This namespace contains astrodynamics and space flight mechanics routines
 * that are related to keplerian motions or models.
 */

namespace kep3
{

namespace detail
{

template <typename T>
struct ensure_64bit_duration_impl {
};

template <typename Rep, typename Period>
    requires std::signed_integral<Rep>
struct ensure_64bit_duration_impl<std::chrono::duration<Rep, Period>> {
    using rep_t = std::conditional_t<(std::numeric_limits<Rep>::digits < std::numeric_limits<std::int64_t>::digits),
                                     std::int64_t, Rep>;
    using type = std::chrono::duration<rep_t, Period>;
};

} // namespace detail

// Small utility that takes the input duration T and checks its Rep type:
// if Rep is at least 64 bits long, then T is returned unchanged, otherwise
// T is widened to use std::int64_t as Rep type.
template <typename T>
using ensure_64bit = typename detail::ensure_64bit_duration_impl<T>::type;

// NOTE: we provide our own definition of microseconds
// in order to portably guarantee that we can represent a time range
// wider than the one guaranteed by the standard's definitions.
// Using at least a 64-bit integer, we are guaranteed to cover a time range
// of at least +-292k years around the reference epoch.
using microseconds = ensure_64bit<std::chrono::microseconds>;

// Similarly, we define our own time point counting ticks via (at least) std::int64_t
// and with microsesecond resolution. We use system_clock as the time point
// clock.
using time_point = std::chrono::time_point<std::chrono::system_clock, microseconds>;

/// epoch class.
/**
 * This class defines and contains a non-Gregorian date (i.e., a date expressed
 * in Julian format). The date is defined in MJD2000 format as a
 * kep_clock::time_point private member. Types of non-Gregorian dates supported:
 *      - Julian Date (JD): the number of days passed since January 1, 4713 BC
 * at 12:00 (noon).
 *      - Modified Julian Date (MJD): the number of days passed since November
 * 17, 1858 at 00:00 (midnight).
 *      - Modified Julian Date 2000 (MJD2000): the number of days passed since
 * Juanuary 1, 2000 at 00:00 (midnight).
 */
class kep3_DLL_PUBLIC epoch
{
    // Offset of 0 MJD2000 wrt Unix time.
    static constexpr auto y2k_offset = microseconds{946684800000000};

    // MJD2000 time point.
    static constexpr auto y2k = time_point{} + y2k_offset;

    // Offset of 0 JD wrt y2k.
    static constexpr auto jd_offset = microseconds{211813444800000000};

    // Offset of 0 MJD wrt y2k.
    static constexpr auto mjd_offset = microseconds{4453401600000000};

    // Used in several places below.
    using seconds_day_ratio = std::ratio<86400>;

    /* Helper functions for constructors */
    static time_point make_tp(int y, unsigned mon, unsigned d, std::int32_t h = 0, std::int32_t min = 0,
                              std::int32_t s = 0, std::int32_t ms = 0, std::int32_t us = 0);

public:
    enum class julian_type { MJD2000, MJD, JD };
    enum class string_format { ISO };

    /** Constructors */
    // Default constructor
    epoch();

    // Constructor from a julian date (as a floating-point value)
    explicit epoch(double epoch_in, julian_type = julian_type::MJD2000);

    // Constructor from string
    explicit epoch(const std::string &, string_format = string_format::ISO);

    // Constructor from time point.
    explicit epoch(time_point);

    /**
     * Constructs an epoch from a std::chrono::duration.
     * The reference point is assumed to be MJD2000.
     * \param[in] time The time as a duration.
     */
    template <typename Rep, typename Period>
    explicit epoch(std::chrono::duration<Rep, Period> duration)
        : m_tp{std::chrono::duration_cast<microseconds>(duration)}
    {
    }

    // Constructor for datetime broken down into its constituents.
    explicit epoch(int y, unsigned mon, unsigned d, std::int32_t h = 0, std::int32_t min = 0, std::int32_t s = 0,
                   std::int32_t ms = 0, std::int32_t us = 0);

    /* Computing non-Gregorian dates */

    /**
     * @return Number of days since 0 JD (including fractional days).
     */
    [[nodiscard]] constexpr double jd() const
    {
        return std::chrono::duration<double, seconds_day_ratio>(m_tp.time_since_epoch() - y2k_offset + jd_offset)
            .count();
    }

    /**
     * @return Number of days since 0 MJD (including fractional days).
     */
    [[nodiscard]] constexpr double mjd() const
    {
        return std::chrono::duration<double, seconds_day_ratio>(m_tp.time_since_epoch() - y2k_offset + mjd_offset)
            .count();
    }

    /**
     * @return Number of days since 0 MJD2000 (including fractional days).
     */
    [[nodiscard]] constexpr double mjd2000() const
    {
        return std::chrono::duration<double, seconds_day_ratio>(m_tp.time_since_epoch() - y2k_offset).count();
    }

    // Conversions
    static time_point tp_from_days(double days);

    // Duration conversions
    static constexpr double as_sec(const microseconds &d)
    {
        return std::chrono::duration<double, std::ratio<1>>(d).count();
    }

    // Printing
    [[nodiscard]] std::string as_utc_string() const;
    static std::string as_utc_string(const time_point &tp);

    /** operator overloads for sum and diff (epoch-days) and comparison
     * operators
     * **/

    kep3_DLL_PUBLIC friend std::ostream &operator<<(std::ostream &s, epoch const &epoch_in);

    template <typename Rep, typename Period>
    epoch &operator+=(const std::chrono::duration<Rep, Period> &d)
    {
        m_tp += std::chrono::duration_cast<microseconds>(d);
        return *this;
    }

    template <typename Rep, typename Period>
    epoch &operator-=(const std::chrono::duration<Rep, Period> &d)
    {
        m_tp -= std::chrono::duration_cast<microseconds>(d);
        return *this;
    }

    kep3_DLL_PUBLIC friend bool operator>(const epoch &, const epoch &);
    kep3_DLL_PUBLIC friend bool operator<(const epoch &, const epoch &);
    kep3_DLL_PUBLIC friend bool operator>=(const epoch &, const epoch &);
    kep3_DLL_PUBLIC friend bool operator<=(const epoch &, const epoch &);
    kep3_DLL_PUBLIC friend bool operator==(const epoch &, const epoch &);
    kep3_DLL_PUBLIC friend bool operator!=(const epoch &, const epoch &);

    template <typename Rep, typename Period>
    epoch operator+(const std::chrono::duration<Rep, Period> &d)
    {
        return epoch(m_tp + std::chrono::duration_cast<microseconds>(d));
    }

    template <typename Rep, typename Period>
    epoch operator-(const std::chrono::duration<Rep, Period> &d)
    {
        return epoch(m_tp - std::chrono::duration_cast<microseconds>(d));
    }

    static constexpr auto days(double value)
    {
        return std::chrono::duration_cast<microseconds>(std::chrono::duration<double, seconds_day_ratio>(value));
    }

    static constexpr auto sec(double value)
    {
        return std::chrono::duration_cast<microseconds>(std::chrono::duration<double, std::ratio<1>>(value));
    }

    kep3_DLL_PUBLIC friend microseconds operator-(const epoch &lhs, const epoch &rhs);

    [[nodiscard]] time_point get_tp() const;

    [[nodiscard]] static epoch now();

private:
    // Serialization code
    friend class boost::serialization::access;

    template <class Archive>
    void save(Archive &ar, const unsigned) const
    {
        const auto count{m_tp.time_since_epoch().count()};
        ar & count;
    }
    template <class Archive>
    void load(Archive &ar, const unsigned)
    {
        microseconds::rep count{0};
        ar & count;
        m_tp = time_point{microseconds(count)};
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int file_version)
    {
        boost::serialization::split_member(ar, *this, file_version);
    }

    time_point m_tp = y2k;
};

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &s, const epoch &epoch_in);

} // end of namespace kep3

template <>
struct fmt::formatter<kep3::epoch> : fmt::ostream_formatter {
};

#endif // kep3_EPOCH_HPP
