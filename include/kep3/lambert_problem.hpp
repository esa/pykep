// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_LAMBERT_PROBLEM_H
#define kep3_LAMBERT_PROBLEM_H

#include <array>
#include <vector>

#include <fmt/ostream.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>

namespace kep3
{

/// Lambert Problem
/**
 * This class represent a Lambert's problem. When instantiated it assumes a
 * prograde orbit (unless otherwise stated) and evaluates all the solutions up
 * to a maximum number of multiple revolutions. After the object is instantiated
 * the solutions can be retreived using the appropriate getters. Note that the
 * number of solutions will be N_max*2 + 1, where N_max is the maximum number of
 * revolutions.
 *
 * NOTE: The class has been tested extensively via monte carlo runs checked with
 * numerical propagation. Compared to the previous Lambert Solver in the
 * keplerian_toolbox it is 1.7 times faster (on average as defined by
 * lambert_test.cpp). With respect to Gooding algorithm it is 1.3 - 1.5 times
 * faster (zero revs - multi revs). The algorithm is described in detail in:
 *
 * Izzo, Dario. "Revisiting Lambert’s problem." Celestial Mechanics and
 * Dynamical Astronomy 121 (2015): 1-15.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class kep3_DLL_PUBLIC lambert_problem
{
    static const std::array<double, 3> default_r0;
    static const std::array<double, 3> default_r1;

public:
    // We choose in this case to friend  the streaming operator as to not expose a number of tedious getters (lambda, chord etc...)
    friend kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const lambert_problem &);
    explicit lambert_problem(const std::array<double, 3> &r0  = default_r0, const std::array<double, 3> &r1 = default_r1,
                             double tof = kep3::pi / 2, double mu = 1., bool cw = false, unsigned multi_revs = 1);
    [[nodiscard]] const std::vector<std::array<double, 3>> &get_v0() const;
    [[nodiscard]] const std::vector<std::array<double, 3>> &get_v1() const;
    [[nodiscard]] const std::array<double, 3> &get_r0() const;
    [[nodiscard]] const std::array<double, 3> &get_r1() const;
    [[nodiscard]] const double &get_tof() const;
    [[nodiscard]] const double &get_mu() const;
    [[nodiscard]] const std::vector<double> &get_x() const;
    [[nodiscard]] const std::vector<unsigned> &get_iters() const;
    [[nodiscard]] unsigned get_Nmax() const;
    [[nodiscard]] const bool &get_cw() const;

private:
    unsigned householder(double, double &, unsigned, double, unsigned) const;
    void dTdx(double &, double &, double &, double, double) const;
    void x2tof(double &tof, double x0, unsigned N) const;
    void x2tof2(double &tof, double x0, unsigned N) const;
    [[nodiscard]] static double hypergeometricF(double z, double tol);
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & m_r0;
        ar & m_r1;
        ar & m_tof;
        ar & m_mu;
        ar & m_v0;
        ar & m_v1;
        ar & m_iters;
        ar & m_x;
        ar & m_s;
        ar & m_c;
        ar & m_lambda;
        ar & m_iters;
        ar & m_Nmax;
        ar & m_has_converged;
        ar & m_multi_revs;
        ar & m_cw;
    }

    std::array<double, 3> m_r0, m_r1;
    double m_tof;
    double m_mu;
    std::vector<std::array<double, 3>> m_v0;
    std::vector<std::array<double, 3>> m_v1;
    std::vector<unsigned> m_iters;
    std::vector<double> m_x;
    double m_s, m_c, m_lambda;
    unsigned m_Nmax;
    bool m_has_converged;
    unsigned m_multi_revs;
    bool m_cw;
};

} // namespace kep3

template <>
struct fmt::formatter<kep3::lambert_problem> : fmt::ostream_formatter {
};

#endif // kep3_LAMBERT_PROBLEM_H