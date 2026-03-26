// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <kep3/core_astro/encodings.hpp>

namespace kep3
{

std::vector<double> alpha2direct(const std::vector<double> &alphas, double tof)
{
    std::vector<double> log_alphas(alphas.size());
    std::transform(alphas.begin(), alphas.end(), log_alphas.begin(), [](double x) { return std::log(x); });
    double sum = std::accumulate(log_alphas.begin(), log_alphas.end(), 0.);
    std::transform(log_alphas.begin(), log_alphas.end(), log_alphas.begin(),
                   [sum, tof](double x) { return x * tof / sum; });
    return log_alphas;
}

std::pair<std::vector<double>, double> direct2alpha(const std::vector<double> &tofs)
{
    std::vector<double> retval(tofs.size());
    double T = std::accumulate(tofs.begin(), tofs.end(), 0.);
    std::transform(tofs.begin(), tofs.end(), retval.begin(), [T](double x) { return std::exp(-x / T); });
    return {retval, T};
}

std::vector<double> eta2direct(const std::vector<double> &etas, double max_tof)
{
    if (etas.empty()) {
        throw std::invalid_argument("etas must be non-empty");
    }

    std::vector<double> tofs(etas.size(), 0.);
    tofs[0] = max_tof * etas[0];
    double cumulative_T = tofs[0];
    for (decltype(etas.size()) i = 1u; i < etas.size(); ++i) {
        tofs[i] = (max_tof - cumulative_T) * etas[i];
        cumulative_T += tofs[i];
    }
    return tofs;
}

std::vector<double> direct2eta(const std::vector<double> &tofs, double max_tof)
{
    if (tofs.empty()) {
        throw std::invalid_argument("tofs must be non-empty");
    }
    auto retval = tofs;
    retval[0] = tofs[0] / max_tof;
    double cumulative_tofs = tofs[0];
    for (decltype(tofs.size()) i = 1u; i < tofs.size(); ++i) {
        retval[i] = tofs[i] / (max_tof - cumulative_tofs);
        cumulative_tofs += tofs[i];
    }
    return retval;
}

} // namespace kep3