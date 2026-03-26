// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/container_hash/hash.hpp>

#include <fmt/core.h>

#include <heyoka/expression.hpp>
#include <heyoka/llvm_state.hpp>
#include <heyoka/model/vsop2013.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/vsop2013.hpp>

namespace kep3::udpla
{

namespace detail
{

namespace
{

// name -> planet index map.
// NOLINTNEXTLINE(cert-err58-cpp)
const std::unordered_map<std::string, std::uint32_t> pl_idx_map
    = {{"mercury", 1u}, {"venus", 2u},  {"earth_moon", 3u}, {"mars", 4u}, {"jupiter", 5u},
       {"saturn", 6u},  {"uranus", 7u}, {"neptune", 8u},    {"pluto", 9u}};

// Hasher for the JIT cache.
struct state_dict_hasher {
    std::size_t operator()(const std::pair<std::uint32_t, double> &p) const noexcept
    {
        std::size_t seed = std::hash<std::uint32_t>{}(p.first);
        boost::hash_combine(seed, std::hash<double>{}(p.second));

        return seed;
    }
};

// The JIT cache.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<std::pair<std::uint32_t, double>, heyoka::llvm_state, state_dict_hasher> state_dict;

// Mutex for thread-safe access to the cache.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex state_dict_mutex;

} // namespace

} // namespace detail

struct vsop2013::impl {
    using fptr_t = void (*)(double *, const double *, const double *, const double *) noexcept;

    heyoka::llvm_state m_state;
    fptr_t eval_f = nullptr;
    double m_mu = 0;
    std::string m_pl_name;
    double m_thresh = 0;

    impl() = default;
    explicit impl(heyoka::llvm_state s, double mu, std::string pl_name, double thresh)
        : m_state(std::move(s)), m_mu(mu), m_pl_name(std::move(pl_name)), m_thresh(thresh)
    {
        eval_f = reinterpret_cast<fptr_t>(m_state.jit_lookup("eval_f"));
    }
    impl(impl &&) noexcept = delete;
    impl(const impl &other)
        : m_state(other.m_state), m_mu(other.m_mu), m_pl_name(other.m_pl_name), m_thresh(other.m_thresh)
    {
        eval_f = reinterpret_cast<fptr_t>(m_state.jit_lookup("eval_f"));
    }
    impl &operator=(impl &&) noexcept = delete;
    impl &operator=(const impl &) = delete;
    ~impl() = default;

    void save(boost::archive::binary_oarchive &ar, unsigned) const
    {
        ar << m_state;
        ar << m_mu;
        ar << m_pl_name;
        ar << m_thresh;
    }
    void load(boost::archive::binary_iarchive &ar, unsigned)
    {
        ar >> m_state;
        ar >> m_mu;
        ar >> m_pl_name;
        ar >> m_thresh;

        eval_f = reinterpret_cast<fptr_t>(m_state.jit_lookup("eval_f"));
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

vsop2013::vsop2013() : vsop2013("mercury") {}

vsop2013::vsop2013(std::string pl_name, double thresh)
{
    if (!std::isfinite(thresh) || thresh < 0) {
        throw std::invalid_argument(
            fmt::format("The threshold value must be finite and non-negative, but it is {} instead", thresh));
    }

    boost::algorithm::to_lower(pl_name);
    const auto pl_it = detail::pl_idx_map.find(pl_name);
    if (pl_it == detail::pl_idx_map.end()) {
        throw std::invalid_argument(fmt::format("The planet name '{}' is invalid", pl_name));
    }

    // Fetch the planet index.
    const auto pl_index = pl_it->second;

    // Fetch the planet mu.
    auto mu = heyoka::model::get_vsop2013_mus()[pl_index];
    // The original mu is in AU**3/day**2, convert it to SI units.
    mu = mu * kep3::AU * kep3::AU * kep3::AU / (kep3::DAY2SEC * kep3::DAY2SEC);

    // Lock down for access to the cache.
    std::lock_guard lock(detail::state_dict_mutex);

    if (auto it = detail::state_dict.find({pl_index, thresh}); it == detail::state_dict.end()) {
        // Cache miss, we need to create a new compiled function.

        // Time variable.
        auto tm = heyoka::expression("tm");

        // Create the expression for the VSOP2013 solution.
        // NOTE: we use (tm - 0.5) / 365250 as time coordinate because:
        // - VSOP2013 counts time from J2000 and not MJD2000 (they differ by 12 hours),
        // - VSOP2013 measures time in thousands of Julian years.
        auto v_ex = heyoka::model::vsop2013_cartesian_icrf(pl_index, heyoka::kw::time_expr = (tm - 0.5) / 365250.,
                                                           heyoka::kw::thresh = thresh);

        // Create the llvm_state and add the compiled function.
        heyoka::llvm_state s;
        heyoka::add_cfunc<double>(s, "eval_f", v_ex, {tm});
        s.compile();

        // Add the compiled state to the cache.
        [[maybe_unused]] auto [_, flag] = detail::state_dict.insert({{pl_index, thresh}, s});
        assert(flag);

        // Build the impl.
        m_impl = std::make_unique<impl>(std::move(s), mu, std::move(pl_name), thresh);
    } else {
        // Cache hit, build the impl from the cache content.
        m_impl = std::make_unique<impl>(it->second, mu, std::move(pl_name), thresh);
    }
}

vsop2013::vsop2013(vsop2013 &&) noexcept = default;

vsop2013::vsop2013(const vsop2013 &other) : m_impl(std::make_unique<impl>(*other.m_impl)) {}

vsop2013 &vsop2013::operator=(vsop2013 &&) noexcept = default;

vsop2013 &vsop2013::operator=(const vsop2013 &other)
{
    if (this != &other) {
        *this = vsop2013(other);
    }

    return *this;
}

void vsop2013::save(boost::archive::binary_oarchive &ar, unsigned) const
{
    ar << m_impl;
}

void vsop2013::load(boost::archive::binary_iarchive &ar, unsigned)
{
    try {
        ar >> m_impl;
        // LCOV_EXCL_START
    } catch (...) {
        *this = vsop2013{};
        throw;
    }
    // LCOV_EXCL_STOP
}

vsop2013::~vsop2013() = default;

std::array<std::array<double, 3>, 2> vsop2013::eph(double mjd2000) const
{
    double out[6]{};
    m_impl->eval_f(out, &mjd2000, nullptr, nullptr);

    return {
        {{{out[0] * kep3::AU, out[1] * kep3::AU, out[2] * kep3::AU}},
         {{out[3] * kep3::AU / kep3::DAY2SEC, out[4] * kep3::AU / kep3::DAY2SEC, out[5] * kep3::AU / kep3::DAY2SEC}}}};
}

std::string vsop2013::get_name() const
{
    return fmt::format("vsop2013 {}, threshold={}", m_impl->m_pl_name, m_impl->m_thresh);
}

} // namespace kep3::udpla

// NOLINTNEXTLINE
KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(kep3::udpla::vsop2013, kep3::detail::planet_iface)
