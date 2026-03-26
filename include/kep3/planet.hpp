// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_PLANET_H
#define kep3_PLANET_H

#include <cmath>
#include <exception>
#include <string>
#include <typeinfo>

#include <boost/core/demangle.hpp>
#include <boost/safe_numerics/safe_integer.hpp>

#include <fmt/core.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>

#include <kep3/detail/tanuki.hpp>

namespace kep3
{

namespace detail
{

// Concepts to detect whether user classes have certain methods implemented.
#define KEP3_UDPLA_CONCEPT_HAS_GET(name, type)                                                                         \
    template <typename T>                                                                                              \
    concept udpla_has_get_##name = requires(const T &p) {                                                              \
        { p.get_##name() } -> std::same_as<type>;                                                                      \
    }

KEP3_UDPLA_CONCEPT_HAS_GET(mu_central_body, double);
KEP3_UDPLA_CONCEPT_HAS_GET(mu_self, double);
KEP3_UDPLA_CONCEPT_HAS_GET(radius, double);
KEP3_UDPLA_CONCEPT_HAS_GET(safe_radius, double);
KEP3_UDPLA_CONCEPT_HAS_GET(name, std::string);
KEP3_UDPLA_CONCEPT_HAS_GET(extra_info, std::string);

#undef KEP3_UDPLA_CONCEPT_HAS_GET

template <typename T>
concept udpla_has_eph = requires(const T &p, double mjd2000) {
    { p.eph(mjd2000) } -> std::same_as<std::array<std::array<double, 3>, 2>>;
};

template <typename T>
concept udpla_has_eph_v = requires(const T &p, const std::vector<double> &mjd2000s) {
    { p.eph_v(mjd2000s) } -> std::same_as<std::vector<double>>;
};

template <typename T>
concept udpla_has_acc = requires(const T &p, double mjd2000) {
    { p.acc(mjd2000) } -> std::same_as<std::array<double, 3>>;
};

template <typename T>
concept udpla_has_acc_v = requires(const T &p, const std::vector<double> &mjd2000s) {
    { p.acc_v(mjd2000s) } -> std::same_as<std::vector<double>>;
};

template <typename T>
concept udpla_has_period = requires(const T &p, double mjd2000) {
    { p.period(mjd2000) } -> std::same_as<double>;
};

template <typename T>
concept udpla_has_elements = requires(const T &p, double mjd2000, kep3::elements_type el_t) {
    { p.elements(mjd2000, el_t) } -> std::same_as<std::array<double, 6>>;
};

// Concept detecting if the type T can be used as a udpla.
template <typename T>
concept any_udpla = std::default_initializable<T> && std::copy_constructible<T> && std::move_constructible<T>
                    && std::destructible<T> && udpla_has_eph<T>;

// Fwd declaration of the interface implementation.
template <typename Base, typename Holder, typename T>
    requires any_udpla<T>
struct planet_iface_impl;

struct kep3_DLL_PUBLIC planet_iface {
    virtual ~planet_iface();

    [[nodiscard]] virtual double get_mu_central_body() const = 0;
    [[nodiscard]] virtual double get_mu_self() const = 0;
    [[nodiscard]] virtual double get_radius() const = 0;
    [[nodiscard]] virtual double get_safe_radius() const = 0;
    [[nodiscard]] virtual std::string get_name() const = 0;
    [[nodiscard]] virtual std::string get_extra_info() const = 0;

    [[nodiscard]] virtual std::array<std::array<double, 3>, 2> eph(double) const = 0;
    [[nodiscard]] virtual std::vector<double> eph_v(const std::vector<double> &) const = 0;
    // If implemented returns the acceleration vector at mjd2000
    [[nodiscard]] virtual std::array<double, 3> acc(double) const = 0;
    [[nodiscard]] virtual std::vector<double> acc_v(const std::vector<double> &) const = 0;


    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] virtual double period(double = 0.) const = 0;
    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] virtual std::array<double, 6> elements(double = 0.,
                                                         kep3::elements_type = kep3::elements_type::KEP_F) const
        = 0;

    // Methods that are non virtual, not implementable by the user in udplas, but visible in kep3::planet
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &) const;
    [[nodiscard]] std::array<double, 3> acc(const epoch &) const;

    [[nodiscard]] std::array<double, 6> elements(const kep3::epoch &,
                                                 kep3::elements_type = kep3::elements_type::KEP_F) const;

    [[nodiscard]] double period(const epoch &) const;

    template <typename Base, typename Holder, typename T>
    using impl = planet_iface_impl<Base, Holder, T>;
};

// Helper macro to implement getters in the planet interface implementation.
#define KEP3_UDPLA_IMPLEMENT_GET(name, type, def_value)                                                                \
    [[nodiscard]] type get_##name() const final                                                                        \
    {                                                                                                                  \
        if constexpr (udpla_has_get_##name<T>) {                                                                       \
            return getval<Holder>(this).get_##name();                                                                  \
        } else {                                                                                                       \
            return def_value;                                                                                          \
        }                                                                                                              \
    }

// NOTE: implement this in the cpp in order to avoid
// instantiating the same code over and over.
kep3_DLL_PUBLIC double period_from_energy(const std::array<double, 3> &, const std::array<double, 3> &, double);
kep3_DLL_PUBLIC std::array<double, 6> elements_from_posvel(const std::array<std::array<double, 3>, 2> &, double,
                                                           kep3::elements_type);
template <typename T>
std::vector<double> default_eph_vectorization(const T *self, const std::vector<double> &mjd2000s)
{
    // We simply call a for loop.
    const auto size = mjd2000s.size();
    using size_type = std::vector<double>::size_type;
    std::vector<double> retval;
    retval.resize(boost::safe_numerics::safe<size_type>(size) * 6);
    for (decltype(mjd2000s.size()) i = 0u; i < size; ++i) {
        auto values = self->eph(mjd2000s[i]);
        retval[6 * i] = values[0][0];
        retval[6 * i + 1] = values[0][1];
        retval[6 * i + 2] = values[0][2];
        retval[6 * i + 3] = values[1][0];
        retval[6 * i + 4] = values[1][1];
        retval[6 * i + 5] = values[1][2];
    }
    return retval;
}

template <typename T>
std::vector<double> default_acc_vectorization(const T *self, const std::vector<double> &mjd2000s)
{
    // We simply call a for loop.
    const auto size = mjd2000s.size();
    std::vector<double> retval(size * 3);
    for (decltype(mjd2000s.size()) i = 0u; i < size; ++i) {
        auto value = self->acc(mjd2000s[i]);
        retval[3*i] = value[0];
        retval[3*i + 1] = value[1];
        retval[3*i + 2] = value[2];

    }
    return retval;
}

// Planet interface implementation.
template <typename Base, typename Holder, typename T>
    requires any_udpla<T>
struct planet_iface_impl : public Base {
    KEP3_UDPLA_IMPLEMENT_GET(mu_central_body, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(mu_self, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(safe_radius, double, -1)
    KEP3_UDPLA_IMPLEMENT_GET(name, std::string, boost::core::demangle(typeid(T).name()))
    KEP3_UDPLA_IMPLEMENT_GET(extra_info, std::string, "")

    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double mjd2000) const final
    {
        return getval<Holder>(this).eph(mjd2000);
    }

    [[nodiscard]] std::vector<double> eph_v(const std::vector<double> &mjd2000s) const final
    {
        if constexpr (udpla_has_eph_v<T>) {
            return getval<Holder>(this).eph_v(mjd2000s);
        } else {
            return default_eph_vectorization(this, mjd2000s);
        }
    }

    [[nodiscard]] std::array<double, 3> acc(double mjd2000) const final
    {
        if constexpr (udpla_has_acc<T>) {
            // User has correctly implemented the method acc, return the user implementation
            return getval<Holder>(this).acc(mjd2000);
        } else {
            // User did not implement the method acc but the udpla has a central body, return the Keplerian acceleration
            // NOTE: this violates the actual ephemerides as it may not be dv/dt, but in most cases will be a good
            // enough approximation, and is actually correct for keplerian motion.
            if constexpr (udpla_has_get_mu_central_body<T>) {
                auto mu = getval<Holder>(this).get_mu_central_body();
                const auto &[pos, vel] = getval<Holder>(this).eph(mjd2000);
                double r3 = std::pow(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2], 1.5);
                return {-mu / r3 * pos[0], -mu / r3 * pos[1], -mu / r3 * pos[2]};
                //  User did not implement the method acc and the udpla does not have a central body ... nothing we can
                //  do
            } else {
                throw not_implemented_error(fmt::format("The acc method has not been implemented correctly (if at all) "
                                                        "in the User Defined planet (UDPLA)"));
            }
        }
    }

    [[nodiscard]] std::vector<double> acc_v(const std::vector<double> &mjd2000s) const final
    {
        if constexpr (udpla_has_acc_v<T>) {
            return getval<Holder>(this).acc_v(mjd2000s);
        } else {
            return default_acc_vectorization(this, mjd2000s);
        }
    }

    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] double period(double mjd2000 = 0.) const final
    {
        // If the user provides an efficient way to compute the period, then use it
        if constexpr (udpla_has_period<T>) {
            return getval<Holder>(this).period(mjd2000);
        } else if constexpr (udpla_has_get_mu_central_body<T>) {
            // If the user provides the central body parameter, then compute the
            // period from the energy at epoch
            auto [r, v] = eph(mjd2000);
            double mu = get_mu_central_body();
            return period_from_energy(r, v, mu);
        } else {
            // There is no way to compute a period for this planet
            throw not_implemented_error(fmt::format("A period nor a central body mu has been declared for '{}', "
                                                    "impossible to provide a default implementation",
                                                    get_name()));
        }
    }

    // NOLINTNEXTLINE(google-default-arguments)
    [[nodiscard]] std::array<double, 6> elements(double mjd2000 = 0.,
                                                 kep3::elements_type el_type = kep3::elements_type::KEP_F) const final
    {
        // If the user provides an efficient way to compute the orbital elements, then use it.
        if constexpr (udpla_has_elements<T>) {
            return getval<Holder>(this).elements(mjd2000, el_type);
        } else if constexpr (udpla_has_get_mu_central_body<T>) {
            // If the user provides the central body parameter, then compute the
            // elements using posvel computed at ep and converted.
            auto pos_vel = eph(mjd2000);
            double mu = get_mu_central_body();
            return elements_from_posvel(pos_vel, mu, el_type);
        } else {
            // There is no way to compute osculating elements for this planet
            throw not_implemented_error(
                fmt::format("A central body mu has not been declared for '{}', "
                            "impossible to provide a default implementation to compute the osculating elements",
                            get_name()));
        }
    }
};

#undef KEP3_UDPLA_IMPLEMENT_GET

// The udpla used in the default constructor of planet.
struct kep3_DLL_PUBLIC null_udpla {
    null_udpla() = default;
    static std::array<std::array<double, 3>, 2> eph(double);

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned) {};
};

// Reference interface for the planet class.
struct kep3_DLL_PUBLIC planet_ref_iface {
    template <typename Wrap>
    struct impl {
        tanuki_kep3_REF_IFACE_MEMFUN(get_mu_central_body) 
        tanuki_kep3_REF_IFACE_MEMFUN(get_mu_self)
        tanuki_kep3_REF_IFACE_MEMFUN(get_radius) 
        tanuki_kep3_REF_IFACE_MEMFUN(get_safe_radius)
        tanuki_kep3_REF_IFACE_MEMFUN(get_name) 
        tanuki_kep3_REF_IFACE_MEMFUN(get_extra_info)
        tanuki_kep3_REF_IFACE_MEMFUN(eph) 
        tanuki_kep3_REF_IFACE_MEMFUN(eph_v)
        tanuki_kep3_REF_IFACE_MEMFUN(acc) 
        tanuki_kep3_REF_IFACE_MEMFUN(acc_v) 
        tanuki_kep3_REF_IFACE_MEMFUN(period)
        tanuki_kep3_REF_IFACE_MEMFUN(elements)

            // Implement the extract functionality.
            template <typename T>
            T *extract()
        {
            return value_ptr<T>(*static_cast<Wrap *>(this));
        }
    };
};

} // namespace detail

#if defined(__GNUC__)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"

#endif

// Definition of the planet class.
class kep3_DLL_PUBLIC planet
{
    using wrap_t
        = tanuki_kep3::wrap<detail::planet_iface, tanuki_kep3::config<detail::null_udpla, detail::planet_ref_iface>{
                                                      .pointer_interface = false}>;

    wrap_t m_wrap;

private:
    friend class boost::serialization::access;

    template <class Archive>
    void save(Archive &ar, const unsigned int /*version*/) const
    {
        ar & m_wrap;
    }

    template <class Archive>
    void load(Archive &ar, const unsigned int /*version*/)
    {
        ar & m_wrap;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    planet() = default;

    template <typename T>
        requires(!std::same_as<planet, std::remove_cvref_t<T>>) && std::constructible_from<wrap_t, T &&>
    explicit planet(T &&p) : m_wrap(std::forward<T>(p))
    {
    }

    planet(const planet &) noexcept = default;
    planet(planet &&) noexcept = default;
    planet &operator=(const planet &) noexcept = default;
    planet &operator=(planet &&) noexcept = default;
    ~planet() = default;

    [[nodiscard]] double get_mu_central_body() const
    {
        return m_wrap.get_mu_central_body();
    }
    [[nodiscard]] double get_mu_self() const
    {
        return m_wrap.get_mu_self();
    }
    [[nodiscard]] double get_radius() const
    {
        return m_wrap.get_radius();
    }
    [[nodiscard]] double get_safe_radius() const
    {
        return m_wrap.get_safe_radius();
    }
    [[nodiscard]] std::string get_name() const
    {
        return m_wrap.get_name();
    }
    [[nodiscard]] std::string get_extra_info() const
    {
        return m_wrap.get_extra_info();
    }
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double mjd2000) const
    {
        return m_wrap.eph(mjd2000);
    }
    [[nodiscard]] std::vector<double> eph_v(const std::vector<double> &mjd2000s) const
    {
        return m_wrap.eph_v(mjd2000s);
    }
    [[nodiscard]] std::array<double, 3> acc(double mjd2000) const
    {
        return m_wrap.acc(mjd2000);
    }
    [[nodiscard]] std::vector<double> acc_v(const std::vector<double> &mjd2000s) const
    {
        return m_wrap.acc_v(mjd2000s);
    }
    [[nodiscard]] double period(double mjd2000 = 0.) const
    {
        return m_wrap.period(mjd2000);
    }
    [[nodiscard]] std::array<double, 6> elements(double mjd2000 = 0.,
                                                 kep3::elements_type el_t = kep3::elements_type::KEP_F) const
    {
        return m_wrap.elements(mjd2000, el_t);
    }

    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(const epoch &ep) const
    {
        return m_wrap.eph(ep);
    }

    [[nodiscard]] std::array<double, 3> acc(const epoch &ep) const
    {
        return m_wrap.acc(ep);
    }

    [[nodiscard]] std::array<double, 6> elements(const kep3::epoch &ep,
                                                 kep3::elements_type el_t = kep3::elements_type::KEP_F) const
    {
        return m_wrap.elements(ep, el_t);
    }
    [[nodiscard]] double period(const epoch &ep) const
    {
        return m_wrap.period(ep);
    }

    template <typename T>
    [[nodiscard]] T *extract()
    {
        return m_wrap.template extract<T>();
    }

    template <typename T>
    [[nodiscard]] const T *extract() const
    {
        return value_ptr<T>(m_wrap);
    }

    [[nodiscard]] std::type_index get_type_index() const
    {
        return value_type_index(m_wrap);
    }
};

kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const kep3::planet &);
} // namespace kep3

// version 2: planet has an acc optional method
BOOST_CLASS_VERSION(kep3::planet, 2)

template <>
struct fmt::formatter<kep3::planet> : fmt::ostream_formatter {
};

// Keep exporting the null udpla under the iface (same as before).
KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(kep3::detail::null_udpla, kep3::detail::planet_iface)

#endif // kep3_PLANET_H
