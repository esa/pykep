// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef KEP3_DETAIL_S11N_HPP
#define KEP3_DETAIL_S11N_HPP

#include <cstddef>
#include <locale>
#include <random>
#include <sstream>
#include <string>
#include <tuple>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

// Include the archives.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace kep3::detail
{

// Implementation of tuple serialization.
template <std::size_t N>
struct tuple_s11n {
    template <class Archive, typename... Args>
    static void serialize(Archive &ar, std::tuple<Args...> &t, unsigned version)
    {
        ar &std::get<N - 1u>(t);
        tuple_s11n<N - 1u>::serialize(ar, t, version);
    }
};

template <>
struct tuple_s11n<0> {
    template <class Archive, typename... Args>
    static void serialize(Archive &, std::tuple<Args...> &, unsigned)
    {
    }
};

} // namespace kep3::detail


// NOTE: These macros implement the “KEY in header / IMPLEMENT in one .cpp” pattern for Boost.Serialization when used
// together with tanuki_kep3 type-erasure wrappers in a shared-library build. The KEY+extern-template macro declares the
// Boost export GUIDs for the wrapped concrete type and, crucially, adds `extern template` declarations for the
// corresponding `tanuki_kep3::v1::holder<ud_type, iface_type, Sem>` specializations so that consuming TUs (tests/users)
// do NOT implicitly instantiate those holders (and thus do not emit their own RTTI/vtables). The IMPLEMENT+instantiate
// macro must be used exactly once in the library: it provides the matching explicit instantiation definitions for the
// holders and emits the export/registration code so that the shared library “owns” a single, consistent type identity
// for the polymorphic erased objects across the dylib/executable boundary (avoiding runtime “unregistered cast/class”
// errors during serialization/deserialization). // [web:60][web:58][web:4]

#define KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(ud_type, iface_type)                         \
    /* Export key (GUID) for both semantics */                                                            \
    tanuki_kep3_S11N_WRAP_EXPORT_KEY(ud_type, iface_type)                                                       \
                                                                                                           \
    /* Prevent consumers from instantiating holder specializations */                                      \
    extern template struct ::tanuki_kep3::holder<ud_type, iface_type, ::tanuki_kep3::wrap_semantics::value>;  \
    extern template struct ::tanuki_kep3::holder<ud_type, iface_type, ::tanuki_kep3::wrap_semantics::reference>;

#define KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(ud_type, iface_type)                       \
    /* Provide the explicit instantiations matching the extern templates */                               \
    template struct ::tanuki_kep3::holder<ud_type, iface_type, ::tanuki_kep3::wrap_semantics::value>;        \
    template struct ::tanuki_kep3::holder<ud_type, iface_type, ::tanuki_kep3::wrap_semantics::reference>;    \
                                                                                                          \
    /* Export implementation (GUID registration) */                                                        \
    tanuki_kep3_S11N_WRAP_EXPORT_IMPLEMENT(ud_type, iface_type)

#define KEP3_S11N_EXPORT_WRAP(ud_type, iface_type)                                      \
    KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(ud_type, iface_type)                      \
    KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(ud_type, iface_type)

#endif
