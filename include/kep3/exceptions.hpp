// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_DETAIL_EXCEPTIONS_HPP
#define kep3_DETAIL_EXCEPTIONS_HPP

#include <stdexcept>

#include <kep3/config.hpp>
#include <kep3/detail/visibility.hpp>

namespace kep3
{
struct kep3_DLL_PUBLIC not_implemented_error final : std::runtime_error {
    using std::runtime_error::runtime_error;
};

} // namespace kep3
#endif