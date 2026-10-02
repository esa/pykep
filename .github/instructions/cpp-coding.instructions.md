---
description: "Use when writing or modifying generic C++ code in this repository. Keep the style minimal and consistent until more domain-specific rules are added."
applyTo: "**/*.cpp,**/*.hpp,**/*.h,**/*.cc,**/*.cxx"
---

- Follow the repository's `.clang-format`: C++20, 4-space indentation, spaces only, 120-column limit, right-aligned pointers, Linux braces, and no one-line control statements. Run the formatter rather than hand-formatting large changes.
- Use the existing file layout: project copyright/SPDX header, standard-library includes, third-party includes, then `kep3` includes, with blank lines between groups. Keep includes alphabetized within each group where practical.
- Keep declarations in the corresponding public header and definitions in the matching `.cpp`. Use `m_` for private data, `get_`/`set_` for accessors, and `is_` for boolean predicates. Preserve the library's snake_case naming.
- Prefer the concrete standard-library types used by this project, especially `std::array` for fixed-size vectors and states, `std::pair`/`std::optional` for multi-value and optional results, and `std::vector` for dynamic sequences.
- Match the local `const` style: use `const` before types for parameters and references, and use `type const` for many local computed values when that matches the surrounding function. Prefer `auto` for structured bindings and obvious iterator/container types, not for obscuring important API types.
- Pass read-only objects as `const &` when appropriate. Use raw pointers only for explicit non-owning interfaces; otherwise prefer values, references, or standard smart pointers with clear ownership.
- Keep validation near the start of public operations. Throw `std::domain_error`, `std::invalid_argument`, or `std::logic_error` consistently with the surrounding code; use assertions only for internal invariants that cannot be caused by user input.
- Put implementation-only helpers in `kep3::detail`, a more specific nested namespace, or an anonymous namespace. Do not expose internal helpers through public headers.
- Use Doxygen-style comments for public classes and functions, including `@brief`, parameter descriptions, and return descriptions where useful. Keep inline comments for intent, invariants, numerical-method choices, and non-obvious workarounds; do not narrate obvious statements.
- Preserve the project's numerical-code conventions: use explicit floating-point literals such as `0.` and `1.`, use `std::` math functions, and retain `// NOLINT` or `// NOLINTNEXTLINE(check-name)` only for documented intentional exceptions to clang-tidy.
- Use `decltype(container.size())` or another matching unsigned/index type for loop counters when the bound comes from a container or library API. Avoid signed/unsigned conversions and unchecked narrowing.
- Keep public APIs and serialization-facing state predictable. Avoid unrelated refactors, new abstractions, and formatting churn when making a focused change.
