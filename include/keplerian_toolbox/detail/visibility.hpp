#ifndef KEP_TOOLBOX_VISIBILITY_HPP
#define KEP_TOOLBOX_VISIBILITY_HPP

// Convenience macros for setting the visibility of entities
// when building/using the shared library. Mostly inspired by:
// https://gcc.gnu.org/wiki/Visibility
// We check first for Windows, where we assume every compiler
// knows dllexport/dllimport. On other platforms, we use the GCC-like
// syntax for GCC, clang and ICC. Otherwise, we leave the definitions
// empty.
#if defined(_WIN32) || defined(__CYGWIN__)

#if defined(keplerian_toolbox_EXPORTS)

#define KEP_TOOLBOX_DLL_PUBLIC __declspec(dllexport)

#else

#define KEP_TOOLBOX_DLL_PUBLIC __declspec(dllimport)

#endif

#define KEP_TOOLBOX_DLL_LOCAL

#elif defined(__clang__) || defined(__GNUC__) || defined(__INTEL_COMPILER)

#define KEP_TOOLBOX_DLL_PUBLIC __attribute__((visibility("default")))
#define KEP_TOOLBOX_DLL_LOCAL __attribute__((visibility("hidden")))

#else

#define KEP_TOOLBOX_DLL_PUBLIC
#define KEP_TOOLBOX_DLL_LOCAL

#endif

// NOTE: it seems like on Windows using dllimport/dllexport on inline classes
// is generally not helpful (and potentially harmful), apart from special use cases:
// https://stackoverflow.com/questions/8876279/c-inline-functions-with-dllimport-dllexport
// https://stackoverflow.com/questions/24511376/how-to-dllexport-a-class-derived-from-stdruntime-error
// https://devblogs.microsoft.com/oldnewthing/20140109-00/?p=2123
// Setting the visibility attribute on GCC-like compilers for inline classes, however, seems to be ok.
// Thus, we use a specialised definition for marking "public"ly visible inline classes.
#if defined(_WIN32) || defined(__CYGWIN__)

#define KEP_TOOLBOX_DLL_PUBLIC_INLINE_CLASS

#else

#define KEP_TOOLBOX_DLL_PUBLIC_INLINE_CLASS KEP_TOOLBOX_DLL_PUBLIC

#endif

#endif
