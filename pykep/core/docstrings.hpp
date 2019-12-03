#ifndef PYKEP_DOCSTRINGS_HPP
#define PYKEP_DOCSTRINGS_HPP

#include <string>

namespace pykep
{
// epoch
std::string epoch_doc();
std::string epoch_from_string_doc();
std::string epoch_from_iso_string_doc();

// Lambert problem
std::string lambert_problem_doc();

// Propagations
std::string propagate_lagrangian_doc();
std::string propagate_taylor_doc();


} // namespace pykep

#endif
