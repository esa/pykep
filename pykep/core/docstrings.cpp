#include <string>

#include "docstrings.hpp"

namespace pykep
{

std::string expression_doc()
{
    return R"(
A precise point in time. Julian dates are supported.
)";
}

} // namespace pykep
