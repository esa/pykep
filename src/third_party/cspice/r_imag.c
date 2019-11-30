#include "f2c.h"

#ifdef KR_headers
double r_imag(z) complex_type *z;
#else
double r_imag(complex_type *z)
#endif
{
return(z->i);
}
