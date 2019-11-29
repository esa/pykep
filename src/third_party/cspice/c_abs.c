#include "f2c.h"

#ifdef KR_headers
extern double f__cabs();

double c_abs(z) complex_type *z;
#else
extern double f__cabs(double, double);

double c_abs(complex_type *z)
#endif
{
return( f__cabs( z->r, z->i ) );
}
