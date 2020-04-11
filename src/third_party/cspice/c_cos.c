#include "f2c.h"

#ifdef KR_headers
extern double sin(), cos(), sinh(), cosh();

VOID c_cos(r, z) complex_type *r, *z;
#else
#undef abs
#include "math.h"

void c_cos(complex_type *r, complex_type *z)
#endif
{
	double zr = z->r;
	r->r =   cos(zr) * cosh(z->i);
	r->i = - sin(zr) * sinh(z->i);
	}
