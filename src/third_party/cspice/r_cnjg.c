#include "f2c.h"

#ifdef KR_headers
VOID r_cnjg(r, z) complex_type *r, *z;
#else
VOID r_cnjg(complex_type *r, complex_type *z)
#endif
{
r->r = z->r;
r->i = - z->i;
}
