#include <stdio.h>
#include <stdlib.h>

#define BOUNDARY_TYPE_FREESURFACE 0
#define BOUNDARY_TYPE_EXPONENT 1

/* ith-difference to 00 position  */
#define Dx00(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith-1)]   - a[iz*nx + ix-(ith  ) ] ) )

#define Dx10(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith-1)]   - a[iz*nx + ix-(ith  ) ] ) )

#define Dx11(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith  )]   - a[iz*nx + ix-(ith-1) ] ) )

#define Dx01(a,fdx,ith,ix,iz)   \
    ( fdx[ith-1] * ( a[iz*nx + ix+(ith  )]   - a[iz*nx + ix-(ith-1) ] ) )

#define Dz00(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+(ith-1))*nx + ix] - a[(iz- ith   )*nx + ix] ) )

#define Dz01(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+(ith-1))*nx + ix] - a[(iz- ith   )*nx + ix] ) )

#define Dz11(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+ ith   )*nx + ix] - a[(iz-(ith-1))*nx + ix] ) )

#define Dz10(a,fdx,ith,ix,iz)   \
    ( fdz[ith-1] * ( a[(iz+ ith   )*nx + ix] - a[(iz-(ith-1))*nx + ix] ) )
