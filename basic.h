/*  basic.h
 *
 *  Copyright 2011  Simon Fristed Eskildsen, Vladimir Fonov,
 *   	      	    Pierrick Coupé, Jose V. Manjon
 *
 *  This file is part of mincbeast.
 *
 *  mincbeast is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mincbeast is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with mincbeast.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For questions and feedback, please contact:
 *  Simon Fristed Eskildsen <eskild@gmail.com> 
 */


#ifndef BASIC_H
#define BASIC_H

#ifdef HAVE_MINC
#include <volume_io.h>
#endif //HAVE_MINC

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef STATUS_OK
#define STATUS_OK 0
#endif

#ifndef STATUS_ERR
#define STATUS_ERR -1
#endif

#ifndef STATUS_GIVE_UP
#define STATUS_GIVE_UP 5
#endif

#ifndef PI
#define PI 3.141592654
#endif

#ifndef EPSILON
#define EPSILON 0.0001
#endif

#define byte unsigned char

#ifndef ABS
#define ABS(a) ((a) > 0.0 ? (a) : -(a))
#endif

#ifndef SQR
#define SQR(a) ((a)*(a))
#endif

#ifdef   MAX
#undef   MAX
#endif
#define  MAX( x, y )  ( ((x) >= (y)) ? (x) : (y) )
#ifdef   MIN
#undef   MIN
#endif
#define  MIN( x, y )  ( ((x) <= (y)) ? (x) : (y) )

#ifdef   MAX3
#undef   MAX3
#endif
#define MAX3(a,b,c) MAX(MAX(a,b),c)
#ifdef   MIN3
#undef   MIN3
#endif
#define MIN3(a,b,c) MIN(MIN(a,b),c)
#ifdef   CLAMP
#undef   CLAMP
#endif
#define CLAMP(a,b,c) ( (a) < (b) ? (b) : ( (a) > (c) ? (c) : (a) ) )

#ifdef ROUND
#undef ROUND
#endif
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))

#ifndef HAVE_MINC
typedef int VIO_BOOL;
typedef double VIO_Real;
#endif

#ifndef _NETCDF_
typedef int nc_type;
#endif

#ifndef HAVE_MINC
#define NC_NAT          0       /* NAT = 'Not A Type' (c.f. NaN) */
#define NC_BYTE         1       /* signed 1 byte integer */
#define NC_CHAR         2       /* ISO/ASCII character */
#define NC_SHORT        3       /* signed 2 byte integer */
#define NC_INT          4       /* signed 4 byte integer */
#define NC_LONG         NC_INT  /* deprecated, but required for backward compatibility. */
#define NC_FLOAT        5       /* single precision floating point number */
#define NC_DOUBLE       6       /* double precision floating point number */
#define NC_UBYTE        7       /* unsigned 1 byte int */
#define NC_USHORT       8       /* unsigned 2-byte int */
#define NC_UINT         9       /* unsigned 4-byte int */
#define NC_INT64        10      /* signed 8-byte int */
#define NC_UINT64       11      /* unsigned 8-byte int */
#define NC_STRING       12      /* string */
#endif

#define NC_DOUBLE_SIZE sizeof(double)

typedef struct {
  float x;
  float y;
  float z;
} point3D;

typedef struct {
  int x;
  int y;
  int z;
} point3D_int;

typedef struct {
  float *start;
  float *step;
  int   *length;
  char  *history;
} image_metadata;




#define FALSE 0
#define TRUE 1
#define STATUS_OK 0
#define STATUS_ERR -1
#define STATUS_GIVE_UP 5
#define byte unsigned char

#define DTOR 0.017453293

#define SQR(a) ((a)*(a))
#define SET_3DPOINT(p,a,b,c) (p).x=(a); (p).y=(b); (p).z=(c)
#define ADD_3DPOINT(p,p1,p2) (p).x=(p1).x+(p2).x; (p).y=(p1).y+(p2).y; (p).z=(p1).z+(p2).z
#define SUB_3DPOINT(p,p1,p2) (p).x=(p1).x-(p2).x; (p).y=(p1).y-(p2).y; (p).z=(p1).z-(p2).z 
#define DIV_3DPOINT(p,p1,a) (p).x=(p1).x/(a); (p).y=(p1).y/(a); (p).z=(p1).z/(a)
#define MUL_3DPOINT(p,p1,a) (p).x=(p1).x*(a); (p).y=(p1).y*(a); (p).z=(p1).z*(a)
#define COPY_3DPOINT(a,b) (a).x=(b).x; (a).y=(b).y; (a).z=(b).z
#define PRINTVEC(a) fprintf(stderr,"X:%3.5f Y:%3.5f Z:%3.5f\n",(a).x,(a).y, (a).z)
#define EUCLID_SQR_DIST(a,b,result) (result)=(SQR((a).x - (b).x) + SQR((a).y - (b).y) + SQR((a).z - (b).z))
#define GROUP_SET_VOXEL(vox,a,b,c) (vox).x=(a); (vox).y=(b); (vox).z=(c);
#define GROUP_GET_VOXEL_VALUE(val,vol,a,b,c,x,y,z) \
                   if(((a)>=0) && ((a)<x) && ((b)>=0) && ((b)<y) && \
                      ((c)>=0) && ((c)<z)) \
                         val=(short)get_volume_voxel_value(vol,a, b, c,0,0); \
                   else val=255;
#define ELEM_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }
#define VEC_TRANS_MULT(t,v) ((t).x*(v).x+(t).y*(v).y+(t).z*(v).z)


#endif
