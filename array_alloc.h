/*  array_alloc.h
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


#ifndef ARRAY_ALLOC_H
#define ARRAY_ALLOC_H

#include "basic.h"

#define alloc_error_check(p) { \
    if ((p) == NULL) { \
        fprintf(stderr, "ERROR! Allocation failure. Probably out of memory in %s:%d.\n",__FILE__,__LINE__); \
        abort(); \
    } \
}

byte **alloc_2d_byte(int n1, int n2);
int **alloc_2d_int(int n1, int n2);
char **alloc_2d_char(int n1, int n2);
float **alloc_2d_float(int n1, int n2);
double **alloc_2d_double(int n1, int n2);

byte ***alloc_3d_byte(int n1, int n2, int n3);
char ***alloc_3d_char(int n1, int n2, int n3);
int ***alloc_3d_int(int n1, int n2, int n3);
float ***alloc_3d_float(int n1, int n2, int n3);
double ***alloc_3d_double(int n1, int n2, int n3);

int ****alloc_4d_int(int n1, int n2, int n3, int n4);
float ****alloc_4d_float(int n1, int n2, int n3, int n4);

void initialize_3d_char(char ***array, int n1, int n2, int n3);

void free_3d_double(double ***ptr);
void free_3d_float(float ***ptr);
void free_2d_double(double **ptr);
void free_2d_float(float **ptr);
void free_3d_char(char ***ptr);
void free_2d_char(char **ptr);
void free_3d_int(int ***ptr);
void free_4d(int ****ptr);


void free_meta(image_metadata* meta);

#endif
