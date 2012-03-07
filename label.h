/*  label.h
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


#ifndef LABEL_H
#define LABEL_H

#include <stdlib.h>
#include "basic.h"

typedef struct Volume_wrap {
  nc_type type;
  byte type_size;
  BOOLEAN sign;
  int sizes[3];
  void *data;
} Volume_wrap;

void ***alloc_data3D(int sizes[3],byte size_element);
void free_wrap(Volume_wrap *wrap);

int getLargestObject_float(float *input, int *sizes, Real lblValue, int object_no);

#endif
