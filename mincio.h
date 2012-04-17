/*  mincio.h
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


#ifndef MINCIO_H
#define MINCIO_H

#include <volume_io.h>
#include "basic.h"

void set_volume(float *data, Volume vol, int *sizes);
void get_volume(float *data, Volume vol, int *sizes);
int write_volume(char *name, Volume vol, float *data);

int write_minc(char *filename, float *image, image_metadata *meta,BOOLEAN binary_mask);
image_metadata * read_minc(char *filename, float **image, int *sizes);

#endif