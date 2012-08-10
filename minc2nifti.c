/*  minc2nifti.c
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


#include "beast.h"

int main(int argc, char  *argv[] )
{
  float *input;
  int sizes[5];
  image_metadata *meta;

  if (argc<3){
    fprintf(stderr,"Usage: minc2nifti <input.mnc> <output.nii>\n");
    return 0;
  }

  meta = read_volume(argv[1], &input, sizes);

  if (meta != NULL)
    write_volume_generic(argv[2], input, meta,FALSE);

  return 0;
}
