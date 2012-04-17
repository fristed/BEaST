/*  mincio.c
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


#include <float.h>
#include "mincio.h"

void set_volume(float *data, Volume vol, int *sizes){
  int i,j,k;

  for (i=0;i<sizes[0];i++)
    for (j=0;j<sizes[1];j++)
      for (k=0;k<sizes[2];k++)
        data[i*sizes[1]*sizes[2]+j*sizes[2]+k] = get_volume_real_value(vol,i,j,k,0,0);
}

void get_volume(float *data, Volume vol, int *sizes){
  int i,j,k;

  for (i=0;i<sizes[0];i++)
    for (j=0;j<sizes[1];j++)
      for (k=0;k<sizes[2];k++)
          set_volume_real_value(vol,i,j,k,0,0,data[i*sizes[1]*sizes[2]+j*sizes[2]+k]);
}

int write_volume(char *name, Volume vol, float *data){
  int i,j,k,index,sizes[5];
  float min=FLT_MAX,max=FLT_MIN;

  fprintf(stderr,"Writing %s\n",name);

  get_volume_sizes(vol,sizes);

  for (i=0;i<sizes[0];i++){
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){	  
	index=i*sizes[2]*sizes[1] + j*sizes[2] + k;
	min=MIN(min,data[index]);
	max=MAX(max,data[index]);
      }
    }
  }
  
  set_volume_real_range(vol,min,max);

  get_volume(data, vol, sizes);

  output_volume( name, NC_FLOAT,FALSE,min,max,vol,NULL, (minc_output_options *)NULL);

  return STATUS_OK;
}

int write_minc(char *filename, float *image, image_metadata *meta,BOOLEAN binary_mask){
  Volume volume;
  int i,j,k,index;
  float min=FLT_MAX,max=FLT_MIN;
  Real dummy[3];

  if(binary_mask)
  {
    volume = create_volume(3,NULL,NC_BYTE,FALSE,0.0,1.0);
    printf("Writing a binary volume...\n");
  }
  else 
    volume = create_volume(3,NULL,NC_FLOAT,FALSE,FLT_MIN,FLT_MAX);
  
  if(!binary_mask)
  {
    for (i=0;i<meta->length[0];i++){
      for (j=0;j<meta->length[1];j++){
        for (k=0;k<meta->length[2];k++){	  
          index=i*meta->length[2]*meta->length[1] + j*meta->length[2] + k;
          min=MIN(min,image[index]);
          max=MAX(max,image[index]);
        }
      }
    }
    set_volume_real_range(volume,min,max);
  } else {
    set_volume_real_range(volume,0.0,1.0);
  }
  
  set_volume_sizes(volume,meta->length);
  dummy[0]=meta->start[0];
  dummy[1]=meta->start[1];
  dummy[2]=meta->start[2];
  set_volume_starts(volume,dummy);
  dummy[0]=meta->step[0];
  dummy[1]=meta->step[1];
  dummy[2]=meta->step[2];
  set_volume_separations(volume,dummy);

  alloc_volume_data(volume);
    
  get_volume(image, volume, meta->length);
    
  if(!binary_mask)
    output_volume( filename, NC_FLOAT,FALSE,min, max,volume,NULL,(minc_output_options *)NULL);
  else
    output_volume( filename, NC_BYTE,FALSE,0, 1.0,volume,NULL,(minc_output_options *)NULL);
    
  delete_volume(volume);

  return 0;
}

image_metadata * read_minc(char *filename, float **image, int *sizes){
  Volume volume;
  Real dummy[3];
  image_metadata *meta;
  
  if( input_volume(filename, 3, NULL, NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE, &volume, (minc_input_options *) NULL ) != OK )
      return( NULL );

    meta = (image_metadata *)calloc( 1 , sizeof(image_metadata) ) ;
    meta->start = calloc(3,sizeof(float));
    meta->step = calloc(3,sizeof(float));
    meta->length = calloc(3,sizeof(int));

    get_volume_sizes(volume,sizes);
    *image=malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(**image));
    set_volume(*image, volume, sizes);

    meta->length[0]=sizes[0];
    meta->length[1]=sizes[1];
    meta->length[2]=sizes[2];    

    get_volume_starts(volume,dummy);
    meta->start[0]=dummy[0];
    meta->start[1]=dummy[1];
    meta->start[2]=dummy[2];

    get_volume_separations(volume,dummy);
    meta->step[0]=dummy[0];
    meta->step[1]=dummy[1];
    meta->step[2]=dummy[2];

    delete_volume(volume);
    
    return meta;
}

