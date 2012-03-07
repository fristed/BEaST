/*  niftiio.c
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


#include <nifti1_io.h>
#include "niftiio.h"

//#define DEBUG

image_metadata * read_nifti(char *filename, float **image, int *sizes){
  nifti_image *nim;
  float *data;
  int x,y,z,index;
  image_metadata *meta;

  meta = (image_metadata *)calloc( 1 , sizeof(image_metadata) ) ;
  meta->start = calloc(3,sizeof(float));
  meta->step = calloc(3,sizeof(float));
  meta->length = calloc(3,sizeof(int));

  nim = nifti_image_read( filename , 1 ) ;
  if( nim == NULL ) return NULL;

#ifdef DEBUG
  fprintf(stderr,"%s\n",nifti_image_to_ascii(nim));
#endif

  data = (float *)nim->data;
  *image = malloc(nim->nvox*sizeof(**image));

  for (z=0;z<nim->nz;z++)
     for (y=0;y<nim->ny;y++)
       for (x=0;x<nim->nx;x++){
	 index = z*(nim->ny)*(nim->nx) + y*(nim->nx) + x;
	 (*image)[index]=data[index]*nim->scl_slope + nim->scl_inter;
       }

  meta->length[0]=sizes[0] = nim->nz;
  meta->length[1]=sizes[1] = nim->ny;
  meta->length[2]=sizes[2] = nim->nx;
   
  meta->step[0]=nim->dz;
  meta->step[1]=nim->dy;
  meta->step[2]=nim->dx;

  if (nim->sform_code){
    meta->start[0]=nim->sto_xyz.m[2][3];
    meta->start[1]=nim->sto_xyz.m[1][3];
    meta->start[2]=nim->sto_xyz.m[0][3];    
  } else {
    meta->start[0]=0;
    meta->start[1]=0;
    meta->start[2]=0;
  }

  nifti_image_free(nim) ;

  return meta;
}


int write_nifti_generic(char *filename, float *image, image_metadata *meta)
{
  nifti_image *nim;
  int ll,ntot;
  //float *data;

#ifdef DEBUG
  fprintf(stderr,"nifti WRITE: Dimension sizes: %d, %d, %d\n",meta->length[0],meta->length[1],meta->length[2]);
  fprintf(stderr,"nifti WRITE: Start coordinates: %f, %f, %f\n",meta->start[0],meta->start[1],meta->start[2]);
  fprintf(stderr,"nifti WRITE: Step values: %f, %f, %f\n",meta->step[0],meta->step[1],meta->step[2]);
#endif

  nim = nifti_simple_init_nim();
  nim->ndim=nim->dim[0]=3;
  nim->nx=nim->dim[1]=meta->length[2];
  nim->ny=nim->dim[2]=meta->length[1];
  nim->nz=nim->dim[3]=meta->length[0];
  nim->dx=meta->step[2];
  nim->dy=meta->step[1];
  nim->dz=meta->step[0];

  nim->sform_code = 1;

  nim->sto_xyz.m[0][0] = 1;
  nim->sto_xyz.m[0][1] = 0;
  nim->sto_xyz.m[0][2] = 0;
  nim->sto_xyz.m[0][3] = meta->start[2];
  nim->sto_xyz.m[1][0] = 0;
  nim->sto_xyz.m[1][1] = 1;
  nim->sto_xyz.m[1][2] = 0;
  nim->sto_xyz.m[1][3] = meta->start[1];
  nim->sto_xyz.m[2][0] = 0;
  nim->sto_xyz.m[2][1] = 0;
  nim->sto_xyz.m[2][2] = 1;
  nim->sto_xyz.m[2][3] = meta->start[0];
  nim->sto_xyz.m[3][0] = 0;
  nim->sto_xyz.m[3][1] = 0;
  nim->sto_xyz.m[3][2] = 0;
  nim->sto_xyz.m[3][3] = 1;

  nim->sto_ijk.m[0][0] = 1;
  nim->sto_ijk.m[0][1] = 0;
  nim->sto_ijk.m[0][2] = 0;
  nim->sto_ijk.m[0][3] = -meta->start[2];
  nim->sto_ijk.m[1][0] = 0;
  nim->sto_ijk.m[1][1] = 1;
  nim->sto_ijk.m[1][2] = 0;
  nim->sto_ijk.m[1][3] = -meta->start[1];
  nim->sto_ijk.m[2][0] = 0;
  nim->sto_ijk.m[2][1] = 0;
  nim->sto_ijk.m[2][2] = 1;
  nim->sto_ijk.m[2][3] = -meta->start[0];
  nim->sto_ijk.m[3][0] = 0;
  nim->sto_ijk.m[3][1] = 0;
  nim->sto_ijk.m[3][2] = 0;
  nim->sto_ijk.m[3][3] = 1;


  nim->nvox=nim->nx*nim->ny*nim->nz;
  nim->nbyper=4;
  nim->datatype=16;
  nim->scl_slope=1;
  nim->scl_inter=0;
  ntot = nifti_get_volsize(nim);

  nim->data = (void *)calloc(1,ntot);  /* create image memory */
  bcopy(image, nim->data, ntot); /* copy data */
  
  ll = strlen(filename) ;
  nim->fname = (char *)calloc(1,ll+1) ; strcpy(nim->fname,filename) ;
  nim->iname = (char *)calloc(1,ll+1) ; strcpy(nim->iname,filename) ;

#ifdef DEBUG
  fprintf(stderr,"%s\n",nifti_image_to_ascii(nim));
#endif

  nifti_image_write(nim) ;
  nifti_image_free(nim) ;

  return 0;
}

int write_nifti(char *filename, nifti_image *nim, float *image)
{
  int ll, outmode=-1, usegzip=0;
  char *tmpstr;
  char fname[255];
  int x,y,z,index;
  float *data;

  data = (float *)nim->data;

  for (z=0;z<nim->nz;z++)
    for (y=0;y<nim->ny;y++)
      for (x=0;x<nim->nx;x++){
	index = z*(nim->ny)*(nim->nx) + y*(nim->nx) + x;
	data[index]=(image[index] - nim->scl_inter)/nim->scl_slope;
      }
  

  /* copy the filename */
  sprintf(fname,"%s",filename);    

  /* check if the file must be zipped */
  if (!strcmp("gz",filename + strlen(filename)-2)){
    /* strip the .gz off the filename */
    fname[strlen(filename)-3] = 0;
    usegzip = 1;
  }

  /* check the format */
  if (!strcmp("hdr",fname + strlen(fname)-3)){
    outmode=0;
  }
  if (!strcmp("img",fname + strlen(fname)-3)){
    outmode=0;
  }
  if (!strcmp("nii",fname + strlen(fname)-3)){
    outmode=1;
  }
  if (!strcmp("nia",fname + strlen(fname)-3)){
    outmode=3;
  }

  if (outmode==-1){
    fprintf(stderr,"Error: Unknown format!\n");
    return 1;
  }

  fname[strlen(fname)-4] = 0;

  nim->nifti_type = outmode ;
  if( nim->fname != NULL ) free(nim->fname) ;
  if( nim->iname != NULL ) free(nim->iname) ;

   ll = strlen(fname) ;
   tmpstr = nifti_makebasename(fname);
   nim->fname = (char *)calloc(1,ll+8) ; strcpy(nim->fname,tmpstr) ;
   nim->iname = (char *)calloc(1,ll+8) ; strcpy(nim->iname,tmpstr) ;
   free(tmpstr);

   if( nim->nifti_type == 1 ){
     strcat(nim->fname,".nii") ;
     strcat(nim->iname,".nii") ;
   } else if ( nim->nifti_type == 3 ){
     strcat(nim->fname,".nia") ;
     strcat(nim->iname,".nia") ;
   } else {
     strcat(nim->fname,".hdr") ;
     strcat(nim->iname,".img") ;
   }
   if (usegzip) {
     strcat(nim->fname,".gz");
     strcat(nim->iname,".gz");
   }

   nifti_image_write( nim ) ;
   nifti_image_free( nim ) ;

  return 0;
}

char * get_datatype(int datatype){
  switch (datatype) {
  case 0: return "UNKNOWN"; break;
  case 1: return "BINARY"; break;
  case 2: return "UINT8"; break;
  case 4: return "INT16"; break;
  case 8: return "INT32"; break;
  case 16: return "FLOAT32"; break;
  case 32: return "COMPLEX64"; break;
  case 64: return "FLOAT64"; break;
  case 128: return "RGB24"; break;
  case 256: return "INT8"; break;
  case 512: return "UINT16"; break;
  case 768: return "UINT32"; break;
  case 1024: return "INT64"; break;
  case 1280: return "UINT64"; break;
  case 1536: return "FLOAT128"; break;
  case 1792: return "COMPLEX128"; break;
  case 2048: return "COMPLEX256"; break;
  case 2304: return "RGBA32"; break;
  }
  
  return "REALLY UNKNOWN";
}
