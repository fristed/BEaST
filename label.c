/*  label.c
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


#include <stdio.h>
#include "label.h"
#include "array_alloc.h"

#define LABEL_PUSH(stack,current,voxel) current++; \
                                        stack[current]=voxel;

#define LABEL_POP(stack,current,voxel) if(current>-1) { \
                                         voxel=stack[current]; \
                                         current--; } \
                                       else \
                                         current=-2;

int label_volume_wrap_real(Volume_wrap *woriginal, Volume_wrap *wlabeled, VIO_Real iso, int **label_sizes, byte connectivity) {

  int sizes[3];
  int x,y,z;
  VIO_Real voxel;
  point3D *stack, current_voxel, new_voxel;
  unsigned short lvoxel,label=1;
  int current, i, marked=0, j,k,count;
  int **mask,label_alloc_size=0;
  short ***sdata;
  VIO_Real ***rdata;

  /* Creating mask corresponding to connectivity */
  if(connectivity==26){
    mask = malloc(connectivity*sizeof(*mask));
    mask[0] = malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    count=0;
    for(i=-1;i<2;i++)
      for(j=-1;j<2;j++)
	for(k=-1;k<2;k++){
	  if(!((i==0)&&(j==0)&&(k==0))){
	    mask[count][0] = i;
	    mask[count][1] = j;
	    mask[count++][2] = k;
	  }
	}
    if(count!=connectivity) fprintf(stderr, "ERROR: error creating mask!\n");
  }
  else if(connectivity==18){
    mask = malloc(connectivity*sizeof(*mask));
    mask[0] = malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    count=0;
    for(i=-1;i<2;i++)
      for(j=-1;j<2;j++)
	for(k=-1;k<2;k++){
	  if(!((ABS(i)==1)&&(ABS(j)==1)&&(ABS(k)==1)) && !((i==0)&&(j==0)&&(k==0))) {
	    mask[count][0] = i;
	    mask[count][1] = j;
	    mask[count++][2] = k;
	  }
	}
    if(count!=connectivity) fprintf(stderr, "ERROR: error creating mask!\n");
  }
  else if(connectivity==6){
    mask = malloc(connectivity*sizeof(*mask));
    mask[0] = malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    mask[0][0]=-1;mask[0][1]=0;mask[0][2]=0;
    mask[1][0]=0;mask[1][1]=-1;mask[1][2]=0;
    mask[2][0]=0;mask[2][1]=+1;mask[2][2]=0;
    mask[3][0]=0;mask[3][1]=0;mask[3][2]=-1;
    mask[4][0]=0;mask[4][1]=0;mask[4][2]=+1;
    mask[5][0]=1;mask[5][1]=0;mask[5][2]=0;
  }
  else {
    fprintf(stderr,"ERROR: the specified connectivity %d is not supported!\n",connectivity);
    return STATUS_ERR;
  }

  sizes[0]=woriginal->sizes[0];  sizes[1]=woriginal->sizes[1];  sizes[2]=woriginal->sizes[2];
  wlabeled->sizes[0] = sizes[0];   wlabeled->sizes[1] = sizes[1];   wlabeled->sizes[2] = sizes[2];
  wlabeled->type = NC_SHORT;
  wlabeled->type_size = sizeof(short);
  wlabeled->sign=TRUE;
  sdata=(short ***)alloc_data3D(sizes,sizeof(short));
  wlabeled->data = sdata;
  rdata = (VIO_Real ***)woriginal->data;

  /* Allocate stack */
  stack = malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(*stack));

  for(x=0;x<sizes[0];x++)
    for(y=0;y<sizes[1];y++)
      for(z=0;z<sizes[2];z++)
	sdata[x][y][z] = 0;
  
  label_alloc_size = sizes[0];
  *label_sizes = malloc(label_alloc_size*sizeof(marked));
#ifdef DEBUG
  fprintf(stderr,"label_alloc_size = %d\n",label_alloc_size);
#endif
  for(x=0;x<sizes[0];x++)
    for(y=0;y<sizes[1];y++)
      for(z=0;z<sizes[2];z++) {
	voxel = rdata[x][y][z];
	lvoxel = sdata[x][y][z];
	if((voxel == iso) && (lvoxel == 0)) {
#ifdef DEBUG
	  fprintf(stderr,"\nFound label at (%d,%d,%d) - filling...",x,y,z);
#endif
	  marked=1;
	  current = -1;
	  SET_3DPOINT(current_voxel,x,y,z);
	  sdata[x][y][z] = label; 
	  do{
	    for(i=0;i<connectivity;i++)
	      if((current_voxel.x+mask[i][0]>-1) && (current_voxel.x+mask[i][0]<sizes[0]) && 
		 (current_voxel.y+mask[i][1]>-1) && (current_voxel.y+mask[i][1]<sizes[1]) && 
		 (current_voxel.z+mask[i][2]>-1) && (current_voxel.z+mask[i][2]<sizes[2])) {
		voxel = rdata[(int)current_voxel.x+mask[i][0]][(int)current_voxel.y+mask[i][1]][(int)current_voxel.z+mask[i][2]];
		lvoxel = sdata[(int)current_voxel.x+mask[i][0]][(int)current_voxel.y+mask[i][1]][(int)current_voxel.z+mask[i][2]];
		if((voxel==iso) && (lvoxel==0)) {
		  SET_3DPOINT(new_voxel,current_voxel.x+mask[i][0],current_voxel.y+mask[i][1],current_voxel.z+mask[i][2]);
		  sdata[(int)current_voxel.x+mask[i][0]][(int)current_voxel.y+mask[i][1]][(int)current_voxel.z+mask[i][2]] = label; 		  
		  LABEL_PUSH(stack,current,new_voxel);
		 
		}
	      }	  
	    
	    LABEL_POP(stack,current,current_voxel);
	    marked++;

	  }while(current!=-2);
#ifdef DEBUG
	  fprintf(stderr,"%d voxels",marked);
#endif

	  if (label >= label_alloc_size){
#ifdef DEBUG
	    fprintf(stderr,"Too many labels - increasing allocation!\n");
#endif
	    label_alloc_size += sizes[0];
	    *label_sizes = realloc(*label_sizes,label_alloc_size*sizeof(marked));
	  }

	  (*label_sizes)[label-1] = marked;

	  label++;
	}
      }

  free(stack);
  free(mask[0]);
  free(mask);

  return label-1;

}

int getLargestObject_float(float *input, int *sizes, VIO_Real lblValue, int object_no){
  int x,y,z,kept,index,i,j,k;
  Volume_wrap wrap, wlabeled;
  short ***sdata_l;  
  double ***ddata;

  wrap.sizes[0] = sizes[0];
  wrap.sizes[1] = sizes[1];
  wrap.sizes[2] = sizes[2];
  wrap.type = NC_DOUBLE;
  wrap.type_size=NC_DOUBLE_SIZE;
  
  ddata = malloc(wrap.sizes[0]*sizeof(*ddata));
  ddata[0] = malloc(wrap.sizes[0]*wrap.sizes[1]*sizeof(**ddata));
  for(i=1;i<wrap.sizes[0];i++)
    ddata[i] = ddata[0] + i*wrap.sizes[1];
  
  ddata[0][0] = malloc(wrap.sizes[0]*wrap.sizes[1]*wrap.sizes[2]*sizeof(***ddata));
  for(i=0;i<wrap.sizes[0];i++)
    for(j=0;j<wrap.sizes[1];j++)
      ddata[i][j] = ddata[0][0] +i*wrap.sizes[1]*wrap.sizes[2] + j*wrap.sizes[2];
  for(i=0;i<wrap.sizes[0];i++)
    for(j=0;j<wrap.sizes[1];j++)
      for(k=0;k<wrap.sizes[2];k++){
	index = i*wrap.sizes[2]*wrap.sizes[1] + j*wrap.sizes[2] + k;
	ddata[i][j][k] = input[index];
      }
  wrap.data = (void *)ddata;


#ifdef DEBUG
  fprintf(stderr,"Dimensions: (%d,%d,%d)\n",wrap.sizes[0],wrap.sizes[1],wrap.sizes[2]);
#endif

  kept=getLargestObject_wrap(&wrap, &wlabeled, sizes, lblValue, object_no);

  sdata_l = (short ***)wlabeled.data;

  /* remove all other objects */
  for (x=0;x<sizes[0];x++){
    for (y=0;y<sizes[1];y++){
      for (z=0;z<sizes[2];z++){	
	if (!sdata_l[x][y][z]){
	  index = x*sizes[2]*sizes[1] + y*sizes[2] + z;
	  input[index] = 0;
	}
      }
    }
  }

  free_wrap(&wrap);
  free_wrap(&wlabeled);
  
  return kept;
  
}

int cmp_int(const void *vp, const void *vq){
  int diff = *(int *)vp - *(int *)vq;
  
  return ((diff>=0) ? ((diff>0) ? +1 : 0) : -1);  
}

int getLargestObject_wrap(Volume_wrap *wvol, Volume_wrap *wlabeled, int *sizes, VIO_Real lblValue, int object_no){
  int x,y,z,i;
  int largest_lbl,*label_sizes=NULL,*array_copy,labels,kept=0;
  short voxel;
  short ***sdata_l;

#ifdef DEBUG
  fprintf(stderr,"Wrap: Min = %f, max = %f\n",get_wrap_min(wvol),get_wrap_max(wvol));
  fprintf(stderr,"Labeling volume...");
  fprintf(stderr,"lblValue = %f\n",lblValue);
#endif

  labels = label_volume_wrap_real(wvol,wlabeled,lblValue,&label_sizes,6);

  array_copy = malloc(labels*sizeof(int));

#ifdef DEBUG
  fprintf(stderr,"Number of labels: %d\n",labels);
#endif

  if (!labels)
    return 0;

  /* make a copy of the array of sizes */
  for (i=0;i<labels;i++){
    array_copy[i] = label_sizes[i];
  }

  /* sort the array copy */
  qsort(array_copy,labels,sizeof(int),cmp_int);

#ifdef DEBUG
  fprintf(stderr,"Largest label: %d voxels\n",array_copy[labels-1]);
#endif

  /* find the appropriate label */
  largest_lbl=0;

#ifdef DEBUG
  fprintf(stderr,"label_sizes[%d] = %d\n",largest_lbl,label_sizes[largest_lbl]);
#endif

  while (array_copy[labels-object_no-1] != label_sizes[largest_lbl]){
#ifdef DEBUG
    fprintf(stderr,"label_sizes[%d] = %d\n",largest_lbl,label_sizes[largest_lbl]);
#endif
    largest_lbl++;
  }
  largest_lbl++;
#ifdef DEBUG    
  fprintf(stderr,"Keeping label %d\n",largest_lbl);
#endif
  sdata_l = (short ***)wlabeled->data;

#ifdef DEBUG
  fprintf(stderr,"Removing objects...\n");
  fprintf(stderr,"Dimensions: (%d,%d,%d)\n",sizes[0],sizes[1],sizes[2]);
#endif


  /* remove all other objects */
  for (x=0;x<sizes[0];x++){
    for (y=0;y<sizes[1];y++){
      for (z=0;z<sizes[2];z++){
	voxel = sdata_l[x][y][z];
	
	if (voxel != largest_lbl){
	  sdata_l[x][y][z] = 0;
	}else{
	  kept++;
	}
      }
    }
  }


  free(label_sizes);
  free(array_copy);

#ifdef DEBUG
  fprintf(stderr,"\tLargest object extracted (%d voxels).\n",kept);
#endif

  return kept;

} // getLargestObject()

void ***alloc_data3D(int sizes[3], byte size_element)
{
  void ***iii, **ii, *i;
  int j,limit;
    
  iii = (void ***) malloc(sizes[0]*sizeof(void **));
  alloc_error_check(iii);
  ii = (void **) malloc(sizes[0] * sizes[1] * sizeof(void *));
  alloc_error_check(ii);
  i = (void *) malloc(sizes[0] * sizes[1] * sizes[2]*size_element);
  alloc_error_check(i);

  iii[0] = ii;
  for (j = 1; j < sizes[0]; j++) {
    iii[j] = iii[0] + j*sizes[1];
  }

  limit = sizes[0] * sizes[1];

  ii[0] = i;
  for (j = 1; j < limit; j++) {
    ii[j] = ii[0] + j*sizes[2]*size_element;
  }

  return iii;
}

void free_wrap(Volume_wrap *wrap) {

  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;

  switch(wrap->type) {
    case NC_CHAR:
    case NC_BYTE:
    case NC_NAT:
      bdata = (byte ***)wrap->data;
      free(bdata[0][0]);
      free(bdata[0]);
      free(bdata);
      break;
    case NC_SHORT:
      sdata = (short ***)wrap->data;
      free(sdata[0][0]);
      free(sdata[0]);
      free(sdata);
      break;
    case NC_INT:
      idata = (int ***)wrap->data;
      free(idata[0][0]);
      free(idata[0]);
      free(idata);
      break;
    case NC_FLOAT:
      fdata = (float ***)wrap->data;
      free(fdata[0][0]);
      free(fdata[0]);
      free(fdata);
      break;
    case NC_DOUBLE:
      ddata = (double ***)wrap->data;
      free(ddata[0][0]);
      free(ddata[0]);
      free(ddata);
      break;
    default:
      fprintf(stderr,"ERROR! free_wrap() is freeing nothing!");
  }
  
  wrap->sizes[0] = 0;
  wrap->sizes[1] = 0;
  wrap->sizes[2] = 0;
}
