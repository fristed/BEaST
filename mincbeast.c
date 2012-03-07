/*  mincbeast.c
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


#include <stdlib.h>
#include "ParseArgv.h"
#include "array_alloc.h"
#include "nlmseg.h"
#include "beast.h"
#include "label.h"

const char LICENSE[]="Copyright (C) 2011\tSimon Fristed Eskildsen, Vladimir Fonov, \n\
\t\t\tPierrick Coupe, Jose V. Manjon\n\n\
This program comes with ABSOLUTELY NO WARRANTY; for details type 'cat COPYING'. \n\
This is free software, and you are welcome to redistribute it under certain\n\
conditions; type 'cat COPYING' for details.\n\
";


int main(int argc, char  *argv[] )
{
  char *input_file,*output_file,*libdir,*mask_file;
  char imagelist[FILENAMELENGTH], masklist[FILENAMELENGTH],meanlist[FILENAMELENGTH],varlist[FILENAMELENGTH];
  char ***images, ***masks,***means,***vars;
  int num_images,i,sizes[3][5],tmpsizes[5],volumesize,*selection,steps=3,filled=0; 
  float *imagedata,*maskdata,*meandata,*vardata,**subject,**mask,**positivemask=NULL,**segsubject,**patchcount,**filtered,max,min,**segmented;
  float *tempdata;
  int scale,scaledvolumesize,scales[3] = {1,2,4};
  int masksize=0,initialscale,targetscale,scalesteps;
  beast_conf configuration[3];
  image_metadata **meta;

  BOOLEAN outputprob = FALSE;
  BOOLEAN flipimages = FALSE;
  BOOLEAN load_moments = FALSE;
  BOOLEAN fill_output = FALSE;
  BOOLEAN verbose = FALSE;
  BOOLEAN medianfilter = FALSE;
  BOOLEAN patchfilter = FALSE;
  BOOLEAN relpath = FALSE;
  int sizepatch = 3;
  int searcharea = 5;
  double alpha = 0.2;
  double beta = 0.25;
  double threshold = 0.97;
  int selectionsize = 20;
  int voxelsize=4;
  int targetvoxelsize=1;
  char *positive_file=NULL;
  char *selection_file=NULL;
  char *count_file=NULL;
  char *conf_file=NULL;

/* Argument table */
ArgvInfo argTable[] = {
  {"-probability", ARGV_CONSTANT, (char *) TRUE, (char *) &outputprob,
     "Output the probability map instead of crisp mask."},
  {"-flip", ARGV_CONSTANT, (char *) TRUE, (char *) &flipimages,
     "Flip images around the mid-sagittal plane to increase patch count."},
  {"-load_moments", ARGV_CONSTANT, (char *) TRUE, (char *) &load_moments,
     "Do not calculate moments instead use precalculated library moments. (for optimization purposes)"},
  {"-fill", ARGV_CONSTANT, (char *) TRUE, (char *) &fill_output,
     "Fill holes in the binary output."},
  {"-median", ARGV_CONSTANT, (char *) TRUE, (char *) &medianfilter,
     "Apply a median filter on the probability map."},
  {"-nlm_filter", ARGV_CONSTANT, (char *) TRUE, (char *) &patchfilter,
     "Apply an NLM filter on the probability map (experimental)."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Enable verbose output."},
  {"-relpath", ARGV_CONSTANT, (char *) TRUE, (char *) &relpath,
     "File paths in the library are relative to library root."},
  {"-selection_num", ARGV_INT, (char *) 1, (char *) &selectionsize,
   "Specify number of selected images."},
  {"-alpha", ARGV_FLOAT, (char *) 1, (char *) &alpha,
   "Specify confidence level Alpha."},
  {"-beta", ARGV_FLOAT, (char *) 1, (char *) &beta,
   "Specify smoothness factor Beta."},
  {"-threshold", ARGV_FLOAT, (char *) 1, (char *) &threshold,
   "Specify threshold for patch selection."},
  {"-patch_size", ARGV_INT, (char *) 1, (char *) &sizepatch,
   "Specify patch size (used for all resolutions)."},
  {"-search_area", ARGV_INT, (char *) 1, (char *) &searcharea,
   "Specify size of search area for voxel size 4 (decreases by 1 for each higher resolution step)."},
  {"-initvoxel_size", ARGV_INT, (char *) 1, (char *) &voxelsize,
   "Specify initial voxel size (4, 2, or 1)."},
  {"-finalvoxel_size", ARGV_INT, (char *) 1, (char *) &targetvoxelsize,
   "Specify final voxel size (4, 2, or 1)."},
  {"-positive", ARGV_STRING, (char *) 1, (char *) &positive_file,
   "Specify mask of positive segmentation (inside mask)."},
  {"-output_selection", ARGV_STRING, (char *) 1, (char *) &selection_file,
   "Specify file to output selected files."},
  {"-count", ARGV_STRING, (char *) 1, (char *) &count_file,
   "Specify file to output the patch count."},
  {"-configuration", ARGV_STRING, (char *) 1, (char *) &conf_file,
   "Specify configuration file."},
  {NULL, ARGV_END, NULL, NULL, NULL}
};        

 fprintf(stderr,"\nmincbeast --\t\tan implementation of BEaST (Brain Extraction\n\t\t\tusing non-local Segmentation Technique) version 0.1\n\n");

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 4)) {
    (void) fprintf(stderr,LICENSE);
    (void) fprintf(stderr, 
		   "\nUsage: %s [options] <library dir> <input> <mask> <output>\n",argv[0]);
    (void) fprintf(stderr,"       %s -help\n\n", argv[0]);
    
    exit(STATUS_ERR);
  }

  libdir = argv[argc-4];
  input_file = argv[argc-3];
  mask_file = argv[argc-2]; 
  output_file = argv[argc-1];

  if (targetvoxelsize>voxelsize){
    fprintf(stderr,"ERROR! Final voxel size must smaller or equal to initial voxel size\n");
    return STATUS_ERR;
  }
  if ((voxelsize>4) || (voxelsize<1) || (voxelsize==3)){
    fprintf(stderr,"ERROR! Initial voxel size must be either 4, 2, or 1\n");
    return STATUS_ERR;
  }
  if ((targetvoxelsize>4) || (targetvoxelsize<1) || (targetvoxelsize==3)){
    fprintf(stderr,"ERROR! Final voxel size must be either 4, 2, or 1\n");
    return STATUS_ERR;
  }
  
  meta = (image_metadata **)malloc(3*sizeof(image_metadata*));

  meta[0] = read_volume(input_file, &tempdata, sizes[0]);
  volumesize=sizes[0][0]*sizes[0][1]*sizes[0][2];

  /* test if there is enough mem */
  if (0){
    if (flipimages)
      imagedata = (float *)malloc(2*4*selectionsize*volumesize*sizeof(*imagedata));
    else
      imagedata = (float *)malloc(4*selectionsize*volumesize*sizeof(*imagedata));
    if (imagedata==NULL){
      fprintf(stderr,"FATAL: Not enough memory available! (need %d MB)\n",4*selectionsize*volumesize*(1+flipimages));
      return STATUS_ERR;
    }
    free(imagedata);
  }
  
  subject = alloc_2d_float(3,volumesize*sizeof(**subject));
  cp_volume(tempdata, subject[0], sizes[0]);
  free(tempdata);

  read_volume(mask_file, &tempdata, tmpsizes);

  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])){
    fprintf(stderr,"ERROR! Mask dimension does not match image dimension!\n");
    return STATUS_ERR;
  }    
  mask = alloc_2d_float(3,volumesize*sizeof(**mask));
  cp_volume(tempdata, mask[0], sizes[0]);
  free(tempdata);

  if (positive_file!=NULL){
    read_volume(positive_file, &tempdata, tmpsizes);
    if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])){
      fprintf(stderr,"ERROR! Positive mask dimension does not match image dimension!\n");
      return STATUS_ERR;
    }    
    positivemask = alloc_2d_float(3,volumesize*sizeof(**mask));
    cp_volume(tempdata, positivemask[0], sizes[0]);
    free(tempdata);

    down_sample(positivemask[0], positivemask[1], 2, sizes[0]);  
    down_sample(positivemask[0], positivemask[2], 4, sizes[0]);  
  }

  segmented = alloc_2d_float(3,volumesize*sizeof(**segmented));

  /* downsample the subject and mask */  
  down_sample(subject[0], subject[1], 2, sizes[0]);
  down_sample(subject[0], subject[2], 4, sizes[0]);  
  down_sample(mask[0], mask[1], 2, sizes[0]);  
  down_sample(mask[0], mask[2], 4, sizes[0]);  

  targetscale=(int)(targetvoxelsize/2);
  initialscale=(int)(voxelsize/2);
  scalesteps=initialscale-targetscale+1;

  fprintf(stderr,"Scale steps: %d\n",scalesteps);

  if (conf_file != NULL){
    steps=read_configuration(conf_file, configuration);
    for (i=0;i<3;i++){
      if (configuration[i].voxelsize != MAX(i*2,1)){
	fprintf(stderr,"Syntax error in configuration file!\n");
	return STATUS_ERR;
      }
    }
  }else{
    for (i=initialscale;i>=targetscale;i--){
      configuration[i].voxelsize = MAX(i*2,1);
      configuration[i].patchsize = sizepatch;
      configuration[i].searcharea = MAX(searcharea - i + targetscale,0);
      configuration[i].alpha = alpha;
      configuration[i].beta = beta;
      configuration[i].threshold = threshold;
      configuration[i].selectionsize = selectionsize;
    }    
    configuration[targetscale].alpha = 0.5;
  }

  for (i=initialscale;i<=targetscale;i++){
    fprintf(stderr,"%d %d %d %4.2lf %4.2lf %4.2lf %d\n",configuration[i].voxelsize,configuration[i].patchsize,configuration[i].searcharea,configuration[i].alpha,configuration[i].beta,configuration[i].threshold,configuration[i].selectionsize);
  }

  images = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  masks = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  means = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  vars = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);

  for (scale=2;scale>=0;scale--){

    sprintf(imagelist,"%s/library.stx.%dmm",libdir,scales[scale]);
    sprintf(masklist,"%s/library.masks.%dmm",libdir,scales[scale]);
    if (load_moments) {
      sprintf(meanlist,"%s/library.means.%dmm",libdir,scales[scale]);
      sprintf(varlist,"%s/library.vars.%dmm",libdir,scales[scale]);
    }

    num_images=read_list(imagelist,images[scale],relpath?libdir:"");
    if (read_list(masklist,masks[scale],relpath?libdir:"")!=num_images){
      fprintf(stderr,"ERROR! Number of images and masks does not match!\n");
      return STATUS_ERR;
    }

    if (num_images<configuration[scale].selectionsize){
      fprintf(stderr,"ERROR! Cannot select more images than in the library!\n\tlibrary images: %d\n\tselection: %d\n",num_images,configuration[scale].selectionsize);
      return STATUS_ERR;
    }
    
    if (load_moments) {
      if (read_list(meanlist,means[scale],relpath?libdir:"")!=num_images){
	fprintf(stderr,"ERROR! Number of images and means does not match!\n");
	return STATUS_ERR;
      }
      if (read_list(varlist,vars[scale],relpath?libdir:"")!=num_images){
	fprintf(stderr,"ERROR! Number of images and vars does not match!\n");
	return STATUS_ERR;
      }
    }
  }
  
  read_volume(images[0][0], &tempdata, tmpsizes);
  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])){
    fprintf(stderr,"ERROR! Image dimension does not match library image dimension!\n");
    return STATUS_ERR;
  }
  free(tempdata);

  meta[1] = read_volume(images[1][0], &tempdata, sizes[1]);
  free(tempdata);

  meta[2] = read_volume(images[2][0], &tempdata, sizes[2]);
  free(tempdata);

  /* make the downsampled masks crisp */
  threshold_data(mask[1],sizes[1],0.5);
  threshold_data(mask[2],sizes[2],0.5);

  segsubject = alloc_2d_float(3,volumesize*sizeof(**segsubject));
  patchcount = alloc_2d_float(3,volumesize*sizeof(**patchcount));
  filtered = alloc_2d_float(3,volumesize*sizeof(**filtered));

  if (verbose) fprintf(stderr,"Initial voxel size: %d\nTarget voxel size: %d\n",scales[initialscale],scales[targetscale]);

  for (scale=initialscale;scale>=targetscale;scale--){
  
    selection = (int *)malloc(configuration[scale].selectionsize*sizeof(*selection));
    pre_selection(subject[scale], mask[scale], images[scale], sizes[scale], num_images, configuration[scale].selectionsize, selection, selection_file,verbose);
    
    fprintf(stderr,"Performing segmentation at %dmm resolution\nReading files ",scales[scale]);
    
    scaledvolumesize = sizes[scale][0]*sizes[scale][1]*sizes[scale][2];

    imagedata = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*imagedata));
    maskdata = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*maskdata));
    meandata = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*meandata));
    vardata = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*vardata));
    
    /* read the libray images, masks, and moments */
    for (i=0;i<configuration[scale].selectionsize;i++){
      fprintf(stderr,".");
      read_volume(images[scale][selection[i]], &tempdata, tmpsizes);     
      cp_volume(tempdata, imagedata+i*scaledvolumesize, tmpsizes);
      free(tempdata);
    }
    fprintf(stderr,"*");    
    for (i=0;i<configuration[scale].selectionsize;i++){
      fprintf(stderr,".");
      read_volume(masks[scale][selection[i]], &tempdata, tmpsizes);     
      cp_volume(tempdata, maskdata+i*scaledvolumesize, tmpsizes);
      free(tempdata);

    }
    fprintf(stderr,"*");

    if (!load_moments){
      /* calculate the mean and variance for the library images */
      /* this must be done if the selected patch size is different from the one used in the precalculation */
      for (i=0;i<configuration[scale].selectionsize;i++){
	fprintf(stderr,"c");
	ComputeFirstMoment(imagedata+i*scaledvolumesize, meandata+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
	ComputeSecondMoment(imagedata+i*scaledvolumesize, meandata+i*scaledvolumesize, vardata+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
      }
    }else{
      for (i=0;i<configuration[scale].selectionsize;i++){
	fprintf(stderr,".");
	read_volume(means[scale][selection[i]], &tempdata, tmpsizes);     
	cp_volume(tempdata, meandata+i*scaledvolumesize, tmpsizes);
	free(tempdata);
      }
      fprintf(stderr,"*");    
      for (i=0;i<configuration[scale].selectionsize;i++){
	fprintf(stderr,".");
	read_volume(vars[scale][selection[i]], &tempdata, tmpsizes);     
	cp_volume(tempdata, vardata+i*scaledvolumesize, tmpsizes);
	free(tempdata);
      }
    }
    fprintf(stderr,"\n");
    /* end of reading files */

    /* remove any disconnected parts */
    masksize = getLargestObject_float(mask[scale], sizes[scale], 1, 0);
    
    fprintf(stderr,"Mask size: %d\nAlpha: %f\n",masksize,configuration[scale].alpha);
    
    /* make sure we starting from a clean slate */
    wipe_data(segsubject[scale],sizes[scale],0.0);

    if (flipimages){
      /* doubling the library selection by flipping images along the mid-sagittal plane */
      imagedata = (float *)realloc(imagedata,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*imagedata));
      maskdata = (float *)realloc(maskdata,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*maskdata));
      meandata = (float *)realloc(meandata,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*meandata));
      vardata = (float *)realloc(vardata,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*vardata));
      
      for (i=0;i<configuration[scale].selectionsize;i++){
	flip_data(imagedata+i*scaledvolumesize, imagedata+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
	flip_data(maskdata+i*scaledvolumesize, maskdata+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
	flip_data(meandata+i*scaledvolumesize, meandata+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
	flip_data(vardata+i*scaledvolumesize, vardata+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
      }

      max = nlmsegFuzzy4D(subject[scale], imagedata, maskdata, meandata, vardata, mask[scale], configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], configuration[scale].selectionsize*2, segsubject[scale], patchcount[scale]);
    }else{
      max = nlmsegFuzzy4D(subject[scale], imagedata, maskdata, meandata, vardata, mask[scale], configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], configuration[scale].selectionsize, segsubject[scale], patchcount[scale]);
    }
    
    free(imagedata);
    free(maskdata);
    free(meandata);
    free(vardata);

    
    if (positive_file!=NULL){
      /* add the certain positive segmentation (inside mask) */
      add_mask_data(segsubject[scale],positivemask[scale],sizes[scale]);  
    }
    
    /* add the certain segmentation from the previous scale */
    add_mask_data(segsubject[scale],segmented[scale],sizes[scale]);  

    if (medianfilter){
      median_filter(segsubject[scale], sizes[scale], 3);
    }

    /* the patch filter is experimental */
    if (patchfilter){
      wipe_data(filtered[scale],sizes[scale],0.0);
      wipe_data(patchcount[scale],sizes[scale],0.0);
      max = nlmfilter(subject[scale], mask[scale], segsubject[scale], 2*configuration[scale].patchsize, 2*configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], filtered[scale], patchcount[scale]);
      combine_maps(segsubject[scale], filtered[scale], mask[scale], sizes[scale]);
    }

    if (scale > targetscale){
      /* if performing a higher resolution step, upsample the result and create new mask */
      resize_trilinear(segsubject[scale], sizes[scale], sizes[scale-1], segsubject[scale-1]);
      masksize=update_mask(segsubject[scale-1], mask[scale-1], segmented[scale-1], sizes[scale-1], configuration[scale].alpha, 1.0-configuration[scale].alpha);
    }
        
    free(selection);
  } // for each scale

  if (!outputprob){
    fprintf(stderr,"Thresholding estimator at %f\n",configuration[targetscale].alpha);
    threshold_data(segsubject[targetscale], sizes[targetscale], configuration[targetscale].alpha);
    getLargestObject_float(segsubject[targetscale], sizes[targetscale], 1, 0);

    if (fill_output){
      wipe_data(mask[targetscale], sizes[targetscale], 1.0);
      filled = flood_fill_float(segsubject[targetscale], mask[targetscale], sizes[targetscale], 0, 0, 0, 0, 6);
      segsubject[targetscale]=mask[targetscale];
    }

  }

  write_volume_generic(output_file, segsubject[targetscale], meta[targetscale]);

  if (count_file!=NULL){
    write_volume_generic(count_file, patchcount[targetscale], meta[targetscale]);
  }

  free_2d_float(subject);  
  if (positive_file!=NULL)
    free_2d_float(positivemask);
  free_2d_float(segmented);
  free_2d_float(segsubject);
  free_2d_float(patchcount);

  free_3d_char(images);
  free_3d_char(masks);
  free_3d_char(means);
  free_3d_char(vars);

  return STATUS_OK;
}
