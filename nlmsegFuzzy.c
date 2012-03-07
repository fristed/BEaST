/*  nlmsegFuzzy.c
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#include "nlmseg.h"

#define MINCOUNT 100

/*Can be use to store the distance and the location*/
/*Some compilation issues can occur here, use typedef strcut with window*/
typedef struct {
    double dist;
    int x;
    int y;
    int z;
    int t;
} data_t;


float nlmsegFuzzy4D(float *subject, float *imagedata, float *maskdata, float *meandata, float *vardata, float *mask, int sizepatch, int searcharea, float beta, float threshold, int dims[3], int librarysize, float *SegSubject, float *PatchCount)
{    
  float *MeansSubj, *VarsSubj, *PatchImg, *PatchTemp, *localmask;
    float w,average,totalweight, d, Mean, TMean, Var, TVar,th,proba, min, max;
    int i,j,k,ii,jj,kk,ni,nj,nk,v,f,ndim;
    data_t storage, *tab;
    int sadims,volumesize,index;
    int count;
    int realcount;
    int p;
    //int lim; /*Number of the closest patches taken into account*/
    int t,mincount=MINCOUNT;
    int notfinished=1;
    double minidist;
    double epsi = 0.0001;
    time_t time1,time2;

    fprintf(stderr,"Patch size: %d\nSearch area: %d\nBeta: %f\nThreshold: %f\nSelection: %d\n",sizepatch,searcharea,beta,threshold,librarysize);
        
    ndim = 3;
    volumesize=dims[0]*dims[1]*dims[2];
        
    /*Patch radius*/
    f = sizepatch;
    
    /*Search Area radius*/
    v = searcharea;
        
    sadims = pow(2*v+1,ndim);
    sadims = sadims * librarysize;
    tab=(data_t*)malloc(sadims*sizeof(*tab));
    
    PatchImg=(float*) malloc((2*f+1)*(2*f+1)*(2*f+1)*sizeof(float));
    PatchTemp=(float*) malloc((2*f+1)*(2*f+1)*(2*f+1)*sizeof(float));
    
    /* allocate memory */
    MeansSubj = (float *)calloc(volumesize,sizeof(*MeansSubj));
    VarsSubj = (float *)calloc(volumesize,sizeof(*VarsSubj));
    localmask = (float *)calloc(volumesize,sizeof(*localmask));

    bcopy(mask,localmask,volumesize*sizeof(*localmask));

    fprintf(stderr,"Dimensions: %d %d %d\n",dims[0],dims[1],dims[2]);
    
    fprintf(stderr,"Computing first moment image...");
    time1=time(NULL);
    ComputeFirstMoment(subject, MeansSubj, dims, f, &min, &max);
    time2=time(NULL);
    fprintf(stderr,"done (%d sec)\nComputing second moment image...",(int)(time2-time1));
    ComputeSecondMoment(subject, MeansSubj, VarsSubj, dims, f, &min, &max);
    fprintf(stderr,"done");
    time1=time(NULL);

    if (1){

    do {
      fprintf(stderr," (%d sec)\nSegmenting          ",(int)(time1-time2));
      time2=time(NULL);
      notfinished=0;
      
      for(i=0;i<dims[0];i++)
	{
	  fprintf(stderr,"\b\b\b\b\b\b\b\b\b%3d / %3d",i+1,dims[0]);
	  for(j=0;j<dims[1];j++)            
	    {
	      for(k=0;k<dims[2];k++)                
		{                
		  index=i*(dims[2]*dims[1])+(j*dims[2])+k;

		  /* mask check */
		  if ( localmask[index]  > 0 )
		    {
		      proba=0.;
		      average=0;
		      totalweight=0;
		      count = 0;
		      realcount=0;
		      minidist = 1000000000.0;
                    
		      ExtractPatch(subject, PatchTemp, i, j, k, f, dims[0], dims[1], dims[2]);

		      TMean = MeansSubj[index];
		      TVar = VarsSubj[index];
      
		      /* go through the search space  */
		      for(ii=-v;ii<=v;ii++)
			{
			  for(jj=-v;jj<=v;jj++)
			    {
			      for(kk=-v;kk<=v;kk++)
				{
				  ni=i+ii;
				  nj=j+jj;
				  nk=k+kk;
                                                                
				  if((ni>=0) && (nj>=0) && (nk>=0) && (ni<dims[0]) && (nj<dims[1]) && (nk<dims[2]))
				    {
                                    
				      for(t=0;t<librarysize;t++)
					{
					  Mean = meandata[t*volumesize+ni*(dims[2]*dims[1])+(nj*dims[2])+nk];
					  Var = vardata[t*volumesize+ni*(dims[2]*dims[1])+(nj*dims[2])+nk];                                  
                                        
					  /*Similar Luminance and contrast -> Cf Wang TIP 2004*/
					  th = ((2 * Mean * TMean + epsi) / ( Mean*Mean + TMean*TMean + epsi))  * ((2 * sqrt(Var) * sqrt(TVar) + epsi) / (Var + TVar + epsi));
					  if(th > threshold)
					    {
					      ExtractPatch4D(imagedata, PatchImg,ni,nj,nk,t,f,dims[0],dims[1],dims[2]);
					      d =  SSDPatch(PatchImg,PatchTemp,f);
					      if (d < minidist) minidist = d;
					      storage.dist  = d;
					      storage.z = ni;
					      storage.y = nj;
					      storage.x = nk;
					      storage.t = t;
					      tab[count] = storage;
					      count ++;
					    }
					  
					}
				      
				    }
				}
			    }
			  
			}

		      /* require a minimum number of selected patches  */
		      if (count >= mincount) {

			/* Sort the distance*/
			/*This can be removed according to the chosen strategy*/
			/*qsort (tab, count, sizeof *tab, cmp);*/
			p = 0;
		      
			/*You can use the closest Patches (i.e. the 'lim' closest ones) or all the preselected patches (count)*/
			//lim = count; /*in this case, you take all the preselected patches into account*/
                    
			if (minidist<=0) minidist = epsi; /*to avoid division by zero*/		   
		      
			while (p < count)
			  {
			    storage = tab[p];
			    w = exp(- ((storage.dist)/(beta*(minidist)) ) ); /*The smoothing parameter is the minimal distance*/
			  
			    if (w>0)  
			      {
				average  = average + maskdata[(storage.t*volumesize)+(storage.z*(dims[2]*dims[1]))+(storage.y*dims[2])+storage.x]*w;			      
				totalweight = totalweight + w;
				realcount++;
			      }

			    p++;                                                                              			  			  
			  } // while

			/* We compute the probability */
			proba = average / totalweight;

			SegSubject[index] = proba;
			PatchCount[index] = realcount;
		      
		      } else {
			/* Not enough similar patches */                    
			notfinished=1;
			SegSubject[index] = -1;
		      }
		      
		    }// mask check                
		  
		} // for k
	    } // for j
	} // for i
      time1=time(NULL);

      if (notfinished){
	threshold=threshold*0.95;
	mincount=mincount*0.95;
	v=v+1;
	count=0;
	for(i=0;i<dims[0];i++)
	  {
	    for(j=0;j<dims[1];j++)	    
	      {
		for(k=0;k<dims[2];k++)                
		{
		  index=i*(dims[2]*dims[1])+(j*dims[2])+k;
		  if (SegSubject[index]<0){
		    localmask[index] = 1;
		    count++;
		  }else{
		    localmask[index] = 0;
		  }
		}
	      }
	  }
	fprintf(stderr," (redoing %d voxels) t=%f, min=%d ",count,threshold,mincount);
	free(tab);
	sadims = pow(2*v+1,ndim);
	sadims = sadims * librarysize;
	tab=(data_t*)malloc(sadims*sizeof(*tab));

      }
      
    }while (notfinished);
    
    }

    fprintf(stderr," done (%d sec, t=%f, min=%d)\n",(int)(time1-time2),threshold,mincount);

    free(tab);
    free(PatchImg);
    free(PatchTemp);
    free(MeansSubj);
    free(VarsSubj);
    free(localmask);

    return max;
}


