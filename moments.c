/*  moments.c
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
#include "basic.h"
#include "nlmseg.h"

void ComputeFirstMoment(float* ima, float* means, const int* dims, int f, float *min, float *max)
{    
    int indice=0;
    int i,j,k;
    int ii,jj,kk;
    int ni,nj,nk;
    float mean = 0.0;

    *min=FLT_MAX;
    *max=FLT_MIN;
       
    for(i=0;i<dims[0];i++)
    {
      for(j=0;j<dims[1];j++)
        {
	  for(k=0;k<dims[2];k++)
            {                                
	      mean=0;
	      indice=0;
              
	      for(ii=-f;ii<=f;ii++)
                {
		  for(jj=-f;jj<=f;jj++)
                    {
		      for(kk=-f;kk<=f;kk++)
                        {
			  ni=i+ii;
			  nj=j+jj;
			  nk=k+kk;
                          
			  if(ni<0) ni=-ni;
			  if(nj<0) nj=-nj;
			  if(nk<0) nk=-nk;
			  if(ni>=dims[0]) ni=2*dims[0]-ni-1;
			  if(nj>=dims[1]) nj=2*dims[1]-nj-1;
			  if(nk>=dims[2]) nk=2*dims[2]-nk-1;
                          
                          
			  mean = mean + ima[ni*(dims[2]*dims[1])+(nj*dims[2])+nk];			    
                          
			  indice=indice+1;
                          
                        }
                    }
                }
	      
	      mean=mean/indice;
	      
	      *min=MIN(*min,mean);
	      *max=MAX(*max,mean);
	      
	      means[i*(dims[2]*dims[1])+(j*dims[2])+k]=mean;                
              
            }
        }
    }
}


void ComputeSecondMoment(float* ima, float* means, float* variance, const int* dims,int f, float *min, float *max)
{
  
  int indice =0;
  int i,j,k;
  int ii,jj,kk;
  int ni,nj,nk;
  float var=0.0;
  
  *min=FLT_MAX;
  *max=FLT_MIN;
  
  
  for(i=0;i<dims[0];i++)
    {
      for(j=0;j<dims[1];j++)
        {
	  for(k=0;k<dims[2];k++)
            {
	      
	      var=0;
	      indice=0;
	      
	      for(ii=-f;ii<=f;ii++)
                {
		  for(jj=-f;jj<=f;jj++)
                    {
		      for(kk=-f;kk<=f;kk++)
                        {
			  ni=i+ii;
			  nj=j+jj;
			  nk=k+kk;
			  if(ni>=0 && nj>=0 && nk>0 && ni<dims[0] && nj<dims[1] && nk<dims[2])
                            {
			      var = var + ((ima[ni*(dims[2]*dims[1])+(nj*dims[2])+nk]-means[i*(dims[2]*dims[1])+(j*dims[2])+k])*(ima[ni*(dims[2]*dims[1])+(nj*dims[2])+nk]-means[i*(dims[2]*dims[1])+(j*dims[2])+k]));
			      indice=indice+1;
                            }
                        }
                    }
                }
	      var=var/(indice-1);
	      
	      *min=MIN(*min,var);
	      *max=MAX(*max,var);                
              
	      variance[i*(dims[2]*dims[1])+(j*dims[2])+k]=var;
              
            }
        }
    }
    
}


void ComputeFirstMoment4D(float* ima,float* atlas, float* means, float* Ameans,const int* dims, int f)
{
    
  int indice=0;
  int i,j,k,t;
  int ii,jj,kk;
  int ni,nj,nk;
  float mean = 0.0;
  float Amean = 0.0;
  
  
  for(k=0;k<dims[2];k++)
    {
      for(i=0;i<dims[1];i++)
        {
	  for(j=0;j<dims[0];j++)
            {
	      
              
	      Amean=0;
	      indice=0;
              
	      for(ii=-f;ii<=f;ii++)
                {
		  for(jj=-f;jj<=f;jj++)
                    {
		      for(kk=-f;kk<=f;kk++)
                        {
			  ni=i+ii;
			  nj=j+jj;
			  nk=k+kk;
                          
			  if(ni<0) ni=-ni;
			  if(nj<0) nj=-nj;
			  if(nk<0) nk=-nk;
			  if(ni>=dims[1]) ni=2*dims[1]-ni-1;
			  if(nj>=dims[0]) nj=2*dims[0]-nj-1;
			  if(nk>=dims[2]) nk=2*dims[2]-nk-1;
                          
                          
			  Amean = Amean + atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
			  indice=indice+1;
                          
                        }
                    }
                }
	      
	      Amean=Amean/indice;
	      Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j]=Amean;
              
            }
        }
    }
  
  for(t=0;t<dims[3];t++)
    {
      for(k=0;k<dims[2];k++)
        {
	  for(i=0;i<dims[1];i++)
            {
	      for(j=0;j<dims[0];j++)
                {		                                      
		  mean=0;
		  indice=0;
                  
		  for(ii=-f;ii<=f;ii++)
                    {
		      for(jj=-f;jj<=f;jj++)
                        {
			  for(kk=-f;kk<=f;kk++)
                            {
			      ni=i+ii;
			      nj=j+jj;
			      nk=k+kk;
                              
			      if(ni<0) ni=-ni;
			      if(nj<0) nj=-nj;
			      if(nk<0) nk=-nk;
			      if(ni>=dims[1]) ni=2*dims[1]-ni-1;
			      if(nj>=dims[0]) nj=2*dims[0]-nj-1;
			      if(nk>=dims[2]) nk=2*dims[2]-nk-1;
                              
			      mean = mean + ima[t*(dims[0]*dims[1]*dims[2])+(nk*(dims[0]*dims[1])+(ni*dims[0])+nj)];
			      indice=indice+1;
                              
                            }
                        }
                    }
		  
		  mean=mean/indice;
		  means[t*(dims[0]*dims[1]*dims[2])+(k*(dims[0]*dims[1])+(i*dims[0])+j)]=mean;
                  
                }
            }
        }           
    }     
}

void ComputeSecondMoment4D(float* ima,float* atlas, float* means, float* Ameans, float* variance, float* Avariance,const int* dims,int f)
{    
  int indice =0;
  int i,j,k,t;
  int ii,jj,kk;
  int ni,nj,nk;
  float var=0.0;
  float Avar=0.0;
  
  for(k=0;k<dims[2];k++)
    {
      for(i=0;i<dims[1];i++)
        {
	  for(j=0;j<dims[0];j++)
            {
	      
              
	      Avar=0;
	      indice=0;
	      for(ii=-f;ii<=f;ii++)
                {
		  for(jj=-f;jj<=f;jj++)
                    {
		      for(kk=-f;kk<=f;kk++)
                        {
			  ni=i+ii;
			  nj=j+jj;
			  nk=k+kk;
			  if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
                            {
			      Avar = Avar + ((atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j])*(atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j]));
			      indice=indice+1;
                            }
                        }
                    }
                }
	      
	      Avar=Avar/(indice-1);
	      Avariance[k*(dims[0]*dims[1])+(i*dims[0])+j]=Avar;
            }
        }
    }
  
  for(t=0;t<dims[3];t++)
    {
      for(k=0;k<dims[2];k++)
        {
	  for(i=0;i<dims[1];i++)
	    {
	      for(j=0;j<dims[0];j++)
                {		                    
		  var=0;
		  indice=0;
                  
		  for(ii=-f;ii<=f;ii++)
                    {
		      for(jj=-f;jj<=f;jj++)
                        {
			  for(kk=-f;kk<=f;kk++)
                            {
			      ni=i+ii;
			      nj=j+jj;
			      nk=k+kk;
			      if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
                                {
				  var = var + ((ima[t*(dims[0]*dims[1]*dims[2])+nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j])*(ima[t*(dims[0]*dims[1]*dims[2])+nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j]));
				  indice=indice+1;
                                }
                            }
                        }
                    }
                    
		  var=var/(indice-1);
		  variance[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j]=var;
                }
            }
        }        
    }  
}
