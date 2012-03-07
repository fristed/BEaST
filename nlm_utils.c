/*  nlm_utils.c
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
#include <stdio.h>

void Gkernel(float* gaussKernel, float sigma, int size)
{
    /*Contruction of a 3D gaussian kernel*/
    float radiusSquared =0.;
    int compt =0;
    int x,y,z;
    
    for(z=1;z<=2*size+1;z++)
    {
        for(x=1;x<=2*size+1;x++)
        {
            for(y=1;y<=2*size+1;y++)
            {
                radiusSquared = (x-(size+1))*(x-(size+1)) + (y-(size+1))*(y-(size+1)) + (z-(size+1))*(z-(size+1));
                gaussKernel[compt] = exp(-radiusSquared/(2*sigma*sigma));
                compt ++;
            }
        }
    }
    
    return;
}


void DisplayPatch(float* Patch, int size)
{
    
    int i,j,k;
    int Psize = 2*size+1;
    
    printf("Patch : \n ");
    for(k=0;k<Psize;k++)
    {
        for(i=0;i<Psize;i++)
        {
            for(j=0;j<Psize;j++)
            {
                printf("%f ", Patch[k*(Psize*Psize)+(i*Psize)+j]);
            }
            
            printf("\n ");
        }
        
        printf("\n ");
    }
    
    
    printf("\n ");
    
}



void ExtractPatch(float* ima, float* Patch, int x, int y, int z, int size, int sx, int sy, int sz)
{    
  int i,j,k;
  int ni1, nj1, nk1;
  int Psize = 2*size +1;
  
  for(i=0;i<(2*size+1)*(2*size+1)*(2*size+1);i++)
    Patch[i]=0.0;
  
  for(i=-size;i<=size;i++)
    {
      for(j=-size;j<=size;j++)
        {
	  for(k=-size;k<=size;k++)
            {
	      ni1=x+i;
	      nj1=y+j;
	      nk1=z+k;
              
	      if(ni1<0) ni1=-ni1;
	      if(nj1<0) nj1=-nj1;
	      if(nk1<0) nk1=-nk1;
              
	      if(ni1>=sx) ni1=2*sx-ni1-1;
	      if(nj1>=sy) nj1=2*sy-nj1-1;
	      if(nk1>=sz) nk1=2*sz-nk1-1;
              
	      Patch[(i+size)*(Psize*Psize)+((j+size)*Psize)+(k+size)] = ima[ni1*(sz*sy)+(nj1*sz)+nk1];
            }
        }
    }  
}


void ExtractPatch4D(float* ima, float* Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz)
{
    
    int i,j,k;
    int ni1, nj1, nk1;
    int Psize =  2*size +1;
    
    for(i=0;i<(2*size+1)*(2*size+1)*(2*size+1);i++)
        Patch[i]=0.0;
    
   
    for(i=-size;i<=size;i++)
      {
        for(j=-size;j<=size;j++)
	  {
            for(k=-size;k<=size;k++)
	      {
                ni1=x+i;
                nj1=y+j;
                nk1=z+k;
                
                if(ni1<0) ni1=-ni1;
                if(nj1<0) nj1=-nj1;
                if(nk1<0) nk1=-nk1;
                
                if(ni1>=sx) ni1=2*sx-ni1-1;
                if(nj1>=sy) nj1=2*sy-nj1-1;
                if(nk1>=sz) nk1=2*sz-nk1-1;
                
                Patch[(i+size)*(Psize*Psize)+((j+size)*Psize)+(k+size)]=ima[t*(sx*sy*sz)+ni1*(sz*sy)+(nj1*sz)+nk1];
            }
        }
    }
    
    
    
    
}


