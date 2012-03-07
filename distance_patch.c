/*  distance_patch.c
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

float SSDPatch(float* PatchImg, float* PatchTemplate, int f)
{
    /*SSD */
    float d;
    int i,indice;
    indice =0;
    d = 0.0;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = d + (PatchImg[i] - PatchTemplate[i]) * (PatchImg[i] - PatchTemplate[i]);
        indice++;
        
    }
    
    d=d/indice;
    
    return d;
    
}


float SSDPatchMean(float* PatchImg, float* PatchTemplate, int f, float Mean, float TMean)
{
    /*SSD */
    float d;
    int i,indice;
    indice =0;
    d = 0.0;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = d + ((PatchImg[i]-Mean) - (PatchTemplate[i]-TMean)) * ((PatchImg[i]-Mean) - (PatchTemplate[i]-TMean));
        indice++;
        
    }
    
    d=d/indice;
    
    return d;
    
}


float WSSDPatch(float* PatchImg, float* PatchTemplate, float* kernel, int f)
{
    /*Weighted SSD */
    float d,distancetotal,weight;
    int i;
    distancetotal=0.0;
    weight=0.0;
    d = 0.0;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = (PatchImg[i] - PatchTemplate[i]) * (PatchImg[i] - PatchTemplate[i]);
        distancetotal = distancetotal + kernel[i]*d;
        weight = weight + kernel[i];
        
    }
    
    d=distancetotal/weight;
    
    return d;
    
}


float NSSDPatch(float* PatchImg, float* PatchTemplate, float Mean, float TMean, float Var, float TVar, int f)
{
    /*Weighted Normalized SSD similar to Weighted CoC */
    float d,distancetotal,weight;
    int i;
    distancetotal=0.0;
    weight=0.0;
    d = 0.0;
    float eps = 0.01;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = ((PatchImg[i]-Mean)/(Var+eps)) - ((PatchTemplate[i]-TMean)/(TVar+eps));
        distancetotal = distancetotal + d*d;
        weight = weight + 1;
        
    }
    
    d=distancetotal/weight;
    
    return d;
    
}

float WNSSDPatch(float* PatchImg, float* PatchTemplate, float* kernel, float Mean, float TMean, float Var, float TVar, int f)
{
    /*Weighted Normalized SSD similar to Weighted CoC */
    float d,distancetotal,weight;
    int i;
    distancetotal=0.0;
    weight=0.0;
    d = 0.0;
    float eps = 0.01;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = ((PatchImg[i]-Mean)/(Var+eps)) - ((PatchTemplate[i]-TMean)/(TVar+eps));
        distancetotal = distancetotal + kernel[i]*d*d;
        weight = weight + kernel[i];
        
    }
    
    d=distancetotal/weight;
    
    return d;
    
}

float WCoCPatch(float* PatchImg, float* PatchTemplate, float* kernel, float Mean, float TMean, float Var, float TVar, int f)
{
     /* http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient*/
    /*Weighted CoC */
    float d,distancetotal,weight;
    int i;
    distancetotal=0.0;
    weight=0.0;
    d = 0.0;
    float eps = 0.1;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        d = d + ((PatchImg[i]-Mean))*((PatchTemplate[i]-TMean))*kernel[i];
        weight = weight + kernel[i];
        
    }
    
    d=d/weight; /*Weighted Covariance x,y*/
    if (d<0) d = 0.; /*inverted intensity*/
    
    d = (d+eps) / (sqrt(Var*TVar)+eps); /*Weighted Correlation*/
    
    return d;
    
}

float SSIMPatch(float* PatchImg, float* PatchTemplate, float Mean, float TMean, float Var, float TVar, int f)
{
    
    float d,weight,CoC,acu;
    int i;
    weight=0.0;
    d = 0.0;
    CoC =0.0;
    acu =0;
    float eps = 0.001;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        CoC = CoC + (PatchImg[i] - Mean) * (PatchTemplate[i] - TMean);
        acu++;
    }
    
    CoC = CoC / acu;
    if (CoC < 0) CoC = 0;
    
    /* Cf Wang TIP 2004*/
    d = (2*Mean*TMean + eps) * (2*CoC + eps);
    d = d / ((Mean*Mean + TMean*TMean +eps) * (Var + TVar + eps));
    
    return d;
    
}

float WSSIMPatch(float* PatchImg, float* PatchTemplate, float* kernel, float Mean, float TMean, float Var, float TVar, int f)
{
    
    float d,weight,CoC;
    int i;
    weight=0.0;
    d = 0.0;
    CoC =0.0;
    
    float eps = 0.000001;
    
    for(i=0; i< (2*f+1)*(2*f+1)*(2*f+1);i++)
    {
        CoC = CoC + (PatchImg[i] - Mean) * (PatchTemplate[i] - TMean)*kernel[i];
        weight = weight+kernel[i];
        
    }
    
    CoC = CoC / weight;
    if (CoC < 0) CoC = 0;
    
    /* Cf Wang TIP 2004*/
    d = (2*Mean*TMean + eps) * (2*CoC + eps);
    d = d / ((Mean*Mean + TMean*TMean +eps) * (Var + TVar + eps));
    
    return d;
    
}


