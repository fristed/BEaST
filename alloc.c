/*  alloc.c
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
#include <stdio.h>
#include <strings.h>
#include "array_alloc.h"
#include "basic.h"

int ***alloc_3d_int(int n1, int n2, int n3)
{
    int ***iii, **ii, *i;
    int j;
    
    iii = (int ***) calloc(n1,sizeof(int **));
    alloc_error_check(iii);
    ii = (int **) calloc(n1 * n2,sizeof(int *));
    alloc_error_check(ii);
    iii[0] = ii;
    for (j = 1; j < n1; j++) {
        iii[j] = iii[j - 1] + n2;
    }
    i = (int *) calloc(n1 * n2 * n3,sizeof(int));
    alloc_error_check(i);
    ii[0] = i;
    for (j = 1; j < n1 * n2; j++) {
        ii[j] = ii[j - 1] + n3;
    }
    return iii;
}

byte ***alloc_3d_byte(int n1, int n2, int n3)
{
    byte ***iii, **ii, *i;
    int j;
    
    iii = (byte ***) calloc(n1,sizeof(byte **));
    alloc_error_check(iii);
    ii = (byte **) calloc(n1 * n2,sizeof(byte *));
    alloc_error_check(ii);
    iii[0] = ii;
    for (j = 1; j < n1; j++) {
        iii[j] = iii[j - 1] + n2;
    }
    i = (byte *) calloc(n1 * n2 * n3,sizeof(byte));
    alloc_error_check(i);
    ii[0] = i;
    for (j = 1; j < n1 * n2; j++) {
        ii[j] = ii[j - 1] + n3;
    }
    return iii;
}

char ***alloc_3d_char(int n1, int n2, int n3)
{
  char ***iii, **ii, *i;
  int j;
    
  iii = (char ***) calloc(n1,sizeof(char **));
  alloc_error_check(iii);
  ii = (char **) calloc(n1 * n2,sizeof(char *));
  alloc_error_check(ii);
  iii[0] = ii;
  for (j = 1; j < n1; j++) {
    iii[j] = iii[j - 1] + n2;
  }
  i = (char *) calloc(n1 * n2 * n3,sizeof(char));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1 * n2; j++) {
    ii[j] = ii[j - 1] + n3;
  }
  return iii;
}

void initialize_3d_char(char ***array, int n1, int n2, int n3){
  bzero(array[0][0],n1*n2*n3*sizeof(char));
}

float ***alloc_3d_float(int n1, int n2, int n3)
{
    float ***iii, **ii, *i;
    int j;
    
    iii = (float ***) calloc(n1,sizeof(float **));
    alloc_error_check(iii);
    ii = (float **) calloc(n1 * n2,sizeof(float *));
    alloc_error_check(ii);
    iii[0] = ii;
    for (j = 1; j < n1; j++) {
        iii[j] = iii[j - 1] + n2;
    }
    i = (float *) calloc(n1 * n2 * n3,sizeof(float));
    alloc_error_check(i);
    ii[0] = i;
    for (j = 1; j < n1 * n2; j++) {
        ii[j] = ii[j - 1] + n3;
    }
    return iii;
}

double ***alloc_3d_double(int n1, int n2, int n3)
{
    double ***iii, **ii, *i;
    int j;
    
    iii = (double ***) calloc(n1,sizeof(double **));
    alloc_error_check(iii);
    ii = (double **) calloc(n1 * n2,sizeof(double *));
    alloc_error_check(ii);
    iii[0] = ii;
    for (j = 1; j < n1; j++) {
        iii[j] = iii[j - 1] + n2;
    }
    i = (double *) calloc(n1 * n2 * n3,sizeof(double));
    alloc_error_check(i);
    ii[0] = i;
    for (j = 1; j < n1 * n2; j++) {
        ii[j] = ii[j - 1] + n3;
    }
    return iii;
}

int **alloc_2d_int(int n1, int n2)
{
  int **ii, *i;
  int j;
    
  ii = (int **) calloc(n1,sizeof(int *));
  alloc_error_check(ii);
  i = (int *) calloc(n1 * n2,sizeof(int));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1; j++) {
    ii[j] = ii[j - 1] + n2;
  }
  return ii;
}

byte **alloc_2d_byte(int n1, int n2)
{
  byte **ii, *i;
  int j;
    
  ii = (byte **) calloc(n1,sizeof(byte *));
  alloc_error_check(ii);
  i = (byte *) calloc(n1 * n2,sizeof(byte));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1; j++) {
    ii[j] = ii[j - 1] + n2;
  }
  return ii;
}

float **alloc_2d_float(int n1, int n2)
{
  float **ii, *i;
  int j;
    
  ii = (float **) calloc(n1,sizeof(float *));
  alloc_error_check(ii);
  i = (float *) calloc(n1 * n2,sizeof(float));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1; j++) {
    ii[j] = ii[j - 1] + n2;
  }
  return ii;
}

double **alloc_2d_double(int n1, int n2)
{
  double **ii, *i;
  int j;
    
  ii = (double **) calloc(n1,sizeof(double *));
  alloc_error_check(ii);
  i = (double *) calloc(n1 * n2,sizeof(double));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1; j++) {
    ii[j] = ii[j - 1] + n2;
  }
  return ii;
}

char **alloc_2d_char(int n1, int n2)
{
  char **ii, *i;
  int j;
    
  ii = (char **) calloc(n1,sizeof(char *));
  alloc_error_check(ii);
  i = (char *) calloc(n1 * n2,sizeof(char));
  alloc_error_check(i);
  ii[0] = i;
  for (j = 1; j < n1; j++) {
    ii[j] = ii[j - 1] + n2;
  }
  return ii;
}

void free_3d_int(int ***ptr){
  free(ptr[0][0]);
  free(ptr[0]);
  free(ptr);
}

void free_3d_char(char ***ptr){

  free(ptr[0][0]);
  free(ptr[0]);
  free(ptr);

}

void free_2d_char(char **ptr){

  free(ptr[0]);
  free(ptr);

}

void free_2d_float(float **ptr){

  free(ptr[0]);
  free(ptr);

}

void free_2d_double(double **ptr){

  free(ptr[0]);
  free(ptr);

}

void free_3d_float(float ***ptr){

  free(ptr[0][0]);
  free(ptr[0]);
  free(ptr);

}

void free_3d_double(double ***ptr){

  free(ptr[0][0]);
  free(ptr[0]);
  free(ptr);

}


int ****alloc_4d_int(int n1, int n2, int n3, int n4)
{
    int ****iiii, ***iii, **ii, *i;
    int j;
    
    iiii = (int ****) calloc(n1,sizeof(int ***));
    alloc_error_check(iiii);
    iii = (int ***) calloc(n1 * n2,sizeof(int **));
    alloc_error_check(iii);
    ii = (int **) calloc(n1 * n2 * n3,sizeof(int *));
    alloc_error_check(ii);
    i = (int *) calloc(n1 * n2 * n3 * n4,sizeof(int));
    alloc_error_check(i);

    iiii[0] = iii;
    for (j = 1; j < n1; j++) {
        iiii[j] = iiii[j - 1] + n2;
    }

    iii[0] = ii;
    for (j = 1; j < n1 * n2; j++) {
        iii[j] = iii[j - 1] + n3;
    }

    ii[0] = i;
    for (j = 1; j < n1 * n2 * n3; j++) {
        ii[j] = ii[j - 1] + n4;
    }

    return iiii;
}

float ****alloc_4d_float(int n1, int n2, int n3, int n4)
{
    float ****iiii, ***iii, **ii, *i;
    int j;
    
    iiii = (float ****) calloc(n1,sizeof(float ***));
    alloc_error_check(iiii);
    iii = (float ***) calloc(n1 * n2,sizeof(float **));
    alloc_error_check(iii);
    ii = (float **) calloc(n1 * n2 * n3,sizeof(float *));
    alloc_error_check(ii);
    i = (float *) calloc(n1 * n2 * n3 * n4,sizeof(float));
    alloc_error_check(i);

    iiii[0] = iii;
    for (j = 1; j < n1; j++) {
        iiii[j] = iiii[j - 1] + n2;
    }

    iii[0] = ii;
    for (j = 1; j < n1 * n2; j++) {
        iii[j] = iii[j - 1] + n3;
    }

    ii[0] = i;
    for (j = 1; j < n1 * n2 * n3; j++) {
        ii[j] = ii[j - 1] + n4;
    }

    return iiii;
}

void free_4d(int ****ptr){

  free(ptr[0][0][0]);
  free(ptr[0][0]);
  free(ptr[0]);
  free(ptr);

}
