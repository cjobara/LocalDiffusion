/*
 *
 *  InferenceMAP v1.0
 *  4/3/2015
 *
 *  Author: Mohamed El Beheiry, Physico-Chimie Curie, Institut Curie
 *  		Jean-Baptiste Masson, Physics of Biological Systems, Institut Pasteur
 *  Contact e-mail: mohamed.elbeheiry@gmail.com
 *  Copyright (c) 2015, Mohamed El Beheiry, Jean-Baptiste Masson, Institut Curie, Institut Pasteur
 *  All rights Reserved.
 *
 *  InferenceMAP is released under an "academic use only" licence.
 *  Details of which are provided in the InferenceMAP_License.doc file.
 *  Usage of InferenceMAP requires acceptance of this license.
 *
 *  User instructions for using InferenceMAP are provided in the InferenceMAP User Manual.
 *
 */

#include "stdafx.h"

#include "globals.h"
#include "clustering.h"

void kmeans_01_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] )

//****************************************************************************80
//
//  Purpose:
//
//   KMEANS_01 applies the K-Means algorithm.
//
//  Discussion:
//
//    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
//    observations are to be allocated to CLUSTER_NUM clusters in such
//    a way that the within-cluster sum of squares is minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2011
//
//  Author:
//
//    Original FORTRAN77 version by David Sparks.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Sparks,
//    Algorithm AS 58:
//    Euclidean Cluster Analysis,
//    Applied Statistics,
//    Volume 22, Number 1, 1973, pages 126-130.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the points.
//
//    Output, int CLUSTER[POINT_NUM], indicates which cluster
//    each point belongs to.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number
//    of points in each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the
//    cluster energies.
//
{
  double dc;
  double de;
  double *f;
  int i;
  int il;
  int ir;
  int j;
  int j2;
  int k;
  double point_energy;
  double point_energy_min;
  int swap;

//  float distanceType = (float)iMAP->irregularMeshingGui->distanceButton->value();

  it_num = 0;
//
//  Idiot checks.
//
  if ( cluster_num < 1 )
  {
	fprintf(stderr,"  CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
	fprintf(stderr,"  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
	fprintf(stderr,"  POINT_NUM < 1.\n");
    exit ( 1 );
  }

//
//  For each observation, calculate the distance from each cluster
//  center, and assign to the nearest.
//
  for ( j = 0; j < point_num; j++ )
  {
    point_energy_min = r8_huge ( );
    cluster[j] = -1;

    for ( k = 0; k < cluster_num; k++ )
    {
      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
          fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
      }

      if ( point_energy < point_energy_min )
      {
        point_energy_min = point_energy;
        cluster[j] = k;
      }
    }
  }
//
//  Determine the cluster population counts.
//
  i4vec_zero ( cluster_num, cluster_population );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    cluster_population[k] = cluster_population[k] + 1;
  }
//
//  Calculate the mean and sum of squares for each cluster.
//
  r8vec_zero ( dim_num * cluster_num, cluster_center );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num]
        + point[i+j*dim_num];
    }
  }

  for ( k = 0; k < cluster_num; k++ )
  {
    if ( 0 < cluster_population[k] )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] /
         ( double ) cluster_population[k];
      }
    }
  }
//
//  Set the point energies.
//
  f = r8vec_zero_new ( point_num );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      f[j] = f[j] + fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
    }
  }
//
//  Set the cluster energies.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    cluster_energy[k] = cluster_energy[k] + f[j];
  }
//
//  Adjust the point energies by a weight factor.
//
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    if ( 1 < cluster_population[k] )
    {
      f[j] = f[j] * ( double ) ( cluster_population[k] )
        / ( double ) ( cluster_population[k] - 1 );
    }
  }
//
//  Examine each observation in turn to see if it should be
//  reassigned to a different cluster.
//
  it_num = 0;
  char label[50];
  while ( it_num < it_max )
  {

	sprintf(label,"Clustering Iteration %i\n",it_num);
	textDisplayUpdate(label);
    it_num = it_num + 1;

    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      il = cluster[j];
      ir = il;

      if ( cluster_population[il] <= 1 )
      {
        continue;
      }

      dc = f[j];

      for ( k = 0; k < cluster_num; k++ )
      {
        if ( k != il )
        {
          de = 0.0;
          for ( i = 0; i < dim_num; i++ )
          {
            de = de + fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
          }
          de = de * ( double ) cluster_population[k]
             / ( double ) ( cluster_population[k] + 1 );

          if ( de < dc )
          {
            dc = de;
            ir = k;
          }
        }
      }
//
//  If the lowest value was obtained by staying in the current cluster,
//  then cycle.
//
      if ( ir == il )
      {
        continue;
      }
//
//  Reassign the point from cluster IL to cluster IR.
//
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+il*dim_num] = ( cluster_center[i+il*dim_num]
          * ( double ) ( cluster_population[il] ) - point[i+j*dim_num] )
          / ( double ) ( cluster_population[il] - 1 );

        cluster_center[i+ir*dim_num] = ( cluster_center[i+ir*dim_num]
          * ( double ) ( cluster_population[ir] ) + point[i+j*dim_num] )
          / ( double ) ( cluster_population[ir] + 1 );
      }
      cluster_energy[il] = cluster_energy[il] - f[j];
      cluster_energy[ir] = cluster_energy[ir] + dc;
      cluster_population[ir] = cluster_population[ir] + 1;
      cluster_population[il] = cluster_population[il] - 1;

      cluster[j] = ir;
//
//  Adjust the value of F for points in clusters IL and IR.
//
      for ( j2 = 0; j2 < point_num; j2++ )
      {
        k = cluster[j2];

        if ( k == il || k == ir )
        {
          f[j2] = 0.0;
          for ( i = 0; i < dim_num; i++ )
          {
            f[j2] = f[j2] + fabs ( point[i+j2*dim_num] - cluster_center[i+k*dim_num] );
          }

          if ( 1 < cluster_population[k] )
          {
            f[j2] = f[j2] * ( double ) ( cluster_population[k] )
              / ( ( double ) ( cluster_population[k] - 1 ) );
          }
        }
      }
      swap = swap + 1;
    }
//
//  Exit if no reassignments were made during this iteration.
//
    if ( swap == 0 )
    {
      break;
    }
  }
//
//  Compute the cluster energies.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_energy = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_energy = point_energy +
        fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
    }
    cluster_energy[k] = cluster_energy[k] + point_energy;
  }

  delete [] f;

  return;
}

void kmeans_01_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] )

//****************************************************************************80
//
//  Purpose:
//
//   KMEANS_01 applies the K-Means algorithm.
//
//  Discussion:
//
//    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
//    observations are to be allocated to CLUSTER_NUM clusters in such
//    a way that the within-cluster sum of squares is minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2011
//
//  Author:
//
//    Original FORTRAN77 version by David Sparks.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Sparks,
//    Algorithm AS 58:
//    Euclidean Cluster Analysis,
//    Applied Statistics,
//    Volume 22, Number 1, 1973, pages 126-130.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the points.
//
//    Output, int CLUSTER[POINT_NUM], indicates which cluster
//    each point belongs to.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number
//    of points in each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the
//    cluster energies.
//
{
  double dc;
  double de;
  double *f;
  int i;
  int il;
  int ir;
  int j;
  int j2;
  int k;
  double point_energy;
  double point_energy_min;
  int swap;

//  float distanceType = (float)iMAP->irregularMeshingGui->distanceButton->value();

  it_num = 0;
//
//  Idiot checks.
//
  if ( cluster_num < 1 )
  {
	fprintf(stderr,"  CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
	fprintf(stderr,"  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
	fprintf(stderr,"  POINT_NUM < 1.\n");
    exit ( 1 );
  }
//
//  For each observation, calculate the distance from each cluster
//  center, and assign to the nearest.
//
  for ( j = 0; j < point_num; j++ )
  {
    point_energy_min = r8_huge ( );
    cluster[j] = -1;

    for ( k = 0; k < cluster_num; k++ )
    {
      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
          pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
      }

      if ( point_energy < point_energy_min )
      {
        point_energy_min = point_energy;
        cluster[j] = k;
      }
    }
  }
//
//  Determine the cluster population counts.
//
  i4vec_zero ( cluster_num, cluster_population );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    cluster_population[k] = cluster_population[k] + 1;
  }
//
//  Calculate the mean and sum of squares for each cluster.
//
  r8vec_zero ( dim_num * cluster_num, cluster_center );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num]
        + point[i+j*dim_num];
    }
  }

  for ( k = 0; k < cluster_num; k++ )
  {
    if ( 0 < cluster_population[k] )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] /
         ( double ) cluster_population[k];
      }
    }
  }
//
//  Set the point energies.
//
  f = r8vec_zero_new ( point_num );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      f[j] = f[j] + pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
  }
//
//  Set the cluster energies.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    cluster_energy[k] = cluster_energy[k] + f[j];
  }
//
//  Adjust the point energies by a weight factor.
//
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    if ( 1 < cluster_population[k] )
    {
      f[j] = f[j] * ( double ) ( cluster_population[k] )
        / ( double ) ( cluster_population[k] - 1 );
    }
  }
//
//  Examine each observation in turn to see if it should be
//  reassigned to a different cluster.
//
  it_num = 0;

  char label[50];
  while ( it_num < it_max )
  {
		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			char progress[60];
			sprintf(progress,"Clustering Iteration %i\n",it_num);
			textDisplayUpdate(progress);
		}
    it_num = it_num + 1;

    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      il = cluster[j];
      ir = il;

      if ( cluster_population[il] <= 1 )
      {
        continue;
      }

      dc = f[j];

      for ( k = 0; k < cluster_num; k++ )
      {
        if ( k != il )
        {
          de = 0.0;
          for ( i = 0; i < dim_num; i++ )
          {
            de = de + pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
          }
          de = de * ( double ) cluster_population[k]
             / ( double ) ( cluster_population[k] + 1 );

          if ( de < dc )
          {
            dc = de;
            ir = k;
          }
        }
      }
//
//  If the lowest value was obtained by staying in the current cluster,
//  then cycle.
//
      if ( ir == il )
      {
        continue;
      }
//
//  Reassign the point from cluster IL to cluster IR.
//
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+il*dim_num] = ( cluster_center[i+il*dim_num]
          * ( double ) ( cluster_population[il] ) - point[i+j*dim_num] )
          / ( double ) ( cluster_population[il] - 1 );

        cluster_center[i+ir*dim_num] = ( cluster_center[i+ir*dim_num]
          * ( double ) ( cluster_population[ir] ) + point[i+j*dim_num] )
          / ( double ) ( cluster_population[ir] + 1 );
      }
      cluster_energy[il] = cluster_energy[il] - f[j];
      cluster_energy[ir] = cluster_energy[ir] + dc;
      cluster_population[ir] = cluster_population[ir] + 1;
      cluster_population[il] = cluster_population[il] - 1;

      cluster[j] = ir;
//
//  Adjust the value of F for points in clusters IL and IR.
//
      for ( j2 = 0; j2 < point_num; j2++ )
      {
        k = cluster[j2];

        if ( k == il || k == ir )
        {
          f[j2] = 0.0;
          for ( i = 0; i < dim_num; i++ )
          {
            f[j2] = f[j2] + pow ( point[i+j2*dim_num] - cluster_center[i+k*dim_num], 2 );
          }

          if ( 1 < cluster_population[k] )
          {
            f[j2] = f[j2] * ( double ) ( cluster_population[k] )
              / ( ( double ) ( cluster_population[k] - 1 ) );
          }
        }
      }
      swap = swap + 1;
    }
//
//  Exit if no reassignments were made during this iteration.
//
    if ( swap == 0 )
    {
      break;
    }
  }
//
//  Compute the cluster energies.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_energy = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_energy = point_energy +
        pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
    cluster_energy[k] = cluster_energy[k] + point_energy;
  }

  delete [] f;

  return;
}

double r8_huge ()

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}

int *i4vec_negone_new ( int n )
//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_NEGONE_NEW creates an I4VEC and sets it to -1.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_NEGONE_NEW[N], a vector of -1's.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = -1;
  }
  return a;
}

double *cluster_initialize_5 ( int dim_num, int point_num, int cluster_num,
  double point[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//   CLUSTER_INITIALIZE_5 initializes the cluster centers to random values.
//
//  Discussion:
//
//    In this case, each cluster center is a random convex combination
//    of the data points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates
//    of the points.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the coordinates of the cluster centers.
//
{
  double *cluster_center;
  double column_sum;
  double *factor;
  int j;
  int k;
//
//  Get a PxC block of random factors.
//
  factor = r8mat_uniform_01_new ( point_num, cluster_num, seed );
//
//  Make each column of factors have unit sum.
//
  for ( k = 0; k < cluster_num; k++ )
  {
    column_sum = 0.0;
    for ( j = 0; j < point_num; j++ )
    {
      column_sum = column_sum + factor[j+k*point_num];
    }
    for ( j = 0; j < point_num; j++ )
    {
      factor[j+k*point_num] = factor[j+k*point_num] / column_sum;
    }
  }
//
//  Set centers = points * factors.
//
  cluster_center = r8mat_mm_new ( dim_num, point_num, cluster_num, point,
    factor );

  delete [] factor;

  return cluster_center;
}

double *cluster_variance_compute ( int dim_num, int point_num, int cluster_num,
  double point[], int cluster[], double cluster_center[] )

//****************************************************************************80
//
//  Purpose:
//
//   CLUSTER_VARIANCE_COMPUTE computes the variance of the clusters.
//
//  Discussion:
//
//    The cluster variance (from the cluster center) is the average of the
//    sum of the squares of the distances of each point in the cluster to the
//    cluster center.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of data points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
//
//    Input, int CLUSTER[POINT_NUM], the cluster to which each
//    data point belongs.
//
//    Input, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the
//    centers associated with the minimal energy clustering.
//
//    Output, double CLUSTER_VARIANCE_COMPUTE[CLUSTER_NUM], the variance
//    associated with each cluster.
//
{
  int *cluster_population;
  double *cluster_variance;
  int i;
  int j;
  int k;
  double point_variance;

  cluster_population = i4vec_zero_new ( cluster_num );
  cluster_variance = r8vec_zero_new ( cluster_num );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_variance = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_variance = point_variance +
        pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
    cluster_variance[k] = cluster_variance[k] + point_variance;
    cluster_population[k] = cluster_population[k] + 1;
  }

  for ( k = 0; k < cluster_num; k++ )
  {
    if ( 0 < cluster_population[k] )
    {
      cluster_variance[k] = cluster_variance[k] / cluster_population[k];
    }
  }

  delete [] cluster_population;

  return cluster_variance;
}

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}

double *r8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}

void hmeans_01_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_01 applies the H-Means algorithm.
//
//  Discussion:
//
//    The data for the H-Means problem is a set of N points X in
//    M-dimensions, and a desired number of clusters K.
//
//    The goal is to determine K points Z, called cluster centers, so that
//    if we associate each point X with its nearest Z value, we minimize
//    the standard deviation or cluster energy.  Writing CLUSTER(I) to
//    indicate the index of the nearest cluster center to point X(I), the
//    energy can be written as:
//
//      Energy = Sum ( 1 <= I <= N ) || X(I) - Z(CLUSTER(I)) ||^2
//
//    where
//
//      || X - Z ||^2 = Sum ( 1 <= J <= M ) ( X(J) - Z(J) )^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Wendy Martinez, Angel Martinez,
//    Computational Statistics Handbook with MATLAB,
//    pages 373-376,
//    Chapman and Hall / CRC, 2002.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of data points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int &IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the
//    centers associated with the minimal energy clustering.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM],
//    the populuation of each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy
//    associated with each cluster.
//
{
  double *centroid;
  bool debug = true;
  int i;
  int j;
  int k;
  int k2;
  int missed;
  double point_energy;
  double point_energy_min;
  int swap;
//
//  Data checks.
//
  if ( cluster_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n  POINT_NUM < 1.\n");
    exit ( 1 );
  }

  if ( it_max < 0 )
  {
    printf("HMEANS_01 - Fatal error!\n  IT_MAX < 0.\n");
    exit ( 1 );
  }
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
  for ( j = 0; j < point_num; j++ )
  {
    if ( cluster[j] < 0 || cluster_num <= cluster[j] )
    {
      point_energy_min = r8_huge ( );

      for ( k = 0; k < cluster_num; k++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
          fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }
      }
    }
  }
  it_num = 0;

  char label[50];
  while ( it_num < it_max )
  {
		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			char progress[60];
			sprintf(progress,"Clustering Iteration %i\n",it_num);
			textDisplayUpdate(progress);
		}
		it_num = it_num + 1;
//
//  #1:
//  Assign each point to the cluster of its nearest center.
//
    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      point_energy_min = r8_huge ( );
      k = cluster[j];

      for ( k2 = 0; k2 < cluster_num; k2++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
          fabs ( point[i+j*dim_num] - cluster_center[i+k2*dim_num] );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k2;
        }
      }

      if ( k != cluster[j] )
      {
        swap = swap + 1;
      }
    }
//
//  Terminate if no points were swapped.
//
    if ( 1 < it_num )
    {
      if ( swap == 0 )
      {
        break;
      }
    }
//
//  #2:
//  Determine the total energy of the new clustering with current centroids.
//
    r8vec_zero ( cluster_num, cluster_energy );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];

      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
        fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
      }
      cluster_energy[k] = cluster_energy[k] + point_energy;
    }
//
//  #3:
//  Determine the centroids of the clusters.
//
    centroid = r8vec_zero_new ( dim_num * cluster_num );
    i4vec_zero ( cluster_num, cluster_population );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];
      cluster_population[k] = cluster_population[k] + 1;
      for ( i = 0; i < dim_num; i++ )
      {
        centroid[i+k*dim_num] = centroid[i+k*dim_num] + point[i+j*dim_num];
      }
    }
//
//  Now divide by the population to get the centroid.
//  But if a center has no population, pick a point at random.
//
    missed = 0;

    for ( k = 0; k < cluster_num; k++ )
    {
      if ( cluster_population[k] != 0 )
      {
        for ( i = 0; i < dim_num; i++ )
        {
          centroid[i+k*dim_num] = centroid[i+k*dim_num]
          / ( double ) ( cluster_population[k] );
        }
      }
      else
      {
        for ( i = 0; i < dim_num; i++ )
        {
          centroid[i+k*dim_num] = point[i+missed*dim_num];
        }
        missed = missed + 1;
      }
    }

    for ( k = 0; k < cluster_num; k++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = centroid[i+k*dim_num];
      }
    }

    delete [] centroid;
//
//  #4:
//  Determine the total energy of the current clustering with new centroids.
//
    r8vec_zero ( cluster_num, cluster_energy );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];

      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
        fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
      }
      cluster_energy[k] = cluster_energy[k] + point_energy;
    }
  }
  return;
}

void hmeans_01_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_01 applies the H-Means algorithm.
//
//  Discussion:
//
//    The data for the H-Means problem is a set of N points X in
//    M-dimensions, and a desired number of clusters K.
//
//    The goal is to determine K points Z, called cluster centers, so that
//    if we associate each point X with its nearest Z value, we minimize
//    the standard deviation or cluster energy.  Writing CLUSTER(I) to
//    indicate the index of the nearest cluster center to point X(I), the
//    energy can be written as:
//
//      Energy = Sum ( 1 <= I <= N ) || X(I) - Z(CLUSTER(I)) ||^2
//
//    where
//
//      || X - Z ||^2 = Sum ( 1 <= J <= M ) ( X(J) - Z(J) )^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Wendy Martinez, Angel Martinez,
//    Computational Statistics Handbook with MATLAB,
//    pages 373-376,
//    Chapman and Hall / CRC, 2002.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of data points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int &IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the
//    centers associated with the minimal energy clustering.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM],
//    the populuation of each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy
//    associated with each cluster.
//
{
  double *centroid;
  bool debug = true;
  int i;
  int j;
  int k;
  int k2;
  int missed;
  double point_energy;
  double point_energy_min;
  int swap;
//
//  Data checks.
//
  if ( cluster_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
    printf("HMEANS_01 - Fatal error!\n  POINT_NUM < 1.\n");
    exit ( 1 );
  }

  if ( it_max < 0 )
  {
    printf("HMEANS_01 - Fatal error!\n  IT_MAX < 0.\n");
    exit ( 1 );
  }
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
  for ( j = 0; j < point_num; j++ )
  {
    if ( cluster[j] < 0 || cluster_num <= cluster[j] )
    {
      point_energy_min = r8_huge ( );

      for ( k = 0; k < cluster_num; k++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
          pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }
      }
    }
  }
  it_num = 0;

  char label[50];
  while ( it_num < it_max )
  {

		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			char progress[60];
			sprintf(progress,"Clustering Iteration %i\n",it_num);
			textDisplayUpdate(progress);
		}
		it_num = it_num + 1;
//
//  #1:
//  Assign each point to the cluster of its nearest center.
//
    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      point_energy_min = r8_huge ( );
      k = cluster[j];

      for ( k2 = 0; k2 < cluster_num; k2++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
          pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k2;
        }
      }

      if ( k != cluster[j] )
      {
        swap = swap + 1;
      }
    }
//
//  Terminate if no points were swapped.
//
    if ( 1 < it_num )
    {
      if ( swap == 0 )
      {
        break;
      }
    }
//
//  #2:
//  Determine the total energy of the new clustering with current centroids.
//
    r8vec_zero ( cluster_num, cluster_energy );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];

      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
        pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
      }
      cluster_energy[k] = cluster_energy[k] + point_energy;
    }
//
//  #3:
//  Determine the centroids of the clusters.
//
    centroid = r8vec_zero_new ( dim_num * cluster_num );
    i4vec_zero ( cluster_num, cluster_population );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];
      cluster_population[k] = cluster_population[k] + 1;
      for ( i = 0; i < dim_num; i++ )
      {
        centroid[i+k*dim_num] = centroid[i+k*dim_num] + point[i+j*dim_num];
      }
    }
//
//  Now divide by the population to get the centroid.
//  But if a center has no population, pick a point at random.
//
    missed = 0;

    for ( k = 0; k < cluster_num; k++ )
    {
      if ( cluster_population[k] != 0 )
      {
        for ( i = 0; i < dim_num; i++ )
        {
          centroid[i+k*dim_num] = centroid[i+k*dim_num]
          / ( double ) ( cluster_population[k] );
        }
      }
      else
      {
        for ( i = 0; i < dim_num; i++ )
        {
          centroid[i+k*dim_num] = point[i+missed*dim_num];
        }
        missed = missed + 1;
      }
    }

    for ( k = 0; k < cluster_num; k++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = centroid[i+k*dim_num];
      }
    }

    delete [] centroid;
//
//  #4:
//  Determine the total energy of the current clustering with new centroids.
//
    r8vec_zero ( cluster_num, cluster_energy );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];

      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        point_energy = point_energy +
        pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
      }
      cluster_energy[k] = cluster_energy[k] + point_energy;
    }
  }
  return;
}

void hmeans_02_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_02 applies the H-Means algorithm.
//
//  Discussion:
//
//    This is a simple routine to group a set of points into K clusters,
//    each with a center point, in such a way that the total cluster
//    energy is minimized.  The total cluster energy is the sum of the
//    squares of the distances of each point to the center of its cluster.
//
//    The algorithm begins with an initial estimate for the cluster centers:
//
//    1. The points are assigned to the nearest cluster centers.
//
//    2. The iteration exit ( 1 );s if the total energy has not changed
//        significantly, or we have reached the maximum number of iterations.
//
//    3. Each cluster center is replaced by the centroid of the points
//       in the cluster.
//
//    4. Return to step 1.
//
//    The algorithm may fail to find the best solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates
//    of the points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the coordinates of the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number of
//    points assigned to each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of
//    the clusters.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
{
  bool debug = false;
  int i;
  int j;
  int k;
  int k2;
  double point_energy;
  double point_energy_min;
  int swap;
//
//  Data checks.
//
  if ( cluster_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  POINT_NUM < 1.\n");
    exit ( 1 );
  }

  if ( it_max < 0 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  IT_MAX < 0.\n");
    exit ( 1 );
  }
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
  for ( j = 0; j < point_num; j++ )
  {
    if ( cluster[j] < 0 || cluster_num <= cluster[j] )
    {
      point_energy_min = r8_huge ( );
      for ( k = 0; k < cluster_num; k++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
            fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
        }
        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }
      }
    }
  }

  it_num = 0;

  char label[50];
  for ( ; ; )
  {
//
//  Given centers, assign points to nearest center.
//
    i4vec_zero ( cluster_num, cluster_population );
    r8vec_zero ( cluster_num, cluster_energy );

    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      point_energy_min = r8_huge ( );
      k = cluster[j];

      for ( k2 = 0; k2 < cluster_num; k2++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
            fabs ( point[i+j*dim_num] - cluster_center[i+k2*dim_num] );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k2;
        }
      }
      if ( k != cluster[j] )
      {
        swap = swap + 1;
      }
      k = cluster[j];
      cluster_energy[k] = cluster_energy[k] + point_energy_min;
      cluster_population[k] = cluster_population[k] + 1;
    }

    if ( debug )
    {
//      cout << "  " << setw(3) << it_num
//           << "  " << setw(14) << r8vec_sum ( cluster_num, cluster_energy ) << "\n";
    }

    if ( 0 < it_num )
    {
      if ( swap == 0 )
      {
        break;
      }
    }

    if ( it_max <= it_num )
    {
      break;
    }

	if (iMAP->updateDisplay) {
		iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
		char progress[60];
		sprintf(progress,"Clustering Iteration %i\n",it_num);
		textDisplayUpdate(progress);
	}
	it_num = it_num + 1;
//
//  Given points in cluster, replace center by centroid.
//
    r8vec_zero ( dim_num * cluster_num, cluster_center );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num]
          + point[i+j*dim_num];
      }
    }

    for ( k = 0; k < cluster_num; k++ )
    {
      if ( cluster_population[k] != 0 )
      {
        for ( i = 0; i < dim_num; i++ )
        {
          cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] /
            ( double ) ( cluster_population[k] );
        }
      }
      else
      {
        j = i4_uniform ( 0, point_num - 1, seed );
        for ( i = 0; i < dim_num; i++ )
        {
          cluster_center[i+k*dim_num] = point[i+j*dim_num];
        }
      }
    }
  }
//
//  Compute the energy based on the final value of the cluster centers.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_energy = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_energy = point_energy +
      fabs ( point[i+j*dim_num] - cluster_center[i+k*dim_num] );
    }
    cluster_energy[k] = cluster_energy[k] + point_energy;
  }
  return;
}

void hmeans_02_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_02 applies the H-Means algorithm.
//
//  Discussion:
//
//    This is a simple routine to group a set of points into K clusters,
//    each with a center point, in such a way that the total cluster
//    energy is minimized.  The total cluster energy is the sum of the
//    squares of the distances of each point to the center of its cluster.
//
//    The algorithm begins with an initial estimate for the cluster centers:
//
//    1. The points are assigned to the nearest cluster centers.
//
//    2. The iteration exit ( 1 );s if the total energy has not changed
//        significantly, or we have reached the maximum number of iterations.
//
//    3. Each cluster center is replaced by the centroid of the points
//       in the cluster.
//
//    4. Return to step 1.
//
//    The algorithm may fail to find the best solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates
//    of the points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the coordinates of the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number of
//    points assigned to each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of
//    the clusters.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
{
  bool debug = false;
  int i;
  int j;
  int k;
  int k2;
  double point_energy;
  double point_energy_min;
  int swap;
//
//  Data checks.
//
  if ( cluster_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  CLUSTER_NUM < 1.\n");
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  DIM_NUM < 1.\n");
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  POINT_NUM < 1.\n");
    exit ( 1 );
  }

  if ( it_max < 0 )
  {
    fprintf(stderr,"HMEANS_02 - Fatal error!\n  IT_MAX < 0.\n");
    exit ( 1 );
  }
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
  for ( j = 0; j < point_num; j++ )
  {
    if ( cluster[j] < 0 || cluster_num <= cluster[j] )
    {
      point_energy_min = r8_huge ( );
      for ( k = 0; k < cluster_num; k++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
            pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
        }
        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }
      }
    }
  }

  it_num = 0;

  char label[50];

  for ( ; ; )
  {
//
//  Given centers, assign points to nearest center.
//
    i4vec_zero ( cluster_num, cluster_population );
    r8vec_zero ( cluster_num, cluster_energy );

    swap = 0;

    for ( j = 0; j < point_num; j++ )
    {
      point_energy_min = r8_huge ( );
      k = cluster[j];

      for ( k2 = 0; k2 < cluster_num; k2++ )
      {
        point_energy = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          point_energy = point_energy +
            pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k2;
        }
      }
      if ( k != cluster[j] )
      {
        swap = swap + 1;
      }
      k = cluster[j];
      cluster_energy[k] = cluster_energy[k] + point_energy_min;
      cluster_population[k] = cluster_population[k] + 1;
    }

    if ( debug )
    {
//      cout << "  " << setw(3) << it_num
//           << "  " << setw(14) << r8vec_sum ( cluster_num, cluster_energy ) << "\n";
    }

    if ( 0 < it_num )
    {
      if ( swap == 0 )
      {
        break;
      }
    }

    if ( it_max <= it_num )
    {
      break;
    }

	if (iMAP->updateDisplay) {
		iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
		char progress[60];
		sprintf(progress,"Clustering Iteration %i\n",it_num);
		textDisplayUpdate(progress);
	}
    it_num = it_num + 1;

	//
//  Given points in cluster, replace center by centroid.
//
    r8vec_zero ( dim_num * cluster_num, cluster_center );

    for ( j = 0; j < point_num; j++ )
    {
      k = cluster[j];
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num]
          + point[i+j*dim_num];
      }
    }

    for ( k = 0; k < cluster_num; k++ )
    {
      if ( cluster_population[k] != 0 )
      {
        for ( i = 0; i < dim_num; i++ )
        {
          cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] /
            ( double ) ( cluster_population[k] );
        }
      }
      else
      {
        j = i4_uniform ( 0, point_num - 1, seed );
        for ( i = 0; i < dim_num; i++ )
        {
          cluster_center[i+k*dim_num] = point[i+j*dim_num];
        }
      }
    }
  }
//
//  Compute the energy based on the final value of the cluster centers.
//
  r8vec_zero ( cluster_num, cluster_energy );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_energy = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_energy = point_energy +
      pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
    cluster_energy[k] = cluster_energy[k] + point_energy;
  }
  return;
}

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf(stderr,"I4_UNIFORM - Fatal error!\n  Input value of SEED = 0.\n");
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
}




int r8vec_min_index ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, int R8VEC_MIN_INDEX, the index of the smallest entry.
//
{
  int i;
  int min_index;

  if ( n <= 0 )
  {
    min_index = -1;
  }
  else
  {
    min_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[i] < a[min_index] )
      {
        min_index = i;
      }
    }
  }

  return min_index;
}

