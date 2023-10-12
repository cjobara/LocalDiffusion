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

#ifndef KMEANS_H_
#define KMEANS_H_

#include "inference.h"

#include "globals.h"

extern Globals *iMAP;

double r8_huge();

void kmeans_01_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
	int &it_num, double point[], int cluster[], double cluster_center[],
	int cluster_population[], double cluster_energy[] );
void kmeans_01_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
	int &it_num, double point[], int cluster[], double cluster_center[],
	int cluster_population[], double cluster_energy[] );


void hmeans_01_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] );
void hmeans_01_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] );
void hmeans_02_l1 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[], int *seed );
void hmeans_02_l2 ( int dim_num, int point_num, int cluster_num, int it_max,
  int &it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[], int *seed );

void i4vec_zero ( int n, int a[] );
void r8vec_zero ( int n, double a[] );
double *r8vec_zero_new ( int n );
int *i4vec_negone_new ( int n );
double *cluster_initialize_5 ( int dim_num, int point_num, int cluster_num, double point[], int *seed );
double *cluster_variance_compute ( int dim_num, int point_num, int cluster_num,
	double point[], int cluster[], double cluster_center[] );
int *i4vec_zero_new ( int n );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double r8vec_sum ( int n, double a[] );
int i4_uniform ( int a, int b, int *seed );
int i4_min ( int i1, int i2 );
int i4_max ( int i1, int i2 );
int r4_nint ( float x );
int r8vec_min_index ( int n, double a[] );
#endif /* KMEANS_H_ */
