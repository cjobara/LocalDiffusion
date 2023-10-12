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

#ifndef INFERENCE_H_
#define INFERENCE_H_

#include <math.h>

#include "inference.h"
#include "file.h"
#include "globals.h"
#include "optimization.h"
#include "trees.h"

// OpenGL Libraries
#ifdef __APPLE__
	#include <omp.h>
#endif

/********* SQUARE MESH **********/

// D Inference
double dPosteriorSquare(double *D);
double dPosteriorSmoothingSquare(double *D);
double dPosteriorSmoothingSquareSelection(double *D);
double dGradDxSquare(double *D, int x, int y);
double dGradDySquare(double *D, int x, int y);
double dPosteriorSmoothingSquareRandomizedOptimization(double *D);
double dPosteriorSmoothingSquareRandomizedOptimizationSelection(double *D);
double dGradDxSquareRandomizedOptimization(double *D, int x, int y);
double dGradDySquareRandomizedOptimization(double *D, int x, int y);

// DF Inference functions
double dfPosteriorSquare(double *FxFyD);
double dfPosteriorSmoothingSquare(double *DFxFy);
double dfPosteriorSmoothingSquareSelection(double *DFxFy);
double dfGradDxSquare(double *DFxFy, int x, int y);
double dfGradDySquare(double *DFxFy, int x, int y);
double dfFxValueSquare(double *V_ij, int x, int y);
double dfFyValueSquare(double *V_ij, int x, int y);
double squareDifferenceSquare(double *V_ij);

double dfPosteriorSmoothingSquareRandomizedOptimization(double *DFxFy);
double dfPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DFxFy);
double dfGradDxSquareRandomizedOptimization(double *DFxFy, int x, int y);
double dfGradDySquareRandomizedOptimization(double *DFxFy, int x, int y);

// DDr Inference
double ddrPosteriorSquare(double *DrxDryD);
double ddrPosteriorSmoothingSquare(double *DDrxDry);
double ddrPosteriorSmoothingSquareSelection(double *DDrxDry);
double ddrGradDxSquare(double *DDrxDry, int x, int y);
double ddrGradDySquare(double *DDrxDry, int x, int y);

double ddrPosteriorSmoothingSquareRandomizedOptimization(double *DDrxDry);
double ddrPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DDrxDry);
double ddrGradDxSquareRandomizedOptimization(double *DDrxDry, int x, int y);
double ddrGradDySquareRandomizedOptimization(double *DDrxDry, int x, int y);

// DV Inference functions
double dvPosteriorSquare(double *DV);
double dvPosteriorSquareSelection(double *DV);
double dvPosteriorSmoothingSquare(double *DV);
double dvPosteriorSmoothingSquareSelection(double *DV);
double dvGradVxSquare(double *DV, int x, int y);
double dvGradVySquare(double *DV, int x, int y);
double dvGradDxSquare(double *DV, int x, int y);
double dvGradDySquare(double *DV, int x, int y);
double dvPosteriorSquareRandomizedOptimization(double *DV);
double dvPosteriorSmoothingSquareRandomizedOptimization(double *DV);
double dvPosteriorSquareRandomizedOptimizationSelection(double *DV);
double dvPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DV);
double dvGradVxSquareRandomizedOptimization(double *DV, int x, int y);
double dvGradVySquareRandomizedOptimization(double *DV, int x, int y);
double dvGradDxSquareRandomizedOptimization(double *DV, int x, int y);
double dvGradDySquareRandomizedOptimization(double *DV, int x, int y);

// Polynomial Inference Functions
double polynomialPosteriorSquare(double *coeff_f);
double polynomialFxValueSquare(double *coeffff, double x, double y);
double polynomialFyValueSquare(double *coefffff, double x, double y);
double polynomialDValueSquare(double *coeffff, int x, int y);
double polynomialVValueSquare(double *coeff, double x, double y );
double polynomialPosteriorSquareSelection(double *coeff_f);
double polynomialFxValueSquareSelection(double *coeffff, double x, double y);
double polynomialFyValueSquareSelection(double *coefffff, double x, double y);
double polynomialDValueSquareSelection(double *coeffff, int x, int y);
double polynomialVValueSquareSelection(double *coeff, double x, double y );

/********** VORONOI MESH **********/

// D Inference
double dPosteriorVoronoi(double *D);
double dPosteriorSmoothingVoronoiSelection(double *D);
double dPosteriorSmoothingVoronoi(double *D);
double dGradDxVoronoi(double *D, int i);
double dGradDyVoronoi(double *D, int i);

double dPosteriorSmoothingVoronoiRandomizedOptimization(double *D);
double dPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *D);
double dGradDxVoronoiRandomizedOptimization(double *D, int i);
double dGradDyVoronoiRandomizedOptimization(double *D, int i);

// DF inference functions
double dfPosteriorVoronoi(double *coeff_f);
double dfPosteriorSmoothingVoronoi(double *DFxFy);
double dfPosteriorSmoothingVoronoiSelection(double *DFxFy);
double dfGradDxVoronoi(double *DFxFy, int i);
double dfGradDyVoronoi(double *DFxFy, int i);
double dfFxValueVoronoi(double *V_ij, int i);
double dfFyValueVoronoi(double *V_ij, int i);
double squareDifferenceVoronoi(double *V_ij);

double dfPosteriorSmoothingVoronoiRandomizedOptimization(double *DFxFy);
double dfPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DFxFy);
double dfGradDxVoronoiRandomizedOptimization(double *DFxFy, int i);
double dfGradDyVoronoiRandomizedOptimization(double *DFxFy, int i);

// DDr Inference
double ddrPosteriorVoronoi(double *DrxDryD);
double ddrPosteriorSmoothingVoronoi(double *DDrxDry);
double ddrPosteriorSmoothingVoronoiSelection(double *DDrxDry);
double ddrGradDxVoronoi(double *DDrxDry, int i);
double ddrGradDyVoronoi(double *DDrxDry, int i);

double ddrPosteriorSmoothingVoronoiRandomizedOptimization(double *DDrxDry);
double ddrPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DDrxDry);
double ddrGradDxVoronoiRandomizedOptimization(double *DDrxDry, int i);
double ddrGradDyVoronoiRandomizedOptimization(double *DDrxDry, int i);

// DV  inference functions
double dvPosteriorVoronoi(double *DV);
double dvPosteriorSmoothingVoronoi(double *DV);
double dvPosteriorVoronoiSelection(double *DV);
double dvPosteriorSmoothingVoronoiSelection(double *DV);
double dvGradVxVoronoi(double *DV, int i);
double dvGradVyVoronoi(double *DV, int i);
double dvGradDxVoronoi(double *DV, int i);
double dvGradDyVoronoi(double *DV, int i);

// DV randomized optimization functions
double dvPosteriorVoronoiRandomizedOptimization(double *DV);
double dvPosteriorSmoothingVoronoiRandomizedOptimization(double *DV);
double dvPosteriorVoronoiRandomizedOptimizationSelection(double *DV);
double dvPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DV);
double dvGradVxVoronoiRandomizedOptimization(double *DV, int c);
double dvGradVyVoronoiRandomizedOptimization(double *DV, int c);
double dvGradDxVoronoiRandomizedOptimization(double *DV, int c);
double dvGradDyVoronoiRandomizedOptimization(double *DV, int c);

// Polynomial Inference Functions
double polynomialPosteriorVoronoi(double *coeff_f);
double polynomialPosteriorVoronoiSelection(double *coeff_f);
double polynomialFxValueVoronoi(double *coeffff, double x, double y);
double polynomialFyValueVoronoi(double *coefffff, double x, double y);
double polynomialDValueVoronoi(double *coeffff, int i);
double polynomialVValueVoronoi(double *coeff, double x, double y );
double polynomialFxValueVoronoiSelection(double *coeffff, double x, double y);
double polynomialFyValueVoronoiSelection(double *coefffff, double x, double y);
double polynomialDValueVoronoiSelection(double *coeffff, int i);
double polynomialVValueVoronoiSelection(double *coeff, double x, double y );

/********** QUAD TREE MESH **********/

// D Inference
double dPosteriorQuadTree(double *D);
void dPosteriorSmoothingSetupQuadTree(double *D, QuadTree *tree);
double dPosteriorSmoothingQuadTree(double *D);
double dPosteriorSmoothingQuadTreeSelection(double *D);
double dGradDxQuadTree(double *D, QuadTree *tree);
double dGradDyQuadTree(double *D, QuadTree *tree);

void dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *D, QuadTree *tree);
double dPosteriorSmoothingQuadTreeRandomizedOptimization(double *D);
double dPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *D);
double dGradDxQuadTreeRandomizedOptimization(double *D, QuadTree *tree);
double dGradDyQuadTreeRandomizedOptimization(double *D, QuadTree *tree);

// DF direct inference functions
double dfPosteriorQuadTree(double *FxFyD);
void dfPosteriorSmoothingSetupQuadTree(double *DFxFy, QuadTree *tree);
double dfPosteriorSmoothingQuadTree(double *DFxFy);
double dfPosteriorSmoothingQuadTreeSelection(double *DFxFy);
double dfGradDxQuadTree(double *DFxFy, QuadTree *tree);
double dfGradDyQuadTree(double *DFxFy, QuadTree *tree);
double dfFxValueQuadTree(double *V_ij, QuadTree *tree);
double dfFyValueQuadTree(double *V_ij, QuadTree *tree);
double squareDifferenceQuadTree(double *V_ij);
void squareDifferenceQuadTree(double *V_ij, QuadTree *tree, double *squareDifference);
void squareDifferenceSetupQuadTree(double *V_ij,QuadTree *tree);
void squareDifferenceQuadTree(double *V_ij, QuadTree *tree);

void dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree);
double dfPosteriorSmoothingQuadTreeRandomizedOptimization(double *DFxFy);
double dfPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DFxFy);
double dfGradDxQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree);
double dfGradDyQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree);

// DDr Inference
double ddrPosteriorQuadTree(double *DrxDryD);
void ddrPosteriorSmoothingSetupQuadTree(double *DDrxDry, QuadTree *tree);
double ddrPosteriorSmoothingQuadTree(double *DDrxDry);
double ddrPosteriorSmoothingQuadTreeSelection(double *DDrxDry);
double ddrGradDxQuadTree(double *DDrxDry, QuadTree *tree);
double ddrGradDyQuadTree(double *DDrxDry, QuadTree *tree);

void ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree);
double ddrPosteriorSmoothingQuadTreeRandomizedOptimization(double *DDrxDry);
double ddrPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DDrxDry);
double ddrGradDxQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree);
double ddrGradDyQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree);

// DV direct inference functions
void dvPosteriorSetupQuadTree(double *DV, QuadTree *tree);
double dvPosteriorQuadTree(double *DV);
double dvPosteriorQuadTreeSelection(double *DV);
void dvPosteriorSmoothingSetupQuadTree(double *DV, QuadTree *tree);
double dvPosteriorSmoothingQuadTree(double *DV);
double dvPosteriorSmoothingQuadTreeSelection(double *DV);
double dvGradVxQuadTree(double *DV, QuadTree *tree);
double dvGradVyQuadTree(double *DV, QuadTree *tree);
double dvGradDxQuadTree(double *DV, QuadTree *tree);
double dvGradDyQuadTree(double *DV, QuadTree *tree);

void dvPosteriorSetupQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);
double dvPosteriorQuadTreeRandomizedOptimization(double *DV);
double dvPosteriorQuadTreeRandomizedOptimizationSelection(double *DV);

void dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);
double dvPosteriorSmoothingQuadTreeRandomizedOptimization(double *DV);
double dvPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DV);

double dvGradVxQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);
double dvGradVyQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);
double dvGradDxQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);
double dvGradDyQuadTreeRandomizedOptimization(double *DV, QuadTree *tree);

// Polynomial Inference Functions
double polynomialPosteriorQuadTree(double *coeff_f);
double polynomialPosteriorQuadTreeSelection(double *coeff_f);
double polynomialFxValueQuadTree(double *coeffff, double x, double y);
double polynomialFyValueQuadTree(double *coefffff, double x, double y);
double polynomialDValueQuadTree(double *coeffff, QuadTree *tree);
double polynomialEpValueQuadTree(double *coeff, double x, double y );

/********** SINGLE TRAJECTORY **********/
double singleTrajectoryPosterior(double *coeff_f);

/********** CUSTOM SELECTION **********/
double customSelectionDFPosterior(double *DFxFy);
double customSelectionDDrPosterior(double *DFxFy);

/********** MISCELLANEOUS ***********/

void textDisplayUpdate(char *update);

#endif /* INFERENCE_H_ */
