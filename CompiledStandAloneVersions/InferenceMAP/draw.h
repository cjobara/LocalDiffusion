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

#ifndef DRAW_H_
#define DRAW_H_

static double maxarg1,maxarg2;
#define FMAX(a,b) ( maxarg1 = (a), maxarg2 = (b) , (maxarg1) > (maxarg2) ? \
(maxarg1) : (maxarg2) )
#define FMIN(a,b) ( maxarg1 = (a), maxarg2 = (b) , (maxarg1) < (maxarg2) ? \
(maxarg1) : (maxarg2) )

//double FMAX(double a, double b) {
//	if (a > b) { return a; }
//	else { return b; }
//}
//
//double FMIN(double a, double b) {
//	if (a < b) { return a; }
//	else { return b; }
//}

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#include <stdio.h>
#include <stddef.h>
#include <iostream>
#include "trees.h"

void updateView();
void updateViewAxis();

void drawFrame();
void drawIntensity();
void drawLocalizations();

void clipPlanes();
void unclipPlanes();

// 2D Trajectories
void drawTrajectories();
void animateTrajectories();

// 3D Trajectories
void drawTrajectories3D();
void animateTrajectories3D();
void drawLocalizations3D();
//void drawRegularMesh3D();
//void drawRegularForces3D();
//void drawArrow3D(float xc, float yc, float zc, float fx, float fy, float fz, float mag);
//void drawCube(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float rgb[], float alpha);
void zSortLocalizations();
void zSortLines();
int compareZ(const void* a, const void* b);

// Square Drawing
void drawSquareMesh();
void drawSquareForces();
void drawSquareLandscape();
void drawSquareRandomizedOptimizationZones();
void drawSquarePosteriorHighlight();

// Voronoi Tessellation Drawing
void drawVoronoiMesh();
void drawVoronoiForces();
void drawVoronoiPosteriorHighlight();
void drawVoronoiLandscape();
void drawVoronoiRandomizedOptimizationZones();

// Quad-Tree Drawing
void drawQuadTreeMesh();
void createQuadTreeOverlay(QuadTree* tree, int *lineProgress, int *quadProgress, int *quadColorProgress);
void drawQuadTreeForces(QuadTree* tree);
void drawQuadtreePosteriorHighlight();
void drawQuadTreeSpotVisualization(QuadTree *tree);
void drawQuadTreeContourPlot();
void drawQuadTreeLabels(QuadTree *tree);
void drawQuadTreeLandscape();
void drawQuadTreeNeighbourConnections(QuadTree *tree);
void drawQuadTreeRandomizedOptimizationZones(QuadTree *tree);
void drawQuadTreeRandomizedOptimizationZones();

void drawVariables();
void billboardBegin();
void billboardEnd();

void drawCircle(float cx, float cy, float r, int num_segments);

#endif /* DRAW_H_ */
