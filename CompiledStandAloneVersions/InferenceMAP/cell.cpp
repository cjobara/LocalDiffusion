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

#include "cell.h"
#include "mesh.h"
#include "globals.h"
#include "draw.h"
#include "file.h"
#include "optimization.h"
#include "inference.h"

SquareCell::SquareCell(SquareMesh *parent) {
	jeffreysPrior = false;
	smoothingPrior = false;
	variance = 0.0;
	this->parent = parent;
	roActivated = false;
	roIteration = 0;
	roIdentifier = -1;
	identifier = -1;
	count = 0;
	row = -1;
	column = -1;
	tMin = xMin = yMin = -9999.0;
	tMax = xMax = yMax = 9999.0;
	activated = true;
	xPointer = NULL;
	yPointer = NULL;
	iPointer = NULL;
	tPointer = NULL;
	indexPointer = NULL;
	xForce = 0.0;
	yForce = 0.0;
	zForce = 0.0;
	forceMagnitude = 0.0;
	diffusion = 0.0;
	potential = 0.0;
	gradVx = 0.0;
	gradVy = 0.0;
	gradVz = 0.0;
	gradDx = 0.0;
	gradDy = 0.0;
	gradDz = 0.0;
	xNeighbourPos = 0;
	yNeighbourPos = 0;
	xNeighbours = 0;
	yNeighbours = 0;
	perimeter = 0.0;
	area = 0.0;
	volume = 0.0;
	xCentroid = 0.0;
	yCentroid = 0.0;
	zAverage = 0.0;
	averageDx = 0.0;
	averageDy = 0.0;
	averageDt = 0.0;
	depth = 0;
	dxPointer = NULL;
	dyPointer = NULL;
	dtPointer = NULL;
	landscape = 0.0;
}

SquareCell::SquareCell() {
	jeffreysPrior = false;
	smoothingPrior = false;
	variance = 0.0;
	this->parent = NULL;
	roActivated = false;
	roIteration = 0;
	roIdentifier = -1;
	identifier = -1;
	count = 0;
	row = -1;
	column = -1;
	tMin = xMin = yMin = -9999.0;
	tMax = xMax = yMax = 9999.0;
	activated = true;
	xPointer = NULL;
	yPointer = NULL;
	iPointer = NULL;
	tPointer = NULL;
	indexPointer = NULL;
	xForce = 0.0;
	yForce = 0.0;
	zForce = 0.0;
	forceMagnitude = 0.0;
	diffusion = 0.0;
	potential = 0.0;
	gradVx = 0.0;
	gradVy = 0.0;
	gradVz = 0.0;
	gradDx = 0.0;
	gradDy = 0.0;
	gradDz = 0.0;
	xNeighbourPos = 0;
	yNeighbourPos = 0;
	xNeighbours = 0;
	yNeighbours = 0;
	perimeter = 0.0;
	area = 0.0;
	volume = 0.0;
	xCentroid = 0.0;
	yCentroid = 0.0;
	zAverage = 0.0;
	averageDx = 0.0;
	averageDy = 0.0;
	averageDt = 0.0;
	depth = 0;
	dxPointer = NULL;
	dyPointer = NULL;
	dtPointer = NULL;
	landscape = 0.0;
}

void SquareCell::reset() {
	variance = 0.0;
	count = 0;
	activated = true;
	row = -1;
	column = -1;
	xPointer = NULL;
	yPointer = NULL;
	iPointer = NULL;
	tPointer = NULL;
	indexPointer = NULL;
	xForce = 0.0;
	yForce = 0.0;
	zForce = 0.0;
	forceMagnitude = 0.0;
	diffusion = 0.0;
	potential = 0.0;
	xNeighbourPos = 0;
	yNeighbourPos = 0;
	xNeighbours = 0;
	yNeighbours = 0;
	averageDx = 0.0;
	averageDy = 0.0;
	averageDt = 0.0;
}

void SquareCell::initialize() {
	xPointer = new double[count];
	yPointer = new double[count];
	iPointer = new double[count];
	tPointer = new double[count];
	dxPointer = new double[count];
	dyPointer = new double[count];
	dtPointer = new double[count];
	indexPointer = new int[count];
	for (int h = 0; h < count; h++) {
		xPointer[h] = 0.0;
		yPointer[h] = 0.0;
		iPointer[h] = 0.0;
		tPointer[h] = 0.0;
		indexPointer[h] = 0;
		dxPointer[h] = 0.0;
		dyPointer[h] = 0.0;
		dtPointer[h] = 0.0;
	}
}


void SquareCell::setDetection(int cellIndex, int fileIndex) {
	xPointer[cellIndex] = file->xPointer[fileIndex];
	yPointer[cellIndex] = file->yPointer[fileIndex];
	iPointer[cellIndex] = file->iPointer[fileIndex];
	tPointer[cellIndex] = file->tPointer[fileIndex];
	indexPointer[cellIndex] = fileIndex;
}

SquareCell::~SquareCell() {
	delete [] xPointer;
	delete [] yPointer;
	delete [] iPointer;
	delete [] tPointer;
	delete [] indexPointer;
}

double SquareCell::getXCentre() {
	if (file->squareMesh->selectionMode()) {
		return file->squareMesh->selection.xMin+(double)(column+0.5)*parent->getDx();
	} else {
		return file->xMin+(double)(column+0.5)*parent->getDx();
	}
}

double SquareCell::getYCentre() {
	if (file->squareMesh->selectionMode()) {
		return file->squareMesh->selection.yMin+(double)(row+0.5)*parent->getDx();
	} else {
		return file->yMin+(double)(row+0.5)*parent->getDx();
	}
}

double SquareCell::getDx() { return parent->getDx(); }
double SquareCell::getDy() { return parent->getDy(); }

VoronoiCell::VoronoiCell() {
	roActivated = false;
	roIteration = 0;
	roIdentifier = -1;
	jeffreysPrior = false;
	smoothingPrior = false;
	nLeftNeighbours = 0;
	nRightNeighbours = 0;
	nTopNeighbours = 0;
	nBottomNeighbours = 0;
	leftNeighbours = NULL;
	rightNeighbours = NULL;
	topNeighbours = NULL;
	bottomNeighbours = NULL;
	leftNeighbourWeights = NULL;
	rightNeighbourWeights = NULL;
	topNeighbourWeights = NULL;
	bottomNeighbourWeights = NULL;
	xArea = 0.0;
	yArea = 0.0;
	variance = 0.0;
	gradDy = 0.0;
	gradDx = 0.0;
	landscape = 0.0;
	this->parent = NULL;
	identifier = -999;
	cellId = 0;
	count = 0;
	activated = false;
	xPointer = NULL;
	yPointer = NULL;
	iPointer = NULL;
	tPointer = NULL;
	indexPointer = NULL;
	xForce = 0.0;
	yForce = 0.0;
	xCentroid = 0.0;
	yCentroid = 0.0;
	forceMagnitude = 0.0;
	diffusion = 0.0;
	potential = 0.0;
	gradVx = 0.0;
	gradVy = 0.0;
	energy = 0.0;
	xCentre = 0.0;
	yCentre = 0.0;
	xMean = 0.0;
	yMean = 0.0;
	variance = 0.0;
	area = 0.0;
	perimeter = 0.0;
	numberOfNeighbours = 0;
	neighbourIndices = NULL;
	vertices = NULL;
	randomValue = 0.0;
	cornerCell = 0;
	nVertices = 0;
	averageDx = 0.0;
	averageDy = 0.0;
	averageDt = 0.0;
}

VoronoiCell::~VoronoiCell() {
	if (xPointer != NULL) { delete [] xPointer; }
	if (yPointer != NULL) { delete [] yPointer; }
	if (iPointer != NULL) { delete [] iPointer; }
	if (tPointer != NULL) { delete [] tPointer; }
	if (neighbourIndices != NULL) { delete [] neighbourIndices; }
	if (indexPointer != NULL) { delete [] indexPointer; }
	delete parent;
}
