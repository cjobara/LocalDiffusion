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

#include "mesh.h"
#include "cell.h"
#include "globals.h"
#include "draw.h"
#include "file.h"
#include "optimization.h"
#include "inference.h"
#include "clustering.h"

#include <FL/fl_ask.H>

#define GTOL 1.e-16
#define K_R 4.e-7

#define PI 3.1415926535897932384626433832795028841971693993751

extern Globals *iMAP;
extern File *file;

// 2D Square Mesh
SquareMesh::SquareMesh(double dx) {

	this->randomizedOptimizationGui = NULL;

	//this->selection = NULL;
	this->maximumNeighbourDistance = 0.5; // 500 nm
	this->selectionEnable = false;
	this->roEnable = false;
	this->landscapeNormals = NULL;
	this->landscapeTriangles = NULL;
	this->costArray = NULL;
	this->optimizationArrayFloat = NULL;
	this->dvGradxArray = NULL;
	this->dvGradyArray = NULL;
	this->identityArray = NULL;
	this->cost = 0.0;

	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->speedMax = -1000000.0;
	this->speedMin = 1000000.0;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dimensions = 0;
	this->xMax = -1000000.0;
	this->xMin = 1000000.0;
	this->yMax = -1000000.0;
	this->yMin = 1000000.0;
	this->activeDetections = 0;
	this->optimizationMode = 0;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->potentialMin = 1000000.0;
	this->activeCells = 0;
	this->numberOfCoefficients = 0;
	this->dtMean = 0.0;
	this->hessian = NULL;

	// set passed variables
	this->dx = dx;
	this->dy = dx;

	// default values
	this->polynomialOrder = 2;
	this->smoothingPriorEnable = false;
	this->jeffreysPriorEnable = true;
	this->sigma = 0.00; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentZoneX = 0;
	this->currentZoneY = 0;
	this->currentZone2X = 0;
	this->currentZone2Y = 0;
	this->beta = 2.0;
	this->zonalPotentialsCalculated = false;

	// calculate number of rows and columns
	this->xCells = (int)ceil(file->xRange/this->dx);
	this->yCells = (int)ceil(file->yRange/this->dy);

	// allocate for 2d Cells array + fillArray for tracking assignment progress
	int** fillArray = new int*[xCells];
	for (int rr = 0; rr < xCells; rr++) {
		fillArray[rr] = new int[yCells];
	}
	cells = new SquareCell*[xCells];
	for (int a = 0; a < xCells; a++) {
		cells[a] = new SquareCell[yCells];
		for (int b = 0; b < yCells; b++) {
			cells[a][b].parent = this;
			cells[a][b].area = this->dx*this->dy;
			cells[a][b].perimeter = 2.0*this->dx + 2.0*this->dy;
			fillArray[a][b] = 0;
			cells[a][b].count = 0;
		}
	}

	// count points in each of the cells
	for (int u = 0; u < file->localizationCount; u++) {
		const int x_pos = (int)floor( (file->xPointer[u] - file->xMin)/this->dx );
		const int y_pos = (int)floor( (file->yPointer[u] - file->yMin)/this->dy );
		cells[x_pos][y_pos].increment();
	}

	// initialize individual cells based on count
	// find cell with maximal points
	int maxCount = -100000;
	int minCount = 100000;
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			// assign indices to cell
			cells[a][b].setColumn(a);
			cells[a][b].setRow(b);
			// deactivate empty cells
			if (cells[a][b].getCount() == 0) {
				cells[a][b].deactivate();
				cells[a][b].setIdentifier(-999); // may need to change for indexing
			}
			else if (cells[a][b].getCount() < 2) {
				cells[a][b].deactivate();
				cells[a][b].setIdentifier(-999); // may need to change for indexing
				cells[a][b].initialize();
			}
			else {
				cells[a][b].activate();
				cells[a][b].initialize();
				if (maxCount < cells[a][b].getCount()) {
					maxCount = cells[a][b].getCount();
					this->maxCell = &cells[a][b];
				}
				if (minCount > cells[a][b].getCount()) {
					minCount = cells[a][b].getCount();
					this->minCell = &cells[a][b];
				}
			}
			if (cells[a][b].getCount() < file->squareMeshGui->minPointsSlider->value()) {
				cells[a][b].deactivate();
			}
			cells[a][b].xCentroid = 0.0;
			cells[a][b].yCentroid = 0.0;
		}
	}

	// total number of occupied cells
	this->totalVariables = counter;

	// assign individual detections in their respective cells
	for (int u = 0; u < file->localizationCount; u++) {
		const int x_pos = (int)floor( (file->xPointer[u] - file->xMin)/this->dx );
		const int y_pos = (int)floor( (file->yPointer[u] - file->yMin)/this->dy );
		cells[x_pos][y_pos].setDetection(fillArray[x_pos][y_pos],u);
		fillArray[x_pos][y_pos]++;
		cells[x_pos][y_pos].xCentroid += file->xPointer[u];
		cells[x_pos][y_pos].yCentroid += file->yPointer[u];
	}

	// assign deltas for inference calculation
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
				for (int k = 0; k < cells[a][b].getCount(); k++) {
					if (cells[a][b].indexPointer[k] < file->localizationCount-1) {
						if (file->nPointer[cells[a][b].indexPointer[k]] == file->nPointer[cells[a][b].indexPointer[k]+1]) {
							cells[a][b].dxPointer[k] = file->xPointer[cells[a][b].indexPointer[k]+1]-file->xPointer[cells[a][b].indexPointer[k]];
							cells[a][b].dyPointer[k] = file->yPointer[cells[a][b].indexPointer[k]+1]-file->yPointer[cells[a][b].indexPointer[k]];
							cells[a][b].dtPointer[k] = file->tPointer[cells[a][b].indexPointer[k]+1]-file->tPointer[cells[a][b].indexPointer[k]];
						} else {
							cells[a][b].dxPointer[k] = -999;
							cells[a][b].dyPointer[k] = -999;
							cells[a][b].dtPointer[k] = -999;
						}
					}
				}
		}
	}


	// find maxima/minima in each cell and for the entire mesh
	double xMinMesh,xMaxMesh,yMinMesh,yMaxMesh,tMinMesh,tMaxMesh;
	xMaxMesh = -1000000.0;
	xMinMesh = 1000000.0;
	yMaxMesh = -1000000.0;
	yMinMesh = 1000000.0;
	tMaxMesh = -1000000.0;
	tMinMesh = 1000000.0;

	int az = 0;
	double xMinCell,xMaxCell,yMinCell,yMaxCell,tMinCell,tMaxCell;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				az++;
				xMaxCell = -1000000.0;
				xMinCell = 1000000.0;
				yMaxCell = -1000000.0;
				yMinCell = 1000000.0;
				tMaxCell = -1000000.0;
				tMinCell = 1000000.0;
				cells[a][b].averageDx = 0.0;
				cells[a][b].averageDy = 0.0;
				cells[a][b].averageDt = 0.0;
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					// find max/min x,y coordinates;
					if (xMinCell > cells[a][b].getX(i)) { xMinCell = cells[a][b].getX(i); }
					if (xMaxCell < cells[a][b].getX(i)) { xMaxCell = cells[a][b].getX(i); }
					if (yMinCell > cells[a][b].getY(i)) { yMinCell = cells[a][b].getY(i); }
					if (yMaxCell < cells[a][b].getY(i)) { yMaxCell = cells[a][b].getY(i); }
					if (tMinCell > cells[a][b].getT(i)) { tMinCell = cells[a][b].getT(i); }
					if (tMaxCell < cells[a][b].getT(i)) { tMaxCell = cells[a][b].getT(i); }
					if (i < cells[a][b].getCount()-1) {
						cells[a][b].averageDx += fabs(cells[a][b].getX(i+1)-cells[a][b].getX(i));
						cells[a][b].averageDy += fabs(cells[a][b].getY(i+1)-cells[a][b].getY(i));
						cells[a][b].averageDt += fabs(cells[a][b].getT(i+1)-cells[a][b].getT(i));
					}
				}
				cells[a][b].averageDx /= cells[a][b].getCount()-1;
				cells[a][b].averageDy /= cells[a][b].getCount()-1;
				cells[a][b].averageDt /= cells[a][b].getCount()-1;
				cells[a][b].setSpeed( sqrt(cells[a][b].averageDx*cells[a][b].averageDx+cells[a][b].averageDy*cells[a][b].averageDy)/cells[a][b].averageDt );
				if (cells[a][b].getSpeed() > speedMax) { speedMax = cells[a][b].getSpeed(); }
				if (cells[a][b].getSpeed() < speedMin) { speedMin = cells[a][b].getSpeed(); }
				if (xMaxMesh < xMaxCell) { xMaxMesh = xMaxCell; }
				if (xMinMesh > xMinCell) { xMinMesh = xMinCell; }
				if (xMaxMesh < xMaxCell) { xMaxMesh = xMaxCell; }
				if (yMinMesh > yMinCell) { yMinMesh = yMinCell; }
				if (yMaxMesh < yMaxCell) { yMaxMesh = yMaxCell; }
				if (tMinMesh > tMinCell) { tMinMesh = tMinCell; }
				if (tMaxMesh < tMaxCell) { tMaxMesh = tMaxCell; }
				cells[a][b].setXMin(xMinCell);
				cells[a][b].setXMax(xMaxCell);
				cells[a][b].setYMin(yMinCell);
				cells[a][b].setYMax(yMaxCell);
				cells[a][b].setTMin(tMinCell);
				cells[a][b].setTMax(tMaxCell);
			}
			// determine centroids
			if (cells[a][b].getCount() > 0) {
				cells[a][b].xCentroid /= cells[a][b].getCount();
				cells[a][b].yCentroid /= cells[a][b].getCount();
			}
			// determine position variance
			cells[a][b].variance = 0.0;
			for (int l = 0; l < cells[a][b].getCount(); l++) {
				cells[a][b].variance += sqrt( (cells[a][b].xCentroid-cells[a][b].getX(l))*(cells[a][b].xCentroid-cells[a][b].getX(l))+
												(cells[a][b].yCentroid-cells[a][b].getY(l))*(cells[a][b].yCentroid-cells[a][b].getY(l)));
			}
			if (cells[a][b].getCount() > 0) { cells[a][b].variance /= cells[a][b].getCount(); }
		}
	}
	this->setXMin(xMinMesh);
	this->setXMax(xMaxMesh);
	this->setYMin(yMinMesh);
	this->setYMax(yMaxMesh);
	this->setTMin(tMinMesh);
	this->setTMax(tMaxMesh);
	
	for (int rr = 0; rr < xCells; rr++) {
		delete [] fillArray[rr];
	}
	delete [] fillArray;

	file->squareMeshGui->updateVariables(az);
}

SquareMesh::SquareMesh(double dx, SelectionCell selection) {

	this->randomizedOptimizationGui = NULL;

	selectionEnable = true;
	this->selection = selection;
//	this->m_sparse_matrix = 0;

	this->maximumNeighbourDistance = 0.5;
	this->roEnable = false;
	this->landscapeNormals = NULL;
	this->landscapeTriangles = NULL;
	this->costArray = NULL;
	this->optimizationArrayFloat = NULL;
	this->dvGradxArray = NULL;
	this->dvGradyArray = NULL;
	this->identityArray = NULL;
	this->cost = 0.0;

	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->speedMax = -1000000.0;
	this->speedMin = 1000000.0;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dimensions = 0;
	this->xMax = -1000000.0;
	this->xMin = 1000000.0;
	this->yMax = -1000000.0;
	this->yMin = 1000000.0;
	this->activeDetections = 0;
	this->optimizationMode = 0;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->potentialMin = 1000000.0;
	this->activeCells = 0;
	this->numberOfCoefficients = 0;
	this->dtMean = 0.0;
	this->hessian = NULL;

	// set passed variables
	this->dx = dx;
	this->dy = dx;

	// default values
	this->polynomialOrder = 2;
	this->smoothingPriorEnable = false;
	this->sigma = 0.00; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentZoneX = 0;
	this->currentZoneY = 0;
	this->beta = 2.0;
	this->zonalPotentialsCalculated = false;

	// calculate number of rows and columns
	this->xCells = (int)ceil(selection.xRange/this->dx);
	this->yCells = (int)ceil(selection.yRange/this->dy);

	// allocate for 2d Cells array + fillArray for tracking assignment progress
	int** fillArray = new int*[xCells];
	for (int rr = 0; rr < xCells; rr++) {
		fillArray[rr] = new int[yCells];
	}
	cells = new SquareCell*[xCells];
	for (int a = 0; a < xCells; a++) {
		cells[a] = new SquareCell[yCells];
		for (int b = 0; b < yCells; b++) {
			cells[a][b].parent = this;
			cells[a][b].area = this->dx*this->dy;
			cells[a][b].perimeter = 2.0*this->dx + 2.0*this->dy;
			fillArray[a][b] = 0;
			cells[a][b].count = 0;
		}
	}

	// count points in each of the cells
	for (int u = 0; u < selection.count; u++) {
		const int x_pos = (int)floor( (selection.xPointer[u] - selection.xMin)/this->dx );
		const int y_pos = (int)floor( (selection.yPointer[u] - selection.yMin)/this->dy );
		cells[x_pos][y_pos].increment();
	}

	// initialize individual cells based on count
	// find cell with maximal points
	int maxCount = -100000;
	int minCount = 100000;
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			// assign indices to cell
			cells[a][b].setColumn(a);
			cells[a][b].setRow(b);
			// deactivate empty cells
			if (cells[a][b].getCount() == 0) {
				cells[a][b].deactivate();
				cells[a][b].setIdentifier(-999); // may need to change for indexing
			}
			else if (cells[a][b].getCount() < 2) {
				cells[a][b].deactivate();
				cells[a][b].setIdentifier(-999); // may need to change for indexing
				cells[a][b].initialize();
			}
			else {
				cells[a][b].activate();
				cells[a][b].initialize();
				if (maxCount < cells[a][b].getCount()) {
					maxCount = cells[a][b].getCount();
					this->maxCell = &cells[a][b];
				}
				if (minCount > cells[a][b].getCount()) {
					minCount = cells[a][b].getCount();
					this->minCell = &cells[a][b];
				}
			}
			if (cells[a][b].getCount() < file->squareMeshGui->minPointsSlider->value()) {
				cells[a][b].deactivate();
			}
			cells[a][b].xCentroid = 0.0;
			cells[a][b].yCentroid = 0.0;
		}
	}

	// total number of occupied cells
	this->totalVariables = counter;

	// assign individual detections in their respective cells
	for (int u = 0; u < selection.count; u++) {
		const int x_pos = (int)floor( (selection.xPointer[u] - selection.xMin)/this->dx );
		const int y_pos = (int)floor( (selection.yPointer[u] - selection.yMin)/this->dy );
		cells[x_pos][y_pos].setDetection(fillArray[x_pos][y_pos],selection.indexArray[u]);
		fillArray[x_pos][y_pos]++;
		cells[x_pos][y_pos].xCentroid += selection.xPointer[u];
		cells[x_pos][y_pos].yCentroid += selection.yPointer[u];
	}

	// assign deltas for inference calculation
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				for (int k = 0; k < cells[a][b].getCount(); k++) {
					if (cells[a][b].indexPointer[k] < file->localizationCount-1) {
						if (file->nPointer[cells[a][b].indexPointer[k]] == file->nPointer[cells[a][b].indexPointer[k]+1]) {
							cells[a][b].dxPointer[k] = file->xPointer[cells[a][b].indexPointer[k]+1]-file->xPointer[cells[a][b].indexPointer[k]];
							cells[a][b].dyPointer[k] = file->yPointer[cells[a][b].indexPointer[k]+1]-file->yPointer[cells[a][b].indexPointer[k]];
							cells[a][b].dtPointer[k] = file->tPointer[cells[a][b].indexPointer[k]+1]-file->tPointer[cells[a][b].indexPointer[k]];
						} else {
							cells[a][b].dxPointer[k] = -999;
							cells[a][b].dyPointer[k] = -999;
							cells[a][b].dtPointer[k] = -999;
						}
					}
				}
			}
		}
	}


	// find maxima/minima in each cell and for the entire mesh
	double xMinMesh,xMaxMesh,yMinMesh,yMaxMesh,tMinMesh,tMaxMesh;
	xMaxMesh = -1000000.0;
	xMinMesh = 1000000.0;
	yMaxMesh = -1000000.0;
	yMinMesh = 1000000.0;
	tMaxMesh = -1000000.0;
	tMinMesh = 1000000.0;

	int az = 0;
	double xMinCell,xMaxCell,yMinCell,yMaxCell,tMinCell,tMaxCell;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				az++;
				xMaxCell = -1000000.0;
				xMinCell = 1000000.0;
				yMaxCell = -1000000.0;
				yMinCell = 1000000.0;
				tMaxCell = -1000000.0;
				tMinCell = 1000000.0;
				cells[a][b].averageDx = 0.0;
				cells[a][b].averageDy = 0.0;
				cells[a][b].averageDt = 0.0;
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					// find max/min x,y coordinates;
					if (xMinCell > cells[a][b].getX(i)) { xMinCell = cells[a][b].getX(i); }
					if (xMaxCell < cells[a][b].getX(i)) { xMaxCell = cells[a][b].getX(i); }
					if (yMinCell > cells[a][b].getY(i)) { yMinCell = cells[a][b].getY(i); }
					if (yMaxCell < cells[a][b].getY(i)) { yMaxCell = cells[a][b].getY(i); }
					if (tMinCell > cells[a][b].getT(i)) { tMinCell = cells[a][b].getT(i); }
					if (tMaxCell < cells[a][b].getT(i)) { tMaxCell = cells[a][b].getT(i); }
					if (i < cells[a][b].getCount()-1) {
						cells[a][b].averageDx += fabs(cells[a][b].getX(i+1)-cells[a][b].getX(i));
						cells[a][b].averageDy += fabs(cells[a][b].getY(i+1)-cells[a][b].getY(i));
						cells[a][b].averageDt += fabs(cells[a][b].getT(i+1)-cells[a][b].getT(i));
					}
				}
				cells[a][b].averageDx /= cells[a][b].getCount()-1;
				cells[a][b].averageDy /= cells[a][b].getCount()-1;
				cells[a][b].averageDt /= cells[a][b].getCount()-1;
				cells[a][b].setSpeed( sqrt(cells[a][b].averageDx*cells[a][b].averageDx+cells[a][b].averageDy*cells[a][b].averageDy)/cells[a][b].averageDt );
				if (cells[a][b].getSpeed() > speedMax) { speedMax = cells[a][b].getSpeed(); }
				if (cells[a][b].getSpeed() < speedMin) { speedMin = cells[a][b].getSpeed(); }
				if (xMinMesh > xMinCell) { xMinMesh = xMinCell; }
				if (xMaxMesh < xMaxCell) { xMaxMesh = xMaxCell; }
				if (yMinMesh > yMinCell) { yMinMesh = yMinCell; }
				if (yMaxMesh < yMaxCell) { yMaxMesh = yMaxCell; }
				if (tMinMesh > tMinCell) { tMinMesh = tMinCell; }
				if (tMaxMesh < tMaxCell) { tMaxMesh = tMaxCell; }
				cells[a][b].setXMin(xMinCell);
				cells[a][b].setXMax(xMaxCell);
				cells[a][b].setYMin(yMinCell);
				cells[a][b].setYMax(yMaxCell);
				cells[a][b].setTMin(tMinCell);
				cells[a][b].setTMax(tMaxCell);
			}
			// determine centroids
			if (cells[a][b].getCount() > 0) {
				cells[a][b].xCentroid /= cells[a][b].getCount();
				cells[a][b].yCentroid /= cells[a][b].getCount();
			}
			// determine position variance
			cells[a][b].variance = 0.0;
			for (int l = 0; l < cells[a][b].getCount(); l++) {
				cells[a][b].variance += sqrt( (cells[a][b].xCentroid-cells[a][b].getX(l))*(cells[a][b].xCentroid-cells[a][b].getX(l))+
												(cells[a][b].yCentroid-cells[a][b].getY(l))*(cells[a][b].yCentroid-cells[a][b].getY(l)));
			}
			if (cells[a][b].getCount() > 0) { cells[a][b].variance /= cells[a][b].getCount(); }
		}
	}
	
	for (int rr = 0; rr < xCells; rr++) {
		delete [] fillArray[rr];
	}
	delete [] fillArray;

	this->setXMin(xMinMesh);
	this->setXMax(xMaxMesh);
	this->setYMin(yMinMesh);
	this->setYMax(yMaxMesh);
	this->setTMin(tMinMesh);
	this->setTMax(tMaxMesh);
	file->squareMeshGui->updateVariables(az);
}

double SquareMesh::getDx(int x, int y, int i) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return -999; }
	return cells[x][y].dxPointer[i];
}

double SquareMesh::getDy(int x, int y, int i) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return -999; }
	return cells[x][y].dyPointer[i];
}

double SquareMesh::getDt(int x, int y, int i) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return -999; }
	return cells[x][y].dtPointer[i];
}

double SquareMesh::getCentroidX(int x,int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return -999; }
	return cells[x][y].getXCentroid();
}

double SquareMesh::getCentroidY(int x,int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return -999; }
	return cells[x][y].getYCentroid();
}

bool SquareMesh::active(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return false; }
	return cells[x][y].active();
}

bool SquareMesh::roActive(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return false; }
	return cells[x][y].roActive();
}

void SquareMesh::deactivateCell(int x, int y) {
	if (x <= xCells && y <= yCells && x >= 0 && y >= 0) {
		cells[x][y].deactivate();
	}
}

void SquareMesh::activateCell(int x, int y) {
	if (x <= xCells && y <= yCells && x >= 0 && y >= 0) {
		// can't activate cell with no points
		if (cells[x][y].getCount() >= 2) {
			cells[x][y].activate();
		}
	}
}

int SquareMesh::getCellCount(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0; }
	return cells[x][y].getCount();
}

int SquareMesh::getIdentifier(int x, int y) {
	return cells[x][y].getIdentifier();
}

int SquareMesh::getRoIdentifier(int x, int y) {
	return cells[x][y].getRoIdentifier();
}

int SquareMesh::getXNeighbours(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0; }
	return cells[x][y].getXNeighbours();
}

int SquareMesh::getYNeighbours(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0; }
	return cells[x][y].getYNeighbours();
}

double SquareMesh::getForceX(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getForceX();
}

double SquareMesh::getForceY(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getForceY();
}

double SquareMesh::getForce(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return sqrt(cells[x][y].getForceY()*cells[x][y].getForceY()+cells[x][y].getForceX()*cells[x][y].getForceX());
}

double SquareMesh::getDiffusion(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getDiffusion();
}

double SquareMesh::getDiffusionLog(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getDiffusionLog();
}

double SquareMesh::getPotential(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getPotential();
}

double SquareMesh::getXMean(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getXCentroid();
}

double SquareMesh::getYMean(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0.0; }
	return cells[x][y].getYCentroid();
}

int SquareMesh::getXNeighbourPos(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0; }
	return cells[x][y].getXNeighbourPos();
}

int SquareMesh::getYNeighbourPos(int x, int y) {
	if (x >= xCells || y >= yCells || x < 0 || y < 0) { return 0; }
	return cells[x][y].getYNeighbourPos();
}

int SquareMesh::getCount(int x, int y) { return cells[x][y].count; }

SquareCell* SquareMesh::getCell(int x, int y) { return &cells[x][y]; }

SquareCell* SquareMesh::getCell(int id) {
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].identifier == id) {
				return &cells[a][b];
			}
		}
	}
	fprintf(stderr,"Cell ID %i not found\n",id);
	return NULL;
}

SquareCell* SquareMesh::getCellRo(int id) {
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].roIdentifier == id) {
				return &cells[a][b];
			}
		}
	}
	fprintf(stderr,"Cell ID %i not found\n",id);
	return NULL;
}

void SquareMesh::clearArrays() {
	if (optimizationArray != NULL) { delete [] optimizationArray; }
	if (potentialArray != NULL) { delete [] potentialArray; }
}

SquareMesh::~SquareMesh() {
	// delete Cells array
	if (cells != NULL) {
		for (int r = 0; r < xCells ; r++) { delete [] cells[r]; }
		delete [] cells;
	}
	delete [] optimizationArray;
	delete [] potentialArray;
	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}
}

void SquareMesh::infer() {
	this->minPointsPerCell = (int)file->squareMeshGui->minPointsSlider->value();
	iMAP->startTime = clock();
	iMAP->stopCalculation = false;

	iMAP->glWindow->deactivate();
	
	file->optimizationMode = file->squareMeshGui->inferenceModeChoice->value();

	if (randomizedOptimizationGui != NULL) { randomizedOptimizationGui->hide(); }

	file->squareMeshGui->overlayDiffusionButton->value(1);
	file->squareMeshGui->overlayDiffusionButton->do_callback();

	switch(file->optimizationMode) {
		case 0: // D Inference

			file->squareMeshGui->pauseButton->deactivate();

			if (file->squareMeshGui->smoothingPriorButton->value()) {
				file->squareMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationD();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDSmoothing();
					}
					inferDSmoothing();
					postInferDSmoothing();
				}
			} else {
				preInferD();
				inferD();
				postInferD();
			}

			file->squareMeshGui->overlayPointNumberButton->value(0);
			file->squareMeshGui->overlayForceArrowsButton->deactivate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			file->squareMeshGui->overlayForceArrowsButton->value(0);
			file->squareMeshGui->overlayDiffusionButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->do_callback();
			file->squareMeshGui->overlayForceMagnitudeButton->deactivate();
			file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);

			file->squareMeshGui->minPosteriorSlider->activate();
			file->squareMeshGui->maxPosteriorSlider->activate();
			file->squareMeshGui->fPosteriorButton->deactivate();
			file->squareMeshGui->fPosteriorButton->copy_label("Force");
			file->squareMeshGui->dPosterioriButton->activate();
			file->squareMeshGui->vPosteriorButton->deactivate();
			file->squareMeshGui->potentialReferenceButton->deactivate();
			file->squareMeshGui->posterioriSampleNumberSlider->activate();

			zonalPotentialsCalculated = false;

			file->squareMeshGui->overlayPotentialButton->deactivate();
			file->squareMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
			file->squareMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->squareMeshGui->deltaTSlider->deactivate();
			file->squareMeshGui->timeStepsSlider->deactivate();
			file->squareMeshGui->generateTrajectoriesButton->deactivate();

			break;
		case 1: // DF Inference
			// no pausing in DF mode
			file->squareMeshGui->pauseButton->deactivate();

			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);

			file->squareMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->squareMeshGui->deltaTSlider->deactivate();
			file->squareMeshGui->timeStepsSlider->deactivate();
			file->squareMeshGui->generateTrajectoriesButton->deactivate();

			if (file->squareMeshGui->smoothingPriorButton->value()) {
				file->squareMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					if (inferRandomizedOptimizationDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->squareMeshGui->overlayPotentialButton->activate();
						file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->squareMeshGui->numberOfTrajectoriesSlider->activate();
						file->squareMeshGui->deltaTSlider->activate();
						file->squareMeshGui->timeStepsSlider->activate();
						file->squareMeshGui->generateTrajectoriesButton->activate();
					}
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDFSmoothing();
					}
					inferDFSmoothing();
					zonalPotentialsCalculated = false;
					if (iMAP->stopCalculation == false) {
						if (postInferDFSmoothing()) {
							iMAP->updateDisplayButton->value(1);
							iMAP->updateDisplayButton->do_callback();
							inferPotentialsDF();
							zonalPotentialsCalculated = true;
							file->squareMeshGui->overlayPotentialButton->activate();
							file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
							file->squareMeshGui->numberOfTrajectoriesSlider->activate();
							file->squareMeshGui->deltaTSlider->activate();
							file->squareMeshGui->timeStepsSlider->activate();
							file->squareMeshGui->generateTrajectoriesButton->activate();
						}
					}
				}
			} else {
				preInferDF();
				inferDF();
				zonalPotentialsCalculated = false;
				if (iMAP->stopCalculation == false) {
					if (postInferDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->squareMeshGui->overlayPotentialButton->activate();
						file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->squareMeshGui->numberOfTrajectoriesSlider->activate();
						file->squareMeshGui->deltaTSlider->activate();
						file->squareMeshGui->timeStepsSlider->activate();
						file->squareMeshGui->generateTrajectoriesButton->activate();
					}
				}
			}

			file->squareMeshGui->overlayDiffusionButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceMagnitudeButton->activate();
			file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->squareMeshGui->minPosteriorSlider->activate();
			file->squareMeshGui->maxPosteriorSlider->activate();
			file->squareMeshGui->fPosteriorButton->activate();
			file->squareMeshGui->fPosteriorButton->copy_label("Force");
			file->squareMeshGui->dPosterioriButton->activate();
			file->squareMeshGui->vPosteriorButton->deactivate();
			file->squareMeshGui->potentialReferenceButton->deactivate();
			file->squareMeshGui->posterioriSampleNumberSlider->activate();

			break;
		case 2: // DDr Inference
			// no pausing in DDr mode
			file->squareMeshGui->pauseButton->deactivate();

			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);

			if (file->squareMeshGui->smoothingPriorButton->value()) {
				file->squareMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationDDr();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDDrSmoothing();
					}
					inferDDrSmoothing();
					postInferDDrSmoothing();
				}
			} else {
				preInferDDr();
				inferDDr();
				postInferDDr();
			}

			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);
			file->squareMeshGui->overlayDiffusionButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceMagnitudeButton->activate();
			file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->squareMeshGui->minPosteriorSlider->activate();
			file->squareMeshGui->maxPosteriorSlider->activate();
			file->squareMeshGui->fPosteriorButton->activate();
			file->squareMeshGui->fPosteriorButton->copy_label("Drift");
			file->squareMeshGui->dPosterioriButton->activate();
			file->squareMeshGui->vPosteriorButton->deactivate();
			file->squareMeshGui->potentialReferenceButton->deactivate();
			file->squareMeshGui->posterioriSampleNumberSlider->activate();

			zonalPotentialsCalculated = false;

			file->squareMeshGui->overlayPotentialButton->deactivate();
			file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->squareMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->squareMeshGui->deltaTSlider->deactivate();
			file->squareMeshGui->timeStepsSlider->deactivate();
			file->squareMeshGui->generateTrajectoriesButton->deactivate();

			break;
		case 3: // DV Inference
			file->squareMeshGui->pauseButton->activate();

			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);

			if (roEnable) {
				iMAP->updateDisplayButton->value(1);
				iMAP->updateDisplayButton->do_callback();
				iMAP->updateDisplayButton->deactivate();
				if (randomizedOptimizationGui != NULL) {
					randomizedOptimizationGui->hide();
					delete randomizedOptimizationGui;
					randomizedOptimizationGui = NULL;
				}
				randomizedOptimizationGui = new RandomizedOptimizationGui();
				inferRandomizedOptimizationDV();
				iMAP->updateDisplayButton->activate();
				randomizedOptimizationGui->finish();
			}
			else {
				if (iMAP->pauseCalculation == false) {
					preInferDV();
				}
				inferDV();
				postInferDV();
			}

			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);
			file->squareMeshGui->overlayDiffusionButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceMagnitudeButton->activate();
			file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayPotentialButton->activate();
			file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);

			file->squareMeshGui->minPosteriorSlider->activate();
			file->squareMeshGui->maxPosteriorSlider->activate();
			file->squareMeshGui->fPosteriorButton->deactivate();
			file->squareMeshGui->fPosteriorButton->copy_label("Force");
			file->squareMeshGui->dPosterioriButton->activate();
			file->squareMeshGui->vPosteriorButton->activate();
			file->squareMeshGui->potentialReferenceButton->activate();
			file->squareMeshGui->posterioriSampleNumberSlider->activate();

			file->squareMeshGui->numberOfTrajectoriesSlider->activate();
			file->squareMeshGui->deltaTSlider->activate();
			file->squareMeshGui->timeStepsSlider->activate();
			file->squareMeshGui->generateTrajectoriesButton->activate();
			break;
		case 4: // Polynomial Potential (i.e. confined trajectory)
			// update data structures for optimization run
			if (iMAP->pauseCalculation == false) {
				preInferPolynomial();
			}
			// find minimum of function
			inferPolynomial();
			// save forces and diffusion values to activated cells;
			postInferPolynomial();
			file->squareMeshGui->overlayForceArrowsButton->activate();
			file->squareMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceArrowsButton->value(1);
			file->squareMeshGui->overlayDiffusionButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->activate();
			file->squareMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayForceMagnitudeButton->activate();
			file->squareMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->squareMeshGui->overlayPotentialButton->activate();
			file->squareMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->squareMeshGui->minPosteriorSlider->activate();
			file->squareMeshGui->maxPosteriorSlider->activate();
			file->squareMeshGui->fPosteriorButton->activate();
			file->squareMeshGui->dPosterioriButton->activate();
			file->squareMeshGui->vPosteriorButton->deactivate();
			file->squareMeshGui->potentialReferenceButton->deactivate();
			file->squareMeshGui->posterioriSampleNumberSlider->activate();
			file->squareMeshGui->numberOfTrajectoriesSlider->activate();
			file->squareMeshGui->deltaTSlider->activate();
			file->squareMeshGui->timeStepsSlider->activate();
			file->squareMeshGui->generateTrajectoriesButton->activate();

			break;
	}

	iMAP->glWindow->activate();

	iMAP->endTime = clock();
	char label[100];
	const float dif = ((float)iMAP->endTime-(float)iMAP->startTime)/CLOCKS_PER_SEC;
	textDisplayUpdate("\n");
	sprintf(label,"Calculation Time: %.2f [s]\n",dif);
	textDisplayUpdate(label);

	if (iMAP->pauseCalculation == true) {
		file->squareMeshGui->pauseButton->activate();
		file->squareMeshGui->stopButton->activate();
		file->squareMeshGui->inferButton->deactivate();
		file->squareMeshGui->resetButton->deactivate();
	}

	// initialize posteriori selections
	int l = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				if (l==0) {
					currentZoneX = a;
					currentZoneY = b;
				} else if (l == 1) {
					currentZone2X = a;
					currentZone2Y = b;
					break;
				}
				l++;
			}
		}
		if (l==1) { break; }
	}

	Fl::redraw();
}

void SquareMesh::inferPolynomial() {
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorSquareSelection, (dfunc));
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorSquare, (dfunc));
	}
}

void SquareMesh::preInferPolynomial() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	polynomialOrder = file->squareMeshGui->getPolynomialOrder();

	// inference parameters
	numberOfCoefficients = (int) floor((double)(polynomialOrder+1)*(double)(polynomialOrder+2)/2.0-1.0);

	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				totalVariables++;
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)] == file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}
				// for global diffusion calculation
				cells[a][b].setIdentifier(counter);
				cells[a][b].resetPriors();
				counter++;
				activeCells++;
			}
			else {
				cells[a][b].setIdentifier(-999);
			}
		}
	}

	xMean /= (double) activeDetections;
	yMean /= (double) activeDetections;
	dtMean /= (double) activeDetections;
	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = totalVariables + numberOfCoefficients; // length of vector to be passed for optimization

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }

	optimizationArray[2] = K_R/(4.e-9);
	optimizationArray[4] = K_R/(4.e-9);

	// initial values for optimization array
	for (int i = numberOfCoefficients; i < dimensions; i++) { optimizationArray[i] = 1.0/2.0*(D_eff_x + D_eff_y); }

}

void SquareMesh::postInferPolynomial() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (selectionMode()) {
				if (this->active(a,b)) {
					cells[a][b].setForceX(polynomialFxValueSquareSelection(optimizationArray,xMin+(double)(a+0.5)*getDx(),yMin+(double)(b+0.5)*getDy()));
					cells[a][b].setForceY(polynomialFyValueSquareSelection(optimizationArray,xMin+(double)(a+0.5)*getDx(),yMin+(double)(b+0.5)*getDy()));
					cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
					cells[a][b].setDiffusion(polynomialDValueSquareSelection(optimizationArray,a,b));
					cells[a][b].setDiffusionLog(log(polynomialDValueSquareSelection(optimizationArray,a,b)));
					cells[a][b].setPotential(polynomialVValueSquareSelection(optimizationArray,xMin+(double)(a+0.5)*getDx(),yMin+(double)(b+0.5)*getDy()));
					if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
					if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
					if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
					if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
					if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
					if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
					if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
					if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
				}
				else {
					cells[a][b].setForceX(0.0);
					cells[a][b].setForceY(0.0);
					cells[a][b].setForceMagnitude(0.0);
					cells[a][b].setDiffusion(0.0);
					cells[a][b].setDiffusionLog(0.0);
					cells[a][b].setPotential(0.0);
				}
			} else {
				if (this->active(a,b)) {
					cells[a][b].setForceX(polynomialFxValueSquare(optimizationArray,cells[a][b].getXCentre(),cells[a][b].getYCentre()));
					cells[a][b].setForceY(polynomialFyValueSquare(optimizationArray,cells[a][b].getXCentre(),cells[a][b].getYCentre()));
					cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
					cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
					cells[a][b].setDiffusion(polynomialDValueSquare(optimizationArray,a,b));
					cells[a][b].setDiffusionLog(log(polynomialDValueSquare(optimizationArray,a,b)));
					cells[a][b].setPotential(polynomialVValueSquare(optimizationArray,cells[a][b].getXCentre(),cells[a][b].getYCentre()));
					if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
					if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
					if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
					if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
					if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
					if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
					if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
					if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
				}
				else {
					cells[a][b].setForceX(0.0);
					cells[a][b].setForceY(0.0);
					cells[a][b].setForceMagnitude(0.0);
					cells[a][b].setDiffusion(0.0);
					cells[a][b].setDiffusionLog(0.0);
					cells[a][b].setPotential(0.0);
				}
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();
}

void SquareMesh::preInferD() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	const float maxDist = getMaximumNeighbourDistance();
	int neighboursTemp = 0;
	float am1dist,ap1dist,bm1dist,bp1dist;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				if (a > 0) {
					am1dist = sqrt( (file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { am1dist = 0.0; }
				if (a < xCells-1) {
					ap1dist = sqrt( (file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { ap1dist = 0.0; }
				if (b > 0) {
					bm1dist = sqrt( (file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
								 	 (file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bm1dist = 0.0; }
				if (b < yCells-1) {
					bp1dist = sqrt( (file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bp1dist = 0.0; }

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

//	// deactivate regions without neighbours in x and y
//	for (int a = 0; a < xCells; a++) {
//		for (int b = 0; b < yCells; b++) {
//			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b)==0 ) {
//				cells[a][b].setIdentifier(-999);
//				cells[a][b].deactivate();
//			}
//		}
//	}

	// assign 1D identifier for active cells
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
				counter++;
			}
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 1; // D_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
		for (int e = 0; e < dimensions; e++) {
			hessian[i][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }

}

void SquareMesh::preInferDSmoothing() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int neighboursTemp = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}


	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
//				cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
				cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
				counter++;
			}
		}
	}

	// length of optimization array
	dimensions = counter; // D_ij,V_ij

	cost = 0.0;

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// for the kernel
	optimizationArrayFloat = new float[dimensions];
	dvGradxArray = new float[dimensions];
	dvGradyArray = new float[dimensions];
	identityArray = new int[xCells*yCells];

	// initialize potentials
	int i = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ((this->active(a,b))) {
				optimizationArray[i] = 0.5*(D_eff_x + D_eff_y);
//				optimizationArray[2*i+1] = cells[a][b].getPotential();
				i++;
				identityArray[a*yCells+b] = cells[a][b].getIdentifier();
			} else {
				identityArray[a*yCells+b] = -999;
			}
		}
	}

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}
}

void SquareMesh::inferD() {
	int c = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				setCurrentZone(a,b);
				c++;

				const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
				const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

				optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);

				dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSquare, (dfunc));

				if (iMAP->updateDisplay) {
					iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
					char progress[100];
					sprintf(progress,"Square Meshing\n1   (D) Inference\n2\n3   Zone %i / %i\n",c,this->getActiveCells());
					textDisplayUpdate(progress);
					findExtremesD();
				}

//				cells[a][b].setXForce(optimizationArray[0]);
//				cells[a][b].setYForce(optimizationArray[1]);
//				cells[a][b].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
				cells[a][b].setDiffusion(optimizationArray[0]);
				cells[a][b].setDiffusionLog(log(optimizationArray[0]));

				if (iMAP->stopCalculation) { return; }

			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}
}

void SquareMesh::inferDSmoothing() {

	char label[50];
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingSquareSelection, (dfunc));
		sprintf(label,"Cost = %f\n",dPosteriorSmoothingSquareSelection(optimizationArray));
		textDisplayUpdate(label);
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingSquare, (dfunc));
		sprintf(label,"Cost = %f\n",dPosteriorSmoothingSquare(optimizationArray));
		textDisplayUpdate(label);
	}

}

int SquareMesh::postInferDSmoothing() {

	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setDiffusion(optimizationArray[getCell(a,b)->getIdentifier()]);
				cells[a][b].setDiffusionLog(log(optimizationArray[getCell(a,b)->getIdentifier()]));
				cells[a][b].setPotential(optimizationArray[getCell(a,b)->getIdentifier()+1]);
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return 1;

}

int SquareMesh::postInferD() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
//				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
//				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
//				cells[a][b].setXForce(0.0);
//				cells[a][b].setYForce(0.0);
//				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return 1;
}

void SquareMesh::preInferDDr() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	const float maxDist = getMaximumNeighbourDistance();
	int neighboursTemp = 0;
	float am1dist,ap1dist,bm1dist,bp1dist;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				if (a > 0) {
					am1dist = sqrt( (file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { am1dist = 0.0; }
				if (a < xCells-1) {
					ap1dist = sqrt( (file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { ap1dist = 0.0; }
				if (b > 0) {
					bm1dist = sqrt( (file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
								 	 (file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bm1dist = 0.0; }
				if (b < yCells-1) {
					bp1dist = sqrt( (file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bp1dist = 0.0; }

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b)==0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}

	// assign 1D identifier for active cells
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
				counter++;
			}
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
		for (int e = 0; e < dimensions; e++) {
			hessian[i][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }

}

void SquareMesh::inferDDr() {
	int c = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				c++;
				setCurrentZone(a,b);

				const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
				const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

				optimizationArray[0] = 0.0;
				optimizationArray[1] = 0.0;
				optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

				dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSquare, (dfunc));

				if (iMAP->updateDisplay) {
					iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
					char progress[100];
					sprintf(progress,"Square Meshing\n1   (D,Drift) Inference\n2\n3   Zone %i / %i\n",c,this->getActiveCells());
					textDisplayUpdate(progress);
					findExtremesDDr();
				}

				cells[a][b].setForceX(optimizationArray[0]);
				cells[a][b].setForceY(optimizationArray[1]);
				cells[a][b].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
				cells[a][b].setDiffusion(optimizationArray[2]);
				cells[a][b].setDiffusionLog(log(optimizationArray[2]));

				if (iMAP->stopCalculation) { return; }

			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}
}

int SquareMesh::postInferDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusion(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return 1;
}

void SquareMesh::preInferDF() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	const float maxDist = getMaximumNeighbourDistance();
	int neighboursTemp = 0;
	float am1dist,ap1dist,bm1dist,bp1dist;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				if (a > 0) {
					am1dist = sqrt( (file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { am1dist = 0.0; }
				if (a < xCells-1) {
					ap1dist = sqrt( (file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { ap1dist = 0.0; }
				if (b > 0) {
					bm1dist = sqrt( (file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
								 	 (file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bm1dist = 0.0; }
				if (b < yCells-1) {
					bp1dist = sqrt( (file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bp1dist = 0.0; }

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b)==0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}

	// assign 1D identifier for active cells
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
				counter++;
			}
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
		for (int e = 0; e < dimensions; e++) {
			hessian[i][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }

}

void SquareMesh::inferDF() {
	int c = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				c++;
				setCurrentZone(a,b);

				const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
				const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

				optimizationArray[0] = 0.0;
				optimizationArray[1] = 0.0;
				optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

				dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSquare, (dfunc));

				if (iMAP->updateDisplay) {
					iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
					char progress[100];
					sprintf(progress,"Square Meshing\n1   (D,F) Inference\n2\n3   Zone %i / %i\n",c,this->getActiveCells());
					textDisplayUpdate(progress);
					findExtremesDF();
				}

				cells[a][b].setForceX(optimizationArray[0]);
				cells[a][b].setForceY(optimizationArray[1]);
				cells[a][b].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
				cells[a][b].setDiffusion(optimizationArray[2]);
				cells[a][b].setDiffusionLog(log(optimizationArray[2]));

				if (iMAP->stopCalculation) { return; }

			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}
}

void SquareMesh::preInferDFSmoothing() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int neighboursTemp = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
//				cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
				cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				counter++;
			}
		}
	}

	// length of optimization array
	dimensions = 3*counter; // D_ij,V_ij

	cost = 0.0;

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// for the kernel
	optimizationArrayFloat = new float[dimensions];
	dvGradxArray = new float[dimensions];
	dvGradyArray = new float[dimensions];
	identityArray = new int[xCells*yCells];

	// initialize potentials
	int i = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ((this->active(a,b))) {
				optimizationArray[3*i] = 0.5*(D_eff_x + D_eff_y);
				optimizationArray[3*i+1] = cells[a][b].getForceX();
				optimizationArray[3*i+2] = cells[a][b].getForceY();
				i++;
				identityArray[a*yCells+b] = cells[a][b].getIdentifier();
			} else {
				identityArray[a*yCells+b] = -999;
			}
		}
	}

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}
}

void SquareMesh::inferDFSmoothing() {

	char label[100];
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingSquareSelection, (dfunc));
		sprintf(label,"Cost =  %f\n",dfPosteriorSmoothingSquareSelection(optimizationArray));
		textDisplayUpdate(label);
	}
	else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingSquare, (dfunc));
		sprintf(label,"Cost =  %f\n",dfPosteriorSmoothingSquare(optimizationArray));
		textDisplayUpdate(label);
	}

}

int SquareMesh::postInferDFSmoothing() {

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceX(optimizationArray[3*getCell(a,b)->getIdentifier()+1]);
				cells[a][b].setForceY(optimizationArray[3*getCell(a,b)->getIdentifier()+2]);
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				cells[a][b].setDiffusion(optimizationArray[3*getCell(a,b)->getIdentifier()]);
				cells[a][b].setDiffusionLog(log(optimizationArray[3*getCell(a,b)->getIdentifier()]));
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void SquareMesh::preInferDDrSmoothing() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int neighboursTemp = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
//				cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
				cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				counter++;
			}
		}
	}

	// length of optimization array
	dimensions = 3*counter; // D_ij,V_ij

	cost = 0.0;

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// for the kernel
	optimizationArrayFloat = new float[dimensions];
	dvGradxArray = new float[dimensions];
	dvGradyArray = new float[dimensions];
	identityArray = new int[xCells*yCells];

	// initialize potentials
	int i = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ((this->active(a,b))) {
				optimizationArray[3*i] = 0.5*(D_eff_x + D_eff_y);
				optimizationArray[3*i+1] = cells[a][b].getForceX();
				optimizationArray[3*i+2] = cells[a][b].getForceY();
				i++;
				identityArray[a*yCells+b] = cells[a][b].getIdentifier();
			} else {
				identityArray[a*yCells+b] = -999;
			}
		}
	}

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}
}

void SquareMesh::inferDDrSmoothing() {

	char label[100];
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingSquareSelection, (dfunc));
		sprintf(label,"Cost =  %f\n",ddrPosteriorSmoothingSquareSelection(optimizationArray));
		textDisplayUpdate(label);
	}
	else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingSquare, (dfunc));
		sprintf(label,"Cost =  %f\n",ddrPosteriorSmoothingSquare(optimizationArray));
		textDisplayUpdate(label);
	}

}

int SquareMesh::postInferDDrSmoothing() {

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceX(optimizationArray[3*getCell(a,b)->getIdentifier()+1]);
				cells[a][b].setForceY(optimizationArray[3*getCell(a,b)->getIdentifier()+2]);
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				cells[a][b].setDiffusion(optimizationArray[3*getCell(a,b)->getIdentifier()]);
				cells[a][b].setDiffusionLog(log(optimizationArray[3*getCell(a,b)->getIdentifier()]));
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return 1;
}

int SquareMesh::postInferDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void SquareMesh::inferPotentialsDF() {
	// do not reinitialize if pause button was pressed
	iMAP->pauseCalculation = false;
//	if (iMAP->stopCalculation == false) {

		int k = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (this->active(a,b)) {
					k++;
				}
			}
		}

		dimensions = totalVariables = k;

		if (potentialArray != NULL) {
			delete [] potentialArray;
			potentialArray = NULL;
		}
		potentialArray = new double[dimensions];

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
			for (int e = 0; e < dimensions; e++) {
				hessian[i][e] = 0.0;
			}
		}

		// initialize potentials
		int i = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ((this->active(a,b))) {
					potentialArray[getIdentifier(a,b)] = -log( (double)cells[a][b].getCount()/(double)maxCell->getCount() );
//					cells[a][b].setPotential(potentialArray[getIdentifier(a,b)]);
					i++;
//					fprintf(stderr,"%i\t[%i][%i]\t%f\t%f\t%f\n",getIdentifier(a,b),a,b,cells[a][b].getForceX(),cells[a][b].getForceY(),cells[a][b].getPotential());
				}
			}
		}
//	}
	iMAP->stopCalculation = false;

//	for (int a = 0; a < xCells; a++) {
//		for (int b = 0; b < yCells; b++) {
//			if ((this->active(a,b))) {
//				fprintf(stderr,"[%i,%i]\t\t%f\t%f\t%f\n",a,b,cells[a][b].getPotential(),cells[a][b].getForceX(),cells[a][b].getForceY());
//			}
//		}
//	}

	// perform optimization of potential array
	iMAP->updatePotentialCalculation = true;
	dfpmin(potentialArray,dimensions,GTOL,iterations,fret,squareDifferenceSquare,dfunc);
	iMAP->updatePotentialCalculation = false;

//	for (int r = 0; r < dimensions; r++) {
//		fprintf(stderr,"%i\t%f\n",r,potentialArray[r]);
//	}

	// set potentials
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	int j = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( (this->active(a,b)) /*&& (this->getXNeighbours(a,b)!=0) && (this->getYNeighbours(a,b)!=0)*/ ) {
				cells[a][b].setPotential(potentialArray[getIdentifier(a,b)]);
				if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
				if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
				j++;
			}
			else {
				cells[a][b].setPotential(1e3);
			}
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	offsetPotentials();

}

void SquareMesh::preInferDV() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int neighboursTemp = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {
				for (int i = 0; i < cells[a][b].getCount(); i++) {
					if (cells[a][b].getIndex(i) < file->localizationCount-1) {
						if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
							xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
							yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
							dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
							activeDetections++;
						}
					}
				}

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active()) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active()) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

				cells[a][b].resetPriors();

			}
			else {
				cells[a][b].setIdentifier(-999);
				cells[a][b].setXNeighbours(0);
				cells[a][b].setYNeighbours(0);
				cells[a][b].deactivate();
			}
		}
	}

	// deactivate regions without neighbours in x and y
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
				cells[a][b].setIdentifier(-999);
				cells[a][b].deactivate();
			}
		}
	}


	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	int counter = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ( cells[a][b].active() ) {
				cells[a][b].setIdentifier(counter);
				cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
				cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
				counter++;
			}
		}
	}

	// length of optimization array
	dimensions = 2*counter; // D_ij,V_ij

	cost = 0.0;

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// for the kernel
	optimizationArrayFloat = new float[dimensions];
	dvGradxArray = new float[dimensions];
	dvGradyArray = new float[dimensions];
	identityArray = new int[xCells*yCells];

	// initialize potentials
	int i = 0;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if ((this->active(a,b))) {
				optimizationArray[2*i] = 0.5*(D_eff_x + D_eff_y);
				optimizationArray[2*i+1] = cells[a][b].getPotential();
				i++;
				identityArray[a*yCells+b] = cells[a][b].getIdentifier();
			} else {
				identityArray[a*yCells+b] = -999;
			}
		}
	}

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}
}

void SquareMesh::inferDV() {

	char label[50];
	if (this->selectionMode()) {
		if (file->squareMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingSquareSelection, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSmoothingSquareSelection(optimizationArray));
			textDisplayUpdate(label);
		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSquareSelection, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSquareSelection(optimizationArray));
			textDisplayUpdate(label);
		}
	}
	else {
		if (file->squareMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingSquare, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSmoothingSquare(optimizationArray));
			textDisplayUpdate(label);

		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSquare, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSquare(optimizationArray));
			textDisplayUpdate(label);
		}
	}

}

int SquareMesh::postInferDV() {

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceX(dvGradVxSquare(optimizationArray,a,b));
				cells[a][b].setForceY(dvGradVySquare(optimizationArray,a,b));
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				cells[a][b].setDiffusion(optimizationArray[2*getCell(a,b)->getIdentifier()]);
				cells[a][b].setDiffusionLog(log(optimizationArray[2*getCell(a,b)->getIdentifier()]));
				cells[a][b].setPotential(optimizationArray[2*getCell(a,b)->getIdentifier()+1]);
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
				if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
				if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return 1;
}

void SquareMesh::inferRandomizedOptimizationDV() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int neighboursTemp = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					for (int i = 0; i < cells[a][b].getCount(); i++) {
						if (cells[a][b].getIndex(i) < file->localizationCount-1) {
							if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
								xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
								yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
								dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
								activeDetections++;
							}
						}
					}

					// count neighbours in x
					neighboursTemp = 0;
					if (a == xCells-1) {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
					}
					else if (a == 0) {
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}
					else {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
					cells[a][b].setXNeighbours(neighboursTemp);

					// count neighbours in y
					neighboursTemp = 0;
					if (b == yCells-1) {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
					}
					else if (b == 0) {
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}
					else {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
					cells[a][b].setYNeighbours(neighboursTemp);

					totalVariables++;

				}
				else {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].setXNeighbours(0);
					cells[a][b].setYNeighbours(0);
					cells[a][b].deactivate();
				}
			}
		}

		// deactivate regions without neighbours in x and y
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].deactivate();
				}
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[2*totalVariables];

		// set default inferred values to cell
		int counter = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( cells[a][b].active() ) {
					cells[a][b].setRoIterations(0);
					cells[a][b].setIdentifier(counter);
					cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
					cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
					intermediateOptimizationArray[2*counter] = cells[a][b].getDiffusion();
					intermediateOptimizationArray[2*counter+1] = cells[a][b].getPotential();
					counter++;
				}
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->squareMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->squareMeshGui->roMaximumIterationsSlider->value();

	int centerId,centerRow,centerCol;
	int loops = 0;

	SquareCell *currentCell;
	roXCoords = roYCoords = NULL;
	const float radius = file->squareMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roXCoords[u]][roYCoords[u]].roDeactivate();
				cells[roXCoords[u]][roYCoords[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = rand() % activeCells;
		centerRow = getCell(centerId)->getRow();
		centerCol = getCell(centerId)->getColumn();

		// count number of cells
		roZones = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) { roZones++; }
				}
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
		dimensions = 2*roZones;
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }
		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roXCoords != NULL) { delete [] roXCoords; delete [] roYCoords; }
		roXCoords = new int[roZones];
		roYCoords = new int[roZones];

		int p = 0;

		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) {
						currentCell = getCell(a,b);
						roXCoords[p] = a;
						roYCoords[p] = b;
						optimizationArray[2*p] = cells[a][b].getDiffusion();
						optimizationArray[2*p+1] = cells[a][b].getPotential();
						cells[a][b].roActivate();
						cells[a][b].setRoIdentifier(p);
						p++;
					}
				}
			}
		}

		if (file->squareMesh->selectionMode()) {
			if (file->squareMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingSquareRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSquareRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
		} else {
			if (file->squareMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingSquareRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSquareRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[2*cells[roXCoords[t]][roYCoords[t]].getIdentifier()] = optimizationArray[2*t];
			intermediateOptimizationArray[2*cells[roXCoords[t]][roYCoords[t]].getIdentifier()+1] = optimizationArray[2*t+1];
		}

		oldCost = cost;

		if (file->squareMesh->selectionMode()) {
			if (file->squareMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingSquareSelection(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorSquareSelection(intermediateOptimizationArray);
			}
		} else {
			if (file->squareMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingSquare(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorSquare(intermediateOptimizationArray);
			}
		}

		for (int t = 0; t < roZones; t++) {
			cells[roXCoords[t]][roYCoords[t]].setForceX(dvGradVxSquareRandomizedOptimization(optimizationArray,roXCoords[t],roYCoords[t]));
			cells[roXCoords[t]][roYCoords[t]].setForceY(dvGradVySquareRandomizedOptimization(optimizationArray,roXCoords[t],roYCoords[t]));
			cells[roXCoords[t]][roYCoords[t]].setForceMagnitude(sqrt(cells[roXCoords[t]][roYCoords[t]].getForceX()*cells[roXCoords[t]][roYCoords[t]].getForceX()+cells[roXCoords[t]][roYCoords[t]].getForceY()*cells[roXCoords[t]][roYCoords[t]].getForceY()));
			cells[roXCoords[t]][roYCoords[t]].setDiffusion(optimizationArray[2*t]);
			cells[roXCoords[t]][roYCoords[t]].setDiffusionLog(log(optimizationArray[2*t]));
			cells[roXCoords[t]][roYCoords[t]].setPotential(optimizationArray[2*t+1]);
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDV();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
				if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
				if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	if (iMAP->pauseCalculation == false) {

		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < roZones; u++) {
		getCellRo(u)->roDeactivate();
		getCellRo(u)->setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();
}

int SquareMesh::inferRandomizedOptimizationDF() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int neighboursTemp = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					for (int i = 0; i < cells[a][b].getCount(); i++) {
						if (cells[a][b].getIndex(i) < file->localizationCount-1) {
							if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
								xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
								yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
								dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
								activeDetections++;
							}
						}
					}

					// count neighbours in x
					neighboursTemp = 0;
					if (a == xCells-1) {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
					}
					else if (a == 0) {
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}
					else {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
					cells[a][b].setXNeighbours(neighboursTemp);

					// count neighbours in y
					neighboursTemp = 0;
					if (b == yCells-1) {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
					}
					else if (b == 0) {
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}
					else {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
					cells[a][b].setYNeighbours(neighboursTemp);

					totalVariables++;

					cells[a][b].resetPriors();


				}
				else {
					cells[a][b].setIdentifier(-999);
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].setXNeighbours(0);
					cells[a][b].setYNeighbours(0);
					cells[a][b].deactivate();
				}
			}
		}

		// deactivate regions without neighbours in x and y
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].setIdentifier(-999);
					cells[a][b].deactivate();
				}
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		// set default inferred values to cell
		int counter = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( cells[a][b].active() ) {
					cells[a][b].setRoIterations(0);
					cells[a][b].setIdentifier(counter);
					cells[a][b].setForceX(0.0);
					cells[a][b].setForceY(0.0);
					cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
					cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
					intermediateOptimizationArray[3*counter] = cells[a][b].getDiffusion();
					intermediateOptimizationArray[3*counter+1] = cells[a][b].getForceX();
					intermediateOptimizationArray[3*counter+2] = cells[a][b].getForceY();
					counter++;
				}
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->squareMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->squareMeshGui->roMaximumIterationsSlider->value();

	int centerId,centerRow,centerCol;
	int loops = 0;

	SquareCell *currentCell;
	roXCoords = roYCoords = NULL;
	const float radius = file->squareMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roXCoords[u]][roYCoords[u]].roDeactivate();
				cells[roXCoords[u]][roYCoords[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = rand() % activeCells;
		centerRow = getCell(centerId)->getRow();
		centerCol = getCell(centerId)->getColumn();

		// count number of cells
		roZones = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) { roZones++; }
				}
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		dimensions = 3*roZones;
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}

		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }
		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roXCoords != NULL) { delete [] roXCoords; delete [] roYCoords; }
		roXCoords = new int[roZones];
		roYCoords = new int[roZones];

		int p = 0;

		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) {
						currentCell = getCell(a,b);
						roXCoords[p] = a;
						roYCoords[p] = b;
						optimizationArray[3*p] = cells[a][b].getDiffusion();
						optimizationArray[3*p+1] = cells[a][b].getForceX();
						optimizationArray[3*p+2] = cells[a][b].getForceY();
						cells[a][b].roActivate();
						cells[a][b].setRoIdentifier(p);
						p++;
					}
				}
			}
		}

		if (file->squareMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingSquareRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingSquareRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()] = optimizationArray[3*t];
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->squareMesh->selectionMode()) {
			cost = dfPosteriorSmoothingSquareSelection(intermediateOptimizationArray);
		} else {
			cost = dfPosteriorSmoothingSquare(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roXCoords[t]][roYCoords[t]].setForceX(optimizationArray[3*t+1]);
			cells[roXCoords[t]][roYCoords[t]].setForceY(optimizationArray[3*t+2]);
			cells[roXCoords[t]][roYCoords[t]].setForceMagnitude(sqrt(cells[roXCoords[t]][roYCoords[t]].getForceX()*cells[roXCoords[t]][roYCoords[t]].getForceX()+cells[roXCoords[t]][roYCoords[t]].getForceY()*cells[roXCoords[t]][roYCoords[t]].getForceY()));
			cells[roXCoords[t]][roYCoords[t]].setDiffusion(optimizationArray[3*t]);
			cells[roXCoords[t]][roYCoords[t]].setDiffusionLog(log(optimizationArray[3*t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDF();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

//	if (iMAP->pauseCalculation == false) {

//		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
//	}

	for (int u = 0; u < roZones; u++) {
		getCellRo(u)->roDeactivate();
		getCellRo(u)->setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);

}

void SquareMesh::inferRandomizedOptimizationDDr() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int neighboursTemp = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					for (int i = 0; i < cells[a][b].getCount(); i++) {
						if (cells[a][b].getIndex(i) < file->localizationCount-1) {
							if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
								xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
								yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
								dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
								activeDetections++;
							}
						}
					}

					// count neighbours in x
					neighboursTemp = 0;
					if (a == xCells-1) {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
					}
					else if (a == 0) {
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}
					else {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
					cells[a][b].setXNeighbours(neighboursTemp);

					// count neighbours in y
					neighboursTemp = 0;
					if (b == yCells-1) {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
					}
					else if (b == 0) {
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}
					else {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
					cells[a][b].setYNeighbours(neighboursTemp);

					totalVariables++;

				}
				else {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].setXNeighbours(0);
					cells[a][b].setYNeighbours(0);
					cells[a][b].deactivate();
				}
			}
		}

		// deactivate regions without neighbours in x and y
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].deactivate();
				}
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		// set default inferred values to cell
		int counter = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( cells[a][b].active() ) {
					cells[a][b].setRoIterations(0);
					cells[a][b].setIdentifier(counter);
					cells[a][b].setForceX(0.0);
					cells[a][b].setForceY(0.0);
					cells[a][b].setPotential( -log((double)cells[a][b].getCount()/(double)maxCell->getCount()) );
					cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
					intermediateOptimizationArray[3*counter] = cells[a][b].getDiffusion();
					intermediateOptimizationArray[3*counter+1] = cells[a][b].getForceX();
					intermediateOptimizationArray[3*counter+2] = cells[a][b].getForceY();
					counter++;
				}
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->squareMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->squareMeshGui->roMaximumIterationsSlider->value();

	int centerId,centerRow,centerCol;
	int loops = 0;

	SquareCell *currentCell;
	roXCoords = roYCoords = NULL;
	const float radius = file->squareMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roXCoords[u]][roYCoords[u]].roDeactivate();
				cells[roXCoords[u]][roYCoords[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = rand() % activeCells;
		centerRow = getCell(centerId)->getRow();
		centerCol = getCell(centerId)->getColumn();

		// count number of cells
		roZones = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) { roZones++; }
				}
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
		dimensions = 3*roZones;
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }
		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roXCoords != NULL) { delete [] roXCoords; delete [] roYCoords; }
		roXCoords = new int[roZones];
		roYCoords = new int[roZones];

		int p = 0;

		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) {
						currentCell = getCell(a,b);
						roXCoords[p] = a;
						roYCoords[p] = b;
						optimizationArray[3*p] = cells[a][b].getDiffusion();
						optimizationArray[3*p+1] = cells[a][b].getForceX();
						optimizationArray[3*p+2] = cells[a][b].getForceY();
						cells[a][b].roActivate();
						cells[a][b].setRoIdentifier(p);
						p++;
					}
				}
			}
		}

		if (file->squareMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingSquareRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingSquareRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()] = optimizationArray[3*t];
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*cells[roXCoords[t]][roYCoords[t]].getIdentifier()+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->squareMesh->selectionMode()) {
			cost = ddrPosteriorSmoothingSquareSelection(intermediateOptimizationArray);
		} else {
			cost = ddrPosteriorSmoothingSquare(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roXCoords[t]][roYCoords[t]].setForceX(optimizationArray[3*t+1]);
			cells[roXCoords[t]][roYCoords[t]].setForceY(optimizationArray[3*t+2]);
			cells[roXCoords[t]][roYCoords[t]].setForceMagnitude(sqrt(cells[roXCoords[t]][roYCoords[t]].getForceX()*cells[roXCoords[t]][roYCoords[t]].getForceX()+cells[roXCoords[t]][roYCoords[t]].getForceY()*cells[roXCoords[t]][roYCoords[t]].getForceY()));
			cells[roXCoords[t]][roYCoords[t]].setDiffusion(optimizationArray[3*t]);
			cells[roXCoords[t]][roYCoords[t]].setDiffusionLog(log(optimizationArray[3*t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDDr();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setForceMagnitude(sqrt(cells[a][b].getForceX()*cells[a][b].getForceX()+cells[a][b].getForceY()*cells[a][b].getForceY()));
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	if (iMAP->pauseCalculation == false) {

		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < roZones; u++) {
		getCellRo(u)->roDeactivate();
		getCellRo(u)->setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();
}

void SquareMesh::inferRandomizedOptimizationD() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->squareMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->squareMeshGui->getBeta();
	sigma = file->squareMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int neighboursTemp = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					for (int i = 0; i < cells[a][b].getCount(); i++) {
						if (cells[a][b].getIndex(i) < file->localizationCount-1) {
							if (file->nPointer[this->getCell(a,b)->getIndex(i)] == file->nPointer[this->getCell(a,b)->getIndex(i)+1]) {
								xMean += fabs(file->xPointer[cells[a][b].getIndex(i)+1]-file->xPointer[cells[a][b].getIndex(i)]);
								yMean += fabs(file->yPointer[cells[a][b].getIndex(i)+1]-file->yPointer[cells[a][b].getIndex(i)]);
								dtMean += file->tPointer[cells[a][b].getIndex(i)+1]-file->tPointer[cells[a][b].getIndex(i)];
								activeDetections++;
							}
						}
					}

					// count neighbours in x
					neighboursTemp = 0;
					if (a == xCells-1) {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
					}
					else if (a == 0) {
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}
					else {
						if (cells[a-1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(-1);
						}
						if (cells[a+1][b].active()) {
							neighboursTemp++;
							cells[a][b].setXNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
					cells[a][b].setXNeighbours(neighboursTemp);

					// count neighbours in y
					neighboursTemp = 0;
					if (b == yCells-1) {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
					}
					else if (b == 0) {
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}
					else {
						if (cells[a][b-1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(-1);
						}
						if (cells[a][b+1].active()) {
							neighboursTemp++;
							cells[a][b].setYNeighbourPos(1);
						}
					}

					if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
					cells[a][b].setYNeighbours(neighboursTemp);

					totalVariables++;

				}
				else {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].setXNeighbours(0);
					cells[a][b].setYNeighbours(0);
					cells[a][b].deactivate();
				}
			}
		}

		// deactivate regions without neighbours in x and y
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( getXNeighbours(a,b)==0 && getYNeighbours(a,b) == 0 ) {
					cells[a][b].setRoIdentifier(-999);
					cells[a][b].deactivate();
				}
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[totalVariables];

		// set default inferred values to cell
		int counter = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if ( cells[a][b].active() ) {
					cells[a][b].setRoIterations(0);
					cells[a][b].setIdentifier(counter);
					cells[a][b].setDiffusion(0.5*(D_eff_x + D_eff_y));
					intermediateOptimizationArray[counter] = cells[a][b].getDiffusion();
					counter++;
				}
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->squareMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->squareMeshGui->roMaximumIterationsSlider->value();

	int centerId,centerRow,centerCol;
	int loops = 0;

	SquareCell *currentCell;
	roXCoords = roYCoords = NULL;
	const float radius = file->squareMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roXCoords[u]][roYCoords[u]].roDeactivate();
				cells[roXCoords[u]][roYCoords[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = rand() % activeCells;
		centerRow = getCell(centerId)->getRow();
		centerCol = getCell(centerId)->getColumn();

		// count number of cells
		roZones = 0;
		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) { roZones++; }
				}
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
		dimensions = roZones;
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }
		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roXCoords != NULL) { delete [] roXCoords; delete [] roYCoords; }
		roXCoords = new int[roZones];
		roYCoords = new int[roZones];

		int p = 0;

		for (int a = 0; a < xCells; a++) {
			for (int b = 0; b < yCells; b++) {
				if (cells[a][b].active()) {
					const float dist = sqrt( (cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid)*(cells[centerCol][centerRow].xCentroid-cells[a][b].xCentroid) +
											 (cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid)*(cells[centerCol][centerRow].yCentroid-cells[a][b].yCentroid) );
					if (dist < radius) {
						currentCell = getCell(a,b);
						roXCoords[p] = a;
						roYCoords[p] = b;
						optimizationArray[p] = cells[a][b].getDiffusion();
						cells[a][b].roActivate();
						cells[a][b].setRoIdentifier(p);
						p++;
					}
				}
			}
		}

		if (file->squareMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingSquareRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingSquareRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[cells[roXCoords[t]][roYCoords[t]].getIdentifier()] = optimizationArray[t];
		}

		oldCost = cost;

		if (file->squareMesh->selectionMode()) {
			cost = dPosteriorSmoothingSquareSelection(intermediateOptimizationArray);
		} else {
			cost = dPosteriorSmoothingSquare(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roXCoords[t]][roYCoords[t]].setDiffusion(optimizationArray[t]);
			cells[roXCoords[t]][roYCoords[t]].setDiffusionLog(log(optimizationArray[t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesD();
		Fl::unlock();
		Fl::check();

		loops++;

	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	if (iMAP->pauseCalculation == false) {

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < roZones; u++) {
		getCellRo(u)->roDeactivate();
		getCellRo(u)->setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->squareMeshGui->saveButton->activate();
}

void SquareMesh::offsetPotentials() {
	// set reference potential to zero
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				cells[a][b].setPotential(cells[a][b].getPotential()-potentialMin);
			}
		}
	}
	potentialMax -= potentialMin;
	potentialMin = 0.0;
}

void SquareMesh::exportMesh() {
	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Square Mesh File");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Mesh File\t*.smesh\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".smesh");
		FILE * writeFile;

		writeFile = fopen(nativeFilename,"wb");

		// File Information
		fprintf(writeFile,"FILE INFORMATION\n");
		fprintf(writeFile,"Filename:\t%s\n",file->fileName);
		fprintf(writeFile,"Total Number of Trajectories:\t%i\n",file->numberOfFiles);
		fprintf(writeFile,"Duration [s]:\t%.3f\n",iMAP->endIntervalSlider->value()-iMAP->startIntervalSlider->value());
		fprintf(writeFile,"Acquisition Time [ms]:\t%.1f\n",1000.0*iMAP->exposureTime);
		if (this->selectionMode()) {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",selection.xMax-selection.xMin,selection.yMax-selection.yMin);
		} else {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",this->xMax-this->xMin,this->yMax-this->yMin);
		}

		// Mesh Information
		fprintf(writeFile,"MESH INFORMATION\n");
		fprintf(writeFile,"Meshing: Square\n");
		fprintf(writeFile,"Selection Mode: %i\n",this->selectionMode());
		fprintf(writeFile,"Side Length [nm]:\t%.1f\n",this->dx*1000.0);
		fprintf(writeFile,"Total Cells:\t%i\n",this->getXCells()*this->getYCells());
		fprintf(writeFile,"Total Localizations:\t%i\n",file->localizationCount);
		fprintf(writeFile,"Active Cells:\t%i\n",this->getActiveCells());
		fprintf(writeFile,"Active Localizations:\t%i\n\n",this->getActiveDetections());

		// Inference Parameters
		fprintf(writeFile,"INFERENCE PARAMETERS\n");
		fprintf(writeFile,"Noise Sigma [nm]:\t%.1f\n",this->getSigma()*1000.0);
		fprintf(writeFile,"Maximum Neighbour Distance [nm]:\t%.1f\n",this->getMaximumNeighbourDistance());
		fprintf(writeFile,"Minimum Points per Cell:\t%i\n",this->getMinPointsPerCell());
		fprintf(writeFile,"Optimization Scheme:\t%i\n",this->getOptimizationMode());
		fprintf(writeFile,"Potential Calculated:\t%i\n",this->zonalPotentialsCalculated);
		fprintf(writeFile,"Prior Enabled:\t%i\n",this->smoothingPriorEnabled());
		if (this->optimizationMode == 4) {
			fprintf(writeFile,"Polynomial Order:\t%i\n",this->getPolynomialOrder());
			fprintf(writeFile,"Polynomial Coefficients:\t");
			for (int q = 0; q < this->getCoefficients()-1; q++) {
				fprintf(writeFile,"%f\t",this->getOptimizationArray(q));
			}
			fprintf(writeFile,"%f\n\n",this->getOptimizationArray(getCoefficients()-1));
		} else {
			fprintf(writeFile,"Polynomial Order:\tN/A\n");
			fprintf(writeFile,"Polynomial Coefficients:\tN/A\n\n");
		}

		// Data
		if (file->optimizationMode == 2) { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tDx\t\tDy\t\tV\n"); }
		else { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tFx\t\tFy\t\tV\n"); }
		for (int a = 0; a < getXCells(); a++) {
			for (int b = 0; b < getYCells(); b++) {
				fprintf(writeFile,"%i\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
						getCell(a,b)->active(),
						getCell(a,b)->getXCentre(),
						getCell(a,b)->getYCentre(),
						getCell(a,b)->getCount(),
						getCell(a,b)->getXCentroid(),
						getCell(a,b)->getYCentroid(),
						getCell(a,b)->getVariance(),
						getCell(a,b)->getDiffusion(),
						getCell(a,b)->getForceX(),
						getCell(a,b)->getForceY(),
						getCell(a,b)->getPotential()
						);
			}
		}

		fclose(writeFile);
	}
}

void SquareMesh::findExtremesD() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	iMAP->overlayAdjustment = true;
}

void SquareMesh::findExtremesDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	iMAP->overlayAdjustment = true;
}

void SquareMesh::findExtremesDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusion(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
			}
			else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
			}
		}
	}

	iMAP->overlayAdjustment = true;
}

void SquareMesh::findExtremesDV() {

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (this->active(a,b)) {
				if (forceMax < cells[a][b].getForceMagnitude()) { forceMax = cells[a][b].getForceMagnitude(); }
				if (forceMin > cells[a][b].getForceMagnitude()) { forceMin = cells[a][b].getForceMagnitude(); }
				if (diffusionMax < cells[a][b].getDiffusion()) { diffusionMax = cells[a][b].getDiffusion(); }
				if (diffusionMin > cells[a][b].getDiffusion()) { diffusionMin = cells[a][b].getDiffusion(); }
				if (diffusionLogMax < cells[a][b].getDiffusionLog()) { diffusionLogMax = cells[a][b].getDiffusionLog(); }
				if (diffusionLogMin > cells[a][b].getDiffusionLog()) { diffusionLogMin = cells[a][b].getDiffusionLog(); }
				if (potentialMax < cells[a][b].getPotential()) { potentialMax = cells[a][b].getPotential(); }
				if (potentialMin > cells[a][b].getPotential()) { potentialMin = cells[a][b].getPotential(); }
			} else {
				cells[a][b].setForceX(0.0);
				cells[a][b].setForceY(0.0);
				cells[a][b].setForceMagnitude(0.0);
				cells[a][b].setDiffusion(0.0);
				cells[a][b].setDiffusionLog(0.0);
				cells[a][b].setPotential(1.0e3);
			}
		}
	}

	iMAP->overlayAdjustment = true;
}

void SquareMesh::updateNeighbours() {

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	const float maxDist = file->squareMesh->getMaximumNeighbourDistance();
	int neighboursTemp = 0;
	float am1dist,ap1dist,bm1dist,bp1dist;
	for (int a = 0; a < xCells; a++) {
		for (int b = 0; b < yCells; b++) {
			if (cells[a][b].active()) {

				if (a > 0) {
					am1dist = sqrt( (file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { am1dist = 0.0; }
				if (a < xCells-1) {
					ap1dist = sqrt( (file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { ap1dist = 0.0; }
				if (b > 0) {
					bm1dist = sqrt( (file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
								 	 (file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bm1dist = 0.0; }
				if (b < yCells-1) {
					bp1dist = sqrt( (file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
									 (file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
				} else { bp1dist = 0.0; }

				// count neighbours in x
				neighboursTemp = 0;
				if (a == xCells-1) {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
				}
				else if (a == 0) {
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}
				else {
					if (cells[a-1][b].active() && am1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(-1);
					}
					if (cells[a+1][b].active() && ap1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setXNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setXNeighbourPos(0); }
				cells[a][b].setXNeighbours(neighboursTemp);

				// count neighbours in y
				neighboursTemp = 0;
				if (b == yCells-1) {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
				}
				else if (b == 0) {
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}
				else {
					if (cells[a][b-1].active() && bm1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(-1);
					}
					if (cells[a][b+1].active() && bp1dist < maxDist) {
						neighboursTemp++;
						cells[a][b].setYNeighbourPos(1);
					}
				}

				if (neighboursTemp == 2) { cells[a][b].setYNeighbourPos(0); }
				cells[a][b].setYNeighbours(neighboursTemp);

				totalVariables++;
				activeCells++;

			}
		}
	}

}

/**************** VORONOI MESHING ******************/

VoronoiMesh::VoronoiMesh() {
	// default values
	this->roEnable = false;
	this->roZones = 0;
	this->randomizedOptimizationGui = NULL;
	this->minVariance = 1000000.0;
	this->maxVariance = -1000000.0;
	this->maximumNeighbourDistance = file->voronoiMeshGui->neighbourDistanceSlider->value()/1000.0;
	this->smoothingPriorEnable = false;
	this->jeffreysPriorEnable = true;
	//this->selection = NULL;
	this->selectionEnable = false;
	this->landscapeTriangles = NULL;
	this->landscapeNormals = NULL;
	this->landscapeColors = NULL;
	this->landscapeVertices = 0;
	this->yMax = file->yMax;
	this->yMin = file->yMin;
	this->xMax = file->xMax;
	this->xMin = file->xMin;
	this->xCentreList = NULL;
	this->yCentreList = NULL;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dtMean = 0.0;
	this->activeDetections = 0;
	this->numberOfIterations = 0;
	this->seed = 123456789;
	this->clusterVariances = NULL;
	this->activeCells = 0;
	this->clusterEnergies = NULL;
	this->clusterCentres = NULL;
	this->maxIterations = 50;
	this->numberOfCoefficients = 0;
	this->minCell = NULL;
	this->maxCell = NULL;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	this->clusterPopulations = NULL;
	this->clusterIndices = NULL;
	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->charDistance = 0.0;
	this->optimizationMode = 0;
	this->potentialMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->point = 0;
	this->cells = NULL;
	this->numberOfClusters = 0;
	this->dimensions = 0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->polynomialOrder = 2;
	this->slidingWindowEnable = false;
	this->sigma = 0.0; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentZone = 0;
	this->currentZone2 = 0;
	this->beta = 2.0;
	this->spatialDimensions = 2;
	this->colorArray = NULL;
	this->distanceType = 2;
	this->clusteringMethod = 0; // kmeans_01 default
	this->zonalPotentialsCalculated = false;
	this->hessian = NULL;
	this->linesOverlay = NULL;
	this->polygonsOverlay = NULL;
}

VoronoiMesh::VoronoiMesh(SelectionCell selection) {
	// default values
	this->roEnable = false;
	this->roZones = 0;
	this->randomizedOptimizationGui = NULL;
	this->maximumNeighbourDistance = file->voronoiMeshGui->neighbourDistanceSlider->value()/1000.0;
	this->minVariance = 1000000.0;
	this->maxVariance = -1000000.0;
	this->selection = selection;
	this->selectionEnable = true;
	this->landscapeTriangles = NULL;
	this->landscapeNormals = NULL;
	this->landscapeColors = NULL;
	this->landscapeVertices = 0;
	this->yMax = selection.yMax;
	this->yMin = selection.yMin;
	this->xMax = selection.xMax;
	this->xMin = selection.xMin;
	this->xCentreList = NULL;
	this->yCentreList = NULL;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dtMean = 0.0;
	this->activeDetections = 0;
	this->numberOfIterations = 0;
	this->seed = 123456789;
	this->clusterVariances = NULL;
	this->activeCells = 0;
	this->clusterEnergies = NULL;
	this->clusterCentres = NULL;
	this->maxIterations = 50;
	this->numberOfCoefficients = 0;
	this->minCell = NULL;
	this->maxCell = NULL;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	this->clusterPopulations = NULL;
	this->clusterIndices = NULL;
	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->charDistance = 0.0;
	this->optimizationMode = 0;
	this->potentialMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->point = 0;
	this->cells = NULL;
	this->numberOfClusters = 0;
	this->dimensions = 0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->polynomialOrder = 2;
	this->slidingWindowEnable = false;
	this->sigma = 0.0; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentZone = 0;
	this->currentZone2 = 0;
	this->beta = 2.0;
	this->spatialDimensions = 2;
	this->colorArray = NULL;
	this->distanceType = 2;
	this->clusteringMethod = 0; // kmeans_01 default
	this->zonalPotentialsCalculated = false;
	this->hessian = NULL;
	this->linesOverlay = NULL;
	this->polygonsOverlay = NULL;
}

bool VoronoiMesh::active(int i) {
	if (i >= numberOfClusters || i < 0) { return false; }
	return cells[i].active();
}

void VoronoiMesh::deactivateCellClick(int i) {
	if (i < numberOfClusters && i >= 0) {
		cells[i].deactivate();
		iMAP->activeZones--;
		file->voronoiMeshGui->updateVariables(iMAP->activeZones);
	}
}

void VoronoiMesh::activateCellClick(int i) {
	if (i < numberOfClusters && i >= 0) {
		// can't activate cell with no points
		if (cells[i].getCount() > 2) {
			cells[i].activate();
			iMAP->activeZones++;
			file->voronoiMeshGui->updateVariables(iMAP->activeZones);
		}
	}
}

int VoronoiMesh::getIdentifier(int i) {
	if (i < numberOfClusters && i >= 0) {
		return cells[i].identifier;
	}
	else {
		return NULL;
	}
}

VoronoiCell* VoronoiMesh::getCell(int i) {
	if (i < numberOfClusters && i >= 0) {
		return &cells[i];
	}
	else {
		return NULL;
	}
}

int VoronoiMesh::getCellCount(int xi) {
	if (xi >= 0 && xi < numberOfClusters) {
		return cells[xi].getCount();
	}
	else {
		return 0;
	}
}

int VoronoiMesh::getNumberOfNeighbours(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].numberOfNeighbours;
	}
	else {
		return 0;
	}
}

double VoronoiMesh::getDiffusion(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].diffusion;
	}
	else {
		return 0.0;
	}
}

double VoronoiMesh::getDiffusionLog(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].diffusionLog;
	}
	else {
		return 0.0;
	}
}

double VoronoiMesh::getPotential(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].potential;
	}
	else {
		return 0.0;
	}
}

double VoronoiMesh::getForceX(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].xForce;
	}
	else {
		return 0.0;
	}
}

double VoronoiMesh::getForceY(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].yForce;
	}
	else {
		return 0.0;
	}
}

double VoronoiMesh::getForce(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return sqrt(cells[i].xForce*cells[i].xForce + cells[i].yForce*cells[i].yForce);
	}
	else {
		return 0.0;
	}
}

int VoronoiMesh::getCount(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].getCount();
	}
	else {
		return 0;
	}
}

double VoronoiMesh::getCentroidX(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].getXCentroid();
	}
	else {
		return 0;
	}
}

double VoronoiMesh::getCentroidY(int i) {
	if (i >= 0 && i < numberOfClusters) {
		return cells[i].getYCentroid();
	}
	else {
		return 0;
	}
}

void VoronoiMesh::createClusters() {

	// create point file
	point = new double[spatialDimensions*file->localizationCount];
	int q = 0;
	for (int a = 0; a < file->localizationCount; a++) {
		// 2d file
		if (spatialDimensions == 2) {
			point[q] = file->xPointer[a];
			q++;
			point[q] = file->yPointer[a];
			q++;
		}
	}

	// initialize kmeans clustering
	numberOfClusters = file->voronoiMeshGui->getCells();
	maxIterations = file->voronoiMeshGui->getMaxIterations();
	seed = 123456789;

	clusterIndices = i4vec_negone_new(file->localizationCount);
	clusterEnergies = new double[numberOfClusters];
	clusterPopulations = new int[numberOfClusters];

	clusterCentres = this->constantPointsInitialization();

	setClusteringMethod(file->voronoiMeshGui->clusteringChoice->value());

	// default to display progress in display window
	// clear text buffer
	iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
	iMAP->updateDisplayButton->value(1);
	iMAP->updateDisplayButton->do_callback();

	int misc = 123456789;

	switch(clusteringMethod) {
		case 0:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				kmeans_01_l1(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
					  clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			case 2:
				kmeans_01_l2(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
					  clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			}
			break;
		case 1:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				hmeans_01_l1(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			case 2:
				hmeans_01_l2(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			}
			break;
		case 2:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				hmeans_02_l1(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies,&misc);
				break;
			case 2:
				hmeans_02_l2(spatialDimensions,file->localizationCount,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies,&misc);
				break;
			}
			break;
	}
	clusterVariances = cluster_variance_compute(spatialDimensions,file->localizationCount,numberOfClusters,point,clusterIndices,clusterCentres);

	// assign clustered values
	cells = new VoronoiCell[numberOfClusters];//(IrregularCell*)calloc(numberOfClusters,sizeof(IrregularCell));
	int *progress = new int[numberOfClusters];

	int maxCount = -1000000;
	int maxIndex = -1;
	int minCount = 1000000;
	int minIndex = -1;
	this->minVariance = 1000000.0;
	this->maxVariance = -1000000.0;

	for (int c = 0; c < numberOfClusters; c++) {
		progress[c] = 0;

		cells[c].parent = this;
		cells[c].count = clusterPopulations[c];
		cells[c].energy = clusterEnergies[c];
		cells[c].variance = clusterVariances[c];
		cells[c].xCentre = clusterCentres[2*c];
		cells[c].yCentre = clusterCentres[2*c+1];

		cells[c].xPointer = new double[clusterPopulations[c]];
		cells[c].yPointer = new double[clusterPopulations[c]];
		cells[c].iPointer = new double[clusterPopulations[c]];
		cells[c].tPointer = new double[clusterPopulations[c]];
		cells[c].indexPointer = new int[clusterPopulations[c]];

		cells[c].cornerCell = -1;
		cells[c].numberOfNeighbours = 0;

		cells[c].setCellId(c);

		if (maxCount < cells[c].count) {
			maxCount = cells[c].count;
			maxIndex = c;
		}
		if (minCount > cells[c].count) {
			minCount = cells[c].count;
			minIndex = c;
		}
		cells[c].variance = 0.0;
		cells[c].randomValue = (float)ran1(&iMAP->seed);//((float)rand()/(float)RAND_MAX);
	}
	this->minCell = &cells[minIndex];
	this->maxCell = &cells[maxIndex];

	// determine cells closest to the 4 corners
	double blMinDistance = 10000000.0;
	int blMinIndex = -1;
	double brMinDistance = 10000000.0;
	int brMinIndex = -1;
	double trMinDistance = 10000000.0;
	int trMinIndex = -1;
	double tlMinDistance = 10000000.0;
	int tlMinIndex = -1;

	for (int q = 0; q < numberOfClusters; q++) {
		// bottom left - xMin,yMin
		if ( blMinDistance > sqrt((cells[q].xCentre-file->xMin)*(cells[q].xCentre-file->xMin)+(cells[q].yCentre-file->yMin)*(cells[q].yCentre-file->yMin)) ) {
			blMinDistance = sqrt((cells[q].xCentre-file->xMin)*(cells[q].xCentre-file->xMin)+(cells[q].yCentre-file->yMin)*(cells[q].yCentre-file->yMin));
			blMinIndex = q;
		}
		// bottom right - xMax,yMin
		if ( brMinDistance > sqrt((cells[q].xCentre-file->xMax)*(cells[q].xCentre-file->xMax)+(cells[q].yCentre-file->yMin)*(cells[q].yCentre-file->yMin)) ) {
			brMinDistance = sqrt((cells[q].xCentre-file->xMax)*(cells[q].xCentre-file->xMax)+(cells[q].yCentre-file->yMin)*(cells[q].yCentre-file->yMin));
			brMinIndex = q;
		}
		// top right - xMax,yMax
		if ( trMinDistance > sqrt((cells[q].xCentre-file->xMax)*(cells[q].xCentre-file->xMax)+(cells[q].yCentre-file->yMax)*(cells[q].yCentre-file->yMax)) ) {
			trMinDistance = sqrt((cells[q].xCentre-file->xMax)*(cells[q].xCentre-file->xMax)+(cells[q].yCentre-file->yMax)*(cells[q].yCentre-file->yMax));
			trMinIndex = q;
		}
		// top left - xMin,yMax
		if ( tlMinDistance > sqrt((cells[q].xCentre-file->xMin)*(cells[q].xCentre-file->xMin)+(cells[q].yCentre-file->yMax)*(cells[q].yCentre-file->yMax)) ) {
			tlMinDistance = sqrt((cells[q].xCentre-file->xMin)*(cells[q].xCentre-file->xMin)+(cells[q].yCentre-file->yMax)*(cells[q].yCentre-file->yMax));
			tlMinIndex = q;
		}
	}
	cornerCellIndices[0] = blMinIndex;
	cornerCellIndices[1] = brMinIndex;
	cornerCellIndices[2] = trMinIndex;
	cornerCellIndices[3] = tlMinIndex;

	cells[blMinIndex].cornerCell = 0;
	cells[brMinIndex].cornerCell = 1;
	cells[trMinIndex].cornerCell = 2;
	cells[tlMinIndex].cornerCell = 3;

	// assign localizations to their appropriate cell
	for (int d = 0; d < file->localizationCount; d++) {
		cells[clusterIndices[d]].xPointer[progress[clusterIndices[d]]] = file->xPointer[d];
		cells[clusterIndices[d]].yPointer[progress[clusterIndices[d]]] = file->yPointer[d];
		cells[clusterIndices[d]].iPointer[progress[clusterIndices[d]]] = file->iPointer[d];
		cells[clusterIndices[d]].tPointer[progress[clusterIndices[d]]] = file->tPointer[d];
		cells[clusterIndices[d]].indexPointer[progress[clusterIndices[d]]] = d;

//		if (cells[clusterIndices[d]].getCount() >= 20) {
//			cells[clusterIndices[d]].activate();
//		}

		cells[clusterIndices[d]].xCentroid += file->xPointer[d];
		cells[clusterIndices[d]].yCentroid += file->yPointer[d];

		progress[clusterIndices[d]]++;
	}

	// find barycentre and dx,dy of cluster localizations
	float maxDistance = -100000.0;
	float minDistance = 100000.0;
	for (int a = 0; a < numberOfClusters; a++) {
		cells[a].xCentroid /= cells[a].count;
		cells[a].yCentroid /= cells[a].count;
		cells[a].averageDx = 0.0;
		cells[a].averageDy = 0.0;
		cells[a].averageDt = 0.0;
		for (int b = 0; b < cells[a].count; b++) {
			if (b < cells[a].count-1) {
			cells[a].averageDx = fabs(cells[a].getX(b+1)-cells[a].getX(b));
			cells[a].averageDy = fabs(cells[a].getY(b+1)-cells[a].getY(b));
			cells[a].averageDt = fabs(cells[a].getT(b+1)-cells[a].getT(b));
			}
			cells[a].variance += sqrt( (cells[a].getX(b)-cells[a].xCentroid)*(cells[a].getX(b)-cells[a].xCentroid)+(cells[a].getY(b)-cells[a].yCentroid)*(cells[a].getY(b)-cells[a].yCentroid));
		}
		cells[a].variance /= cells[a].count;
		if (maxVariance < cells[a].variance) { maxVariance = cells[a].variance; }
		if (minVariance > cells[a].variance) { minVariance = cells[a].variance; }

//		for (int l = 0; l < cells[a].nVertices; l++) {
//			cells[a].distanceToCentroid[l] /= maxDistance;
//		}
	}

	// for rendering cluster points
	if (this->colorArray != NULL) { delete [] this->colorArray; }
	this->colorArray = new float[4*file->localizationCount];
	for (int e = 0; e < file->localizationCount; e++) {
		const float *rgb = colormap(file->voronoiMeshGui->getColormap(),cells[clusterIndices[e]].randomValue,file->cMin,file->cMax,false);
		colorArray[4*e] = rgb[0];
		colorArray[4*e+1] = rgb[1];
		colorArray[4*e+2] = rgb[2];
		colorArray[4*e+3] = 0.5;
	}

	// create voronoi diagram
	xCentreList = new float[numberOfClusters];
	yCentreList = new float[numberOfClusters];
	for (int f = 0; f < numberOfClusters; f++) {
		xCentreList[f] = (float)cells[f].xCentre;
		yCentreList[f] = (float)cells[f].yCentre;
	}

	voronoiDiagram.generateVoronoi(xCentreList,yCentreList,numberOfClusters,file->xMin,file->xMax,file->yMin,file->yMax,0.0,false);
	voronoiDiagram.resetIterator();

	float xTotal,yTotal,perimeter,area;
	areaMax = -10000000.0;
	areaMin = 10000000.0;
	perimeterMax = -10000000.0;
	perimeterMin = 10000000.0;
	int m = 0;
	for (int k = 0; k < numberOfClusters; k++) {
		// copy vertices to IrregularCell members
		voronoiDiagram.cells[k].copy(&cells[k]);

		// determine centre-of-mass for each cell
		xTotal = yTotal = 0.0;
		for (int j = 0; j < cells[k].nVertices; j++) {
			xTotal += cells[k].vertices[2*j];
			yTotal += cells[k].vertices[2*j+1];
		}
		cells[k].xMean = xTotal/cells[k].nVertices;
		cells[k].yMean = yTotal/cells[k].nVertices;

		// calculate perimeter and area of each cell
		perimeter = 0.0;
		area = 0.0;
		for (int w = 0; w < cells[k].nVertices; w++) {
			if (w < cells[k].nVertices-1) {
				perimeter += sqrt( (cells[k].vertices[2*w]-cells[k].vertices[2*(w+1)])*(cells[k].vertices[2*w]-cells[k].vertices[2*(w+1)]) +
								   (cells[k].vertices[2*w+1]-cells[k].vertices[2*(w+1)+1])*(cells[k].vertices[2*w+1]-cells[k].vertices[2*(w+1)+1]) );
			} else {
				perimeter += sqrt( (cells[k].vertices[2*w]-cells[k].vertices[0])*(cells[k].vertices[2*w]-cells[k].vertices[0]) +
								   (cells[k].vertices[2*w+1]-cells[k].vertices[1])*(cells[k].vertices[2*w+1]-cells[k].vertices[1]) );
			}
			cells[k].perimeter = perimeter;

			m = cells[k].nVertices-1;
			for (int i=0; i<cells[k].nVertices ; i++){
				area += (cells[k].vertices[2*m] + cells[k].vertices[2*i])*(cells[k].vertices[2*m+1] - cells[k].vertices[2*i+1]);
				m = i;
			}
			area *= -0.5;
			cells[k].area = area;
		}

		if (areaMax < cells[k].area) { areaMax = cells[k].area;	}
		if (areaMin > cells[k].area) { areaMin = cells[k].area; }

		if (perimeterMax < cells[k].perimeter) { perimeterMax = cells[k].perimeter;	}
		if (perimeterMin > cells[k].perimeter) { perimeterMin = cells[k].perimeter; }

		for (int f = 0; f < cells[k].nVertices; f++) {
			cells[k].distanceToCentroid[f] = sqrtf( (cells[k].vertices[2*f]-(float)cells[k].xCentroid)*(cells[k].vertices[2*f]-(float)cells[k].xCentroid) +
												    (cells[k].vertices[2*f+1]-(float)cells[k].yCentroid)*(cells[k].vertices[2*f+1]-(float)cells[k].yCentroid) );

			if (minDistance > cells[k].distanceToCentroid[f]) { minDistance = cells[k].distanceToCentroid[f]; }
			if (maxDistance < cells[k].distanceToCentroid[f]) { maxDistance = cells[k].distanceToCentroid[f]; }
		}
	}


	totalVariables = 0;
	// find neighbour of each voronoi cell
	for (int k = 0; k < numberOfClusters; k++) {
		findNeighbours(k);
		// deactivate any cells with less than 2 points
		if (cells[k].count < 20) { cells[k].deactivate(); }
		else {
			cells[k].activate();
			totalVariables++;
		}
 	}

	charDistance = sqrt(file->xRange*file->yRange/(float)file->voronoiMesh->getNumberOfClusters());

	zonalPotentialsCalculated = false;

	delete [] point;
	delete [] progress;
}

void VoronoiMesh::createClustersSelection() {

	// create point file
	point = new double[spatialDimensions*selection.count];
	int q = 0;
	for (int a = 0; a < selection.count; a++) {
		// 2d file
		if (spatialDimensions == 2) {
			point[q] = selection.xPointer[a];
			q++;
			point[q] = selection.yPointer[a];
			q++;
		}
	}

	// initialize kmeans clustering
	numberOfClusters = file->voronoiMeshGui->getCells();
	maxIterations = file->voronoiMeshGui->getMaxIterations();
	seed = 123456789;

	clusterIndices = i4vec_negone_new(selection.count);
	clusterEnergies = new double[numberOfClusters];
	clusterPopulations = new int[numberOfClusters];

	clusterCentres = this->constantPointsInitializationSelection();

	setClusteringMethod(file->voronoiMeshGui->clusteringChoice->value());

	// clear text buffer
	iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
	iMAP->updateDisplayButton->value(1);
	iMAP->updateDisplayButton->do_callback();

	int misc = 123456789;

	switch(clusteringMethod) {
		case 0:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				kmeans_01_l1(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
					  clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			case 2:
				kmeans_01_l2(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
					  clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			}
			break;
		case 1:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				hmeans_01_l1(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			case 2:
				hmeans_01_l2(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies);
				break;
			}
			break;
		case 2:
			switch (file->voronoiMeshGui->distanceButton->value()) {
			case 1:
				hmeans_02_l1(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies,&misc);
				break;
			case 2:
				hmeans_02_l2(spatialDimensions,selection.count,numberOfClusters,maxIterations,numberOfIterations,point,
						clusterIndices,clusterCentres,clusterPopulations,clusterEnergies,&misc);
				break;
			}
			break;
	}
	clusterVariances = cluster_variance_compute(spatialDimensions,selection.count,numberOfClusters,point,clusterIndices,clusterCentres);

	// assign clustered values
	cells = new VoronoiCell[numberOfClusters];//(IrregularCell*)calloc(numberOfClusters,sizeof(IrregularCell));
	int *progress = new int[numberOfClusters];

	int maxCount = -1000000;
	int maxIndex = -1;
	int minCount = 1000000;
	int minIndex = -1;

	for (int c = 0; c < numberOfClusters; c++) {
		progress[c] = 0;

		cells[c].parent = this;
		cells[c].count = clusterPopulations[c];
		cells[c].energy = clusterEnergies[c];
		cells[c].variance = clusterVariances[c];
		cells[c].xCentre = clusterCentres[2*c];
		cells[c].yCentre = clusterCentres[2*c+1];

		cells[c].xPointer = new double[clusterPopulations[c]];
		cells[c].yPointer = new double[clusterPopulations[c]];
		cells[c].iPointer = new double[clusterPopulations[c]];
		cells[c].tPointer = new double[clusterPopulations[c]];
		cells[c].indexPointer = new int[clusterPopulations[c]];

		cells[c].cornerCell = -1;
		cells[c].numberOfNeighbours = 0;

		cells[c].setCellId(c);

		if (maxCount < cells[c].count) {
			maxCount = cells[c].count;
			maxIndex = c;
		}
		if (minCount > cells[c].count) {
			minCount = cells[c].count;
			minIndex = c;
		}
		cells[c].variance = 0.0;
		cells[c].randomValue = (float)ran1(&iMAP->seed);//((float)rand()/(float)RAND_MAX);
	}
	this->minCell = &cells[minIndex];
	this->maxCell = &cells[maxIndex];

	// determine cells closest to the 4 corners
	double blMinDistance = 10000000.0;
	int blMinIndex = -1;
	double brMinDistance = 10000000.0;
	int brMinIndex = -1;
	double trMinDistance = 10000000.0;
	int trMinIndex = -1;
	double tlMinDistance = 10000000.0;
	int tlMinIndex = -1;

	for (int q = 0; q < numberOfClusters; q++) {
		// bottom left - xMin,yMin
		if ( blMinDistance > sqrt((cells[q].xCentre-selection.xMin)*(cells[q].xCentre-selection.xMin)+(cells[q].yCentre-selection.yMin)*(cells[q].yCentre-selection.yMin)) ) {
			blMinDistance = sqrt((cells[q].xCentre-selection.xMin)*(cells[q].xCentre-selection.xMin)+(cells[q].yCentre-selection.yMin)*(cells[q].yCentre-selection.yMin));
			blMinIndex = q;
		}
		// bottom right - xMax,yMin
		if ( brMinDistance > sqrt((cells[q].xCentre-selection.xMax)*(cells[q].xCentre-selection.xMax)+(cells[q].yCentre-selection.yMin)*(cells[q].yCentre-selection.yMin)) ) {
			brMinDistance = sqrt((cells[q].xCentre-selection.xMax)*(cells[q].xCentre-selection.xMax)+(cells[q].yCentre-selection.yMin)*(cells[q].yCentre-selection.yMin));
			brMinIndex = q;
		}
		// top right - xMax,yMax
		if ( trMinDistance > sqrt((cells[q].xCentre-selection.xMax)*(cells[q].xCentre-selection.xMax)+(cells[q].yCentre-selection.yMax)*(cells[q].yCentre-selection.yMax)) ) {
			trMinDistance = sqrt((cells[q].xCentre-selection.xMax)*(cells[q].xCentre-selection.xMax)+(cells[q].yCentre-selection.yMax)*(cells[q].yCentre-selection.yMax));
			trMinIndex = q;
		}
		// top left - xMin,yMax
		if ( tlMinDistance > sqrt((cells[q].xCentre-selection.xMin)*(cells[q].xCentre-selection.xMin)+(cells[q].yCentre-selection.yMax)*(cells[q].yCentre-selection.yMax)) ) {
			tlMinDistance = sqrt((cells[q].xCentre-selection.xMin)*(cells[q].xCentre-selection.xMin)+(cells[q].yCentre-selection.yMax)*(cells[q].yCentre-selection.yMax));
			tlMinIndex = q;
		}
	}
	cornerCellIndices[0] = blMinIndex;
	cornerCellIndices[1] = brMinIndex;
	cornerCellIndices[2] = trMinIndex;
	cornerCellIndices[3] = tlMinIndex;

	cells[blMinIndex].cornerCell = 0;
	cells[brMinIndex].cornerCell = 1;
	cells[trMinIndex].cornerCell = 2;
	cells[tlMinIndex].cornerCell = 3;

	// assign localizations to their appropriate cell
	for (int d = 0; d < selection.count; d++) {
		cells[clusterIndices[d]].xPointer[progress[clusterIndices[d]]] = selection.xPointer[d];
		cells[clusterIndices[d]].yPointer[progress[clusterIndices[d]]] = selection.yPointer[d];
		cells[clusterIndices[d]].iPointer[progress[clusterIndices[d]]] = selection.iPointer[d];
		cells[clusterIndices[d]].tPointer[progress[clusterIndices[d]]] = selection.tPointer[d];
		cells[clusterIndices[d]].indexPointer[progress[clusterIndices[d]]] = selection.indexArray[d];

//		if (cells[clusterIndices[d]].getCount() >= 20) {
			cells[clusterIndices[d]].activate();
//		}

		cells[clusterIndices[d]].xCentroid += selection.xPointer[d];
		cells[clusterIndices[d]].yCentroid += selection.yPointer[d];

		progress[clusterIndices[d]]++;
	}

	// find barycentre and dx,dy of cluster localizations
	for (int a = 0; a < numberOfClusters; a++) {
		cells[a].xCentroid /= cells[a].count;
		cells[a].yCentroid /= cells[a].count;
		cells[a].averageDx = 0.0;
		cells[a].averageDy = 0.0;
		cells[a].averageDt = 0.0;
		for (int b = 0; b < cells[a].count; b++) {
			if (b < cells[a].count-1) {
			cells[a].averageDx = fabs(cells[a].getX(b+1)-cells[a].getX(b));
			cells[a].averageDy = fabs(cells[a].getY(b+1)-cells[a].getY(b));
			cells[a].averageDt = fabs(cells[a].getT(b+1)-cells[a].getT(b));
			}
			cells[a].variance += sqrt( (cells[a].getX(b)-cells[a].xCentroid)*(cells[a].getX(b)-cells[a].xCentroid)+(cells[a].getY(b)-cells[a].yCentroid)*(cells[a].getY(b)-cells[a].yCentroid));
		}
		cells[a].variance /= cells[a].count;
		if (maxVariance < cells[a].variance) { maxVariance = cells[a].variance; }
		if (minVariance > cells[a].variance) { minVariance = cells[a].variance; }
	}

	// for rendering cluster points
	if (this->colorArray != NULL) { delete [] this->colorArray; }
	this->colorArray = new float[4*selection.count];
	for (int e = 0; e < selection.count; e++) {
		const float *rgb = colormap(file->voronoiMeshGui->getColormap(),cells[clusterIndices[e]].randomValue,file->cMin,file->cMax,false);
		colorArray[4*e] = rgb[0];
		colorArray[4*e+1] = rgb[1];
		colorArray[4*e+2] = rgb[2];
		colorArray[4*e+3] = 0.5;
	}

	// create voronoi diagram
	xCentreList = new float[numberOfClusters];
	yCentreList = new float[numberOfClusters];
	for (int f = 0; f < numberOfClusters; f++) {
		xCentreList[f] = (float)cells[f].xCentre;
		yCentreList[f] = (float)cells[f].yCentre;
	}

	voronoiDiagram.generateVoronoi(xCentreList,yCentreList,numberOfClusters,selection.xMin,selection.xMax,selection.yMin,selection.yMax,0.0,false);
	voronoiDiagram.resetIterator();

	float xTotal,yTotal,perimeter,area;
	areaMax = -10000000.0;
	areaMin = 10000000.0;
	perimeterMax = -10000000.0;
	perimeterMin = 10000000.0;
	int m = 0;
	for (int k = 0; k < numberOfClusters; k++) {
		// copy vertices to IrregularCell members
		voronoiDiagram.cells[k].copy(&cells[k]);

		// determine centre-of-mass for each cell
		xTotal = yTotal = 0.0;
		for (int j = 0; j < cells[k].nVertices; j++) {
			xTotal += cells[k].vertices[2*j];
			yTotal += cells[k].vertices[2*j+1];
		}
		cells[k].xMean = xTotal/cells[k].nVertices;
		cells[k].yMean = yTotal/cells[k].nVertices;

		// calculate perimeter and area of each cell
		perimeter = 0.0;
		area = 0.0;
		for (int w = 0; w < cells[k].nVertices; w++) {
			if (w < cells[k].nVertices-1) {
				perimeter += sqrt( (cells[k].vertices[2*w]-cells[k].vertices[2*(w+1)])*(cells[k].vertices[2*w]-cells[k].vertices[2*(w+1)]) +
								   (cells[k].vertices[2*w+1]-cells[k].vertices[2*(w+1)+1])*(cells[k].vertices[2*w+1]-cells[k].vertices[2*(w+1)+1]) );
			} else {
				perimeter += sqrt( (cells[k].vertices[2*w]-cells[k].vertices[0])*(cells[k].vertices[2*w]-cells[k].vertices[0]) +
								   (cells[k].vertices[2*w+1]-cells[k].vertices[1])*(cells[k].vertices[2*w+1]-cells[k].vertices[1]) );
			}
			cells[k].perimeter = perimeter;

			m = cells[k].nVertices-1;
			for (int i=0; i<cells[k].nVertices ; i++){
				area += (cells[k].vertices[2*m] + cells[k].vertices[2*i])*(cells[k].vertices[2*m+1] - cells[k].vertices[2*i+1]);
				m = i;
			}
			area *= -0.5;
			cells[k].area = area;
		}

//		fprintf(stderr,"%i:\t%f\t%f\n",cells[k].count,cells[k].area,PI*cells[k].variance*cells[k].variance);

		if (areaMax < cells[k].area) { areaMax = cells[k].area;	}
		if (areaMin > cells[k].area) { areaMin = cells[k].area; }

		if (perimeterMax < cells[k].perimeter) { perimeterMax = cells[k].perimeter;	}
		if (perimeterMin > cells[k].perimeter) { perimeterMin = cells[k].perimeter; }

	}

	totalVariables = 0;
	// find neighbour of each voronoi cell
	for (int k = 0; k < numberOfClusters; k++) {
		findNeighbours(k);
		// deactivate any cells with less than 2 points
		if (cells[k].count < 20) { cells[k].deactivate(); }
		else {
			cells[k].activate();
			totalVariables++;
		}
 	}

	charDistance = sqrt(selection.xRange*selection.yRange/(float)file->voronoiMesh->getNumberOfClusters());

	zonalPotentialsCalculated = false;

	delete [] point;
	delete [] progress;
}

void VoronoiMesh::infer() {

	this->minPointsPerCell = (int)file->voronoiMeshGui->minPointsSlider->value();

	iMAP->startTime = clock();
	iMAP->stopCalculation = false;

	iMAP->glWindow->deactivate();

	file->optimizationMode = file->voronoiMeshGui->inferenceModeChoice->value();

	if (randomizedOptimizationGui != NULL) { randomizedOptimizationGui->hide(); }

	file->voronoiMeshGui->overlayDiffusionButton->value(1);
	file->voronoiMeshGui->overlayDiffusionButton->do_callback();

	switch(file->optimizationMode) {
		case 0: // D Inference
			// no pausing in DF mode
			file->voronoiMeshGui->pauseButton->deactivate();

			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				file->voronoiMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationD();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDSmoothing();
					}
					inferDSmoothing();
					postInferDSmoothing();
				}
			} else {
				preInferD();
				inferD();
				postInferD();
			}

			// activate overlay buttons
			file->voronoiMeshGui->overlayPointNumberButton->value(0);
			file->voronoiMeshGui->overlayForceArrowsButton->deactivate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			file->voronoiMeshGui->overlayForceArrowsButton->value(0);
			file->voronoiMeshGui->overlayDiffusionButton->value(1);
			file->voronoiMeshGui->overlayDiffusionButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionButton->do_callback();
			file->voronoiMeshGui->overlayForceMagnitudeButton->deactivate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);

			file->voronoiMeshGui->minPosteriorSlider->activate();
			file->voronoiMeshGui->maxPosteriorSlider->activate();
			file->voronoiMeshGui->fPosteriorButton->deactivate();
			file->voronoiMeshGui->fPosteriorButton->copy_label("Force");
			file->voronoiMeshGui->dPosterioriButton->activate();
			file->voronoiMeshGui->vPosteriorButton->deactivate();
			file->voronoiMeshGui->potentialReferenceButton->deactivate();
			file->voronoiMeshGui->posterioriSampleNumberSlider->activate();

			zonalPotentialsCalculated = false;

			file->voronoiMeshGui->overlayPotentialButton->deactivate();
			file->voronoiMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
			file->voronoiMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->voronoiMeshGui->deltaTSlider->deactivate();
			file->voronoiMeshGui->timeStepsSlider->deactivate();
			file->voronoiMeshGui->generateTrajectoriesButton->deactivate();

			break;
		case 1: // DF Inference
			// no pausing in DF mode
			file->voronoiMeshGui->pauseButton->deactivate();
			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);

			file->voronoiMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->voronoiMeshGui->deltaTSlider->deactivate();
			file->voronoiMeshGui->timeStepsSlider->deactivate();
			file->voronoiMeshGui->generateTrajectoriesButton->deactivate();

			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				file->voronoiMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					if (inferRandomizedOptimizationDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->voronoiMeshGui->overlayPotentialButton->activate();
						file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->voronoiMeshGui->numberOfTrajectoriesSlider->activate();
						file->voronoiMeshGui->deltaTSlider->activate();
						file->voronoiMeshGui->timeStepsSlider->activate();
						file->voronoiMeshGui->generateTrajectoriesButton->activate();
					}
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDFSmoothing();
					}
					inferDFSmoothing();
					zonalPotentialsCalculated = false;
					if (iMAP->stopCalculation == false) {
						if (postInferDFSmoothing()) {
							iMAP->updateDisplayButton->value(1);
							iMAP->updateDisplayButton->do_callback();
							inferPotentialsDF();
							zonalPotentialsCalculated = true;
							file->voronoiMeshGui->overlayPotentialButton->activate();
							file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
							file->voronoiMeshGui->numberOfTrajectoriesSlider->activate();
							file->voronoiMeshGui->deltaTSlider->activate();
							file->voronoiMeshGui->timeStepsSlider->activate();
							file->voronoiMeshGui->generateTrajectoriesButton->activate();
						}
					}
				}
			} else {
				preInferDF();
				inferDF();
				zonalPotentialsCalculated = false;
				if (iMAP->stopCalculation == false) {
					if (postInferDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->voronoiMeshGui->overlayPotentialButton->activate();
						file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->voronoiMeshGui->numberOfTrajectoriesSlider->activate();
						file->voronoiMeshGui->deltaTSlider->activate();
						file->voronoiMeshGui->timeStepsSlider->activate();
						file->voronoiMeshGui->generateTrajectoriesButton->activate();
					}
				}
			}

			// activate overlay buttons
			file->voronoiMeshGui->overlayDiffusionButton->activate();
			file->voronoiMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->activate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->voronoiMeshGui->minPosteriorSlider->activate();
			file->voronoiMeshGui->maxPosteriorSlider->activate();
			file->voronoiMeshGui->fPosteriorButton->activate();
			file->voronoiMeshGui->fPosteriorButton->copy_label("Force");
			file->voronoiMeshGui->dPosterioriButton->activate();
			file->voronoiMeshGui->vPosteriorButton->deactivate();
			file->voronoiMeshGui->potentialReferenceButton->deactivate();
			file->voronoiMeshGui->posterioriSampleNumberSlider->activate();

			break;
		case 2: // DDr Inference
			// no pausing in DF mode
			file->voronoiMeshGui->pauseButton->deactivate();

			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);

			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				file->voronoiMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationDDr();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDDrSmoothing();
					}
					inferDDrSmoothing();
					postInferDDrSmoothing();
				}
			} else {
				preInferDDr();
				inferDDr();
				postInferDDr();
			}

			// activate overlay buttons
			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);
			file->voronoiMeshGui->overlayDiffusionButton->activate();
			file->voronoiMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceMagnitudeButton->activate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->voronoiMeshGui->minPosteriorSlider->activate();
			file->voronoiMeshGui->maxPosteriorSlider->activate();
			file->voronoiMeshGui->fPosteriorButton->activate();
			file->voronoiMeshGui->fPosteriorButton->copy_label("Drift");
			file->voronoiMeshGui->dPosterioriButton->activate();
			file->voronoiMeshGui->vPosteriorButton->deactivate();
			file->voronoiMeshGui->potentialReferenceButton->deactivate();
			file->voronoiMeshGui->posterioriSampleNumberSlider->activate();

			zonalPotentialsCalculated = false;

			file->voronoiMeshGui->overlayPotentialButton->deactivate();
			file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->voronoiMeshGui->deltaTSlider->deactivate();
			file->voronoiMeshGui->timeStepsSlider->deactivate();
			file->voronoiMeshGui->generateTrajectoriesButton->deactivate();

			break;
		case 3: // DV Inference
			file->voronoiMeshGui->pauseButton->activate();

			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);

			if (roEnable) {
				iMAP->updateDisplayButton->value(1);
				iMAP->updateDisplayButton->do_callback();
				iMAP->updateDisplayButton->deactivate();
				if (randomizedOptimizationGui != NULL) {
					randomizedOptimizationGui->hide();
					delete randomizedOptimizationGui;
					randomizedOptimizationGui = NULL;
				}
				randomizedOptimizationGui = new RandomizedOptimizationGui();
				inferRandomizedOptimizationDV();
				iMAP->updateDisplayButton->activate();
				randomizedOptimizationGui->finish();
			}
			else {
				if (iMAP->pauseCalculation == false) {
					preInferDV();
				}
				inferDV();
				postInferDV();
			}

			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);
			file->voronoiMeshGui->overlayDiffusionButton->activate();
			file->voronoiMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceMagnitudeButton->activate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayPotentialButton->activate();
			file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);

			file->voronoiMeshGui->minPosteriorSlider->activate();
			file->voronoiMeshGui->maxPosteriorSlider->activate();
			file->voronoiMeshGui->fPosteriorButton->deactivate();
			file->voronoiMeshGui->fPosteriorButton->copy_label("Force");
			file->voronoiMeshGui->dPosterioriButton->activate();
			file->voronoiMeshGui->vPosteriorButton->activate();
			file->voronoiMeshGui->potentialReferenceButton->activate();
			file->voronoiMeshGui->posterioriSampleNumberSlider->activate();

			file->voronoiMeshGui->numberOfTrajectoriesSlider->activate();
			file->voronoiMeshGui->deltaTSlider->activate();
			file->voronoiMeshGui->timeStepsSlider->activate();
			file->voronoiMeshGui->generateTrajectoriesButton->activate();
			break;
		case 4: // Polynomial Potential
			if (iMAP->pauseCalculation == false) {
				preInferPolynomial();
			}
			inferPolynomial();

			// save forces and diffusion values to activated cells;
			postInferPolynomial();
			// active overlay buttons
			file->voronoiMeshGui->overlayForceArrowsButton->activate();
			file->voronoiMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceArrowsButton->value(1);
			file->voronoiMeshGui->overlayDiffusionButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->activate();
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayForceMagnitudeButton->activate();
			file->voronoiMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->overlayPotentialButton->activate();
			file->voronoiMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->voronoiMeshGui->numberOfTrajectoriesSlider->activate();
			file->voronoiMeshGui->deltaTSlider->activate();
			file->voronoiMeshGui->timeStepsSlider->activate();
			file->voronoiMeshGui->generateTrajectoriesButton->activate();
			break;
	}

	iMAP->endTime = clock();
	char label[100];
	const float dif = ((float)iMAP->endTime-(float)iMAP->startTime)/CLOCKS_PER_SEC;
	textDisplayUpdate("\n");
	sprintf(label,"Calculation Time: %.2f [s]\n",dif);
	textDisplayUpdate(label);

	iMAP->glWindow->activate();

	if (iMAP->pauseCalculation == true) {
		file->voronoiMeshGui->pauseButton->activate();
		file->voronoiMeshGui->stopButton->activate();
		file->voronoiMeshGui->inferButton->deactivate();
		file->voronoiMeshGui->resetButton->deactivate();
	}

	// initialize posteriori selections
	int l = 0;
	for (int a = 0; a < numberOfClusters; a++) {
		if (cells[a].active()) {
			if (l==0) {
				currentZone = a;
			} else if (l == 1) {
				currentZone2 = a;
				break;
			}
			l++;
		}
	}

	Fl::redraw();
}

void VoronoiMesh::inferPolynomial() {
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorVoronoiSelection, (dfunc));
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorVoronoi, (dfunc));
	}
}

void VoronoiMesh::preInferPolynomial() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get polynomial order from slider
	polynomialOrder = file->voronoiMeshGui->getPolynomialOrder();

	// inference parameters
	numberOfCoefficients = (int) floor((double)(polynomialOrder+1)*(double)(polynomialOrder+2)/2.0-1.0);

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)];
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) activeDetections;
	yMean /= (double) activeDetections;
	dtMean /= (double) activeDetections;
	const double D_eff_x = file->averageDx*file->averageDx/dtMean;
	const double D_eff_y = file->averageDy*file->averageDy/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = totalVariables + numberOfCoefficients; // length of vector to be passed for optimization

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete []optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < numberOfCoefficients; i++) { optimizationArray[i] = 0.0; }

	optimizationArray[2] = K_R/(4.e-9);
	optimizationArray[4] = K_R/(4.e-9);

	// initial values for optimization array
	for (int i = numberOfCoefficients; i < dimensions; i++) { optimizationArray[i] = 1.0/2.0*(D_eff_x + D_eff_y); }


}

void VoronoiMesh::postInferPolynomial() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces and diffusion for each of the activated zones
	if (selectionMode()) {
		for (int a = 0; a < numberOfClusters; a++) {
			if (this->active(a)) {
				cells[a].setForceX(polynomialFxValueVoronoiSelection(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				cells[a].setForceY(polynomialFyValueVoronoiSelection(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				cells[a].setForceMagnitude(sqrt(cells[a].getForceX()*cells[a].getForceX()+cells[a].getForceY()*cells[a].getForceY()));
				cells[a].setDiffusion(polynomialDValueVoronoiSelection(optimizationArray,a));
				cells[a].setDiffusionLog(log(polynomialDValueVoronoiSelection(optimizationArray,a)));
				cells[a].setPotential(polynomialVValueVoronoiSelection(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
				if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
				if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
				if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
				if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
				if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
				if (potentialMax < cells[a].getPotential()) { potentialMax = cells[a].getPotential(); }
				if (potentialMin > cells[a].getPotential()) { potentialMin = cells[a].getPotential(); }
			}
			else {
				cells[a].setForceX(0.0);
				cells[a].setForceY(0.0);
				cells[a].setForceMagnitude(0.0);
				cells[a].setDiffusion(0.0);
				cells[a].setDiffusionLog(0.0);
				cells[a].setPotential(0.0);
			}
		}
	} else {
		for (int a = 0; a < numberOfClusters; a++) {
			if (this->active(a)) {
				cells[a].setForceX(polynomialFxValueVoronoi(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				cells[a].setForceY(polynomialFyValueVoronoi(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				cells[a].setForceMagnitude(sqrt(cells[a].getForceX()*cells[a].getForceX()+cells[a].getForceY()*cells[a].getForceY()));
				cells[a].setDiffusion(polynomialDValueVoronoi(optimizationArray,a));
				cells[a].setDiffusionLog(log(polynomialDValueVoronoi(optimizationArray,a)));
				cells[a].setPotential(polynomialVValueVoronoi(optimizationArray,cells[a].getXCentre(),cells[a].getYCentre()));
				if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
				if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
				if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
				if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
				if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
				if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
				if (potentialMax < cells[a].getPotential()) { potentialMax = cells[a].getPotential(); }
				if (potentialMin > cells[a].getPotential()) { potentialMin = cells[a].getPotential(); }
			}
			else {
				cells[a].setForceX(0.0);
				cells[a].setForceY(0.0);
				cells[a].setForceMagnitude(0.0);
				cells[a].setDiffusion(0.0);
				cells[a].setDiffusionLog(0.0);
				cells[a].setPotential(0.0);
			}
		}
	}

	// set reference potential to zero
	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();
}

int VoronoiMesh::getCellSelection(double x, double y) {
	double minDistance = 1000000.0;
	int minIndex = -1;
	for (int r = 0; r < numberOfClusters; r++) {
		const double md = (cells[r].xCentre - x)*(cells[r].xCentre - x) + (cells[r].yCentre - y)*(cells[r].yCentre - y);
		if (minDistance > md) {
			minDistance = md;
			minIndex = r;
		}
	}
	return minIndex;
}

void VoronoiMesh::findNeighbours(int c) {
	float xTemp,yTemp;
	cells[c].numberOfNeighbours = 0;
	int neighboursTemp[24]; // arbitrarily high number
	bool duplicate = false;

	// count number of neighbours
	for (int v = 0; v < cells[c].nVertices; v++) {
		xTemp = cells[c].vertices[2*v];
		yTemp = cells[c].vertices[2*v+1];

		for (int d = 0; d < numberOfClusters; d++) {
			if (d != c) {
				for (int e = 0; e < cells[d].nVertices; e++) {

					if (fabsf(xTemp-cells[d].vertices[2*e]) < 0.000001f && fabsf(yTemp-cells[d].vertices[2*e+1]) < 0.000001f) {
						duplicate = false;
						for (int g = 0; g < cells[c].numberOfNeighbours; g++) {
							if (neighboursTemp[g] == d) { duplicate = true; }
						}
						if (!duplicate) {
							neighboursTemp[cells[c].numberOfNeighbours] = d;
							cells[c].numberOfNeighbours++;
						}
						break;
					}

				}
			}
		}
	}
	cells[c].neighbourIndices = new int[cells[c].numberOfNeighbours];
	cells[c].nRightNeighbours = 0;
	cells[c].nLeftNeighbours = 0;
	cells[c].nTopNeighbours = 0;
	cells[c].nBottomNeighbours = 0;
	// count number of neighbours
	for (int h = 0; h < cells[c].numberOfNeighbours; h++) {
		cells[c].neighbourIndices[h] = neighboursTemp[h];
		const double xDiff = fabs(cells[cells[c].neighbourIndices[h]].xCentre-cells[c].xCentre);
		const double yDiff = fabs(cells[cells[c].neighbourIndices[h]].yCentre-cells[c].yCentre);
		const double hyp2 = xDiff*xDiff+yDiff*yDiff;
		// set centroid distance-based condition for neighbouring zones (zones far away from each other are not neighbours!)
		if (sqrt(hyp2) < this->maximumNeighbourDistance) {
			if (xDiff*xDiff/hyp2 > 0.5) {
				if (cells[cells[c].neighbourIndices[h]].xCentre > cells[c].xCentre) {
					cells[c].nRightNeighbours++;
				} else {
					cells[c].nLeftNeighbours++;
				}
			} else {
				if (cells[cells[c].neighbourIndices[h]].yCentre > cells[c].yCentre) {
					cells[c].nTopNeighbours++;
				} else {
					cells[c].nBottomNeighbours++;
				}
			}
		}
	}
	// assign neighbours
	cells[c].rightNeighbours = new VoronoiCell*[cells[c].nRightNeighbours];
	cells[c].rightNeighbourWeights = new double[cells[c].nRightNeighbours];
	cells[c].leftNeighbours = new VoronoiCell*[cells[c].nLeftNeighbours];
	cells[c].leftNeighbourWeights = new double[cells[c].nLeftNeighbours];
	cells[c].topNeighbours = new VoronoiCell*[cells[c].nTopNeighbours];
	cells[c].topNeighbourWeights = new double[cells[c].nTopNeighbours];
	cells[c].bottomNeighbours = new VoronoiCell*[cells[c].nBottomNeighbours];
	cells[c].bottomNeighbourWeights = new double[cells[c].nBottomNeighbours];
	int rn = 0;
	int ln = 0;
	int tn = 0;
	int bn = 0;
	cells[c].xArea = 0.0;
	cells[c].yArea = 0.0;
	for (int h = 0; h < cells[c].numberOfNeighbours; h++) {
		const double xDiff = fabs(cells[cells[c].neighbourIndices[h]].xCentre-cells[c].xCentre);
		const double yDiff = fabs(cells[cells[c].neighbourIndices[h]].yCentre-cells[c].yCentre);
		const double hyp2 = xDiff*xDiff+yDiff*yDiff;
		// set centroid distance-based condition for neighbouring zones (zones far away from each other are not neighbours!)
		if (sqrt(hyp2) < this->maximumNeighbourDistance) {
			if (xDiff*xDiff/hyp2 > 0.5) {
				if (cells[cells[c].neighbourIndices[h]].xCentre > cells[c].xCentre) {
					cells[c].rightNeighbours[rn] = &cells[cells[c].neighbourIndices[h]];
					cells[c].rightNeighbourWeights[rn] = xDiff*xDiff/hyp2;
					cells[c].xArea += cells[c].rightNeighbours[rn]->area;
					rn++;
				} else {
					cells[c].leftNeighbours[ln] = &cells[cells[c].neighbourIndices[h]];
					cells[c].leftNeighbourWeights[ln] = xDiff*xDiff/hyp2;
					cells[c].xArea += cells[c].leftNeighbours[ln]->area;
					ln++;
				}
			} else {
				if (cells[cells[c].neighbourIndices[h]].yCentre > cells[c].yCentre) {
					cells[c].topNeighbours[tn] = &cells[cells[c].neighbourIndices[h]];
					cells[c].topNeighbourWeights[tn] = yDiff*yDiff/hyp2;
					cells[c].yArea += cells[c].topNeighbours[tn]->area;
					tn++;
				} else {
					cells[c].bottomNeighbours[bn] = &cells[cells[c].neighbourIndices[h]];
					cells[c].bottomNeighbourWeights[bn] = yDiff*yDiff/hyp2;
					cells[c].yArea += cells[c].bottomNeighbours[bn]->area;
					bn++;
				}
			}
		}
	}
}

void VoronoiMesh::findExtremesD() {
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}
	iMAP->overlayAdjustment = true;
}

void VoronoiMesh::findExtremesDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
			if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}
	iMAP->overlayAdjustment = true;
}

void VoronoiMesh::findExtremesDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
			if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}
	iMAP->overlayAdjustment = true;
}

void VoronoiMesh::findExtremesDV() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
			if (potentialMax < cells[i].getPotential()) { potentialMax = cells[i].getPotential(); }
			if (potentialMin > cells[i].getPotential()) { potentialMin = cells[i].getPotential(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}
	iMAP->overlayAdjustment = true;
}

void VoronoiMesh::updateNeighbours() {
	for (int c = 0; c < this->numberOfClusters; c++) {
		cells[c].nRightNeighbours = 0;
		cells[c].nLeftNeighbours = 0;
		cells[c].nTopNeighbours = 0;
		cells[c].nBottomNeighbours = 0;
		// count number of neighbours
		for (int h = 0; h < cells[c].numberOfNeighbours; h++) {
			const double xDiff = fabs(cells[cells[c].neighbourIndices[h]].xCentre-cells[c].xCentre);
			const double yDiff = fabs(cells[cells[c].neighbourIndices[h]].yCentre-cells[c].yCentre);
			const double hyp2 = xDiff*xDiff+yDiff*yDiff;
			// set centroid distance-based condition for neighbouring zones (zones far away from each other are not neighbours!)
			if (sqrt(hyp2) < this->maximumNeighbourDistance) {
				if (xDiff*xDiff/hyp2 > 0.5) {
					if (cells[cells[c].neighbourIndices[h]].xCentre > cells[c].xCentre) {
						cells[c].nRightNeighbours++;
					} else {
						cells[c].nLeftNeighbours++;
					}
				} else {
					if (cells[cells[c].neighbourIndices[h]].yCentre > cells[c].yCentre) {
						cells[c].nTopNeighbours++;
					} else {
						cells[c].nBottomNeighbours++;
					}
				}
			}
		}
		// assign neighbours
		if (cells[c].rightNeighbours != NULL) { delete [] cells[c].rightNeighbours; }
		if (cells[c].rightNeighbourWeights != NULL) { delete [] cells[c].rightNeighbourWeights; }
		if (cells[c].leftNeighbours != NULL) { delete [] cells[c].leftNeighbours; }
		if (cells[c].leftNeighbourWeights != NULL) { delete [] cells[c].leftNeighbourWeights; }
		if (cells[c].topNeighbours != NULL) { delete [] cells[c].topNeighbours; }
		if (cells[c].topNeighbourWeights != NULL) { delete [] cells[c].topNeighbourWeights; }
		if (cells[c].bottomNeighbours != NULL) { delete [] cells[c].bottomNeighbours; }
		if (cells[c].bottomNeighbourWeights != NULL) { delete [] cells[c].bottomNeighbourWeights; }

		cells[c].rightNeighbours = new VoronoiCell*[cells[c].nRightNeighbours];
		cells[c].rightNeighbourWeights = new double[cells[c].nRightNeighbours];
		cells[c].leftNeighbours = new VoronoiCell*[cells[c].nLeftNeighbours];
		cells[c].leftNeighbourWeights = new double[cells[c].nLeftNeighbours];
		cells[c].topNeighbours = new VoronoiCell*[cells[c].nTopNeighbours];
		cells[c].topNeighbourWeights = new double[cells[c].nTopNeighbours];
		cells[c].bottomNeighbours = new VoronoiCell*[cells[c].nBottomNeighbours];
		cells[c].bottomNeighbourWeights = new double[cells[c].nBottomNeighbours];

		int rn = 0;
		int ln = 0;
		int tn = 0;
		int bn = 0;
		cells[c].xArea = 0.0;
		cells[c].yArea = 0.0;
		for (int h = 0; h < cells[c].numberOfNeighbours; h++) {
			const double xDiff = fabs(cells[cells[c].neighbourIndices[h]].xCentre-cells[c].xCentre);
			const double yDiff = fabs(cells[cells[c].neighbourIndices[h]].yCentre-cells[c].yCentre);
			const double hyp2 = xDiff*xDiff+yDiff*yDiff;
			// set centroid distance-based condition for neighbouring zones (zones far away from each other are not neighbours!)
			if (sqrt(hyp2) < this->maximumNeighbourDistance) {
				if (xDiff*xDiff/hyp2 > 0.5) {
					if (cells[cells[c].neighbourIndices[h]].xCentre > cells[c].xCentre) {
						cells[c].rightNeighbours[rn] = &cells[cells[c].neighbourIndices[h]];
						cells[c].rightNeighbourWeights[rn] = xDiff*xDiff/hyp2;
						cells[c].xArea += cells[c].rightNeighbours[rn]->area;
						rn++;
					} else {
						cells[c].leftNeighbours[ln] = &cells[cells[c].neighbourIndices[h]];
						cells[c].leftNeighbourWeights[ln] = xDiff*xDiff/hyp2;
						cells[c].xArea += cells[c].leftNeighbours[ln]->area;
						ln++;
					}
				} else {
					if (cells[cells[c].neighbourIndices[h]].yCentre > cells[c].yCentre) {
						cells[c].topNeighbours[tn] = &cells[cells[c].neighbourIndices[h]];
						cells[c].topNeighbourWeights[tn] = yDiff*yDiff/hyp2;
						cells[c].yArea += cells[c].topNeighbours[tn]->area;
						tn++;
					} else {
						cells[c].bottomNeighbours[bn] = &cells[cells[c].neighbourIndices[h]];
						cells[c].bottomNeighbourWeights[bn] = yDiff*yDiff/hyp2;
						cells[c].yArea += cells[c].bottomNeighbours[bn]->area;
						bn++;
					}
				}
			}
		}
	}
}

void VoronoiMesh::preInferD() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) activeDetections;
	yMean /= (double) activeDetections;
	dtMean /= (double) activeDetections;

	// length of optimization array
	dimensions = 1; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void VoronoiMesh::inferD() {
	for (int c = 0; c < numberOfClusters; c++) {
		if (active(c)) {
			setCurrentZone(c);
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorVoronoi, (dfunc));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Voronoi Tessellation\n1   (D) Inference\n2\n3   Zone %i / %i\n",c+1,this->getActiveCells());
				textDisplayUpdate(progress);
				findExtremesD();
			}

//			cells[c].setXForce(optimizationArray[0]);
//			cells[c].setYForce(optimizationArray[1]);
//			cells[c].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
			cells[c].setDiffusion(optimizationArray[0]);
			cells[c].setDiffusionLog(log(optimizationArray[0]));

			if (iMAP->stopCalculation) { return; }
		} else {
			cells[c].setForceX(0.0);
			cells[c].setForceY(0.0);
			cells[c].setForceMagnitude(0.0);
			cells[c].setDiffusion(0.0);
			cells[c].setDiffusionLog(0.0);
		}
	}
}

int VoronoiMesh::postInferD() {
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return 0;

//	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void VoronoiMesh::preInferDSmoothing() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ( cells[i].active() ) {
			cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
		}
	}

	// length of optimization array
	dimensions = counter; // D_ij,V_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// initialize potentials
	int g = 0;
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			optimizationArray[g] = 0.5*(D_eff_x + D_eff_y);
			g++;
		}
	}
}

void VoronoiMesh::inferDSmoothing() {
	char label[50];
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingVoronoiSelection, (dfunc));
		sprintf(label,"Cost = %f\n",dPosteriorSmoothingVoronoiSelection(optimizationArray));
		textDisplayUpdate(label);
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingVoronoi, (dfunc));
		sprintf(label,"Cost = %f\n",dPosteriorSmoothingVoronoi(optimizationArray));
		textDisplayUpdate(label);
	}
}

int VoronoiMesh::postInferDSmoothing() {
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setDiffusion(optimizationArray[getCell(i)->getIdentifier()]);
			cells[i].setDiffusionLog(log(optimizationArray[getCell(i)->getIdentifier()]));
			cells[i].setPotential(optimizationArray[getCell(i)->getIdentifier()+1]);
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	return 1;
}

void VoronoiMesh::preInferDDr() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) activeDetections;
	yMean /= (double) activeDetections;
	dtMean /= (double) activeDetections;

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void VoronoiMesh::inferDDr() {
	for (int c = 0; c < numberOfClusters; c++) {
		if (active(c)) {

			setCurrentZone(c);
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.0;
			optimizationArray[1] = 0.0;
			optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorVoronoi, (dfunc));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Voronoi Tessellation\n1   (D,Drift) Inference\n2\n3   Zone %i / %i\n",c+1,this->getActiveCells());
				textDisplayUpdate(progress);
				findExtremesDDr();
			}

			cells[c].setForceX(optimizationArray[0]);
			cells[c].setForceY(optimizationArray[1]);
			cells[c].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
			cells[c].setDiffusion(optimizationArray[2]);
			cells[c].setDiffusionLog(log(optimizationArray[2]));

			if (iMAP->stopCalculation) { return; }
		} else {
			cells[c].setForceX(0.0);
			cells[c].setForceY(0.0);
			cells[c].setForceMagnitude(0.0);
			cells[c].setDiffusion(0.0);
			cells[c].setDiffusionLog(0.0);
		}
	}
}

int VoronoiMesh::postInferDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
			if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return 0;

//	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void VoronoiMesh::preInferDF() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) activeDetections;
	yMean /= (double) activeDetections;
	dtMean /= (double) activeDetections;

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void VoronoiMesh::inferDF() {
	for (int c = 0; c < numberOfClusters; c++) {
		if (active(c)) {
			setCurrentZone(c);
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.0;
			optimizationArray[1] = 0.0;
			optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorVoronoi, (dfunc));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Voronoi Tessellation\n1   (D,F) Inference\n2\n3   Zone %i / %i\n",c+1,this->getActiveCells());
				textDisplayUpdate(progress);
				findExtremesDF();
			}

			cells[c].setForceX(optimizationArray[0]);
			cells[c].setForceY(optimizationArray[1]);
			cells[c].setForceMagnitude(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));
			cells[c].setDiffusion(optimizationArray[2]);
			cells[c].setDiffusionLog(log(optimizationArray[2]));

			if (iMAP->stopCalculation) { return; }
		} else {
			cells[c].setForceX(0.0);
			cells[c].setForceY(0.0);
			cells[c].setForceMagnitude(0.0);
			cells[c].setDiffusion(0.0);
			cells[c].setDiffusionLog(0.0);
		}
	}
}

int VoronoiMesh::postInferDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	// save x,y-forces and diffusion for each of the activated zones
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			if (forceMax < cells[a].getForceMagnitude()) { forceMax = cells[a].getForceMagnitude(); }
			if (forceMin > cells[a].getForceMagnitude()) { forceMin = cells[a].getForceMagnitude(); }
			if (diffusionMax < cells[a].getDiffusion()) { diffusionMax = cells[a].getDiffusion(); }
			if (diffusionMin > cells[a].getDiffusion()) { diffusionMin = cells[a].getDiffusion(); }
			if (diffusionLogMax < cells[a].getDiffusionLog()) { diffusionLogMax = cells[a].getDiffusionLog(); }
			if (diffusionLogMin > cells[a].getDiffusionLog()) { diffusionLogMin = cells[a].getDiffusionLog(); }
		}
		else {
			cells[a].setForceX(0.0);
			cells[a].setForceY(0.0);
			cells[a].setForceMagnitude(0.0);
			cells[a].setDiffusion(0.0);
			cells[a].setDiffusionLog(0.0);
			// potential applied afterwards for zonal (if requested)
		}
	}

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void VoronoiMesh::preInferDFSmoothing() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ( cells[i].active() ) {
//			cells[i].setPotential( -log((double)cells[i].getCount()/(double)maxCell->getCount()) );
			cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
		}
	}

	// length of optimization array
	dimensions = 3*counter;

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// initialize potentials
	int g = 0;
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			optimizationArray[3*g] = 0.5*(D_eff_x + D_eff_y);
			optimizationArray[3*g+1] = 0.0;
			optimizationArray[3*g+2] = 0.0;
			g++;
		}
	}

}

void VoronoiMesh::preInferDDrSmoothing() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ( cells[i].active() ) {
//			cells[i].setPotential( -log((double)cells[i].getCount()/(double)maxCell->getCount()) );
			cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
		}
	}

	// length of optimization array
	dimensions = 3*counter;

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// initialize potentials
	int g = 0;
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			optimizationArray[3*g] = 0.5*(D_eff_x + D_eff_y);
			optimizationArray[3*g+1] = 0.0;
			optimizationArray[3*g+2] = 0.0;
			g++;
		}
	}

}

void VoronoiMesh::inferDFSmoothing() {
	if (file->voronoiMesh->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingVoronoiSelection, (dfunc));
	}
	else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingVoronoi, (dfunc));
	}
}

void VoronoiMesh::inferDDrSmoothing() {
	if (file->voronoiMesh->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingVoronoiSelection, (dfunc));
	}
	else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingVoronoi, (dfunc));
	}
}

int VoronoiMesh::postInferDFSmoothing() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceX(optimizationArray[3*getCell(i)->getIdentifier()+1]);
			cells[i].setForceY(optimizationArray[3*getCell(i)->getIdentifier()+2]);
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			cells[i].setDiffusion(optimizationArray[3*getCell(i)->getIdentifier()]);
			cells[i].setDiffusionLog(log(optimizationArray[3*getCell(i)->getIdentifier()]));
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

int VoronoiMesh::postInferDDrSmoothing() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceX(optimizationArray[3*getCell(i)->getIdentifier()+1]);
			cells[i].setForceY(optimizationArray[3*getCell(i)->getIdentifier()+2]);
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			cells[i].setDiffusion(optimizationArray[3*getCell(i)->getIdentifier()]);
			cells[i].setDiffusionLog(log(optimizationArray[3*getCell(i)->getIdentifier()]));
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	return 1;
}

void VoronoiMesh::inferPotentialsDF() {

	iMAP->pauseCalculation = false;

	// do not reinitialize if pause button was pressed
	if (iMAP->stopCalculation == false) {
		int k = 0;
		for (int b = 0; b < numberOfClusters; b++) {
			if ( this->active(b) ) {
				k++;
			}
		}

		dimensions = totalVariables = k;

		if (potentialArray != NULL) {
			delete [] potentialArray;
			potentialArray = NULL;
		}
		potentialArray = new double[dimensions];

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// initialize potentials
		int i = 0;
		for (int c = 0; c < numberOfClusters; c++) {
			if ((this->active(c))) {
				potentialArray[getIdentifier(c)] = -log( (double)cells[c].getCount()/(double)maxCell->getCount() );
				i++;
			}
		}
	}
	iMAP->stopCalculation = false;

	// perform optimization of potential array
	iMAP->updatePotentialCalculation = true;
	dfpmin(potentialArray,dimensions,GTOL,iterations,fret,squareDifferenceVoronoi,dfunc);
	iMAP->updatePotentialCalculation = false;

	// set potentials
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	int j = 0;
	for (int c = 0; c < numberOfClusters; c++) {
		if ( (this->active(c)) /*&& (this->getXNeighbours(a,b)!=0) && (this->getYNeighbours(a,b)!=0)*/ ) {
			cells[c].setPotential(potentialArray[getIdentifier(c)]);
			if (potentialMax < cells[c].getPotential()) { potentialMax = cells[c].getPotential(); }
			if (potentialMin > cells[c].getPotential()) { potentialMin = cells[c].getPotential(); }
			j++;
		}
		else {
			cells[c].setPotential(1e3);
		}
	}

	offsetPotentials();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}
}

// constant point per zone initialization (before feeding into kmeans/hmeans)
double* VoronoiMesh::constantPointsInitialization() {

	const int N = file->localizationCount;
	const int pointsPerZone = (int)ceil((double)N/(double)this->numberOfClusters);

	double radius = 0.0;
	const double inc = 0.005;
	int totalCount = 0;
	int zoneCount = 0;
	int zone = 0;

	int *check = new int[N];

	// find left-most point
	double minLeftValue = 1000000.0;
	int mi;
	for (int h = 0; h < N; h++) {
		if (minLeftValue > file->xPointer[h]) {
			minLeftValue = file->xPointer[h];
			mi = h;
		}
		// initialize check array to -1
		check[h] = -1;
	}

	char label[60];
	while (totalCount < N) {
		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			sprintf(label,"Mesh Initialization: %2.0f%%\n\n",100.0*(float)totalCount/(float)N);

			textDisplayUpdate(label);
		}
		else { fprintf(stderr,"Mesh Initialization: %2.0f%%\n\n",100.0*(float)totalCount/(float)N); }

		zoneCount = 0;
		radius = 0.0;
		while (zoneCount < pointsPerZone && totalCount + zoneCount < N) {
			radius += inc;
			for (int q = 0; q < N; q++) {
				if (sqrt( (file->xPointer[q]-file->xPointer[mi])*(file->xPointer[q]-file->xPointer[mi])+(file->yPointer[q]-file->yPointer[mi])*(file->yPointer[q]-file->yPointer[mi]) ) <= radius) {
					if (check[q] == -1) {
						check[q] = zone;
						zoneCount++;

						if (zoneCount == pointsPerZone) {
							// find new left-most point
							minLeftValue = 1000000.0;
							for (int h = 0; h < N; h++) {
								if (check[h] == -1) {
									if (minLeftValue > file->xPointer[h]) {
										minLeftValue = file->xPointer[h];
										mi = h;
									}
								}
							}
							break;
						}

					}
				}
			}
		}
//		printf("zone %i : %i points\n",zone,zoneCount);
		zone++;
		totalCount += zoneCount;
	}

	double *centres = new double[2*numberOfClusters];
	int *num = new int[numberOfClusters];
	for (int j = 0; j < numberOfClusters; j++) {
		centres[2*j] = 0.0;
		centres[2*j+1] = 0.0;
		num[j] = 0;
	}

	// get barycentres of zones
	for (int p = 0; p < N; p++) {
		centres[2*check[p]] += file->xPointer[p];
		centres[2*check[p]+1] += file->yPointer[p];
		num[check[p]]++;
	}

	// average totals
	for (int e = 0; e < numberOfClusters; e++) {
		centres[2*e] /= (double)num[e];
		centres[2*e+1] /= (double)num[e];
//		printf("%i : [%f,%f]\n",e,centres[2*e],centres[2*e+1]);
	}

	delete [] check;
	delete [] num;

	return centres;
}

double* VoronoiMesh::constantPointsInitializationSelection() {

	const int N = selection.count;
	const int pointsPerZone = (int)ceil((double)N/(double)this->numberOfClusters);

	double radius = 0.0;
	const double inc = 0.005;
	int totalCount = 0;
	int zoneCount = 0;
	int zone = 0;

	int *check = new int[N];

	// find left-most point
	double minLeftValue = 1000000.0;
	int mi;
	for (int h = 0; h < N; h++) {
		if (minLeftValue > selection.xPointer[h]) {
			minLeftValue = selection.xPointer[h];
			mi = h;
		}
		// initialize check array to -1
		check[h] = -1;
	}
	char label[60];
	while (totalCount < N) {
		if (iMAP->updateDisplay) {
			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
			sprintf(label,"Mesh Initialization: %2.0f%%\n\n",100.0*(float)totalCount/(float)N);

			textDisplayUpdate(label);
		}
		else { fprintf(stderr,"Mesh Initialization: %2.0f%%\n\n",100.0*(float)totalCount/(float)N); }

		zoneCount = 0;
		radius = 0.0;
		while (zoneCount < pointsPerZone && totalCount + zoneCount < N) {
			radius += inc;
			for (int q = 0; q < N; q++) {
				if (sqrt( (selection.xPointer[q]-selection.xPointer[mi])*(selection.xPointer[q]-selection.xPointer[mi])+(selection.yPointer[q]-selection.yPointer[mi])*(selection.yPointer[q]-selection.yPointer[mi]) ) <= radius) {
					if (check[q] == -1) {
						check[q] = zone;
						zoneCount++;

						if (zoneCount == pointsPerZone) {
							// find new left-most point
							minLeftValue = 1000000.0;
							for (int h = 0; h < N; h++) {
								if (check[h] == -1) {
									if (minLeftValue > selection.xPointer[h]) {
										minLeftValue = selection.xPointer[h];
										mi = h;
									}
								}
							}
							break;
						}

					}
				}
			}
		}
//		printf("zone %i : %i points\n",zone,zoneCount);
		zone++;
		totalCount += zoneCount;
	}

	double *centres = new double[2*numberOfClusters];
	int *num = new int[numberOfClusters];
	for (int j = 0; j < numberOfClusters; j++) {
		centres[2*j] = 0.0;
		centres[2*j+1] = 0.0;
		num[j] = 0;
	}

	// get barycentres of zones
	for (int p = 0; p < N; p++) {
		centres[2*check[p]] += selection.xPointer[p];
		centres[2*check[p]+1] += selection.yPointer[p];
		num[check[p]]++;
	}

	// average totals
	for (int e = 0; e < numberOfClusters; e++) {
		centres[2*e] /= (double)num[e];
		centres[2*e+1] /= (double)num[e];
//		printf("%i : [%f,%f]\n",e,centres[2*e],centres[2*e+1]);
	}

	delete [] check;
	delete [] num;

	return centres;
}

void VoronoiMesh::preInferDV() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	// get beta from slider
	beta = file->voronoiMeshGui->getBeta();

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;
	int counter = 0;

	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if (cells[i].active()) {
			totalVariables++;
			for (int k = 0; k < cells[i].getCount(); k++) {
				if (cells[i].getIndex(k) != file->localizationCount-1) {
					if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
						xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
						yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
						dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
						activeDetections++;
					}
				}
			}
			// for global diffusion calculation
			cells[i].setIdentifier(counter);
			counter++;
			activeCells++;
			cells[i].resetPriors();

		}
		else {
			cells[i].setIdentifier(-999);
		}
	}

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// set default inferred values to cell
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ( cells[i].active() ) {
			cells[i].setPotential( -log((double)cells[i].getCount()/(double)maxCell->getCount()) );
			cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
		}
	}

	// length of optimization array
	dimensions = 2*counter; // D_ij,V_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];

	// initialize potentials
	int g = 0;
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			optimizationArray[2*g] = 0.5*(D_eff_x + D_eff_y);
			optimizationArray[2*g+1] = cells[i].getPotential();
			g++;
		}
	}

}

void VoronoiMesh::inferDV() {
	if (file->voronoiMesh->selectionMode()) {
		if (file->voronoiMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingVoronoiSelection, (dfunc));
		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorVoronoiSelection, (dfunc));
		}
	}
	else {
		if (file->voronoiMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingVoronoi, (dfunc));
		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorVoronoi, (dfunc));
		}
	}
}

int VoronoiMesh::postInferDV() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceX(dvGradVxVoronoi(optimizationArray,i));
			cells[i].setForceY(dvGradVyVoronoi(optimizationArray,i));
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			cells[i].setDiffusion(optimizationArray[2*getCell(i)->getIdentifier()]);
			cells[i].setDiffusionLog(log(optimizationArray[2*getCell(i)->getIdentifier()]));
			cells[i].setPotential(optimizationArray[2*getCell(i)->getIdentifier()+1]);
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
			if (potentialMax < cells[i].getPotential()) { potentialMax = cells[i].getPotential(); }
			if (potentialMin > cells[i].getPotential()) { potentialMin = cells[i].getPotential(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	return 1;
}

bool VoronoiMesh::roActive(int c) {
	return cells[c].roActive();
}

int VoronoiMesh::getActiveCellId(int counter) {
	for (int g = 0; g < numberOfClusters; g++) {
		if (cells[g].identifier == counter) {
			return g;
		}
	}
	fprintf(stderr,"cluster not found\n");
	return -1;
}

void VoronoiMesh::inferRandomizedOptimizationDV() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int counter = 0;

		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
				totalVariables++;
				for (int k = 0; k < cells[i].getCount(); k++) {
					if (cells[i].getIndex(k) != file->localizationCount-1) {
						if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
							xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
							yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
							dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
							activeDetections++;
						}
					}
				}
				// for global diffusion calculation
				cells[i].setIdentifier(counter);
				counter++;
				activeCells++;
			}
			else {
				cells[i].setIdentifier(-999);
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[2*totalVariables];

		// set default inferred values to cell
		counter = 0;
		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
//				cells[i].setRoIdentifier(-1);
				cells[i].setRoIterations(0);
				cells[i].roDeactivate();
				cells[i].setIdentifier(counter);
				cells[i].setPotential( -log((double)cells[i].getCount()/(double)maxCell->getCount()) );
				cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
				intermediateOptimizationArray[2*counter] = cells[i].getDiffusion();
				intermediateOptimizationArray[2*counter+1] = cells[i].getPotential();
				counter++;
			}
		}

		activeCells = counter;
	}

	const double tolerance = (double)file->voronoiMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->voronoiMeshGui->roMaximumIterationsSlider->value();

	int loops = 0;

	int centerId;

	hessian = NULL;
	optimizationArray = NULL;
	roClusters = NULL;

	const float radius = file->voronoiMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roClusters[u]].roDeactivate();
				cells[roClusters[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = getActiveCellId(rand() % activeCells);

		// count number of cells
		roZones = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) { roZones++; }
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }

		dimensions = 2*roZones;

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roClusters != NULL) { delete [] roClusters; }
		roClusters = new int[roZones];

		int p = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) {
					roClusters[p] = a;
					optimizationArray[2*p] = cells[a].getDiffusion();
					optimizationArray[2*p+1] = cells[a].getPotential();
					cells[a].roActivate();
					cells[a].setRoIdentifier(p);
					p++;
				}
			}
		}

		if (file->voronoiMesh->selectionMode()) {
			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingVoronoiRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorVoronoiRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
		} else {
			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingVoronoiRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorVoronoiRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[2*cells[roClusters[t]].getIdentifier()] = optimizationArray[2*t];
			intermediateOptimizationArray[2*cells[roClusters[t]].getIdentifier()+1] = optimizationArray[2*t+1];
		}

		oldCost = cost;

		if (file->voronoiMesh->selectionMode()) {
			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingVoronoiSelection(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorVoronoiSelection(intermediateOptimizationArray);
			}
		} else {
			if (file->voronoiMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingVoronoi(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorVoronoi(intermediateOptimizationArray);
			}
		}

		for (int t = 0; t < roZones; t++) {
			cells[roClusters[t]].setForceX(dvGradVxVoronoiRandomizedOptimization(optimizationArray,roClusters[t]));
			cells[roClusters[t]].setForceY(dvGradVyVoronoiRandomizedOptimization(optimizationArray,roClusters[t]));
			cells[roClusters[t]].setForceMagnitude( sqrt(cells[roClusters[t]].getForceX()*cells[roClusters[t]].getForceX() + cells[roClusters[t]].getForceY()*cells[roClusters[t]].getForceY()) );
			cells[roClusters[t]].setDiffusion(optimizationArray[2*t]);
			cells[roClusters[t]].setDiffusionLog(log(optimizationArray[2*t]));
			cells[roClusters[t]].setPotential(optimizationArray[2*t+1]);
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDV();
		Fl::unlock();
		Fl::check();

		loops++;

	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
			if (potentialMax < cells[i].getPotential()) { potentialMax = cells[i].getPotential(); }
			if (potentialMin > cells[i].getPotential()) { potentialMin = cells[i].getPotential(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < getNumberOfClusters(); u++) {
		cells[u].roDeactivate();
		cells[u].setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();
}

void VoronoiMesh::inferRandomizedOptimizationD() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int counter = 0;

		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
				totalVariables++;
				for (int k = 0; k < cells[i].getCount(); k++) {
					if (cells[i].getIndex(k) != file->localizationCount-1) {
						if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
							xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
							yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
							dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
							activeDetections++;
						}
					}
				}
				// for global diffusion calculation
				cells[i].setIdentifier(counter);
				counter++;
				activeCells++;
			}
			else {
				cells[i].setIdentifier(-999);
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[totalVariables];

		// set default inferred values to cell
		counter = 0;
		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
//				cells[i].setRoIdentifier(-1);
				cells[i].setRoIterations(0);
				cells[i].roDeactivate();
				cells[i].setIdentifier(counter);
				cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
				intermediateOptimizationArray[counter] = cells[i].getDiffusion();
				counter++;
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->voronoiMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->voronoiMeshGui->roMaximumIterationsSlider->value();

	int loops = 0;

	int centerId;

	hessian = NULL;
	optimizationArray = NULL;
	roClusters = NULL;

	const float radius = file->voronoiMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roClusters[u]].roDeactivate();
				cells[roClusters[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = getActiveCellId(rand() % activeCells);

		// count number of cells
		roZones = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) { roZones++; }
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }

		dimensions = roZones;

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roClusters != NULL) { delete [] roClusters; }
		roClusters = new int[roZones];

		int p = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) {
					roClusters[p] = a;
					optimizationArray[p] = cells[a].getDiffusion();
					cells[a].roActivate();
					cells[a].setRoIdentifier(p);
					p++;
				}
			}
		}

		if (file->voronoiMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingVoronoiRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingVoronoiRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[cells[roClusters[t]].getIdentifier()] = optimizationArray[t];
		}

		oldCost = cost;

		if (file->voronoiMesh->selectionMode()) {
			cost = dPosteriorSmoothingVoronoiSelection(intermediateOptimizationArray);
		} else {
			cost = dPosteriorSmoothingVoronoi(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roClusters[t]].setDiffusion(optimizationArray[t]);
			cells[roClusters[t]].setDiffusionLog(log(optimizationArray[t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesD();
		Fl::unlock();
		Fl::check();

		loops++;

	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < getNumberOfClusters(); u++) {
		cells[u].roDeactivate();
		cells[u].setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();
}

int VoronoiMesh::inferRandomizedOptimizationDF() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int counter = 0;

		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
				totalVariables++;
				for (int k = 0; k < cells[i].getCount(); k++) {
					if (cells[i].getIndex(k) != file->localizationCount-1) {
						if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
							xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
							yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
							dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
							activeDetections++;
						}
					}
				}
				// for global diffusion calculation
				cells[i].setIdentifier(counter);
				counter++;
				activeCells++;
			}
			else {
				cells[i].setIdentifier(-999);
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		// set default inferred values to cell
		counter = 0;
		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
//				cells[i].setRoIdentifier(-1);
				cells[i].setRoIterations(0);
				cells[i].roDeactivate();
				cells[i].setIdentifier(counter);
				cells[i].setForceX(0.0);
				cells[i].setForceY(0.0);
				cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
				intermediateOptimizationArray[3*counter] = cells[i].getDiffusion();
				intermediateOptimizationArray[3*counter+1] = cells[i].getForceX();
				intermediateOptimizationArray[3*counter+2] = cells[i].getForceY();
				counter++;
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->voronoiMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->voronoiMeshGui->roMaximumIterationsSlider->value();

	int loops = 0;

	int centerId;

	hessian = NULL;
	optimizationArray = NULL;
	roClusters = NULL;

	const float radius = file->voronoiMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roClusters[u]].roDeactivate();
				cells[roClusters[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = getActiveCellId(rand() % activeCells);

		// count number of cells
		roZones = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) { roZones++; }
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }

		dimensions = 3*roZones;

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roClusters != NULL) { delete [] roClusters; }
		roClusters = new int[roZones];

		int p = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) {
					roClusters[p] = a;
					optimizationArray[3*p] = cells[a].getDiffusion();
					optimizationArray[3*p+1] = cells[a].getForceX();
					optimizationArray[3*p+2] = cells[a].getForceY();
					cells[a].roActivate();
					cells[a].setRoIdentifier(p);
					p++;
				}
			}
		}

		if (file->voronoiMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingVoronoiRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingVoronoiRandomizedOptimization, (dfunc), maxIterations, loops+1);

		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()] = optimizationArray[3*t];
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->voronoiMesh->selectionMode()) {
			cost = dfPosteriorSmoothingVoronoiSelection(intermediateOptimizationArray);
		} else {
			cost = dfPosteriorSmoothingVoronoi(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roClusters[t]].setForceX(optimizationArray[3*t+1]);
			cells[roClusters[t]].setForceY(optimizationArray[3*t+2]);
			cells[roClusters[t]].setForceMagnitude( sqrt(cells[roClusters[t]].getForceX()*cells[roClusters[t]].getForceX() + cells[roClusters[t]].getForceY()*cells[roClusters[t]].getForceY()) );
			cells[roClusters[t]].setDiffusion(optimizationArray[3*t]);
			cells[roClusters[t]].setDiffusionLog(log(optimizationArray[3*t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDF();
		Fl::unlock();
		Fl::check();

		loops++;

	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < getNumberOfClusters(); u++) {
		cells[u].roDeactivate();
		cells[u].setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void VoronoiMesh::inferRandomizedOptimizationDDr() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->voronoiMeshGui->saveButton->deactivate();
	sigma = file->voronoiMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		totalVariables = 0;
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;
		int counter = 0;

		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
				totalVariables++;
				for (int k = 0; k < cells[i].getCount(); k++) {
					if (cells[i].getIndex(k) != file->localizationCount-1) {
						if (file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)] == file->nPointer[file->voronoiMesh->getCell(i)->getIndex(k)+1]) {
							xMean += fabs(file->xPointer[cells[i].getIndex(k)+1]-file->xPointer[cells[i].getIndex(k)]);
							yMean += fabs(file->yPointer[cells[i].getIndex(k)+1]-file->yPointer[cells[i].getIndex(k)]);
							dtMean += fabs(file->tPointer[cells[i].getIndex(k)+1]-file->tPointer[cells[i].getIndex(k)]);
							activeDetections++;
						}
					}
				}
				// for global diffusion calculation
				cells[i].setIdentifier(counter);
				counter++;
				activeCells++;
			}
			else {
				cells[i].setIdentifier(-999);
			}
		}

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		// set default inferred values to cell
		counter = 0;
		for (int i = 0; i < this->getNumberOfClusters(); i++) {
			if (cells[i].active()) {
//				cells[i].setRoIdentifier(-1);
				cells[i].setRoIterations(0);
				cells[i].roDeactivate();
				cells[i].setIdentifier(counter);
				cells[i].setForceX(0.0);
				cells[i].setForceY(0.0);
				cells[i].setDiffusion(0.5*(D_eff_x + D_eff_y));
				intermediateOptimizationArray[3*counter] = cells[i].getDiffusion();
				intermediateOptimizationArray[3*counter+1] = cells[i].getForceX();
				intermediateOptimizationArray[3*counter+2] = cells[i].getForceY();
				counter++;
			}
		}
		activeCells = counter;
	}

	const double tolerance = (double)file->voronoiMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->voronoiMeshGui->roMaximumIterationsSlider->value();

	int loops = 0;

	int centerId;

	hessian = NULL;
	optimizationArray = NULL;
	roClusters = NULL;

	const float radius = file->voronoiMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		if (loops > 0) {
			for (int u = 0; u < roZones; u++) {
				cells[roClusters[u]].roDeactivate();
				cells[roClusters[u]].setRoIdentifier(-1);
			}
		}

		// select center of block randomly
		centerId = getActiveCellId(rand() % activeCells);

		// count number of cells
		roZones = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) { roZones++; }
			}
		}

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		// clear optimization array if re-run
		if (optimizationArray != NULL) { delete [] optimizationArray; }

		dimensions = 3*roZones;

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		optimizationArray = new double[dimensions];

		// select at random the zones to optimize
		if (roClusters != NULL) { delete [] roClusters; }
		roClusters = new int[roZones];

		int p = 0;
		for (int a = 0; a < getNumberOfClusters(); a++) {
			if (cells[a].active()) {
				const float dist = sqrt( (cells[centerId].xCentroid-cells[a].xCentroid)*(cells[centerId].xCentroid-cells[a].xCentroid) +
										 (cells[centerId].yCentroid-cells[a].yCentroid)*(cells[centerId].yCentroid-cells[a].yCentroid) );
				if (dist < radius) {
					roClusters[p] = a;
					optimizationArray[3*p] = cells[a].getDiffusion();
					optimizationArray[3*p+1] = cells[a].getForceX();
					optimizationArray[3*p+2] = cells[a].getForceY();
					cells[a].roActivate();
					cells[a].setRoIdentifier(p);
					p++;
				}
			}
		}

		if (file->voronoiMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingVoronoiRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingVoronoiRandomizedOptimization, (dfunc), maxIterations, loops+1);

		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()] = optimizationArray[3*t];
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*cells[roClusters[t]].getIdentifier()+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->voronoiMesh->selectionMode()) {
			cost = ddrPosteriorSmoothingVoronoiSelection(intermediateOptimizationArray);
		} else {
			cost = ddrPosteriorSmoothingVoronoi(intermediateOptimizationArray);
		}

		for (int t = 0; t < roZones; t++) {
			cells[roClusters[t]].setForceX(optimizationArray[3*t+1]);
			cells[roClusters[t]].setForceY(optimizationArray[3*t+2]);
			cells[roClusters[t]].setForceMagnitude( sqrt(cells[roClusters[t]].getForceX()*cells[roClusters[t]].getForceX() + cells[roClusters[t]].getForceY()*cells[roClusters[t]].getForceY()) );
			cells[roClusters[t]].setDiffusion(optimizationArray[3*t]);
			cells[roClusters[t]].setDiffusionLog(log(optimizationArray[3*t]));
		}

		Fl::lock();
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		Fl::unlock();
		Fl::check();

		loops++;

	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// save x,y-forces, diffusion, and potential for each of the activated zones
	for (int i = 0; i < this->getNumberOfClusters(); i++) {
		if ((this->active(i))) {
			cells[i].setForceMagnitude(sqrt(cells[i].getForceX()*cells[i].getForceX()+cells[i].getForceY()*cells[i].getForceY()));
			if (forceMax < cells[i].getForceMagnitude()) { forceMax = cells[i].getForceMagnitude(); }
			if (forceMin > cells[i].getForceMagnitude()) { forceMin = cells[i].getForceMagnitude(); }
			if (diffusionMax < cells[i].getDiffusion()) { diffusionMax = cells[i].getDiffusion(); }
			if (diffusionMin > cells[i].getDiffusion()) { diffusionMin = cells[i].getDiffusion(); }
			if (diffusionLogMax < cells[i].getDiffusionLog()) { diffusionLogMax = cells[i].getDiffusionLog(); }
			if (diffusionLogMin > cells[i].getDiffusionLog()) { diffusionLogMin = cells[i].getDiffusionLog(); }
		} else {
			cells[i].setForceX(0.0);
			cells[i].setForceY(0.0);
			cells[i].setForceMagnitude(0.0);
			cells[i].setDiffusion(0.0);
			cells[i].setDiffusionLog(0.0);
			cells[i].setPotential(1.0e3);
		}
	}

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;
	}

	for (int u = 0; u < getNumberOfClusters(); u++) {
		cells[u].roDeactivate();
		cells[u].setRoIdentifier(-1);
	}

	// enable inference variable overlays
	file->inferred = true;
	file->voronoiMeshGui->saveButton->activate();
}

void VoronoiMesh::offsetPotentials() {
	// set reference potential to zero
	for (int a = 0; a < numberOfClusters; a++) {
		if (this->active(a)) {
			cells[a].setPotential(cells[a].getPotential()-potentialMin);
		}
	}
	potentialMax -= potentialMin;
	potentialMin = 0.0;
}

void VoronoiMesh::exportMesh() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Voronoi Mesh File");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Mesh File\t*.vmesh\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".vmesh");
		FILE * writeFile;

		writeFile = fopen(nativeFilename,"wb");

		// File Information
		fprintf(writeFile,"FILE INFORMATION\n");
		fprintf(writeFile,"Filename:\t%s\n",file->fileName);
		fprintf(writeFile,"Total Number of Trajectories:\t%i\n",file->numberOfFiles);
		fprintf(writeFile,"Duration [s]:\t%.3f\n",iMAP->endIntervalSlider->value()-iMAP->startIntervalSlider->value());
		fprintf(writeFile,"Acquisition Time [ms]:\t%.1f\n",1000.0*iMAP->exposureTime);
		if (this->selectionMode()) {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",selection.xMax-selection.xMin,selection.yMax-selection.yMin);
		} else {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",this->xMax-this->xMin,this->yMax-this->yMin);
		}

		// Mesh Information
		fprintf(writeFile,"MESH INFORMATION\n");
		fprintf(writeFile,"Meshing: Voronoi Tessellation\n");
		fprintf(writeFile,"Selection Mode: %i\n",this->selectionMode());
		fprintf(writeFile,"Number of Clusters:\t%i\n",getNumberOfClusters());
		fprintf(writeFile,"Total Cells:\t%i\n",getNumberOfClusters());
		fprintf(writeFile,"Total Localizations:\t%i\n",file->localizationCount);
		fprintf(writeFile,"Active Cells:\t%i\n",this->getActiveCells());
		fprintf(writeFile,"Active Localizations:\t%i\n\n",this->getActiveDetections());

		// Inference Parameters
		fprintf(writeFile,"INFERENCE PARAMETERS\n");
		fprintf(writeFile,"Noise Sigma [nm]:\t%.1f\n",this->getSigma()*1000.0);
		fprintf(writeFile,"Maximum Neighbour Distance [nm]:\t%.1f\n",this->getMaximumNeighbourDistance()*1000.0);
		fprintf(writeFile,"Minimum Points per Cell:\t%i\n",this->getMinPointsPerCell());
		fprintf(writeFile,"Optimization Scheme:\t%i\n",this->getOptimizationMode());
		fprintf(writeFile,"Potential Calculated:\t%i\n",this->zonalPotentialsCalculated);
		fprintf(writeFile,"Prior Enabled:\t%i\n",this->smoothingPriorEnabled());
		if (this->optimizationMode == 4) {
			fprintf(writeFile,"Polynomial Order:\t%i\n",this->getPolynomialOrder());
			fprintf(writeFile,"Polynomial Coefficients:\t");
			for (int q = 0; q < this->getCoefficients()-1; q++) {
				fprintf(writeFile,"%f\t",this->getOptimizationArray(q));
			}
			fprintf(writeFile,"%f\n\n",this->getOptimizationArray(file->squareMesh->getCoefficients()-1));
		} else {
			fprintf(writeFile,"Polynomial Order:\tN/A\n");
			fprintf(writeFile,"Polynomial Coefficients:\tN/A\n\n");
		}

		// Data
		if (file->optimizationFunction == 2) { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tDx\t\tDy\t\tV\n"); }
		else { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tFx\t\tFy\t\tV\n"); }
		for (int a = 0; a < getNumberOfClusters(); a++) {
				fprintf(writeFile,"%i\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
						getCell(a)->active(),
						getCell(a)->getXCentre(),
						getCell(a)->getYCentre(),
						getCell(a)->getCount(),
						getCell(a)->getXCentroid(),
						getCell(a)->getYCentroid(),
						getCell(a)->getVariance(),
						getCell(a)->getDiffusion(),
						getCell(a)->getForceX(),
						getCell(a)->getForceY(),
						getCell(a)->getPotential()
						);

		}

		fclose(writeFile);
	}

}

void VoronoiMesh::exportMeshBatch() {

	char nativeFilename[FILENAME_MAX];
	char fileNameTemp[FILENAME_MAX];

	const int fileLength = strlen(file->fileName);
	for (int g = 0; g < fileLength-6; g++) {
		fileNameTemp[g] = file->fileName[g];
	}
	fileNameTemp[fileLength-6] = '\0';

	strcpy(nativeFilename,file->filePath);
	strcat(nativeFilename,fileNameTemp);

	// store filename
	sprintf(nativeFilename,"%s%s",nativeFilename,".vmesh");
	FILE * writeFile;

	writeFile = fopen(nativeFilename,"wb");

	// File Information
	fprintf(writeFile,"FILE INFORMATION\n");
	fprintf(writeFile,"Filename:\t%s\n",file->fileName);
	fprintf(writeFile,"Number of Trajectories:\t%i\n",file->numberOfFiles);
	fprintf(writeFile,"Duration [s]:\t%.3f\n",file->tMax-file->tMin);
	fprintf(writeFile,"Acquisition Time [ms]:\t%.1f\n",1000.0*iMAP->exposureTime);
	fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",this->xMax-this->xMin,this->yMax-this->yMin);

	// Mesh Information
	fprintf(writeFile,"MESH INFORMATION\n");
	fprintf(writeFile,"Meshing: Voronoi Tessellation\n");
	fprintf(writeFile,"Selection Mode: %i\n",this->selectionMode());
	fprintf(writeFile,"Number of Clusters:\t%i\n",getNumberOfClusters());
	fprintf(writeFile,"Total Cells:\t%i\n",getNumberOfClusters());
	fprintf(writeFile,"Total Localizations:\t%i\n",file->localizationCount);
	fprintf(writeFile,"Active Cells:\t%i\n",this->getActiveCells());
	fprintf(writeFile,"Active Localizations:\t%i\n\n",this->getActiveDetections());

	// Inference Parameters
	fprintf(writeFile,"INFERENCE PARAMETERS\n");
	fprintf(writeFile,"Noise Sigma [nm]:\t%.1f\n",this->getSigma()*1000.0);
	fprintf(writeFile,"Maximum Neighbour Distance [nm]:\t%.1f\n",this->getMaximumNeighbourDistance()*1000.0);
	fprintf(writeFile,"Minimum Points per Cell:\t%i\n",this->getMinPointsPerCell());
	fprintf(writeFile,"Optimization Scheme:\t%i\n",this->getOptimizationMode());
	fprintf(writeFile,"Potential Calculated:\t%i\n",this->zonalPotentialsCalculated);
	fprintf(writeFile,"Prior Enabled:\t%i\n",this->smoothingPriorEnabled());
	if (this->optimizationMode == 4) {
		fprintf(writeFile,"Polynomial Order:\t%i\n",this->getPolynomialOrder());
		fprintf(writeFile,"Polynomial Coefficients:\t");
		for (int q = 0; q < this->getCoefficients()-1; q++) {
			fprintf(writeFile,"%f\t",this->getOptimizationArray(q));
		}
		fprintf(writeFile,"%f\n\n",this->getOptimizationArray(file->squareMesh->getCoefficients()-1));
	} else {
		fprintf(writeFile,"Polynomial Order:\tN/A\n");
		fprintf(writeFile,"Polynomial Coefficients:\tN/A\n\n");
	}

	// Data
	fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tFx\t\tFy\t\tV\n");
	for (int a = 0; a < getNumberOfClusters(); a++) {
			fprintf(writeFile,"%i\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					getCell(a)->active(),
					getCell(a)->getXCentre(),
					getCell(a)->getYCentre(),
					getCell(a)->getCount(),
					getCell(a)->getXCentroid(),
					getCell(a)->getYCentroid(),
					getCell(a)->getVariance(),
					getCell(a)->getDiffusion(),
					getCell(a)->getForceX(),
					getCell(a)->getForceY(),
					getCell(a)->getPotential()
					);

	}

	fclose(writeFile);

}

VoronoiMesh::~VoronoiMesh() {
	// delete Cells array
	if (file->inferred) {
		delete [] optimizationArray;
		delete [] potentialArray;
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}
	if (file->voronoiMeshOverlay) {
		delete [] cells;
		delete [] colorArray;
		delete [] xCentreList;
		delete [] yCentreList;
	}
}

/************* QUAD TREE MESHING *************/
TreeMesh::TreeMesh() {
	roQuadTrees = NULL;
	roIdentifiers = NULL;
	roZones = 0;
	roEnable = false;
	randomizedOptimizationGui = NULL;

	quadTreeLocalizationPointer = NULL;
	selectionEnable = false;
	//selection = NULL;
	smoothingPriorEnable = false;
	jeffreysPriorEnable = true;

	minCapacity = file->treeMeshGui->minCapacitySlider->value();

	// 3D view variables
	landscapeBorderTriangles = NULL;
	landscapeBorderNormals = NULL;
	landscapeBorderColors = NULL;
	landscapeBorderVertices = 0;
	landscapeTriangles = NULL;
	landscapeNormals = NULL;
	landscapeColors = NULL;
	landscapeVertices = 0;
	linesOverlay = NULL;
	quadsOverlay = NULL;
	quadsColor = NULL;
	leafSelected = false;
	leaf2Selected = false;
	identifierCounter = 0;
	maxQuadTreePower = -1;
	minCount = 1000000;
	maxCount = -1000000;
	quadTree = NULL;
	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dimensions = 0;
	this->xMax = file->xMax;
	this->xMin = file->xMin;
	this->yMax = file->yMax;
	this->yMin = file->yMin;
	this->activeDetections = 0;
	this->optimizationMode = 0;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	this->diffusionLogMax = -1000000.0;
	this->diffusionLogMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->potentialMin = 1000000.0;
	this->activeCells = 0;
	this->numberOfCoefficients = 0;
	this->dtMean = 0.0;
	this->hessian = NULL;
	this->allVariables = 0;

	// default values
	this->polynomialOrder = 2;
	this->slidingWindowEnable = false;
	this->sigma = 0.00; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentQuadZone = NULL;
	this->beta = 2.0;
	this->zonalPotentialsCalculated = false;
}

TreeMesh::TreeMesh(SelectionCell selection) {

	roQuadTrees = NULL;
	roIdentifiers = NULL;
	roZones = 0;
	roEnable = false;
	randomizedOptimizationGui = NULL;

	quadTreeLocalizationPointer = NULL;
	landscapeVertices = 0;
	selectionEnable = true;
	this->selection = selection;
	minCapacity = file->treeMeshGui->minCapacitySlider->value();

	// 3D view variables
	landscapeTriangles = NULL;
	landscapeNormals = NULL;
	landscapeColors = NULL;
	linesOverlay = NULL;
	quadsOverlay = NULL;
	quadsColor = NULL;
	leafSelected = false;
	identifierCounter = 0;
	maxQuadTreePower = -1;
	minCount = 1000000;
	maxCount = -1000000;
	quadTree = NULL;
	this->tMax = -1000000.0;
	this->tMin = 1000000.0;
	this->forceMax = -1000000.0;
	this->forceMin = 1000000.0;
	this->xMean = 0.0;
	this->yMean = 0.0;
	this->dimensions = 0;
	this->xMax = file->xMax;
	this->xMin = file->xMin;
	this->yMax = file->yMax;
	this->yMin = file->yMin;
	this->activeDetections = 0;
	this->optimizationMode = 0;
	this->diffusionMax = -1000000.0;
	this->diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	this->potentialMax = -1000000.0;
	this->potentialMin = 1000000.0;
	this->activeCells = 0;
	this->numberOfCoefficients = 0;
	this->dtMean = 0.0;
	this->hessian = NULL;
	this->allVariables = 0;

	// default values
	this->polynomialOrder = 2;
	this->slidingWindowEnable = false;
	this->sigma = 0.00; // 20 nm
	this->totalVariables = 0;
	this->optimizationArray = NULL;
	this->potentialArray = NULL;
	this->currentQuadZone = NULL;
	this->beta = 2.0;
	this->zonalPotentialsCalculated = false;

}

TreeMesh::~TreeMesh() {
	if (file->inferred) {
		delete [] optimizationArray;
		delete [] potentialArray;
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
		}
	}
}

void TreeMesh::infer() {

	iMAP->glWindow->deactivate();

	iMAP->startTime = clock();
	iMAP->stopCalculation = false;

	// assign a QuadTree pointer to each localization (to reduce calculation time)
	if (file->treeMesh->quadTreeLocalizationPointer != NULL) {
		delete [] quadTreeLocalizationPointer;
		quadTreeLocalizationPointer = NULL;
	}
	if (file->treeMesh->selectionMode()) {
		file->treeMesh->quadTreeLocalizationPointer = new QuadTree*[selection.count];
		for (int l = 0; l < selection.count; l++) {
			quadTreeLocalizationPointer[l] = NULL;
		}
		for (int l = 0; l < selection.count; l++) {
			bool found = false;
			getLeafSelection(selection.xPointer[l],file->treeMesh->selection.yPointer[l],quadTree,&found);
			quadTreeLocalizationPointer[l] = selectedQuadLeaf;
		}
	} else {
		file->treeMesh->quadTreeLocalizationPointer = new QuadTree*[file->localizationCount];
		for (int l = 0; l < file->localizationCount; l++) {
			quadTreeLocalizationPointer[l] = NULL;
		}
		for (int l = 0; l < file->localizationCount; l++) {
			bool found = false;
			getLeafSelection(file->xPointer[l],file->yPointer[l],quadTree,&found);
			quadTreeLocalizationPointer[l] = selectedQuadLeaf;
		}
	}

	this->minPointsPerCell = (int)file->treeMeshGui->minPointsSlider->value();
	
	file->optimizationMode = file->treeMeshGui->inferenceModeChoice->value();

	if (randomizedOptimizationGui != NULL) { randomizedOptimizationGui->hide(); }

	file->treeMeshGui->overlayDiffusionButton->value(1);
	file->treeMeshGui->overlayDiffusionButton->do_callback();

	switch(file->optimizationMode) {
		case 0: // D Inference
		{
			// no pausing in DF mode
			file->treeMeshGui->pauseButton->deactivate();

			if (file->treeMeshGui->smoothingPriorButton->value()) {
				file->treeMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationD();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDSmoothing();
					}
					inferDSmoothing(quadTree);
					postInferDSmoothing();
				}
			} else {
				preInferD();
				int num = 0;
				inferD(quadTree,&num);
				postInferD();
			}

			file->treeMeshGui->overlayPointNumberButton->value(0);
			file->treeMeshGui->overlayForceArrowsButton->deactivate();
			file->treeMeshGui->overlayForceArrowsButton->labelcolor(iMAP->bgColor);
			file->treeMeshGui->overlayForceArrowsButton->value(0);
			file->treeMeshGui->overlayDiffusionButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceMagnitudeButton->deactivate();
			file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(iMAP->bgColor);

			file->treeMeshGui->minPosteriorSlider->activate();
			file->treeMeshGui->maxPosteriorSlider->activate();
			file->treeMeshGui->fPosteriorButton->deactivate();
			file->treeMeshGui->fPosteriorButton->copy_label("Force");
			file->treeMeshGui->dPosteriorButton->activate();
			file->treeMeshGui->vPosteriorButton->deactivate();
			file->treeMeshGui->potentialReferenceButton->deactivate();
			file->treeMeshGui->posterioriSampleNumberSlider->activate();

			zonalPotentialsCalculated = false;

			file->treeMeshGui->overlayPotentialButton->deactivate();
			file->treeMeshGui->overlayPotentialButton->labelcolor(iMAP->bgColor);
			file->treeMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->treeMeshGui->deltaTSlider->deactivate();
			file->treeMeshGui->timeStepsSlider->deactivate();
			file->treeMeshGui->generateTrajectoriesButton->deactivate();

		}
			break;
		case 1: // DF Inference
		{
			// no pausing in DF mode
			file->treeMeshGui->pauseButton->deactivate();

			file->treeMeshGui->overlayForceArrowsButton->activate();
			file->treeMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceArrowsButton->value(1);

			file->treeMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->treeMeshGui->deltaTSlider->deactivate();
			file->treeMeshGui->timeStepsSlider->deactivate();
			file->treeMeshGui->generateTrajectoriesButton->deactivate();

			if (file->treeMeshGui->smoothingPriorButton->value()) {
				file->treeMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					if (inferRandomizedOptimizationDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->treeMeshGui->overlayPotentialButton->activate();
						file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->treeMeshGui->numberOfTrajectoriesSlider->activate();
						file->treeMeshGui->deltaTSlider->activate();
						file->treeMeshGui->timeStepsSlider->activate();
						file->treeMeshGui->generateTrajectoriesButton->activate();
					}
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDFSmoothing();
					}
					inferDFSmoothing(quadTree);
					zonalPotentialsCalculated = false;
					if (iMAP->stopCalculation == false) {
						if (postInferDFSmoothing()) {
							iMAP->updateDisplayButton->value(1);
							iMAP->updateDisplayButton->do_callback();
							inferPotentialsDF();
							zonalPotentialsCalculated = true;
							file->treeMeshGui->overlayPotentialButton->activate();
							file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
							file->treeMeshGui->numberOfTrajectoriesSlider->activate();
							file->treeMeshGui->deltaTSlider->activate();
							file->treeMeshGui->timeStepsSlider->activate();
							file->treeMeshGui->generateTrajectoriesButton->activate();
						}
					}
				}
			} else {
				preInferDF();
				int num = 0;
				inferDF(quadTree,&num);
				zonalPotentialsCalculated = false;
				if (iMAP->stopCalculation == false) {
					if (postInferDF()) {
						iMAP->updateDisplayButton->value(1);
						iMAP->updateDisplayButton->do_callback();
						inferPotentialsDF();
						zonalPotentialsCalculated = true;
						file->treeMeshGui->overlayPotentialButton->activate();
						file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
						file->treeMeshGui->numberOfTrajectoriesSlider->activate();
						file->treeMeshGui->deltaTSlider->activate();
						file->treeMeshGui->timeStepsSlider->activate();
						file->treeMeshGui->generateTrajectoriesButton->activate();
					}
				}
			}

			file->treeMeshGui->overlayDiffusionButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceMagnitudeButton->activate();
			file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->treeMeshGui->minPosteriorSlider->activate();
			file->treeMeshGui->maxPosteriorSlider->activate();
			file->treeMeshGui->fPosteriorButton->activate();
			file->treeMeshGui->fPosteriorButton->copy_label("Force");
			file->treeMeshGui->dPosteriorButton->activate();
			file->treeMeshGui->vPosteriorButton->deactivate();
			file->treeMeshGui->potentialReferenceButton->deactivate();
			file->treeMeshGui->posterioriSampleNumberSlider->activate();

		}
			break;
		case 2: // DDr Inference
		{
			// no pausing in DF mode
			file->treeMeshGui->pauseButton->deactivate();

			if (file->treeMeshGui->smoothingPriorButton->value()) {
				file->treeMeshGui->pauseButton->activate();
				if (roEnable) {
					iMAP->updateDisplayButton->value(1);
					iMAP->updateDisplayButton->do_callback();
					iMAP->updateDisplayButton->deactivate();
					if (randomizedOptimizationGui != NULL) {
						randomizedOptimizationGui->hide();
						delete randomizedOptimizationGui;
						randomizedOptimizationGui = NULL;
					}
					randomizedOptimizationGui = new RandomizedOptimizationGui();
					inferRandomizedOptimizationDDr();
					iMAP->updateDisplayButton->activate();
					randomizedOptimizationGui->finish();
				} else {
					if (iMAP->pauseCalculation == false) {
						preInferDDrSmoothing();
					}
					inferDDrSmoothing(quadTree);
					postInferDDrSmoothing();
				}
			} else {
				preInferDDr();
				int num = 0;
				inferDDr(quadTree,&num);
				postInferDDr();
			}

			file->treeMeshGui->overlayForceArrowsButton->activate();
			file->treeMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceArrowsButton->value(1);
			file->treeMeshGui->overlayDiffusionButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceMagnitudeButton->activate();
			file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);

			file->treeMeshGui->minPosteriorSlider->activate();
			file->treeMeshGui->maxPosteriorSlider->activate();
			file->treeMeshGui->fPosteriorButton->activate();
			file->treeMeshGui->fPosteriorButton->copy_label("Drift");
			file->treeMeshGui->dPosteriorButton->activate();
			file->treeMeshGui->vPosteriorButton->deactivate();
			file->treeMeshGui->potentialReferenceButton->deactivate();
			file->treeMeshGui->posterioriSampleNumberSlider->activate();

			// activate potential overlay buttons
			file->treeMeshGui->overlayPotentialButton->deactivate();
			file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->treeMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->treeMeshGui->deltaTSlider->deactivate();
			file->treeMeshGui->timeStepsSlider->deactivate();
			file->treeMeshGui->generateTrajectoriesButton->deactivate();
		}
			break;
		case 3: // DV Inference
			file->treeMeshGui->pauseButton->activate();

			if (roEnable) {
				iMAP->updateDisplayButton->value(1);
				iMAP->updateDisplayButton->do_callback();
				iMAP->updateDisplayButton->deactivate();
				if (randomizedOptimizationGui != NULL) {
					randomizedOptimizationGui->hide();
					delete randomizedOptimizationGui;
					randomizedOptimizationGui = NULL;
				}
				randomizedOptimizationGui = new RandomizedOptimizationGui();
				inferRandomizedOptimizationDV();
				iMAP->updateDisplayButton->activate();
				randomizedOptimizationGui->finish();
			} else {
				if (iMAP->pauseCalculation == false) {
					preInferDV();
				}
				inferDV(quadTree);
				postInferDV();
			}

			file->treeMeshGui->overlayForceArrowsButton->activate();
			file->treeMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceArrowsButton->value(1);
			file->treeMeshGui->overlayDiffusionButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceMagnitudeButton->activate();
			file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayPotentialButton->activate();
			file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);

			file->treeMeshGui->minPosteriorSlider->activate();
			file->treeMeshGui->maxPosteriorSlider->activate();
			file->treeMeshGui->fPosteriorButton->deactivate();
			file->treeMeshGui->fPosteriorButton->copy_label("Force");
			file->treeMeshGui->dPosteriorButton->activate();
			file->treeMeshGui->vPosteriorButton->activate();
			file->treeMeshGui->potentialReferenceButton->activate();
			file->treeMeshGui->posterioriSampleNumberSlider->activate();

			file->treeMeshGui->numberOfTrajectoriesSlider->deactivate();
			file->treeMeshGui->deltaTSlider->deactivate();
			file->treeMeshGui->timeStepsSlider->deactivate();
			file->treeMeshGui->generateTrajectoriesButton->deactivate();
			break;

		case 4: // Polynomial Potential
			// update data structures for optimization run
			if (iMAP->pauseCalculation == false) {
				preInferPolynomial();
			}
			// find minimum of function
			inferPolynomial();
			// save forces and diffusion values to activated cells;
			postInferPolynomial();
			file->treeMeshGui->overlayForceArrowsButton->activate();
			file->treeMeshGui->overlayForceArrowsButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceArrowsButton->value(1);
			file->treeMeshGui->overlayDiffusionButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->activate();
			file->treeMeshGui->overlayDiffusionLogButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayDiffusionButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayForceMagnitudeButton->activate();
			file->treeMeshGui->overlayForceMagnitudeButton->labelcolor(FL_WHITE);
			file->treeMeshGui->overlayPotentialButton->activate();
			file->treeMeshGui->overlayPotentialButton->labelcolor(FL_WHITE);
			file->treeMeshGui->numberOfTrajectoriesSlider->activate();
			file->treeMeshGui->deltaTSlider->activate();
			file->treeMeshGui->timeStepsSlider->activate();
			file->treeMeshGui->generateTrajectoriesButton->activate();
			file->treeMeshGui->fPosteriorButton->activate();
			file->treeMeshGui->fPosteriorButton->copy_label("Force");
			file->treeMeshGui->dPosteriorButton->activate();
			file->treeMeshGui->vPosteriorButton->deactivate();
			file->treeMeshGui->potentialReferenceButton->deactivate();
			break;
	}

	iMAP->endTime = clock();
	char label[100];
	const float dif = ((float)iMAP->endTime-(float)iMAP->startTime)/CLOCKS_PER_SEC;
	textDisplayUpdate("\n");
	sprintf(label,"Calculation Time: %.2f [s]\n",dif);
	textDisplayUpdate(label);

	iMAP->glWindow->activate();

	if (iMAP->pauseCalculation == true) {
		file->treeMeshGui->pauseButton->activate();
		file->treeMeshGui->stopButton->activate();
		file->treeMeshGui->inferButton->deactivate();
		file->treeMeshGui->resetButton->deactivate();
	}

	// initialize posteriori selections
	int l = 0;
	initializePosterioriSelectionZones(quadTree,&l);

	Fl::redraw();
}

void TreeMesh::initializePosterioriSelectionZones(QuadTree* tree, int *sel) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		for (int i = 0; i < tree->getCount(); i++) {
			if (tree->getIndex(i) < file->localizationCount-1) {
				if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
					xMean += fabs(file->xPointer[tree->getIndex(i)+1]-file->xPointer[tree->getIndex(i)]);
					yMean += fabs(file->yPointer[tree->getIndex(i)+1]-file->yPointer[tree->getIndex(i)]);
					dtMean += file->tPointer[tree->getIndex(i)+1]-file->tPointer[tree->getIndex(i)];
					activeDetections++;
				}
			}
		}
		return;
	}
	initializeTreeMesh(tree->nw);
	initializeTreeMesh(tree->ne);
	initializeTreeMesh(tree->sw);
	initializeTreeMesh(tree->se);
}

void TreeMesh::inferPolynomial() {
	if (this->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorQuadTreeSelection, (dfunc));
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, polynomialPosteriorQuadTree, (dfunc));
	}
}

void TreeMesh::preInferDF() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void TreeMesh::preInferPolynomial() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	polynomialOrder = file->treeMeshGui->getPolynomialOrder();
	numberOfCoefficients = (int) floor((double)(polynomialOrder+1)*(double)(polynomialOrder+2)/2.0-1.0);

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	totalVariables = 0;
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	identifierCounter = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);
	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = totalVariables + numberOfCoefficients; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }

	optimizationArray[2] = K_R/(4.e-9);
	optimizationArray[4] = K_R/(4.e-9);

	// initial values for optimization array
	for (int i = numberOfCoefficients; i < dimensions; i++) { optimizationArray[i] = 1.0/2.0*(D_eff_x + D_eff_y); }
}

void TreeMesh::postInferPolynomial() {

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	findExtremesPolynomial(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}
	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

}

void TreeMesh::offsetPotentials() {
	adjustPotentials(file->treeMesh->quadTree);
	potentialMax -= potentialMin;
	potentialMin = 0.0;
}

void TreeMesh::findExtremesD() {
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	findExtremesD(quadTree);
	iMAP->overlayAdjustment = true;
}

void TreeMesh::findExtremesDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	findExtremesDF(quadTree);
	iMAP->overlayAdjustment = true;
}

void TreeMesh::findExtremesDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	findExtremesDF(quadTree);
	iMAP->overlayAdjustment = true;
}

void TreeMesh::findExtremesDV() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;
	findExtremesDV(quadTree);
	iMAP->overlayAdjustment = true;
}

void TreeMesh::findExtremesPolynomial(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->xForce = polynomialFxValueQuadTree(optimizationArray,
							tree->getXCentre(),
							tree->getYCentre() );
			tree->yForce = polynomialFyValueQuadTree(optimizationArray,
							tree->getXCentre(),
							tree->getYCentre() );
			tree->forceMagnitude = sqrt(tree->xForce*tree->xForce+tree->yForce*tree->yForce);
			tree->diffusion = polynomialDValueQuadTree(optimizationArray,tree);
			tree->diffusionLog = log(polynomialDValueQuadTree(optimizationArray,tree));
			tree->potential = polynomialEpValueQuadTree(optimizationArray,
							tree->getXCentre(),
							tree->getYCentre() );
			if (forceMax < tree->forceMagnitude) { forceMax = tree->forceMagnitude; }
			if (forceMin > tree->forceMagnitude) { forceMin = tree->forceMagnitude; }
			if (diffusionMax < tree->diffusion) { diffusionMax = tree->diffusion; }
			if (diffusionMin > tree->diffusion) { diffusionMin = tree->diffusion; }
			if (diffusionLogMax < tree->diffusionLog) { diffusionLogMax = tree->diffusionLog; }
			if (diffusionLogMin > tree->diffusionLog) { diffusionLogMin = tree->diffusionLog; }
			if (potentialMax < tree->potential) { potentialMax = tree->potential; }
			if (potentialMin > tree->potential) { potentialMin = tree->potential; }
		} else {
			tree->xForce = 0.0;
			tree->yForce = 0.0;
			tree->forceMagnitude = 0.0;
			tree->diffusion = 0.0;
			tree->diffusionLog = 0.0;
			tree->potential = 0.0;
		}
		return;
	}
	findExtremesPolynomial(tree->nw);
	findExtremesPolynomial(tree->ne);
	findExtremesPolynomial(tree->sw);
	findExtremesPolynomial(tree->se);
	return;
}

void TreeMesh::initializeTreeMesh(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		for (int i = 0; i < tree->getCount(); i++) {
			if (tree->getIndex(i) < file->localizationCount-1) {
				if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
					xMean += fabs(file->xPointer[tree->getIndex(i)+1]-file->xPointer[tree->getIndex(i)]);
					yMean += fabs(file->yPointer[tree->getIndex(i)+1]-file->yPointer[tree->getIndex(i)]);
					dtMean += file->tPointer[tree->getIndex(i)+1]-file->tPointer[tree->getIndex(i)];
					activeDetections++;
					tree->resetPriors();

				}
			}
		}
		return;
	}
	initializeTreeMesh(tree->nw);
	initializeTreeMesh(tree->ne);
	initializeTreeMesh(tree->sw);
	initializeTreeMesh(tree->se);
}

void TreeMesh::initializeTreeMesh(QuadTree* tree,SelectionCell selection) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		for (int i = 0; i < tree->getCount(); i++) {
			if (tree->getIndex(i) < file->localizationCount-1) {
				if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
					xMean += fabs(file->treeMesh->selection.xPointer[tree->getIndex(i)+1]-file->treeMesh->selection.xPointer[tree->getIndex(i)]);
					yMean += fabs(file->treeMesh->selection.yPointer[tree->getIndex(i)+1]-file->treeMesh->selection.yPointer[tree->getIndex(i)]);
					dtMean += file->treeMesh->selection.tPointer[tree->getIndex(i)+1]-file->treeMesh->selection.tPointer[tree->getIndex(i)];
					activeDetections++;
				}
			}
		}
		return;
	}
	initializeTreeMesh(tree->nw);
	initializeTreeMesh(tree->ne);
	initializeTreeMesh(tree->sw);
	initializeTreeMesh(tree->se);
}


void TreeMesh::preInferD() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 1; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void TreeMesh::inferD(QuadTree* tree, int *clusterNumber) {
	// only infer in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			// set current zone for posteriori function
			setCurrentZone(tree);
			// average effective diffusion
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorQuadTree, (dfunc));

			*clusterNumber = *clusterNumber + 1;

			tree->setD(optimizationArray[0]);
			tree->setDlog(log(optimizationArray[0]));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Quad-Tree Meshing\n1   (D) Inference\n2\n3   Zone %i / %i\n",*clusterNumber,this->getVariables());
				textDisplayUpdate(progress);
				findExtremesD();
			}

			if (iMAP->stopCalculation) { return; }

		} else {
			tree->setFx(0.0);
			tree->setFy(0.0);
			tree->setD(0.0);
			tree->setDlog(0.0);
			tree->setFmag(0.0);
		}
		return;
	}

	if (iMAP->stopCalculation) { return; }

	inferD(tree->nw,clusterNumber);
	inferD(tree->ne,clusterNumber);
	inferD(tree->sw,clusterNumber);
	inferD(tree->se,clusterNumber);

	return;
}

int TreeMesh::postInferD() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// call recursive function to traverse through tree
	findExtremesD(quadTree);

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return 1;
}

void TreeMesh::preInferDSmoothing() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// find neighbours
	resetNeighbourCount(quadTree);
	findNeighbours(quadTree);
	assignNeighbours(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	identifierCounter = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = totalVariables; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) {
		optimizationArray[i] = 1.0/2.0*(D_eff_x + D_eff_y);
		getTree(this->quadTree,i); // assigns leaf with current identity to identifierLeaf
//		optimizationArray[2*i+1] = ( -log((double)identifiedQuadLeaf->getCount()/(double)maxCount) );
	}

	// assign areas for smoothing prior
	if (file->treeMeshGui->smoothingPriorButton->value()) {
		assignAreas(quadTree);
	}
}

void TreeMesh::inferDSmoothing(QuadTree* tree) {
	char label[50];
	if (file->treeMesh->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingQuadTreeSelection, (dfunc));
		sprintf(label,"Cost =  %f\n",dPosteriorSmoothingQuadTreeSelection(optimizationArray));
		textDisplayUpdate(label);
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingQuadTree, (dfunc));
		sprintf(label,"Cost =  %f\n",dPosteriorSmoothingQuadTree(optimizationArray));
		textDisplayUpdate(label);
	}
}

int TreeMesh::postInferDSmoothing() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	findExtremesDSmoothing(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	return 1;
}

void TreeMesh::preInferDDr() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	// length of optimization array
	dimensions = 3; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) { optimizationArray[i] = 0.0; }
}

void TreeMesh::inferDDr(QuadTree* tree,int *clusterNumber) {
	// only infer in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			// set current zone for posteriori function
			setCurrentZone(tree);
			// average effective diffusion
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.0;
			optimizationArray[1] = 0.0;
			optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorQuadTree, (dfunc));

			*clusterNumber = *clusterNumber + 1;

			tree->setFx(optimizationArray[0]);
			tree->setFy(optimizationArray[1]);
			tree->setD(optimizationArray[2]);
			tree->setDlog(log(optimizationArray[2]));
			tree->setFmag(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Quad-Tree Meshing\n1   (D,Drift) Inference\n2\n3   Zone %i / %i\n",*clusterNumber,this->getVariables());
				textDisplayUpdate(progress);
				findExtremesDDr();
			}

			if (iMAP->stopCalculation) { return; }

		} else {
			tree->setFx(0.0);
			tree->setFy(0.0);
			tree->setD(0.0);
			tree->setDlog(0.0);
			tree->setFmag(0.0);
		}
		return;
	}

	if (iMAP->stopCalculation) { return; }

	inferDDr(tree->nw,clusterNumber);
	inferDDr(tree->ne,clusterNumber);
	inferDDr(tree->sw,clusterNumber);
	inferDDr(tree->se,clusterNumber);

	return;
}

int TreeMesh::postInferDDr() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// call recursive function to traverse through tree
	findExtremesDF(quadTree);

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return 1;
}

void TreeMesh::inferDF(QuadTree* tree,int *clusterNumber) {
	// only infer in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			// set current zone for posteriori function
			setCurrentZone(tree);
			// average effective diffusion
			const double D_eff_x = file->averageDx*file->averageDx/file->averageDt;
			const double D_eff_y = file->averageDy*file->averageDy/file->averageDt;

			optimizationArray[0] = 0.0;
			optimizationArray[1] = 0.0;
			optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorQuadTree, (dfunc));

			*clusterNumber = *clusterNumber + 1;

			tree->setFx(optimizationArray[0]);
			tree->setFy(optimizationArray[1]);
			tree->setD(optimizationArray[2]);
			tree->setDlog(log(optimizationArray[2]));
			tree->setFmag(sqrt(optimizationArray[0]*optimizationArray[0]+optimizationArray[1]*optimizationArray[1]));

			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[100];
				sprintf(progress,"Quad-Tree Meshing\n1   (D,F) Inference\n2\n3   Zone %i / %i\n",*clusterNumber,this->getVariables());
				textDisplayUpdate(progress);
				findExtremesDF();
			}

			if (iMAP->stopCalculation) { return; }

		} else {
			tree->setFx(0.0);
			tree->setFy(0.0);
			tree->setD(0.0);
			tree->setDlog(0.0);
			tree->setFmag(0.0);
		}
		return;
	}

	if (iMAP->stopCalculation) { return; }

	inferDF(tree->nw,clusterNumber);
	inferDF(tree->ne,clusterNumber);
	inferDF(tree->sw,clusterNumber);
	inferDF(tree->se,clusterNumber);

	return;
}

int TreeMesh::postInferDF() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	// call recursive function to traverse through tree
	findExtremesDF(quadTree);

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (hessian != NULL) {
		for (int d = 0; d < dimensions; d++) {
			delete [] hessian[d];
		}
		delete [] hessian;
		hessian = NULL;
	}

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

void TreeMesh::preInferDFSmoothing() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// find neighbours
	resetNeighbourCount(quadTree);
	findNeighbours(quadTree);
	assignNeighbours(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	identifierCounter = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = 3*totalVariables; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) {
		optimizationArray[3*i] = 0.5*(D_eff_x + D_eff_y);
		optimizationArray[3*i+1] = 0.0;
		optimizationArray[3*i+2] = 0.0;
		getTree(this->quadTree,i); // assigns leaf with current identity to identifierLeaf
//		optimizationArray[2*i+1] = ( -log((double)identifiedQuadLeaf->getCount()/(double)maxCount) );
	}

	// assign areas for smoothing prior
	if (file->treeMeshGui->smoothingPriorButton->value()) {
		assignAreas(quadTree);
	}
}

void TreeMesh::preInferDDrSmoothing() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// find neighbours
	resetNeighbourCount(quadTree);
	findNeighbours(quadTree);
	assignNeighbours(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	identifierCounter = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = 3*totalVariables; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions; i++) {
		optimizationArray[3*i] = 0.5*(D_eff_x + D_eff_y);
		optimizationArray[3*i+1] = 0.0;
		optimizationArray[3*i+2] = 0.0;
		getTree(this->quadTree,i); // assigns leaf with current identity to identifierLeaf
//		optimizationArray[2*i+1] = ( -log((double)identifiedQuadLeaf->getCount()/(double)maxCount) );
	}

	// assign areas for smoothing prior
	if (file->treeMeshGui->smoothingPriorButton->value()) {
		assignAreas(quadTree);
	}
}

void TreeMesh::inferDFSmoothing(QuadTree* tree) {
	char label[50];
	if (file->treeMesh->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingQuadTreeSelection, (dfunc));
		sprintf(label,"Cost = %f\n",dfPosteriorSmoothingQuadTreeSelection(optimizationArray));
		textDisplayUpdate(label);
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingQuadTree, (dfunc));
		sprintf(label,"Cost = %f\n",dfPosteriorSmoothingQuadTree(optimizationArray));
		textDisplayUpdate(label);
	}
}

void TreeMesh::inferDDrSmoothing(QuadTree* tree) {
	char label[50];
	if (file->treeMesh->selectionMode()) {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingQuadTreeSelection, (dfunc));
		sprintf(label,"Cost = %f\n",ddrPosteriorSmoothingQuadTreeSelection(optimizationArray));
		textDisplayUpdate(label);
	} else {
		dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingQuadTree, (dfunc));
		sprintf(label,"Cost = %f\n",ddrPosteriorSmoothingQuadTree(optimizationArray));
		textDisplayUpdate(label);
	}
}

int TreeMesh::postInferDFSmoothing() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	findExtremesDFSmoothing(quadTree);

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	return fl_choice("Compute potentials?","No","Yes",NULL);
}

int TreeMesh::postInferDDrSmoothing() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	findExtremesDDrSmoothing(quadTree);

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	return 1;
}

void TreeMesh::inferPotentialsDF() {

	iMAP->pauseCalculation = false;

	// do not reinitialize if pause button was pressed
	if (iMAP->stopCalculation == false) {
		dimensions = totalVariables;
		if (potentialArray != NULL) {
			delete [] potentialArray;
			potentialArray = NULL;
		}
		potentialArray = new double[dimensions];

//		fprintf(stder,"total variables = %i\n",totalVariables);

		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// find neighbours
		resetNeighbourCount(quadTree);
		findNeighbours(quadTree);
		assignNeighbours(quadTree);

		// initalize potential values in each leaf
		identifierCounter = 0;
		initializePotentials(quadTree);
	}
	iMAP->stopCalculation = false;

	// perform optimization of potential array

	iMAP->updatePotentialCalculation = true;
	dfpmin(potentialArray,dimensions,GTOL,iterations,fret,squareDifferenceQuadTree,dfunc);
	iMAP->updatePotentialCalculation = false;

	potentialMax = -1000000.0;
	potentialMin = 1000000.0;

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// assign potential values to each leaf
	assignPotentials(quadTree);

	offsetPotentials();
}

void TreeMesh::initializePotentials(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->setEp( -log((double)tree->getCount()/(double)maxCount) );
			tree->identifier = identifierCounter;
			this->potentialArray[identifierCounter] = tree->potential;
//			printf("%i\t%f\n",identifierCounter,tree->potential);
			identifierCounter++;
		} else {
			tree->setEp(0.0);
		}

		return;
	}
	initializePotentials(tree->nw);
	initializePotentials(tree->ne);
	initializePotentials(tree->sw);
	initializePotentials(tree->se);

	return;
}

void TreeMesh::resetNeighbourCount(QuadTree* tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		tree->nTopNeighbours = 0;
		tree->nLeftNeighbours = 0;
		tree->nBottomNeighbours = 0;
		tree->nRightNeighbours = 0;
		tree->nTopNeighboursLandscape = 0;
		tree->nLeftNeighboursLandscape = 0;
		tree->nBottomNeighboursLandscape = 0;
		tree->nRightNeighboursLandscape = 0;
		tree->topProgress = 0;
		tree->bottomProgress = 0;
		tree->rightProgress = 0;
		tree->leftProgress = 0;
		tree->topProgressLandscape = 0;
		tree->bottomProgressLandscape = 0;
		tree->rightProgressLandscape = 0;
		tree->leftProgressLandscape = 0;
		if (tree->topNeighbours != NULL) {
			delete [] tree->topNeighbours;
			tree->topNeighbours = NULL;
		}
		if (tree->bottomNeighbours != NULL) {
			delete [] tree->bottomNeighbours;
			tree->bottomNeighbours = NULL;
		}
		if (tree->leftNeighbours != NULL) {
			delete [] tree->leftNeighbours;
			tree->leftNeighbours = NULL;
		}
		if (tree->rightNeighbours != NULL) {
			delete [] tree->rightNeighbours;
			tree->rightNeighbours = NULL;
		}
		if (tree->topNeighboursLandscape != NULL) {
			delete [] tree->topNeighboursLandscape;
			tree->topNeighboursLandscape = NULL;
		}
		if (tree->bottomNeighboursLandscape != NULL) {
			delete [] tree->bottomNeighboursLandscape;
			tree->bottomNeighboursLandscape = NULL;
		}
		if (tree->leftNeighboursLandscape != NULL) {
			delete [] tree->leftNeighboursLandscape;
			tree->leftNeighboursLandscape = NULL;
		}
		if (tree->rightNeighboursLandscape != NULL) {
			delete [] tree->rightNeighboursLandscape;
			tree->rightNeighboursLandscape = NULL;
		}
		return;
	}
	resetNeighbourCount(tree->nw);
	resetNeighbourCount(tree->ne);
	resetNeighbourCount(tree->sw);
	resetNeighbourCount(tree->se);
	return;
}

void TreeMesh::findNeighbours(QuadTree* tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		findNeighbours(quadTree,tree);
		return;
	}
	findNeighbours(tree->nw);
	findNeighbours(tree->ne);
	findNeighbours(tree->sw);
	findNeighbours(tree->se);
	return;
}

void TreeMesh::findNeighbours(QuadTree* tree, QuadTree* root) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		const double distance = sqrt((tree->xAverage-root->xAverage)*(tree->xAverage-root->xAverage)+(tree->yAverage-root->yAverage)*(tree->yAverage-root->yAverage));
		const double maxDistance = file->treeMeshGui->neighbourDistanceSlider->value()/1000.0;
		// find top neighbours
		if (fabs(tree->getYMin() - root->getYMax()) < 0.00001) {
			if (root->power > tree->power) {
				if (root->getXCentre() >= tree->getXMin() && root->getXCentre() <= tree->getXMax()) {
					if (distance < maxDistance && tree->active()) { root->nTopNeighbours++; }
					root->nTopNeighboursLandscape++;
				}
			} else if (root->power < tree->power) {
				if (tree->getXCentre() >= root->getXMin() && tree->getXCentre() <= root->getXMax()) {
					if (distance < maxDistance && tree->active()) { root->nTopNeighbours++; }
					root->nTopNeighboursLandscape++;
				}
			} else if (root->power == tree->power) {
				if (fabs(tree->getXMin() - root->getXMin()) < 0.00001 && fabs(tree->getXMax() - root->getXMax()) < 0.00001) {
					if (distance < maxDistance && tree->active()) { root->nTopNeighbours++; }
					root->nTopNeighboursLandscape++;
				}
			}
		}
		// find bottom neighbours
		if (fabs(tree->getYMax() - root->getYMin()) < 0.00001) {
			if (root->power > tree->power) {
				if (root->getXCentre() >= tree->getXMin() && root->getXCentre() <= tree->getXMax()) {
					if (distance < maxDistance && tree->active()) { root->nBottomNeighbours++; }
					root->nBottomNeighboursLandscape++;
				}
			} else if (root->power < tree->power) {
				if (tree->getXCentre() >= root->getXMin() && tree->getXCentre() <= root->getXMax()) {
					if (distance < maxDistance && tree->active()) { root->nBottomNeighbours++; }
					root->nBottomNeighboursLandscape++;
				}
			} else if (root->power == tree->power) {
				if (fabs(tree->getXMin() - root->getXMin()) < 0.00001 && fabs(tree->getXMax() - root->getXMax()) < 0.00001) {
					if (distance < maxDistance && tree->active()) { root->nBottomNeighbours++; }
					root->nBottomNeighboursLandscape++;
				}
			}
		}
		// find right neighbours
		if (fabs(tree->getXMin() - root->getXMax()) < 0.00001) {
			if (root->power > tree->power) {
				if (root->getYCentre() >= tree->getYMin() && root->getYCentre() <= tree->getYMax()) {
					if (distance < maxDistance && tree->active()) { root->nRightNeighbours++; }
					root->nRightNeighboursLandscape++;
				}
			} else if (root->power < tree->power) {
				if (tree->getYCentre() >= root->getYMin() && tree->getYCentre() <= root->getYMax()) {
					if (distance < maxDistance && tree->active()) { root->nRightNeighbours++; }
					root->nRightNeighboursLandscape++;
				}
			} else if (root->power == tree->power) {
				if (fabs(tree->getYMin() - root->getYMin()) < 0.00001 && fabs(tree->getYMax() - root->getYMax()) < 0.00001) {
					if (distance < maxDistance && tree->active()) { root->nRightNeighbours++; }
					root->nRightNeighboursLandscape++;
				}
			}
		}
		// find left neighbours
		if (fabs(tree->getXMax() - root->getXMin()) < 0.00001) {
			if (root->power > tree->power) {
				if (root->getYCentre() >= tree->getYMin() && root->getYCentre() <= tree->getYMax()) {
					if (distance < maxDistance && tree->active()) { root->nLeftNeighbours++; }
					root->nLeftNeighboursLandscape++;
				}
			} else if (root->power < tree->power) {
				if (tree->getYCentre() >= root->getYMin() && tree->getYCentre() <= root->getYMax()) {
					if (distance < maxDistance && tree->active()) { root->nLeftNeighbours++; }
					root->nLeftNeighboursLandscape++;
				}
			} else if (root->power == tree->power) {
				if (fabs(tree->getYMin() - root->getYMin()) < 0.00001 && fabs(tree->getYMax() - root->getYMax()) < 0.00001) {
					if (distance < maxDistance && tree->active()) { root->nLeftNeighbours++; }
					root->nLeftNeighboursLandscape++;
				}
			}
		}
		return;
	}
	findNeighbours(tree->nw,root);
	findNeighbours(tree->ne,root);
	findNeighbours(tree->sw,root);
	findNeighbours(tree->se,root);
	return;
}


void TreeMesh::assignNeighbours(QuadTree* tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->nTopNeighbours > 0) {
			tree->topProgress = 0;
			tree->topNeighbours = new QuadTree*[tree->nTopNeighbours];
		}
		if (tree->nBottomNeighbours > 0) {
			tree->bottomProgress = 0;
			tree->bottomNeighbours = new QuadTree*[tree->nBottomNeighbours];
		}
		if (tree->nLeftNeighbours > 0) {
			tree->leftProgress = 0;
			tree->leftNeighbours = new QuadTree*[tree->nLeftNeighbours];
		}
		if (tree->nRightNeighbours > 0) {
			tree->rightProgress = 0;
			tree->rightNeighbours = new QuadTree*[tree->nRightNeighbours];
		}
		if (tree->nTopNeighboursLandscape > 0) {
			tree->topProgressLandscape = 0;
			tree->topNeighboursLandscape = new QuadTree*[tree->nTopNeighboursLandscape];
		}
		if (tree->nBottomNeighboursLandscape > 0) {
			tree->bottomProgressLandscape = 0;
			tree->bottomNeighboursLandscape = new QuadTree*[tree->nBottomNeighboursLandscape];
		}
		if (tree->nLeftNeighboursLandscape > 0) {
			tree->leftProgressLandscape = 0;
			tree->leftNeighboursLandscape = new QuadTree*[tree->nLeftNeighboursLandscape];
		}
		if (tree->nRightNeighboursLandscape > 0) {
			tree->rightProgressLandscape = 0;
			tree->rightNeighboursLandscape = new QuadTree*[tree->nRightNeighboursLandscape];
		}
		assignNeighbours(quadTree,tree);
		return;
	}
	assignNeighbours(tree->nw);
	assignNeighbours(tree->ne);
	assignNeighbours(tree->sw);
	assignNeighbours(tree->se);
	return;
}

void TreeMesh::assignAreas(QuadTree* tree) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->xArea = 0.0;
			tree->yArea = 0.0;
			for (int l = 0; l < tree->nLeftNeighboursLandscape; l++) {
				tree->xArea += tree->getLeftNeighbourLandscape(l)->area;
			}
			for (int r = 0; r < tree->nRightNeighboursLandscape; r++) {
				tree->xArea += tree->getRightNeighbourLandscape(r)->area;
			}
			for (int t = 0; t < tree->nTopNeighboursLandscape; t++) {
				tree->yArea += tree->getTopNeighbourLandscape(t)->area;
			}
			for (int b = 0; b < tree->nBottomNeighboursLandscape; b++) {
				tree->yArea += tree->getBottomNeighbourLandscape(b)->area;
			}
		}
		return;
	}
	assignAreas(tree->nw);
	assignAreas(tree->ne);
	assignAreas(tree->sw);
	assignAreas(tree->se);
	return;

}

void TreeMesh::roDeactivateZones(QuadTree* tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
//		if (tree->sgdActive()) { sgdTest++; }
		tree->roDeactivate();
		tree->setRoIdentifier(-1);
		return;
	}
	roDeactivateZones(tree->nw);
	roDeactivateZones(tree->ne);
	roDeactivateZones(tree->sw);
	roDeactivateZones(tree->se);
	return;
}


void TreeMesh::assignNeighbours(QuadTree* tree, QuadTree* root) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
			const double distance = sqrt((tree->xAverage-root->xAverage)*(tree->xAverage-root->xAverage)+(tree->yAverage-root->yAverage)*(tree->yAverage-root->yAverage));
			const double maxDistance = file->treeMeshGui->neighbourDistanceSlider->value()/1000.0;
			// find top neighbours
			if (fabs(tree->getYMin() - root->getYMax()) < 0.00001) {
				if (root->power > tree->power) {
					if (root->getXCentre() >= tree->getXMin() && root->getXCentre() <= tree->getXMax()) {
						if (distance < maxDistance && tree->active()) {
							root->topNeighbours[root->topProgress] = tree;
							root->topProgress++;
						}
						root->topNeighboursLandscape[root->topProgressLandscape] = tree;
						root->topProgressLandscape++;
					}
				} else if (root->power < tree->power) {
					if (tree->getXCentre() >= root->getXMin() && tree->getXCentre() <= root->getXMax()) {
						if (distance < maxDistance && tree->active()) {
							root->topNeighbours[root->topProgress] = tree;
							root->topProgress++;
						}
						root->topNeighboursLandscape[root->topProgressLandscape] = tree;
						root->topProgressLandscape++;
					}
				} else if (root->power == tree->power) {
					if (fabs(tree->getXMin() - root->getXMin()) < 0.00001 && fabs(tree->getXMax() - root->getXMax()) < 0.00001) {
						if (distance < maxDistance && tree->active()) {
							root->topNeighbours[root->topProgress] = tree;
							root->topProgress++;
						}
						root->topNeighboursLandscape[root->topProgressLandscape] = tree;
						root->topProgressLandscape++;
					}
				}
			}
			// find bottom neighbours
			if (fabs(tree->getYMax() - root->getYMin()) < 0.00001) {
				if (root->power > tree->power) {
					if (root->getXCentre() >= tree->getXMin() && root->getXCentre() <= tree->getXMax()) {
						if (distance < maxDistance && tree->active()) {
							root->bottomNeighbours[root->bottomProgress] = tree;
							root->bottomProgress++;
						}
						root->bottomNeighboursLandscape[root->bottomProgressLandscape] = tree;
						root->bottomProgressLandscape++;
					}
				} else if (root->power < tree->power) {
					if (tree->getXCentre() >= root->getXMin() && tree->getXCentre() <= root->getXMax()) {
						if (distance < maxDistance && tree->active()) {
							root->bottomNeighbours[root->bottomProgress] = tree;
							root->bottomProgress++;
						}
						root->bottomNeighboursLandscape[root->bottomProgressLandscape] = tree;
						root->bottomProgressLandscape++;
					}
				} else if (root->power == tree->power) {
					if (fabs(tree->getXMin() - root->getXMin()) < 0.00001 && fabs(tree->getXMax() - root->getXMax()) < 0.00001) {
						if (distance < maxDistance && tree->active()) {
							root->bottomNeighbours[root->bottomProgress] = tree;
							root->bottomProgress++;
						}
						root->bottomNeighboursLandscape[root->bottomProgressLandscape] = tree;
						root->bottomProgressLandscape++;
					}
				}
			}
			// find right neighbours
			if (fabs(tree->getXMin() - root->getXMax()) < 0.00001) {
				if (root->power > tree->power) {
					if (root->getYCentre() >= tree->getYMin() && root->getYCentre() <= tree->getYMax()) {
						if (distance < maxDistance && tree->active()) {
							root->rightNeighbours[root->rightProgress] = tree;
							root->rightProgress++;
						}
						root->rightNeighboursLandscape[root->rightProgressLandscape] = tree;
						root->rightProgressLandscape++;
					}
				} else if (root->power < tree->power) {
					if (tree->getYCentre() >= root->getYMin() && tree->getYCentre() <= root->getYMax()) {
						if (distance < maxDistance && tree->active()) {
							root->rightNeighbours[root->rightProgress] = tree;
							root->rightProgress++;
						}
						root->rightNeighboursLandscape[root->rightProgressLandscape] = tree;
						root->rightProgressLandscape++;
					}
				} else if (root->power == tree->power) {
					if (fabs(tree->getYMin() - root->getYMin()) < 0.00001 && fabs(tree->getYMax() - root->getYMax()) < 0.00001) {
						if (distance < maxDistance && tree->active()) {
							root->rightNeighbours[root->rightProgress] = tree;
							root->rightProgress++;
						}
						root->rightNeighboursLandscape[root->rightProgressLandscape] = tree;
						root->rightProgressLandscape++;
					}
				}
			}
			// find left neighbours
			if (fabs(tree->getXMax() - root->getXMin()) < 0.00001) {
				if (root->power > tree->power) {
					if (root->getYCentre() >= tree->getYMin() && root->getYCentre() <= tree->getYMax()) {
						if (distance < maxDistance && tree->active()) {
							root->leftNeighbours[root->leftProgress] = tree;
							root->leftProgress++;
						}
						root->leftNeighboursLandscape[root->leftProgressLandscape] = tree;
						root->leftProgressLandscape++;
					}
				} else if (root->power < tree->power) {
					if (tree->getYCentre() >= root->getYMin() && tree->getYCentre() <= root->getYMax()) {
						if (distance < maxDistance && tree->active()) {
							root->leftNeighbours[root->leftProgress] = tree;
							root->leftProgress++;
						}
						root->leftNeighboursLandscape[root->leftProgressLandscape] = tree;
						root->leftProgressLandscape++;
					}
				} else if (root->power == tree->power) {
					if (fabs(tree->getYMin() - root->getYMin()) < 0.00001 && fabs(tree->getYMax() - root->getYMax()) < 0.00001) {
						if (distance < maxDistance && tree->active()) {
							root->leftNeighbours[root->leftProgress] = tree;
							root->leftProgress++;
						}
						root->leftNeighboursLandscape[root->leftProgressLandscape] = tree;
						root->leftProgressLandscape++;
					}
				}
			}

		return;
	}
	assignNeighbours(tree->nw,root);
	assignNeighbours(tree->ne,root);
	assignNeighbours(tree->sw,root);
	assignNeighbours(tree->se,root);

	return;
}

void TreeMesh::assignPotentials(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->setEp( potentialArray[tree->identifier] );
			if (potentialMax < tree->potential) { potentialMax = tree->potential; }
			if (potentialMin > tree->potential) { potentialMin = tree->potential; }
		} else {
			tree->setEp(0.0);
		}

		return;
	}
	assignPotentials(tree->nw);
	assignPotentials(tree->ne);
	assignPotentials(tree->sw);
	assignPotentials(tree->se);

	return;
}

void TreeMesh::findExtremesD(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			if (diffusionMax < tree->getD()) { diffusionMax = tree->getD(); }
			if (diffusionMin > tree->getD()) { diffusionMin = tree->getD(); }
			if (diffusionLogMax < tree->getDlog()) { diffusionLogMax = tree->getDlog(); }
			if (diffusionLogMin > tree->getDlog()) { diffusionLogMin = tree->getDlog(); }
		} else {
			tree->setFx(0.0);
			tree->setFy(0.0);
			tree->setD(0.0);
			tree->setDlog(0.0);
			tree->setFmag(0.0);
		}
		return;
	}
	findExtremesD(tree->nw);
	findExtremesD(tree->ne);
	findExtremesD(tree->sw);
	findExtremesD(tree->se);

	return;
}

void TreeMesh::findExtremesDF(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			if (forceMax < tree->getFmag()) { forceMax = tree->getFmag(); }
			if (forceMin > tree->getFmag()) { forceMin = tree->getFmag(); }
			if (diffusionMax < tree->getD()) { diffusionMax = tree->getD(); }
			if (diffusionMin > tree->getD()) { diffusionMin = tree->getD(); }
			if (diffusionLogMax < tree->getDlog()) { diffusionLogMax = tree->getDlog(); }
			if (diffusionLogMin > tree->getDlog()) { diffusionLogMin = tree->getDlog(); }
		} else {
			tree->setFx(0.0);
			tree->setFy(0.0);
			tree->setD(0.0);
			tree->setDlog(0.0);
			tree->setFmag(0.0);
		}
		return;
	}
	findExtremesDF(tree->nw);
	findExtremesDF(tree->ne);
	findExtremesDF(tree->sw);
	findExtremesDF(tree->se);

	return;
}

void TreeMesh::applyTreeMesh(QuadTree *tree) {

	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (tree->getCount() >= file->treeMeshGui->getMinPoints()) {
			tree->activate();
			if (file->treeMesh->maxCount < tree->getCount()) {
				file->treeMesh->maxCount = tree->getCount();
				file->treeMesh->setMaxTree(tree);
			}
			if (file->treeMesh->minCount > tree->getCount()) {
				file->treeMesh->minCount = tree->getCount();
				file->treeMesh->setMinTree(tree);
			}
			totalVariables++;
		} else {
			tree->deactivate();
		}
//		tree->area = tree->bounds.w*tree->bounds.h;
//		fprintf(stderr,"%i:\t%f\t%f\t%f\n",tree->count,tree->area,PI*tree->variance*tree->variance,9.0*PI*tree->variance*tree->variance);
		tree->area = PI*tree->variance*tree->variance;

		allVariables++;

		return;
	}
	applyTreeMesh(tree->nw);
	applyTreeMesh(tree->ne);
	applyTreeMesh(tree->sw);
	applyTreeMesh(tree->se);
	return;

}

void TreeMesh::updateTreeMesh(QuadTree *tree) {

	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (tree->active()) {
			if (file->treeMesh->maxCount < tree->getCount()) {
				file->treeMesh->maxCount = tree->getCount();
				file->treeMesh->setMaxTree(tree);
			}
			if (file->treeMesh->minCount > tree->getCount()) {
				file->treeMesh->minCount = tree->getCount();
				file->treeMesh->setMinTree(tree);
			}
			totalVariables++;
			tree->identifier = identifierCounter;
			identifierCounter++;
		} else {
			tree->deactivate();
		}

		return;
	}
	updateTreeMesh(tree->nw);
	updateTreeMesh(tree->ne);
	updateTreeMesh(tree->sw);
	updateTreeMesh(tree->se);
	return;

}

void TreeMesh::getLeafSelection(double x, double y, QuadTree *tree, bool *found) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (x >= tree->bounds.x && x <= tree->bounds.x+tree->bounds.w &&
			y >= tree->bounds.y && y <= tree->bounds.y+tree->bounds.h) {
			this->selectedQuadLeaf = tree;
			*found = true;
		}

		return;
	}

	if (*found != true) {
		getLeafSelection(x,y,tree->nw,found);
		getLeafSelection(x,y,tree->ne,found);
		getLeafSelection(x,y,tree->sw,found);
		getLeafSelection(x,y,tree->se,found);
	}
	return;
}

void TreeMesh::getLeafSelection2(double x, double y, QuadTree *tree, bool *found) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {

		if (x >= tree->bounds.x && x <= tree->bounds.x+tree->bounds.w &&
			y >= tree->bounds.y && y <= tree->bounds.y+tree->bounds.h) {
			this->selectedQuadLeaf2 = tree;
			*found = true;
		}
		return;
	}

	if (*found != true) {
		getLeafSelection2(x,y,tree->nw,found);
		getLeafSelection2(x,y,tree->ne,found);
		getLeafSelection2(x,y,tree->sw,found);
		getLeafSelection2(x,y,tree->se,found);
	}
	return;
}

void TreeMesh::preInferDV() {

	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	// count number of variables
	// get estimate of diffusion coefficient in x and y
	// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
	activeDetections = 0;
	activeCells = 0;
	xMean = 0.0;
	yMean = 0.0;
	dtMean = 0.0;

	// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
	initializeTreeMesh(quadTree);

	// find neighbours
	resetNeighbourCount(quadTree);
	findNeighbours(quadTree);
	assignNeighbours(quadTree);

	// determine whether leaf should be active or not
	minCount = 1000000;
	maxCount = -10000000;
	totalVariables = 0;
	identifierCounter = 0;
	updateTreeMesh(quadTree);

	xMean /= (double) (activeDetections-1);
	yMean /= (double) (activeDetections-1);
	dtMean /= (double) (activeDetections-1);

	const double D_eff_x = xMean*xMean/dtMean;
	const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

	// length of optimization array
	dimensions = 2*totalVariables; // D_ij,Fx_ij,Fy_ij

	hessian = new double*[dimensions];
	for(int i = 0; i < dimensions; i++) {
		hessian[i] = new double[dimensions];
	}
	for (int d = 0; d < dimensions; d++) {
		for (int e = 0; e < dimensions; e++) {
			hessian[d][e] = 0.0;
		}
	}

	// clear optimization array if re-run
	if (optimizationArray != NULL) { delete [] optimizationArray; }

	optimizationArray = new double[dimensions];
	for (int i = 0; i < dimensions/2.0; i++) {
		optimizationArray[2*i] = 1.0/2.0*(D_eff_x + D_eff_y);
		getTree(this->quadTree,i); // assigns leaf with current identity to identifierLeaf
		optimizationArray[2*i+1] = ( -log((double)identifiedQuadLeaf->getCount()/(double)maxCount) );
	}

	// assign areas for smoothing prior
	if (file->treeMeshGui->smoothingPriorButton->value()) {
		assignAreas(quadTree);
	}
}

void TreeMesh::inferDV(QuadTree* tree) {
	char label[50];
	if (file->treeMesh->selectionMode()) {
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingQuadTreeSelection, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSmoothingQuadTreeSelection(optimizationArray));
			textDisplayUpdate(label);
		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorQuadTreeSelection, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorQuadTreeSelection(optimizationArray));
			textDisplayUpdate(label);
		}
	}
	else {
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingQuadTree, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorSmoothingQuadTree(optimizationArray));
			textDisplayUpdate(label);
		} else {
			dfpmin(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorQuadTree, (dfunc));
			sprintf(label,"Cost =  %f\n",dvPosteriorQuadTree(optimizationArray));
			textDisplayUpdate(label);
		}
	}
}

int TreeMesh::postInferDV() {
	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	findExtremesDV(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {
		offsetPotentials();

		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}
	}

	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	return 1;
}

void TreeMesh::inferRandomizedOptimizationDV() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;

		// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
		initializeTreeMesh(quadTree);

		// find neighbours
		resetNeighbourCount(quadTree);
		findNeighbours(quadTree);
		assignNeighbours(quadTree);

		// determine whether leaf should be active or not
		minCount = 1000000;
		maxCount = -10000000;
		totalVariables = 0;
		identifierCounter = 0;
		updateTreeMesh(quadTree);

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[2*totalVariables];

		activeCells = identifierCounter;

		for (int i = 0; i < activeCells; i++) {
			getTree(this->quadTree,i);
			identifiedQuadLeaf->setDiffusion(1.0/2.0*(D_eff_x + D_eff_y));
			identifiedQuadLeaf->setPotential(-log((double)identifiedQuadLeaf->getCount()/(double)maxCount));
			intermediateOptimizationArray[2*i] = 0.5*(D_eff_x + D_eff_y);
			intermediateOptimizationArray[2*i+1] = -log((double)identifiedQuadLeaf->getCount()/(double)maxCount);
		}

		// assign areas for smoothing prior
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			assignAreas(quadTree);
		}
	}

	const double tolerance = (double)file->treeMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->treeMeshGui->roMaximumIterationsSlider->value();

	roQuadTrees = NULL;
	roIdentifiers = NULL;
	optimizationArray = NULL;
	QuadTree *currentTree = NULL;

	int neighbourCount = 0;
	int loops = 0;
	const float radius = file->treeMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		// SGD deactivate all zones
		if (loops == 0) { this->roDeactivateZones(quadTree); }

		if (loops > 0) {
			for (int d = 0; d < roZones; d++) {
				roQuadTrees[d]->roDeactivate();
				roQuadTrees[d]->setRoIdentifier(-1);
			}
		}

		// select leaf at random
		const int currentCell = rand() % activeCells;
		getTree(quadTree,currentCell);
		currentTree = identifiedQuadLeaf;

		neighbourCount = 0;
		// count number of cells in defined radius (based on barycentre coordinate)
		countWithinCircle(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&neighbourCount);

		roZones = neighbourCount;

		if (roQuadTrees != NULL) {
			delete [] roQuadTrees;
			delete [] roIdentifiers;
			roQuadTrees = NULL;
			roIdentifiers = NULL;
		}

		roQuadTrees = new QuadTree*[roZones];
		roIdentifiers = new int[roZones];

		// define optimization array
		dimensions = 2*roZones;
		if (optimizationArray != NULL) {
			delete [] optimizationArray;
			optimizationArray = NULL;
		}

		optimizationArray = new double[dimensions];

		int progress = 0;
		assignWithinCircleDV(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&progress);

		if (hessian != NULL) {
			delete [] hessian;
			hessian = NULL;
		}
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// perform inference calculation
		if (file->treeMesh->selectionMode()) {
			if (file->treeMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingQuadTreeRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorQuadTreeRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
			}
		} else {
			if (file->treeMeshGui->smoothingPriorButton->value()) {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorSmoothingQuadTreeRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
			else {
				dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dvPosteriorQuadTreeRandomizedOptimization, (dfunc), maxIterations, loops+1);
			}
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			getTree(this->quadTree,roIdentifiers[t]);
			intermediateOptimizationArray[2*roIdentifiers[t]] = optimizationArray[2*t];
			intermediateOptimizationArray[2*roIdentifiers[t]+1] = optimizationArray[2*t+1];
		}

		oldCost = cost;

		if (file->treeMesh->selectionMode()) {
			if (file->treeMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingQuadTreeSelection(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorQuadTreeSelection(intermediateOptimizationArray);
			}
		} else {
			if (file->treeMeshGui->smoothingPriorButton->value()) {
				cost = dvPosteriorSmoothingQuadTree(intermediateOptimizationArray);
			}
			else {
				cost = dvPosteriorQuadTree(intermediateOptimizationArray);
			}
		}

		if (iMAP->pauseCalculation == false) {
			for (int t = 0; t < roZones; t++) {
				roQuadTrees[t]->setForceX(dvGradVxQuadTreeRandomizedOptimization(optimizationArray,roQuadTrees[t]));
				roQuadTrees[t]->setForceY(dvGradVyQuadTreeRandomizedOptimization(optimizationArray,roQuadTrees[t]));
				roQuadTrees[t]->setForceMagnitude(sqrt(roQuadTrees[t]->getForceX()*roQuadTrees[t]->getForceX()+roQuadTrees[t]->getForceY()*roQuadTrees[t]->getForceY()));
				roQuadTrees[t]->setDiffusion(optimizationArray[2*t]);
				roQuadTrees[t]->setDiffusionLog(log(optimizationArray[2*t]));
				roQuadTrees[t]->setPotential(optimizationArray[2*t+1]);
			}
		}

		Fl::lock();
		// update progress window
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDV();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	delete [] optimizationArray;
	optimizationArray = NULL;
	optimizationArray = intermediateOptimizationArray;

	findExtremesDV(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {

		offsetPotentials();

		// clean up
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;

	}
	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (roQuadTrees != NULL) {
		delete [] roQuadTrees;
		delete [] roIdentifiers;
		roQuadTrees = NULL;
		roIdentifiers = NULL;
	}
}

int TreeMesh::inferRandomizedOptimizationDF() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;

		// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
		initializeTreeMesh(quadTree);

		// find neighbours
		resetNeighbourCount(quadTree);
		findNeighbours(quadTree);
		assignNeighbours(quadTree);

		// determine whether leaf should be active or not
		minCount = 1000000;
		maxCount = -10000000;
		totalVariables = 0;
		identifierCounter = 0;
		updateTreeMesh(quadTree);

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		activeCells = identifierCounter;

		for (int i = 0; i < activeCells; i++) {
			getTree(this->quadTree,i);
			identifiedQuadLeaf->setDiffusion(0.5*(D_eff_x + D_eff_y));
			identifiedQuadLeaf->setForceX(0.0);
			identifiedQuadLeaf->setForceY(0.0);
			intermediateOptimizationArray[3*i] = identifiedQuadLeaf->getDiffusion();
			intermediateOptimizationArray[3*i+1] = identifiedQuadLeaf->getForceX();
			intermediateOptimizationArray[3*i+2] = identifiedQuadLeaf->getForceY();
		}

		// assign areas for smoothing prior
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			assignAreas(quadTree);
		}
	}

	const double tolerance = (double)file->treeMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->treeMeshGui->roMaximumIterationsSlider->value();

	roQuadTrees = NULL;
	roIdentifiers = NULL;
	optimizationArray = NULL;
	QuadTree *currentTree = NULL;

	int neighbourCount = 0;
	int loops = 0;
	const float radius = file->treeMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		// SGD deactivate all zones
		if (loops == 0) { this->roDeactivateZones(quadTree); }

		if (loops > 0) {
			for (int d = 0; d < roZones; d++) {
				roQuadTrees[d]->roDeactivate();
				roQuadTrees[d]->setRoIdentifier(-1);
			}
		}

		// select leaf at random
		const int currentCell = rand() % activeCells;
		getTree(quadTree,currentCell);
		currentTree = identifiedQuadLeaf;

		neighbourCount = 0;
		// count number of cells in defined radius (based on barycentre coordinate)
		countWithinCircle(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&neighbourCount);

		roZones = neighbourCount;

		if (roQuadTrees != NULL) {
			delete [] roQuadTrees;
			delete [] roIdentifiers;
			roQuadTrees = NULL;
			roIdentifiers = NULL;
		}

		roQuadTrees = new QuadTree*[roZones];
		roIdentifiers = new int[roZones];

		// define optimization array
		dimensions = 3*roZones;
		if (optimizationArray != NULL) {
			delete [] optimizationArray;
			optimizationArray = NULL;
		}

		optimizationArray = new double[dimensions];

		int progress = 0;
		assignWithinCircleDF(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&progress);

		if (hessian != NULL) {
			delete [] hessian;
			hessian = NULL;
		}
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// perform inference calculation
		if (file->treeMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingQuadTreeRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dfPosteriorSmoothingQuadTreeRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			getTree(this->quadTree,roIdentifiers[t]);
			intermediateOptimizationArray[3*roIdentifiers[t]] = optimizationArray[3*t];
			intermediateOptimizationArray[3*roIdentifiers[t]+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*roIdentifiers[t]+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->treeMesh->selectionMode()) {
			cost = dfPosteriorSmoothingQuadTreeSelection(intermediateOptimizationArray);
		} else {
			cost = dfPosteriorSmoothingQuadTree(intermediateOptimizationArray);
		}

		if (iMAP->pauseCalculation == false) {
			for (int t = 0; t < roZones; t++) {
				roQuadTrees[t]->setForceX(optimizationArray[3*t+1]);
				roQuadTrees[t]->setForceY(optimizationArray[3*t+2]);
				roQuadTrees[t]->setForceMagnitude(sqrt(roQuadTrees[t]->getForceX()*roQuadTrees[t]->getForceX()+roQuadTrees[t]->getForceY()*roQuadTrees[t]->getForceY()));
				roQuadTrees[t]->setDiffusion(optimizationArray[3*t]);
				roQuadTrees[t]->setDiffusionLog(log(optimizationArray[3*t]));
			}
		}

		Fl::lock();
		// update progress window
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesDF();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	delete [] optimizationArray;
	optimizationArray = NULL;
	optimizationArray = intermediateOptimizationArray;

	findExtremesDFSmoothing(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {

		// clean up
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;

	}
	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (roQuadTrees != NULL) {
		delete [] roQuadTrees;
		delete [] roIdentifiers;
		roQuadTrees = NULL;
		roIdentifiers = NULL;
	}

	return fl_choice("Compute potentials?","No","Yes",NULL);

}

void TreeMesh::inferRandomizedOptimizationDDr() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;

		// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
		initializeTreeMesh(quadTree);

		// find neighbours
		resetNeighbourCount(quadTree);
		findNeighbours(quadTree);
		assignNeighbours(quadTree);

		// determine whether leaf should be active or not
		minCount = 1000000;
		maxCount = -10000000;
		totalVariables = 0;
		identifierCounter = 0;
		updateTreeMesh(quadTree);

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[3*totalVariables];

		activeCells = identifierCounter;

		for (int i = 0; i < activeCells; i++) {
			getTree(this->quadTree,i);
			identifiedQuadLeaf->setDiffusion(1.0/2.0*(D_eff_x + D_eff_y));
			identifiedQuadLeaf->setForceX(0.0);
			identifiedQuadLeaf->setForceY(0.0);
			intermediateOptimizationArray[3*i] = identifiedQuadLeaf->getDiffusion();
			intermediateOptimizationArray[3*i+1] = identifiedQuadLeaf->getForceX();
			intermediateOptimizationArray[3*i+2] = identifiedQuadLeaf->getForceY();
		}

		// assign areas for smoothing prior
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			assignAreas(quadTree);
		}
	}

	const double tolerance = (double)file->treeMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->treeMeshGui->roMaximumIterationsSlider->value();

	roQuadTrees = NULL;
	roIdentifiers = NULL;
	optimizationArray = NULL;
	QuadTree *currentTree = NULL;

	int neighbourCount = 0;
	int loops = 0;
	const float radius = file->treeMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		// SGD deactivate all zones
		if (loops == 0) { this->roDeactivateZones(quadTree); }

		if (loops > 0) {
			for (int d = 0; d < roZones; d++) {
				roQuadTrees[d]->roDeactivate();
				roQuadTrees[d]->setRoIdentifier(-1);
			}
		}

		// select leaf at random
		const int currentCell = rand() % activeCells;
		getTree(quadTree,currentCell);
		currentTree = identifiedQuadLeaf;

		neighbourCount = 0;
		// count number of cells in defined radius (based on barycentre coordinate)
		countWithinCircle(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&neighbourCount);

		roZones = neighbourCount;

		if (roQuadTrees != NULL) {
			delete [] roQuadTrees;
			delete [] roIdentifiers;
			roQuadTrees = NULL;
			roIdentifiers = NULL;
		}

		roQuadTrees = new QuadTree*[roZones];
		roIdentifiers = new int[roZones];

		// define optimization array
		dimensions = 3*roZones;
		if (optimizationArray != NULL) {
			delete [] optimizationArray;
			optimizationArray = NULL;
		}

		optimizationArray = new double[dimensions];

		int progress = 0;
		assignWithinCircleDF(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&progress);

		if (hessian != NULL) {
			delete [] hessian;
			hessian = NULL;
		}
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// perform inference calculation
		if (file->treeMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingQuadTreeRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, ddrPosteriorSmoothingQuadTreeRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			getTree(this->quadTree,roIdentifiers[t]);
			intermediateOptimizationArray[3*roIdentifiers[t]] = optimizationArray[3*t];
			intermediateOptimizationArray[3*roIdentifiers[t]+1] = optimizationArray[3*t+1];
			intermediateOptimizationArray[3*roIdentifiers[t]+2] = optimizationArray[3*t+2];
		}

		oldCost = cost;

		if (file->treeMesh->selectionMode()) {
			cost = ddrPosteriorSmoothingQuadTreeSelection(intermediateOptimizationArray);
		} else {
			cost = ddrPosteriorSmoothingQuadTree(intermediateOptimizationArray);
		}

		if (iMAP->pauseCalculation == false) {
			for (int t = 0; t < roZones; t++) {
				roQuadTrees[t]->setForceX(optimizationArray[3*t+1]);
				roQuadTrees[t]->setForceY(optimizationArray[3*t+2]);
				roQuadTrees[t]->setForceMagnitude(sqrt(roQuadTrees[t]->getForceX()*roQuadTrees[t]->getForceX()+roQuadTrees[t]->getForceY()*roQuadTrees[t]->getForceY()));
				roQuadTrees[t]->setDiffusion(optimizationArray[3*t]);
				roQuadTrees[t]->setDiffusionLog(log(optimizationArray[3*t]));
			}
		}

		Fl::lock();
		// update progress window
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;

		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	delete [] optimizationArray;
	optimizationArray = NULL;
	optimizationArray = intermediateOptimizationArray;

	findExtremesDDrSmoothing(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {

		// clean up
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;

	}
	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (roQuadTrees != NULL) {
		delete [] roQuadTrees;
		delete [] roIdentifiers;
		roQuadTrees = NULL;
		roIdentifiers = NULL;
	}

}

void TreeMesh::inferRandomizedOptimizationD() {
	// turn off overlays if inference re-run
	file->inferred = false;
	file->treeMeshGui->saveButton->deactivate();

	// get polynomial order from slider
	beta = file->treeMeshGui->getBeta();
	sigma = file->treeMeshGui->localizationPrecisionSlider->value()/1000.0;

	if (iMAP->pauseCalculation == false) {
		// count number of variables
		// get estimate of diffusion coefficient in x and y
		// set identifier for cell (variable) such that it corresponds to an entry in the optimization array
		activeDetections = 0;
		activeCells = 0;
		xMean = 0.0;
		yMean = 0.0;
		dtMean = 0.0;

		// calculate xMean,yMean,dtMean,activeDetections,totalVariables,activeCells
		initializeTreeMesh(quadTree);

		// find neighbours
		resetNeighbourCount(quadTree);
		findNeighbours(quadTree);
		assignNeighbours(quadTree);

		// determine whether leaf should be active or not
		minCount = 1000000;
		maxCount = -10000000;
		totalVariables = 0;
		identifierCounter = 0;
		updateTreeMesh(quadTree);

		xMean /= (double) (activeDetections-1);
		yMean /= (double) (activeDetections-1);
		dtMean /= (double) (activeDetections-1);

		const double D_eff_x = xMean*xMean/dtMean;
		const double D_eff_y = yMean*yMean/dtMean; // valeur initial de diffusion

		intermediateOptimizationArray = new double[totalVariables];

		activeCells = identifierCounter;

		for (int i = 0; i < activeCells; i++) {
			getTree(this->quadTree,i);
			identifiedQuadLeaf->setDiffusion(0.5*(D_eff_x + D_eff_y));
			intermediateOptimizationArray[i] = 0.5*(D_eff_x + D_eff_y);
		}

		// assign areas for smoothing prior
		if (file->treeMeshGui->smoothingPriorButton->value()) {
			assignAreas(quadTree);
		}
	}

	const double tolerance = (double)file->treeMeshGui->roToleranceSlider->value()/100.0;
	const int maxIterations = (int)file->treeMeshGui->roMaximumIterationsSlider->value();

	roQuadTrees = NULL;
	roIdentifiers = NULL;
	optimizationArray = NULL;
	QuadTree *currentTree = NULL;

	int neighbourCount = 0;
	int loops = 0;
	const float radius = file->treeMeshGui->roRadiusSlider->value()/1000.0;

	iMAP->pauseCalculation = false;
	iMAP->stopCalculation = false;

	/* for cost evaluation */
	double cost = 10000.0;
	double oldCost = 10000000.0;

	while ( iMAP->pauseCalculation == false && iMAP->stopCalculation == false && fabs( (oldCost-cost)/cost ) > tolerance ) {

		// SGD deactivate all zones
		if (loops == 0) { this->roDeactivateZones(quadTree); }

		if (loops > 0) {
			for (int d = 0; d < roZones; d++) {
				roQuadTrees[d]->roDeactivate();
				roQuadTrees[d]->setRoIdentifier(-1);
			}
		}

		// select leaf at random
		const int currentCell = rand() % activeCells;
		getTree(quadTree,currentCell);
		currentTree = identifiedQuadLeaf;

		neighbourCount = 0;
		// count number of cells in defined radius (based on barycentre coordinate)
		countWithinCircle(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&neighbourCount);

		roZones = neighbourCount;

		if (roQuadTrees != NULL) {
			delete [] roQuadTrees;
			delete [] roIdentifiers;
			roQuadTrees = NULL;
			roIdentifiers = NULL;
		}

		roQuadTrees = new QuadTree*[roZones];
		roIdentifiers = new int[roZones];

		// define optimization array
		dimensions = roZones;
		if (optimizationArray != NULL) {
			delete [] optimizationArray;
			optimizationArray = NULL;
		}

		optimizationArray = new double[dimensions];

		int progress = 0;
		assignWithinCircleD(quadTree,currentTree->getXCentroid(),currentTree->getYCentroid(),radius,&progress);

		if (hessian != NULL) {
			delete [] hessian;
			hessian = NULL;
		}
		hessian = new double*[dimensions];
		for(int i = 0; i < dimensions; i++) {
			hessian[i] = new double[dimensions];
		}
		for (int d = 0; d < dimensions; d++) {
			for (int e = 0; e < dimensions; e++) {
				hessian[d][e] = 0.0;
			}
		}

		// perform inference calculation
		if (file->treeMesh->selectionMode()) {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingQuadTreeRandomizedOptimizationSelection, (dfunc), maxIterations, loops+1);
		} else {
			dfpminRO(optimizationArray, dimensions, GTOL, iterations, fret, dPosteriorSmoothingQuadTreeRandomizedOptimization, (dfunc), maxIterations, loops+1);
		}

		// assign values to optimization array to perform "global" cost calculation
		for (int t = 0; t < roZones; t++) {
			getTree(this->quadTree,roIdentifiers[t]);
			intermediateOptimizationArray[roIdentifiers[t]] = optimizationArray[t];
		}

		oldCost = cost;

		if (file->treeMesh->selectionMode()) {
			cost = dPosteriorSmoothingQuadTreeSelection(intermediateOptimizationArray);
		} else {
			cost = dPosteriorSmoothingQuadTree(intermediateOptimizationArray);
		}

		if (iMAP->pauseCalculation == false) {
			for (int t = 0; t < roZones; t++) {
				roQuadTrees[t]->setDiffusion(optimizationArray[t]);
				roQuadTrees[t]->setDiffusionLog(log(optimizationArray[t]));
			}
		}

		Fl::lock();
		// update progress window
		randomizedOptimizationGui->update(cost,roZones);
		file->drawRandomizedOptimizationZones = true;
		findExtremesD();
		Fl::unlock();
		Fl::check();

		loops++;
	}

	file->drawRandomizedOptimizationZones = false;

	forceMax = -1000000.0;
	forceMin = 1000000.0;
	potentialMax = -1000000.0;
	potentialMin = 1000000.0;
	diffusionMax = -1000000.0;
	diffusionMin = 1000000.0;
	diffusionLogMax = -1000000.0;
	diffusionLogMin = 1000000.0;

	delete [] optimizationArray;
	optimizationArray = NULL;
	optimizationArray = intermediateOptimizationArray;

	findExtremesD(file->treeMesh->quadTree);

	if (iMAP->pauseCalculation == false) {

		// clean up
		if (hessian != NULL) {
			for (int d = 0; d < dimensions; d++) {
				delete [] hessian[d];
			}
			delete [] hessian;
			hessian = NULL;
		}

		delete [] intermediateOptimizationArray;

	}
	// enable inference variable overlays
	file->inferred = true;
	file->treeMeshGui->saveButton->activate();

	if (roQuadTrees != NULL) {
		delete [] roQuadTrees;
		delete [] roIdentifiers;
		roQuadTrees = NULL;
		roIdentifiers = NULL;
	}
}

void TreeMesh::countWithinCircle(QuadTree *tree, float xCentre, float yCentre, float radius, int *count) {

	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			const float dist = sqrt( (xCentre-tree->getXCentroid())*(xCentre-tree->getXCentroid()) + (yCentre-tree->getYCentroid())*(yCentre-tree->getYCentroid()) );
			if (dist < radius) {
				*count = *count + 1;
			}
		}
		return;
	}
	countWithinCircle(tree->nw,xCentre,yCentre,radius,count);
	countWithinCircle(tree->ne,xCentre,yCentre,radius,count);
	countWithinCircle(tree->sw,xCentre,yCentre,radius,count);
	countWithinCircle(tree->se,xCentre,yCentre,radius,count);
	return;

}

void TreeMesh::assignWithinCircleDV(QuadTree *tree, float xCentre, float yCentre, float radius, int *count) {

	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			const float dist = sqrt( (xCentre-tree->getXCentroid())*(xCentre-tree->getXCentroid()) + (yCentre-tree->getYCentroid())*(yCentre-tree->getYCentroid()) );
			if (dist < radius) {
				roQuadTrees[*count] = tree;
				roIdentifiers[*count] = tree->getIdentifier();

				optimizationArray[2*(*count)] = tree->getDiffusion();
				optimizationArray[2*(*count)+1] = tree->getPotential();
				tree->roActivate();
				tree->setRoIdentifier(*count);

				*count = *count + 1;
			}
		}
		return;
	}
	assignWithinCircleDV(tree->nw,xCentre,yCentre,radius,count);
	assignWithinCircleDV(tree->ne,xCentre,yCentre,radius,count);
	assignWithinCircleDV(tree->sw,xCentre,yCentre,radius,count);
	assignWithinCircleDV(tree->se,xCentre,yCentre,radius,count);
	return;

}

void TreeMesh::assignWithinCircleDF(QuadTree *tree, float xCentre, float yCentre, float radius, int *count) {

	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			const float dist = sqrt( (xCentre-tree->getXCentroid())*(xCentre-tree->getXCentroid()) + (yCentre-tree->getYCentroid())*(yCentre-tree->getYCentroid()) );
			if (dist < radius) {
				roQuadTrees[*count] = tree;
				roIdentifiers[*count] = tree->getIdentifier();

				optimizationArray[3*(*count)] = tree->getDiffusion();
				optimizationArray[3*(*count)+1] = tree->getForceX();
				optimizationArray[3*(*count)+2] = tree->getForceY();
				tree->roActivate();
				tree->setRoIdentifier(*count);

				*count = *count + 1;
			}
		}
		return;
	}
	assignWithinCircleDF(tree->nw,xCentre,yCentre,radius,count);
	assignWithinCircleDF(tree->ne,xCentre,yCentre,radius,count);
	assignWithinCircleDF(tree->sw,xCentre,yCentre,radius,count);
	assignWithinCircleDF(tree->se,xCentre,yCentre,radius,count);
	return;

}

void TreeMesh::assignWithinCircleDDr(QuadTree *tree, float xCentre, float yCentre, float radius, int *count) {

	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			const float dist = sqrt( (xCentre-tree->getXCentroid())*(xCentre-tree->getXCentroid()) + (yCentre-tree->getYCentroid())*(yCentre-tree->getYCentroid()) );
			if (dist < radius) {
				roQuadTrees[*count] = tree;
				roIdentifiers[*count] = tree->getIdentifier();

				optimizationArray[3*(*count)] = tree->getDiffusion();
				optimizationArray[3*(*count)+1] = tree->getForceX();
				optimizationArray[3*(*count)+2] = tree->getForceY();
				tree->roActivate();
				tree->setRoIdentifier(*count);

				*count = *count + 1;
			}
		}
		return;
	}
	assignWithinCircleDDr(tree->nw,xCentre,yCentre,radius,count);
	assignWithinCircleDDr(tree->ne,xCentre,yCentre,radius,count);
	assignWithinCircleDDr(tree->sw,xCentre,yCentre,radius,count);
	assignWithinCircleDDr(tree->se,xCentre,yCentre,radius,count);
	return;

}

void TreeMesh::assignWithinCircleD(QuadTree *tree, float xCentre, float yCentre, float radius, int *count) {

	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			const float dist = sqrt( (xCentre-tree->getXCentroid())*(xCentre-tree->getXCentroid()) + (yCentre-tree->getYCentroid())*(yCentre-tree->getYCentroid()) );
			if (dist < radius) {
				roQuadTrees[*count] = tree;
				roIdentifiers[*count] = tree->getIdentifier();

				optimizationArray[(*count)] = tree->getDiffusion();
				tree->roActivate();
				tree->setRoIdentifier(*count);

				*count = *count + 1;
			}
		}
		return;
	}
	assignWithinCircleD(tree->nw,xCentre,yCentre,radius,count);
	assignWithinCircleD(tree->ne,xCentre,yCentre,radius,count);
	assignWithinCircleD(tree->sw,xCentre,yCentre,radius,count);
	assignWithinCircleD(tree->se,xCentre,yCentre,radius,count);
	return;

}

void TreeMesh::findExtremesDSmoothing(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
//			fprintf(stderr,"%i\t%f,\t%f\n",tree->identifier,tree->xForce,tree->yForce);
			tree->diffusion = optimizationArray[tree->identifier];
			tree->diffusionLog = log(optimizationArray[tree->identifier]);
			if (diffusionMax < tree->diffusion) { diffusionMax = tree->diffusion; }
			if (diffusionMin > tree->diffusion) { diffusionMin = tree->diffusion; }
			if (diffusionLogMax < tree->diffusionLog) { diffusionLogMax = tree->diffusionLog; }
			if (diffusionLogMin > tree->diffusionLog) { diffusionLogMin = tree->diffusionLog; }
		} else {
			tree->xForce = 0.0;
			tree->yForce = 0.0;
			tree->forceMagnitude = 0.0;
			tree->diffusion = 0.0;
			tree->diffusionLog = 0.0;
			tree->potential = 0.0;
		}
		return;
	}
	findExtremesDSmoothing(tree->nw);
	findExtremesDSmoothing(tree->ne);
	findExtremesDSmoothing(tree->sw);
	findExtremesDSmoothing(tree->se);
	return;
}

void TreeMesh::findExtremesDV(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			if (roEnable == false) {
				tree->xForce = dvGradVxQuadTree(optimizationArray,tree);
				tree->yForce = dvGradVyQuadTree(optimizationArray,tree);
	//			fprintf(stderr,"%i\t%f,\t%f\n",tree->identifier,tree->xForce,tree->yForce);
				tree->forceMagnitude = sqrt(tree->xForce*tree->xForce+tree->yForce*tree->yForce);
				tree->diffusion = optimizationArray[2*tree->identifier];
				tree->diffusionLog = log(optimizationArray[2*tree->identifier]);
				tree->potential = optimizationArray[2*tree->identifier+1];
			}
			if (forceMax < tree->forceMagnitude) { forceMax = tree->forceMagnitude; }
			if (forceMin > tree->forceMagnitude) { forceMin = tree->forceMagnitude; }
			if (diffusionMax < tree->diffusion) { diffusionMax = tree->diffusion; }
			if (diffusionMin > tree->diffusion) { diffusionMin = tree->diffusion; }
			if (diffusionLogMax < tree->diffusionLog) { diffusionLogMax = tree->diffusionLog; }
			if (diffusionLogMin > tree->diffusionLog) { diffusionLogMin = tree->diffusionLog; }
			if (potentialMax < tree->potential) { potentialMax = tree->potential; }
			if (potentialMin > tree->potential) { potentialMin = tree->potential; }
		} else {
			tree->xForce = 0.0;
			tree->yForce = 0.0;
			tree->forceMagnitude = 0.0;
			tree->diffusion = 0.0;
			tree->diffusionLog = 0.0;
			tree->potential = 0.0;
		}
		return;
	}
	findExtremesDV(tree->nw);
	findExtremesDV(tree->ne);
	findExtremesDV(tree->sw);
	findExtremesDV(tree->se);
	return;
}

void TreeMesh::findExtremesDFSmoothing(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->xForce = optimizationArray[3*tree->identifier+1];
			tree->yForce = optimizationArray[3*tree->identifier+2];
//			fprintf(stderr,"%i\t%f,\t%f\n",tree->identifier,tree->xForce,tree->yForce);
			tree->forceMagnitude = sqrt(tree->xForce*tree->xForce+tree->yForce*tree->yForce);
			tree->diffusion = optimizationArray[3*tree->identifier];
			tree->diffusionLog = log(optimizationArray[3*tree->identifier]);
//			tree->potential = optimizationArray[2*tree->identifier+1];
			if (forceMax < tree->forceMagnitude) { forceMax = tree->forceMagnitude; }
			if (forceMin > tree->forceMagnitude) { forceMin = tree->forceMagnitude; }
			if (diffusionMax < tree->diffusion) { diffusionMax = tree->diffusion; }
			if (diffusionMin > tree->diffusion) { diffusionMin = tree->diffusion; }
			if (diffusionLogMax < tree->diffusionLog) { diffusionLogMax = tree->diffusionLog; }
			if (diffusionLogMin > tree->diffusionLog) { diffusionLogMin = tree->diffusionLog; }
		} else {
			tree->xForce = 0.0;
			tree->yForce = 0.0;
			tree->forceMagnitude = 0.0;
			tree->diffusion = 0.0;
			tree->diffusionLog = 0.0;
			tree->potential = 0.0;
		}
		return;
	}
	findExtremesDFSmoothing(tree->nw);
	findExtremesDFSmoothing(tree->ne);
	findExtremesDFSmoothing(tree->sw);
	findExtremesDFSmoothing(tree->se);
	return;
}

void TreeMesh::findExtremesDDrSmoothing(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->xForce = optimizationArray[3*tree->identifier+1];
			tree->yForce = optimizationArray[3*tree->identifier+2];
//			fprintf(stderr,"%i\t%f,\t%f\n",tree->identifier,tree->xForce,tree->yForce);
			tree->forceMagnitude = sqrt(tree->xForce*tree->xForce+tree->yForce*tree->yForce);
			tree->diffusion = optimizationArray[3*tree->identifier];
			tree->diffusionLog = log(optimizationArray[3*tree->identifier]);
//			tree->potential = optimizationArray[2*tree->identifier+1];
			if (forceMax < tree->forceMagnitude) { forceMax = tree->forceMagnitude; }
			if (forceMin > tree->forceMagnitude) { forceMin = tree->forceMagnitude; }
			if (diffusionMax < tree->diffusion) { diffusionMax = tree->diffusion; }
			if (diffusionMin > tree->diffusion) { diffusionMin = tree->diffusion; }
			if (diffusionLogMax < tree->diffusionLog) { diffusionLogMax = tree->diffusionLog; }
			if (diffusionLogMin > tree->diffusionLog) { diffusionLogMin = tree->diffusionLog; }
		} else {
			tree->xForce = 0.0;
			tree->yForce = 0.0;
			tree->forceMagnitude = 0.0;
			tree->diffusion = 0.0;
			tree->diffusionLog = 0.0;
			tree->potential = 0.0;
		}
		return;
	}
	findExtremesDDrSmoothing(tree->nw);
	findExtremesDDrSmoothing(tree->ne);
	findExtremesDDrSmoothing(tree->sw);
	findExtremesDDrSmoothing(tree->se);
	return;
}

void TreeMesh::adjustPotentials(QuadTree* tree) {
	// only initialize in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()) {
			tree->potential = tree->potential - potentialMin;
		}
		return;
	}
	adjustPotentials(tree->nw);
	adjustPotentials(tree->ne);
	adjustPotentials(tree->sw);
	adjustPotentials(tree->se);
	return;
}

void TreeMesh::getTree(QuadTree *tree,int identifier) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->identifier == identifier) { this->identifiedQuadLeaf = tree; }
		return;
	}
	getTree(tree->nw,identifier);
	getTree(tree->ne,identifier);
	getTree(tree->sw,identifier);
	getTree(tree->se,identifier);
	return;
}

void TreeMesh::getTreeRo(QuadTree *tree,int identifier) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roIdentifier == identifier) { this->identifiedQuadLeaf = tree; }
		return;
	}
	getTree(tree->nw,identifier);
	getTree(tree->ne,identifier);
	getTree(tree->sw,identifier);
	getTree(tree->se,identifier);
	return;
}

void TreeMesh::updateNeighbours(QuadTree* tree) {
	file->treeMesh->resetNeighbourCount(tree);
	file->treeMesh->findNeighbours(tree);
	file->treeMesh->assignNeighbours(tree);
}

void TreeMesh::exportMesh() {

	Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
	char nativeFilename[FILENAME_MAX];

	native->title("Save Quad-Tree Mesh File");
	native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native->filter("Mesh File\t*.qmesh\n");
    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
    native->preset_file(file->fileName);

	if (native->show() == 0) {
		// store filename
		sprintf(nativeFilename,"%s%s",native->filename(),".qmesh");
		FILE * writeFile;

		writeFile = fopen(nativeFilename,"wb");

		// File Information
		fprintf(writeFile,"FILE INFORMATION\n");
		fprintf(writeFile,"Filename:\t%s\n",file->fileName);
		fprintf(writeFile,"Total Number of Trajectories:\t%i\n",file->numberOfFiles);
		fprintf(writeFile,"Duration [s]:\t%.3f\n",iMAP->endIntervalSlider->value()-iMAP->startIntervalSlider->value());
		fprintf(writeFile,"Acquisition Time [ms]:\t%.1f\n",1000.0*iMAP->exposureTime);
		if (this->selectionMode()) {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",selection.xMax-selection.xMin,selection.yMax-selection.yMin);
		} else {
			fprintf(writeFile,"Bounds [um]:\t[%f,%f]\n\n",this->xMax-this->xMin,this->yMax-this->yMin);
		}

		// Mesh Information
		fprintf(writeFile,"MESH INFORMATION\n");
		fprintf(writeFile,"Meshing: Quad-Tree\n");
		fprintf(writeFile,"Selection Mode: %i\n",this->selectionMode());
		fprintf(writeFile,"Minimum Leaf Power:\t%i\n",this->getMinLeafPower());
		fprintf(writeFile,"Minimum Leaf Capacity:\t%i\n",this->getMinCapacity());
		fprintf(writeFile,"Total Cells:\t%i\n",this->totalVariables);
		fprintf(writeFile,"Total Localizations:\t%i\n",file->localizationCount);
		fprintf(writeFile,"Active Cells:\t%i\n",this->getActiveCells());
		fprintf(writeFile,"Active Localizations:\t%i\n\n",this->getActiveDetections());

		// Inference Parameters
		fprintf(writeFile,"INFERENCE PARAMETERS\n");
		fprintf(writeFile,"Noise Sigma [nm]:\t%.1f\n",this->getSigma()*1000.0);
		fprintf(writeFile,"Maximum Neighbour Distance [nm]:\t%.1f\n",this->getMaximumNeighbourDistance());
		fprintf(writeFile,"Minimum Points per Cell:\t%i\n",this->getMinNumberOfPointsPerCell());
		fprintf(writeFile,"Optimization Scheme:\t%i\n",this->getOptimizationMode());
		fprintf(writeFile,"Potential Calculated:\t%i\n",this->zonalPotentialsCalculated);
		fprintf(writeFile,"Prior Enabled:\t%i\n",this->smoothingPriorEnabled());
		if (this->optimizationMode == 4) {
			fprintf(writeFile,"Polynomial Order:\t%i\n",this->getPolynomialOrder());
			fprintf(writeFile,"Polynomial Coefficients:\t");
			for (int q = 0; q < this->getCoefficients()-1; q++) {
				fprintf(writeFile,"%f\t",this->getOptimizationArray(q));
			}
			fprintf(writeFile,"%f\n\n",this->getOptimizationArray(getCoefficients()-1));
		} else {
			fprintf(writeFile,"Polynomial Order:\tN/A\n");
			fprintf(writeFile,"Polynomial Coefficients:\tN/A\n\n");
		}

		// Data
		if (file->optimizationMode == 2) { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tDx\t\tDy\t\tV\t\tPower\tWidth\n"); }
		else { fprintf(writeFile,"Active\tx-Center\ty-Center\tPoints\tx-Centroid\ty-Centroid\tVariance\tD\t\tFx\t\tFy\t\tV\t\tPower\tWidth\n"); }
		exportMesh(this->quadTree,writeFile);
		fclose(writeFile);
	}
}

void TreeMesh::exportMesh(QuadTree* tree,FILE * writeFile) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		fprintf(writeFile,"%i\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%f\n",
				tree->active(),
				tree->getXCentre(),
				tree->getYCentre(),
				tree->getCount(),
				tree->getXCentroid(),
				tree->getYCentroid(),
				tree->getVariance(),
				tree->getDiffusion(),
				tree->getForceX(),
				tree->getForceY(),
				tree->getPotential(),
				tree->getPower(),
				tree->getWidth()
		);
		return;
	}
	exportMesh(tree->nw,writeFile);
	exportMesh(tree->ne,writeFile);
	exportMesh(tree->sw,writeFile);
	exportMesh(tree->se,writeFile);
	return;
}


void randomizedOptimizationSaveCallback(Fl_Button*w,int*v) {

	char meshTypeString [2];
	sprintf(meshTypeString,"%i",v);
	const int meshType = atoi(meshTypeString);

	switch (meshType) {
	case 0: // regular mesh
		file->squareMesh->randomizedOptimizationGui->save();
		break;
	case 1: // voronoi tessellation
		file->voronoiMesh->randomizedOptimizationGui->save();
		break;
	case 2: // quad-tree
		file->treeMesh->randomizedOptimizationGui->save();
		break;
	}

}

void RandomizedOptimizationGui::update(double cost, int zones) {
	const int w = 500;
	const int tableWidth = w-20;
	const int h = 250;
	const int tableHeight = h-20;

	window->make_current();

	addNode(cost);

	if (n > 1) {

		fl_line_style(FL_SOLID);

		fl_color(FL_WHITE); fl_rectf(10,10,tableWidth,tableHeight);
		fl_color(FL_BLACK); fl_rect(10,10,tableWidth,tableHeight);

		fl_line(75,tableHeight-25,tableWidth-20,tableHeight-25);
		fl_line(75,tableHeight-25,75,30);

		fl_color(FL_BLUE);

		findExtremes();

		char xLabel[10];
		char yLabel[10];

		const double spacing = ((double)(tableWidth-90)/(double)n);
		DataNode *oldNode;
		DataNode *newNode = head;

		double oldVal,newVal;
		for (int g = 0; g < n-1; g++) {
			oldNode = newNode;
			newNode = newNode->next;
			oldVal = (oldNode->data-minCost)/(maxCost-minCost);
			newVal = (newNode->data-minCost)/(maxCost-minCost);
			const int intOldVal = (int)(oldVal*(double)(tableHeight-60));
			const int intNewVal = (int)(newVal*(double)(tableHeight-60));
			fl_line(76+(g)*spacing,(tableHeight-25)-intOldVal,76+(g+1)*spacing,(tableHeight-25)-intNewVal);
		}
		// draw x-axis labels
		fl_color(FL_BLACK);
		fl_font(iMAP->boldFont,10);
		for (int t = 0; t < 8; t++) {
			sprintf(xLabel,"%i",(int)(t*(n/8.0)));
			fl_draw(xLabel,72+t*(tableWidth-75)/8,tableHeight-15);
			fl_line_style(FL_DASH);
			fl_line(75+t*(tableWidth-75)/8,tableHeight-25,75+t*(tableWidth-75)/8,30);
		}
		for (int t = 0; t < 5; t++) {
			sprintf(yLabel,"%.3le",minCost+t*((maxCost-minCost)/5.0));
			fl_draw(yLabel,15,tableHeight-23-t*(tableHeight-50)/5);
			fl_line_style(FL_DASH);
			fl_line(75,tableHeight-25-(t)*(tableHeight-50)/5,tableWidth-20,tableHeight-25-(t)*(tableHeight-50)/5);
		}

		// plot annotations
		fl_font(iMAP->boldFont,12);
		fl_draw("Cost\n[a.u.]",20,20,50,30,FL_ALIGN_INSIDE|FL_ALIGN_CENTER,NULL,0);
		fl_draw("Iterations",tableWidth/2-20,h-20);

		fl_line_style(FL_SOLID);
		fl_draw_box(FL_BORDER_BOX,tableWidth-90,15,95,40,FL_GRAY);
		fl_color(FL_BLACK);
		char label[100];
		fl_font(iMAP->boldFont,10);
		sprintf(label,"%i Zones\n%.1f%% Tolerance\n%.1f",zones,file->randomizedOptimizationTolerance*100.0,cost);
		fl_draw(label,tableWidth-90,15,90,40,FL_ALIGN_RIGHT,NULL,0);
	}
}

RandomizedOptimizationGui::RandomizedOptimizationGui() {
	// initialize first node
	head = new DataNode;
	head->next = NULL;
	tail = NULL;
	minCost = 100000000.0;
	maxCost = -100000000.0;
	n = 0;
	const int w = 500;
	const int h = 250;

	window = new Fl_Window(w,h,"Randomized Optimization");
//			window->callback(treeMeshWindowCallback);
	window->begin();
	window->color(iMAP->bgColor);
	window->set_non_modal();

	saveButton = new Fl_Button(w/2-30,20,60,25,"Save");
	saveButton->labelsize(12);
	saveButton->labelfont(iMAP->normalFont);
	saveButton->callback((Fl_Callback*)randomizedOptimizationSaveCallback,(int*)file->meshType);
	saveButton->deactivate();

//			totalZonesBox = new Fl_Box(4*w/5,10,w/5-10,20,"Total Zones:");
//			totalZonesBox->labelsize(12);
//			totalZonesBox->labelfont(1);
//			totalZonesBox->align(FL_ALIGN_LEFT);
//			totalZonesBox->show();
//
//			optimizationZonesBox = new Fl_Box(4*w/5,30,w/5-10,20,"Optimization Zones:");
//			optimizationZonesBox->labelsize(12);
//			optimizationZonesBox->labelfont(1);
//			optimizationZonesBox->align(FL_ALIGN_LEFT);
//			optimizationZonesBox->show();
//
//			toleranceTargetBox = new Fl_Box(4*w/5,50,w/5-10,20,"Tolerance Target [%]:");
//			toleranceTargetBox->labelsize(12);
//			toleranceTargetBox->labelfont(1);
//			toleranceTargetBox->align(FL_ALIGN_LEFT);
//			toleranceTargetBox->show();

	fl_color(FL_WHITE); fl_rectf(10,10,w-w/5-10,h-20);
	fl_color(FL_BLACK); fl_rect(10,10,w-w/5-10,h-20);

	window->end();
	window->show();
}
