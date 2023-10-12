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

#ifndef CELL_H_
#define CELL_H_

#include <stdio.h>
#include <iostream>
#include <stddef.h>
#include <math.h>

class SquareMesh;
class VoronoiMesh;

class SquareCell {

	SquareMesh *parent;

	float landscape;
	int identifier;
	int roIdentifier;
	int count;
	int totalCount;
	int row;
	int column;
	int depth;
	bool activated;
	double *xPointer;
	double *yPointer;
	double *dxPointer;
	double *dyPointer;
	double *dtPointer;
	double *iPointer;
	double *tPointer;
	int *indexPointer;
	double xMin,xMax;
	double yMin,yMax;
	double tMin,tMax;
	double xCentroid;
	double yCentroid;
	double zAverage;
	double xForce;
	double yForce;
	double zForce;
	double forceMagnitude;
	double diffusion;
	double diffusionLog;
	double potential;
	double speed;
	double gradVx;
	double gradVy;
	double gradVz;
	double gradDx;
	double gradDy;
	double gradDz;
	float area;
	float perimeter;
	float volume;
	double averageDx;
	double averageDy;
	double averageDt;
	bool roActivated;
	int roIteration;
	double variance;
	bool jeffreysPrior;
	bool smoothingPrior;

	int xNeighbourPos;
	int yNeighbourPos;
	int xNeighbours;
	int yNeighbours;

	public:
		SquareCell();
		SquareCell(SquareMesh *parent);

		// methods
		void initialize();

		// accessors
		bool active() { return activated; }
		bool roActive() { return roActivated; }
		int getRoIterations() { return roIteration; }
		int getCount() { return count; }
		int getTotalCount() { return totalCount; }
		int getIdentifier() { return identifier; }
		int getRoIdentifier() { return roIdentifier; }
		double getX(int i) { return xPointer[i]; }
		double getY(int i) { return yPointer[i]; }
		double getI(int i) { return iPointer[i]; }
		double getT(int i) { return tPointer[i]; }
		int getIndex(int i) { return indexPointer[i]; }
		double getDx();
		double getDy();
		double getDz();
		double getAverageDx() { return averageDx; }
		double getAverageDy() { return averageDy; }
		double getAverageDt() { return averageDt; }
		int getRow() { return row; }
		int getColumn() { return column; }
		double getXMin() { return xMin; }
		double getXMax() { return xMax; }
		double getXMean() { return (xMin+xMax)/2.0; }
		double getXCentre();
		double getYMin() { return yMin; }
		double getYMax() { return yMax; }
		double getYMean() { return (yMin+yMax)/2.0; }
		double getYCentre();
		double getForceX() { return xForce; }
		double getForceY() { return yForce; }
		double getForceMagnitude() { return forceMagnitude; }
		SquareMesh* getParent() { return parent; }
		double getDiffusion() { return diffusion; }
		double getDiffusionLog() { return diffusionLog; }
		double getPotential() { return potential; }
		double getGradVx() { return gradVx; }
		double getGradVy() { return gradVy; }
		double getGradDx() { return gradDx; }
		double getGradDy() { return gradDy; }
		int getXNeighbours() { return xNeighbours; }
		int getYNeighbours() { return yNeighbours; }
		int getXNeighbourPos() { return xNeighbourPos; }
		int getYNeighbourPos() { return yNeighbourPos; }
		float getArea() { return area; }
		float getPerimeter() { return perimeter; }
		float getVolume() { return volume; }
		double getXCentroid() { return xCentroid; }
		double getYCentroid() { return yCentroid; }
		double getZCentroid() { return zAverage; }
		float getLandscape() { return landscape; }
		double getSpeed() { return speed; }
		double getVariance() { return variance; }

		void applyPriors() { jeffreysPrior = true; smoothingPrior = true; }
		bool priorsApplied() { return (jeffreysPrior || smoothingPrior); }
		void unapplyJeffreysPrior() { jeffreysPrior = false; }
		void applyJeffreysPrior() { jeffreysPrior = true; }
		void applySmoothingPrior() { smoothingPrior = true; }
		void unapplySmoothingPrior() { smoothingPrior = false; }
		void resetPriors() { jeffreysPrior = false; smoothingPrior = false; }

		// mutators
		void roActivate() { roActivated = true; }
		void roDeactivate() { roActivated = false; }
		void setRoIterations(int it) { roIteration = it; }
		void setIdentifier(int counter) { identifier = counter; }
		void setRoIdentifier(int id) { roIdentifier = id; }
		void setXMin(double value) { xMin = value; }
		void setXMax(double value) { xMax = value; }
		void setYMin(double value) { yMin = value; }
		void setYMax(double value) { yMax = value; }
		void setTMin(double value) { tMin = value; }
		void setTMax(double value) { tMax = value; }
		void setRow(int r) { row = r; }
		void setColumn(int c) { column = c; }
		void setDepth(int d) { depth = d; }
		void setDetection(int cellIndex, int fileIndex);
		void activate() { activated = true; }
		void deactivate() { activated = false; }
		void increment() { count++; }
		void decrement() { count--; }
		void reset();
		void setForceX(double force) { xForce = force; }
		void setForceY(double force) { yForce = force; }
		void setZForce(double force) { zForce = force; }
		void setForceMagnitude(double force) { forceMagnitude = force; }
		void setDiffusion(double d) { diffusion = d; }
		void setDiffusionLog(double d) { diffusionLog = d; }
		void setPotential(double v) { potential = v; }
		void setGradVx(double gvx) { gradVx = gvx; }
		void setGradVy(double gvy) { gradVy = gvy; }
		void setGradVz(double gvz) { gradVy = gvz; }
		void setGradDx(double gdx) { gradDx = gdx; }
		void setGradDy(double gdy) { gradDy = gdy; }
		void setGradDz(double gdz) { gradDz = gdz; }
		void setXNeighbours(int n) { xNeighbours = n; }
		void setYNeighbours(int n) { yNeighbours = n; }
		void setXNeighbourPos(int v) { xNeighbourPos = v; }
		void setYNeighbourPos(int v) { yNeighbourPos = v; }
		void setLandscape(float l) { landscape = l; }
		void setSpeed(double s) { speed = s; }
		void setVariance(double g) { variance = g; }

		virtual ~SquareCell();

		friend class SquareMesh;
		friend class Simulation;

};

class VoronoiCell {
	VoronoiMesh *parent;

	int identifier;
	int cellId;
	int count;
	bool activated;
	double *xPointer;
	double *yPointer;
	double *iPointer;
	double *tPointer;
	int *indexPointer;
	double xForce;
	double yForce;
	double forceMagnitude;
	double diffusion;
	double diffusionLog;
	double potential;
	double gradVx;
	double gradVy;
	double gradDx;
	double gradDy;
	double energy;
	double xCentre;
	double yCentre;
	float xMean;
	float yMean;
	double xCentroid;
	double yCentroid;
	double variance;
	float area;
	float perimeter;
	int numberOfNeighbours;
	int *neighbourIndices;
	int nLeftNeighbours;
	int nRightNeighbours;
	int nTopNeighbours;
	int nBottomNeighbours;
	VoronoiCell **leftNeighbours;
	double *leftNeighbourWeights;
	VoronoiCell **rightNeighbours;
	double *rightNeighbourWeights;
	VoronoiCell **topNeighbours;
	double *topNeighbourWeights;
	VoronoiCell **bottomNeighbours;
	double *bottomNeighbourWeights;
	float randomValue;
	double averageDx;
	double averageDy;
	double averageDt;
	float landscape;
	bool jeffreysPrior;
	bool smoothingPrior;
	bool roActivated;
	int roIteration;
	int roIdentifier;

	// for smoothing prior
	double xArea;
	double yArea;

	// 0:bl ; 1:br ; 2:tr ; 3:tl
	int cornerCell;

	public:
		VoronoiCell();
		VoronoiCell(VoronoiMesh *parent);

		// drawing
		int nVertices;
		float *vertices;
		float *distanceToCentroid;

		// accessors
		VoronoiCell* getLeftNeighbours(int l) { return leftNeighbours[l]; }
		int getNLeftNeighbours() { return nLeftNeighbours; }
		double getLeftNeighbourWeight(int l) { return leftNeighbourWeights[l]; }
		VoronoiCell* getRightNeighbours(int r) { return rightNeighbours[r]; }
		int getNRightNeighbours() { return nRightNeighbours; }
		double getRightNeighbourWeight(int r) { return rightNeighbourWeights[r]; }
		VoronoiCell* getTopNeighbours(int t) { return topNeighbours[t]; }
		int getNTopNeighbours() { return nTopNeighbours; }
		double getTopNeighbourWeight(int t) { return topNeighbourWeights[t]; }
		VoronoiCell* getBottomNeighbours(int b) { return bottomNeighbours[b]; }
		int getNBottomNeighbours() { return nBottomNeighbours; }
		double getBottomNeighbourWeight(int b) { return bottomNeighbourWeights[b]; }

		bool active() { return activated; }
		bool roActive() { return roActivated; }
		double getVariance() { return variance; }
		float getLandscape() { return landscape; }
		int getCount() { return count; }
		int getIdentifier() { return identifier; }
		int getRoIdentifier() { return roIdentifier; }
		float getRandomValue() { return randomValue; }
		int getCellId() { return cellId; }
		double getX(int i) { return xPointer[i]; }
		double getY(int i) { return yPointer[i]; }
		double getI(int i) { return iPointer[i]; }
		double getT(int i) { return tPointer[i]; }
		int getIndex(int i) { return indexPointer[i]; }
		double getXCentre() { return xCentre; }
		double getYCentre() { return yCentre; }
		float getXMean() { return xMean; }
		float getYMean() { return yMean; }
		double getXCentroid() { return xCentroid; }
		double getYCentroid() { return yCentroid; }
		double getForceX() { return xForce; }
		double getForceY() { return yForce; }
		double getForceMagnitude() { return forceMagnitude; }
		VoronoiMesh* getParent() { return parent; }
		double getDiffusion() { return diffusion; }
		double getDiffusionLog() { return diffusionLog; }
		double getPotential() { return potential; }
		double getGradVx() { return gradVx; }
		double getGradVy() { return gradVy; }
		double getGradDx() { return gradDx; }
		double getGradDy() { return gradDy; }
		int getCornerCell() { return cornerCell; }
		float getArea() { return area; }
		float getPerimeter() { return perimeter; }
		int getNumberOfNeighbours() { return numberOfNeighbours; }
		int getNeighbour(int i) { return neighbourIndices[i]; }
		double getAreaX() { return xArea; }
		double getAreaY() { return yArea; }

		bool jeffreysPriorApplied() { return jeffreysPrior; }
		bool smoothingPriorApplied() { return smoothingPrior; }
		void applyPriors() { jeffreysPrior = true; smoothingPrior = true; }
		bool priorsApplied() { return (jeffreysPrior || smoothingPrior); }
		void unapplyJeffreysPrior() { jeffreysPrior = false; }
		void applyJeffreysPrior() { jeffreysPrior = true; }
		void applySmoothingPrior() { smoothingPrior = true; }
		void unapplySmoothingPrior() { smoothingPrior = false; }
		void resetPriors() { jeffreysPrior = false; smoothingPrior = false; }

		void roActivate() { roActivated = true; }
		void roDeactivate() { roActivated = false; }
		void setRoIterations(int it) { roIteration = it; }
		void setRoIdentifier(int id) { roIdentifier = id; }

		// mutators
		void setLandscape(float l) { landscape = l; }
		void setRandomValue(float r) { randomValue = r; }
		void setIdentifier(int counter) { identifier = counter; }
		void setCellId(int id) { cellId = id; }
		void setDetection(int cellIndex, int fileIndex);
		void activate() { activated = true; }
		void deactivate() { activated = false; }
		void increment() { count++; }
		void decrement() { count--; }
		void reset();
		void setForceX(double force) { xForce = force; }
		void setForceY(double force) { yForce = force; }
		void setForceMagnitude(double force) { forceMagnitude = force; }
		void setDiffusion(double d) { diffusion = d; }
		void setDiffusionLog(double d) { diffusionLog = d; }
		void setPotential(double v) { potential = v; }
		void setGradVx(double gvx) { gradVx = gvx; }
		void setGradVy(double gvy) { gradVy = gvy; }
		void setGradDx(double gdx) { gradDx = gdx; }
		void setGradDy(double gdy) { gradDy = gdy; }
		void setVariance(double v) { variance = v; }
		void setNVertices(int n) {
			nVertices = n;
			vertices = new float[2*nVertices];
			distanceToCentroid = new float[nVertices];
		}
		void setVertices(int index, float x, float y) {
			if (index >= 0 && index < nVertices) {
				vertices[2*index] = x;
				vertices[2*index+1] = y;
			}
		}
		void setCornerCell(int c) { cornerCell = c; }

		virtual ~VoronoiCell();

		friend class VoronoiMesh;
		friend class Simulation;
};

#endif /* CELL_H_ */
