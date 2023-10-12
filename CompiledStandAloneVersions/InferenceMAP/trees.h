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

#ifndef QUADTREE_H_
#define QUADTREE_H_

#pragma once
#include <vector>

class QuadTree;

void generateTree(QuadTree *tree);
void generateSelectionTree(QuadTree *tree);
void tabulatePoints(QuadTree *tree);
void tabulateSelectionPoints(QuadTree *tree);
void recalibrateTree(QuadTree *tree);
void assignTree(QuadTree *tree);
void assignSelectionTree(QuadTree *tree);
void updateTree(QuadTree *tree);

struct XYI
{
	double x, y;
	int i;
	XYI(double X, double Y, int I) : x(X), y(Y), i(I) { }
	XYI() : x(0.0), y(0.0), i(0) { }
};

struct Rect
{
	double x, y, w, h;
	Rect() : x(0.0), y(0.0), w(0.0), h(0.0) { }
	Rect(double X, double Y, double W, double H) : x(X), y(Y), w(W), h(H) { }
//	Rect(SDL_Rect& R) : x(R.x), y(R.y), w(R.w), h(R.h) { }

	inline bool Contains(XYI& p) { return (p.x >= x && p.y >= y && p.x < x + w && p.y < y + h); }
	inline bool Intersects(Rect& r) { return !(r.x > (x + w) || (r.x + r.w) < x || r.y > (y + h) || (r.y + r.h) < y); }
};

// A single layer of a quad tree
class QuadTree {
	public:
		// Arbitrary constant to indicate how many elements can be stored in this quad tree node
		int capacity;

		// Axis-aligned bounding box stored as a center with half-dimensions
		// to represent the boundaries of this quad tree
		Rect bounds;

		//data inside
		XYI *points;
		int power;
		int identifier;
		double minSide;
		double xAverage,yAverage;

		QuadTree* nw;
		QuadTree* ne;
		QuadTree* sw;
		QuadTree* se;

		//Create a new quadtree
		QuadTree(Rect bounds, int count, double minSide);

		bool Insert(XYI p);
		bool Subdivide();
		void DelChildren();
		void Resize(Rect NewBounds);
		bool lastLevel();
		bool validSubdivision();

		// accessors
		int getCapacity() { return capacity; }
		int getCount() { return count; }
		int getIndex(int i) { return indexPointer[i]; }
		double getFx() { return xForce; }
		double getFy() { return yForce; }
		double getD() { return diffusion; }
		double getDlog() { return diffusionLog; }
		double getEp() { return potential; }
		double getFmag() { return forceMagnitude; }
		bool active() { return activated; }
		double getRandomValue() { return randomValue; }

		// mutators
		void setFx(double fx) { xForce = fx; }
		void setFy(double fy) { yForce = fy; }
		void setD(double d) { diffusion = d; }
		void setDlog(double d) { diffusionLog = d; }
		void setEp(double v) { potential = v; }
		void setFmag(double fm) { forceMagnitude = fm; }
		void setLocalization(double x, double y, int i, int p) {
			xPointer[p] = x;
			yPointer[p] = y;
			indexPointer[p] = i;
		}
		void activate() { activated = true; }
		void deactivate() { activated = false; }

		// Inference-related parameters
		double *xPointer;
		double *yPointer;
		int *indexPointer;
		int count;
		float rgb[3];
		double yForce;
		double xForce;
		double forceMagnitude;
		double potential;
		double diffusion;
		double diffusionLog;
		double gradVx;
		double gradVy;
		double gradDx;
		double gradDy;
		bool activated;
		double randomValue;
		float landscape;
		int vertexIndex;
		bool jeffreysPriorEnable;
		bool smoothingPriorEnable;

		int nLeftNeighbours;
		int nRightNeighbours;
		int nTopNeighbours;
		int nBottomNeighbours;

		// for smoothing prior
		double area;
		double xArea;
		double yArea;
		double variance;

		// for randomized optmization
		bool roActivated;
		int roIdentifier;
		bool roActive() { return roActivated; }
		void roActivate() { roActivated = true; }
		void roDeactivate() { roActivated = false; }
		void setRoIdentifier(int id) { roIdentifier = id; }
		int getIdentifier() { return identifier; }
		int getRoIdentifier() { return roIdentifier; }

		int getNumberOfNeighbours() { return nLeftNeighbours +
											 nRightNeighbours +
											 nTopNeighbours +
											 nBottomNeighbours;
		}
		int getNumberOfLeftNeighbours() { return nLeftNeighbours; }
		int getNumberOfRightNeighbours() { return nRightNeighbours; }
		int getNumberOfTopNeighbours() { return nTopNeighbours; }
		int getNumberOfBottomNeighbours() { return nBottomNeighbours; }
		float getLandscape() { return landscape; }
		int getVertexIndex() { return vertexIndex; }
		void setVertexIndex(int v) { vertexIndex = v; }

		void setLandscape(float l) { landscape = l; }
		//void setLandscapeVertices(int v) { landscape = v; }
		QuadTree **leftNeighbours;
		QuadTree **rightNeighbours;
		QuadTree **topNeighbours;
		QuadTree **bottomNeighbours;
		QuadTree* getLeftNeighbour(int q) { return leftNeighbours[q]; }
		QuadTree* getRightNeighbour(int q) { return rightNeighbours[q]; }
		QuadTree* getTopNeighbour(int q) { return topNeighbours[q]; }
		QuadTree* getBottomNeighbour(int q) { return bottomNeighbours[q]; }
		int leftProgress;
		int rightProgress;
		int topProgress;
		int bottomProgress;

		int nLeftNeighboursLandscape;
		int nRightNeighboursLandscape;
		int nTopNeighboursLandscape;
		int nBottomNeighboursLandscape;
		QuadTree **leftNeighboursLandscape;
		QuadTree **rightNeighboursLandscape;
		QuadTree **topNeighboursLandscape;
		QuadTree **bottomNeighboursLandscape;
		int leftProgressLandscape;
		int rightProgressLandscape;
		int topProgressLandscape;
		int bottomProgressLandscape;
		int getNumberOfLeftNeighboursLandscape() { return nLeftNeighboursLandscape; }
		int getNumberOfRightNeighboursLandscape() { return nRightNeighboursLandscape; }
		int getNumberOfTopNeighboursLandscape() { return nTopNeighboursLandscape; }
		int getNumberOfBottomNeighboursLandscape() { return nBottomNeighboursLandscape; }
		QuadTree* getLeftNeighbourLandscape(int q) { return leftNeighboursLandscape[q]; }
		QuadTree* getRightNeighbourLandscape(int q) { return rightNeighboursLandscape[q]; }
		QuadTree* getTopNeighbourLandscape(int q) { return topNeighboursLandscape[q]; }
		QuadTree* getBottomNeighbourLandscape(int q) { return bottomNeighboursLandscape[q]; }

		double getXMin() { return bounds.x; }
		double getXMax() { return bounds.x+bounds.w; }
		double getYMin() { return bounds.y; }
		double getYMax() { return bounds.y+bounds.h; }
		double getXCentre() { return bounds.x+bounds.w/2.0; }
		double getYCentre() { return bounds.y+bounds.h/2.0; }
		double getForceX() { return xForce; }
		double getForceY() { return yForce; }
		double getForce() { return forceMagnitude; }
		double getDiffusion() { return diffusion; }
		double getDiffusionLog() { return diffusionLog; }
		double getPotential() { return potential; }
		double getGradVx() { return gradVx; }
		double getGradVy() { return gradVy; }
		double getGradDx() { return gradDx; }
		double getGradDy() { return gradDy; }
		double getXCentroid() { return xAverage; }
		double getYCentroid() { return yAverage; }
		double getVariance() { return variance; }
		int getPower() { return power; }
		double getWidth() { return bounds.w; }
		double getHeight() { return bounds.h; }

		void applyPriors() { jeffreysPriorEnable = true; smoothingPriorEnable = true; }
		bool priorsApplied() { return (jeffreysPriorEnable || smoothingPriorEnable); }
		void applySmoothingPrior() { smoothingPriorEnable = true; }
//		void disableSmoothingPrior() { smoothingPriorEnable = false; }
		void applyJeffreysPrior() { jeffreysPriorEnable = true; }
//		void disableJeffreysPrior() { jeffreysPriorEnable = false; }
		void resetPriors() { jeffreysPriorEnable = false; smoothingPriorEnable = false; }

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

	private:
};

#endif /* QUADTREE_H_ */
