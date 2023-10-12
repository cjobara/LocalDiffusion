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

#include "trees.h"

#include "file.h"
#include "globals.h"
#include "annotations.h"
#include "draw.h"
#include "text.h"

extern File *file;
extern Globals *iMAP;

QuadTree::QuadTree(Rect Bounds, int n, double minSide) {
	roActivated = false;
	roIdentifier = -1;
	jeffreysPriorEnable = false;
	smoothingPriorEnable = false;
	leftProgressLandscape = 0;
	rightProgressLandscape = 0;
	topProgressLandscape = 0;
	bottomProgressLandscape = 0;
	nLeftNeighboursLandscape = 0;
	nRightNeighboursLandscape = 0;
	nTopNeighboursLandscape = 0;
	nBottomNeighboursLandscape = 0;
	leftNeighboursLandscape = NULL;
	rightNeighboursLandscape = NULL;
	topNeighboursLandscape = NULL;
	bottomNeighboursLandscape = NULL;
	variance = 0.0;
	area = 0.0;
	xArea = 0.0;
	yArea = 0.0;
	vertexIndex = 0;
	gradDx = 0.0;
	gradDy = 0.0;
	landscape = 0.0;
	xAverage = 0.0;
	yAverage = 0.0;
	this->minSide = minSide;
	identifier = -999;
	nLeftNeighbours = 0;
	nRightNeighbours = 0;
	nTopNeighbours = 0;
	nBottomNeighbours = 0;
	leftNeighbours = NULL;
	rightNeighbours = NULL;
	bottomNeighbours = NULL;
	topNeighbours = NULL;
	leftProgress = 0;
	rightProgress = 0;
	topProgress = 0;
	bottomProgress = 0;
	randomValue = (double)rand()/(double)RAND_MAX;
	activated = false;
	count = 0;
	forceMagnitude = 0.0;
	xForce = 0.0;
	yForce = 0.0;
	diffusion = 0.0;
	potential = 0.0;
	gradVx = 0.0;
	gradVy = 0.0;
	indexPointer = NULL;
	yPointer = NULL;
	xPointer = NULL;
	points = new XYI[n];
	capacity = n;
	bounds = Bounds;
	power = (int)(-logf(bounds.w/file->maxRange)/logf(2.0));
//	fprintf(stderr,"power = %i\n",power);
	for (int i = 0; i < n; i++) { points[i] = XYI(); }
//	sz = 0;
	ne = 0;
	se = 0;
	nw = 0;
	sw = 0;
	if (file->treeMesh->maxQuadTreePower < power) { file->treeMesh->maxQuadTreePower = power; }
}

bool QuadTree::Insert(XYI p) {
	//Ignore objects which are outside
	if (!bounds.Contains(p)) {
		return false;
	}

//	// if below capacity, add
//	if (count < capacity) {
//		// prevent subdivisions smaller than average step in trajectory
//		if (bounds.w > minSide) {
//			points[count] = p;
//			count++;
//		}
//		return true;
//	}

	// if below capacity, add
	if (count < capacity || bounds.w/2 < minSide) {

		count++;

		// have to resize point array
		if (count >= capacity) {
			// create temporary array
			XYI *pointsTemp = new XYI[count];
			// copy old data into temporary array
			for (int t = 0; t < count-1; t++) {
				pointsTemp[t] = points[t];
			}
			pointsTemp[count-1] = p;
			delete [] points;
			points = pointsTemp;
		} else {
			// prevent subdivisions smaller than average step in trajectory
			points[count-1] = p;
		}
		return true;
	}

	// Otherwise, we need to subdivide then add the point to whichever node will accept it
	if (nw == 0 && nw == 0 && sw == 0 && se == 0 && count >= capacity) {
		Subdivide();
	}

	if (nw->Insert(p) || ne->Insert(p) || sw->Insert(p) || se->Insert(p)) {
		return true;
	}

	// Otherwise, the point cannot be inserted for some unknown reason (which should never happen)
	return false;
}

bool QuadTree::Subdivide() {

	ne = new QuadTree(Rect(bounds.x + bounds.w/2.0, bounds.y, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	se = new QuadTree(Rect(bounds.x + bounds.w/2.0, bounds.y + bounds.h/2.0, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	nw = new QuadTree(Rect(bounds.x, bounds.y, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	sw = new QuadTree(Rect(bounds.x, bounds.y + bounds.h/2.0, bounds.w/2.0, bounds.h/2.0),capacity,minSide);

	for (int b = 0; b < count; b++) {
		if (ne->bounds.Contains(points[b])) {
			ne->points[ne->count] = points[b];
			ne->count++;
		}
		else if (nw->bounds.Contains(points[b])) {
			nw->points[nw->count] = points[b];
			nw->count++;
		}
		else if (se->bounds.Contains(points[b])) {
			se->points[se->count] = points[b];
			se->count++;
		}
		else if (sw->bounds.Contains(points[b])) {
			sw->points[sw->count] = points[b];
			sw->count++;
		}
	}

	return true;
}

void QuadTree::DelChildren() {
	if (nw != 0)
	{
		nw->DelChildren();
		delete nw;
		nw = 0;
	}
	if (ne != 0)
	{
		ne->DelChildren();
		delete ne;
		ne = 0;
	}
	if (sw != 0)
	{
		sw->DelChildren();
		delete sw;
		sw = 0;
	}
	if (se != 0)
	{
		se->DelChildren();
		delete se;
		se = 0;
	}
}

void QuadTree::Resize(Rect New) {
	bounds = New;
	ne = new QuadTree(Rect(bounds.x + bounds.w/2.0, bounds.y, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	se = new QuadTree(Rect(bounds.x + bounds.w/2.0, bounds.y + bounds.h/2.0, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	nw = new QuadTree(Rect(bounds.x, bounds.y, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
	sw = new QuadTree(Rect(bounds.x, bounds.y + bounds.h/2.0, bounds.w/2.0, bounds.h/2.0),capacity,minSide);
}

bool QuadTree::lastLevel() {
	return (nw == 0 && ne == 0 && sw == 0 && se == 0);
}

void generateTree(QuadTree *tree) {
	tree->DelChildren();

	for (int i = 0; i < file->localizationCount; i++) {
		tree->Insert(XYI(file->xPointer[i],file->yPointer[i],i));
	}

	recalibrateTree(tree);
	tabulatePoints(tree);
	assignTree(tree);
}

void generateSelectionTree(QuadTree *tree) {
	tree->DelChildren();

	for (int i = 0; i < file->treeMesh->selection.count; i++) {
		tree->Insert(XYI(file->treeMesh->selection.xPointer[i],
						 file->treeMesh->selection.yPointer[i],
						 file->treeMesh->selection.indexArray[i]));
	}
	recalibrateTree(tree);
	tabulateSelectionPoints(tree);

	assignSelectionTree(tree);

	// tabulate points
//	tabulatePoints(tree);
	// recalibrate to ensure all lowest-level leafs have greater than minimum points
//	recalibrateQuadTree(tree);
}

void tabulatePoints(QuadTree *tree) {
	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		tree->count = 0;
		XYI p;
		tree->xAverage = 0.0;
		tree->yAverage = 0.0;
		for (int i = 0; i < file->localizationCount; i++) {
			p.x = file->xPointer[i];
			p.y = file->yPointer[i];
			p.i = i;
			if ( tree->bounds.Contains(p) ) {
				tree->count++;
				tree->xAverage += p.x;
				tree->yAverage += p.y;
			}
		}
		if (tree->count > 0) {
			tree->xAverage /= tree->count;
			tree->yAverage /= tree->count;
		} else {
			tree->xAverage = tree->getXCentre();
			tree->yAverage = tree->getYCentre();
		}

		return;
	}
	tabulatePoints(tree->nw);
	tabulatePoints(tree->ne);
	tabulatePoints(tree->sw);
	tabulatePoints(tree->se);
	return;
}

void tabulateSelectionPoints(QuadTree *tree) {
	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		tree->count = 0;
		XYI p;
		tree->xAverage = 0.0;
		tree->yAverage = 0.0;
		for (int i = 0; i < file->treeMesh->selection.count; i++) {
			p.x = file->treeMesh->selection.xPointer[i];
			p.y = file->treeMesh->selection.yPointer[i];
			p.i = file->treeMesh->selection.indexArray[i];
			if ( tree->bounds.Contains(p) ) {
				tree->count++;
				tree->xAverage += p.x;
				tree->yAverage += p.y;
			}
		}
		if (tree->count > 0) {
			tree->xAverage /= tree->count;
			tree->yAverage /= tree->count;
		} else {
			tree->xAverage = tree->getXCentre();
			tree->yAverage = tree->getYCentre();
		}

		return;
	}
	tabulateSelectionPoints(tree->nw);
	tabulateSelectionPoints(tree->ne);
	tabulateSelectionPoints(tree->sw);
	tabulateSelectionPoints(tree->se);
	return;
}

void recalibrateTree(QuadTree *tree) {
	// recalibrate to ensure all quads have greater than minimum points

	if (tree->nw != 0 && tree->ne != 0 && tree->sw != 0 && tree->se != 0) {
		// ensure second-to-last level

		if (tree->nw->lastLevel() && tree->ne->lastLevel() && tree->sw->lastLevel() && tree->se->lastLevel()) {
			// if any children have less than minimum points, delete parent branch
			// delete any leafs which have widths smaller than the average trajectory step size
//			if (tree->nw->count <= file->treeMeshGui->getMinPoints() ||
//				tree->ne->count <= file->treeMeshGui->getMinPoints() ||
//				tree->sw->count <= file->treeMeshGui->getMinPoints() ||
//				tree->se->count <= file->treeMeshGui->getMinPoints() /*||
//				tree->nw->bounds.w <= tree->minSide ||
//				tree->ne->bounds.w <= tree->minSide ||
//				tree->sw->bounds.w <= tree->minSide ||
//				tree->se->bounds.w <= tree->minSide*/
//				) {
//
//				tree->DelChildren();
//
//				return;
//
//			}
			if (tree->nw->count <= 20 ||
				tree->ne->count <= 20 ||
				tree->sw->count <= 20 ||
				tree->se->count <= 20 /*||
				tree->nw->bounds.w <= tree->minSide ||
				tree->ne->bounds.w <= tree->minSide ||
				tree->sw->bounds.w <= tree->minSide ||
				tree->se->bounds.w <= tree->minSide*/
				) {

				tree->DelChildren();

				return;

			}
		}
		recalibrateTree(tree->nw);
		recalibrateTree(tree->ne);
		recalibrateTree(tree->sw);
		recalibrateTree(tree->se);
	}
	return;
}

void assignTree(QuadTree *tree) {
	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->getCount() > 0) {
//			fprintf(stderr,"count = %i\n",tree->getCount());
			tree->xPointer = new double[tree->getCount()];
			tree->yPointer = new double[tree->getCount()];
			tree->indexPointer = new int[tree->getCount()];

			XYI p;
			int progress = 0;
			tree->variance = 0.0;
			for (int i = 0; i < file->localizationCount; i++) {
				p.x = file->xPointer[i];
				p.y = file->yPointer[i];
				p.i = i;
				if ( tree->bounds.Contains(p) ) {
					tree->setLocalization(p.x,p.y,p.i,progress);
					progress++;
					// for variance calculation
					tree->variance += sqrt( (p.x-tree->xAverage)*(p.x-tree->xAverage)+(p.y-tree->yAverage)*(p.y-tree->yAverage) );
				}
			}
			tree->variance /= tree->count;
		}

		// calculate variance
//			for (int v = 0; v < tree->count; v++) {
//
//			}
		return;
	}
	assignTree(tree->nw);
	assignTree(tree->ne);
	assignTree(tree->sw);
	assignTree(tree->se);
	return;
}

void assignSelectionTree(QuadTree *tree) {
	// only count points in terminal leaf node
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->getCount() > 0) {
//			fprintf(stderr,"count = %i\n",tree->getCount());
			tree->xPointer = new double[tree->getCount()];
			tree->yPointer = new double[tree->getCount()];
			tree->indexPointer = new int[tree->getCount()];

			XYI p;
			int progress = 0;
			for (int i = 0; i < file->treeMesh->selection.count; i++) {
				p.x = file->treeMesh->selection.xPointer[i];
				p.y = file->treeMesh->selection.yPointer[i];
				p.i = file->treeMesh->selection.indexArray[i];
				if ( tree->bounds.Contains(p) ) {
					tree->setLocalization(p.x,p.y,p.i,progress);
					progress++;
				}
			}
		}
		return;
	}
	assignSelectionTree(tree->nw);
	assignSelectionTree(tree->ne);
	assignSelectionTree(tree->sw);
	assignSelectionTree(tree->se);
	return;
}
