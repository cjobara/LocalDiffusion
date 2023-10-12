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

#include <math.h>

#include "inference.h"
#include "optimization.h"
#include "file.h"

extern Globals *iMAP;
extern File *file;

#define PI 3.1415926535897932384626433832795028841971693993751

/********** SQUARE MESH **********/

double polynomialPosteriorSquare(double *coeff_f) {

	int i;
	int xPos, yPos;
	double x_debut, y_debut, D_loc, fx, fy;
	double result, dt, deltaa_x, deltaa_y, D_bruit;

	result = 0.0;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				deltaa_x = file->xPointer[i+1] - x_debut;
				deltaa_y = file->yPointer[i+1] - y_debut;

				fx = polynomialFxValueSquare(coeff_f,
					file->squareMesh->getCell(xPos,yPos)->getXCentre(),
					file->squareMesh->getCell(xPos,yPos)->getYCentre() );
				fy = polynomialFyValueSquare(coeff_f,
					file->squareMesh->getCell(xPos,yPos)->getXCentre(),
					file->squareMesh->getCell(xPos,yPos)->getYCentre() );

				D_loc  = polynomialDValueSquare(coeff_f, xPos , yPos);

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Jeffrey's Prior
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double polynomialPosteriorSquareSelection(double *coeff_f) {

	int i;
	int xPos, yPos;
	double x_debut, y_debut, D_loc, fx, fy;
	double result, dt, deltaa_x, deltaa_y, D_bruit;

	result = 0.0;

	#pragma omp for
	for (i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				deltaa_x = file->xPointer[ii+1] - x_debut;
				deltaa_y = file->yPointer[ii+1] - y_debut;
				fx = polynomialFxValueSquareSelection( coeff_f,file->squareMesh->selection.xMin+(double)(xPos+0.5)*file->squareMesh->getDx(),file->squareMesh->selection.yMin+(double)(yPos+0.5)*file->squareMesh->getDy() );
				fy = polynomialFyValueSquareSelection( coeff_f,file->squareMesh->selection.xMin+(double)(xPos+0.5)*file->squareMesh->getDx(),file->squareMesh->selection.yMin+(double)(yPos+0.5)*file->squareMesh->getDy() );

				D_loc  = polynomialDValueSquareSelection(coeff_f,xPos,yPos);

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
//				res -= log(4.0*PI*(D_loc + D_bruit)*dt);
//				res -= (pow((deltaa_x-D_loc*fx*dt), 2.0) - pow((deltaa_y-D_loc*fy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Jeffrey's Prior
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dPosteriorSquare(double *D) {

	double result = 0.0;

	const int a = file->squareMesh->getCurrentZoneX();
	const int b = file->squareMesh->getCurrentZoneY();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->getCellCount(a,b); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (file->squareMesh->getCell(a,b)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)] == file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)+1]) {

				const double dx = file->squareMesh->getDx(a,b,i);
				const double dy = file->squareMesh->getDy(a,b,i);
				const double dt = file->squareMesh->getDt(a,b,i);

				const double D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;

				result += - log(4.0*PI*(D[0]+D_bruit)*dt) - (dx*dx)/(4.0*(D[0]+D_bruit)*dt) - (dy*dy)/(4.0*(D[0]+D_bruit)*dt);

			}
		}
	}

	// Jeffrey's Prior
	if (file->squareMesh->getCell(a,b)->priorsApplied() == false) {
		if (file->squareMesh->jeffreysPriorEnabled()) {
			result += -log(D[0]*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
			file->squareMesh->getCell(a,b)->applyJeffreysPrior();
		}
	}

	return -result;
}

double dPosteriorSmoothingSquareSelection(double *D) {

	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dGradDxSquare(D,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dGradDySquare(D,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = D[file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;

}

double dPosteriorSmoothingSquare(double *D) {

	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dGradDxSquare(D,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dGradDySquare(D,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = D[file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}
	return -result;
}

double dGradDxSquare(double *D, int x, int y) {

	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = D[file->squareMesh->getIdentifier(x-1,y)];
		q[1] = D[file->squareMesh->getIdentifier(x,y)];
		q[2] = D[file->squareMesh->getIdentifier(x+1,y)];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = D[file->squareMesh->getIdentifier(x,y)];
			q[2] = D[file->squareMesh->getIdentifier(x+1,y)];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = D[file->squareMesh->getIdentifier(x-1,y)];
			q[1] = D[file->squareMesh->getIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dGradDySquare(double *D, int x, int y) {

	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = D[file->squareMesh->getIdentifier(x,y-1)];
		q[1] = D[file->squareMesh->getIdentifier(x,y)];
		q[2] = D[file->squareMesh->getIdentifier(x,y+1)];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = D[file->squareMesh->getIdentifier(x,y)];
			q[2] = D[file->squareMesh->getIdentifier(x,y+1)];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = D[file->squareMesh->getIdentifier(x,y-1)];
			q[1] = D[file->squareMesh->getIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;

}

double dPosteriorSmoothingSquareRandomizedOptimization(double *D) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dGradDxSquareRandomizedOptimization(D,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dGradDySquareRandomizedOptimization(D,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = D[file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dPosteriorSmoothingSquareRandomizedOptimizationSelection(double *D) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dGradDxSquareRandomizedOptimization(D,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dGradDySquareRandomizedOptimization(D,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = D[file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dGradDxSquareRandomizedOptimization(double *D, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? D[file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
		q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? D[file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? D[file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? D[file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
			q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dGradDySquareRandomizedOptimization(double *D, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? D[file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
		q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? D[file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? D[file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? D[file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
			q[1] = D[file->squareMesh->getRoIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double ddrPosteriorSquare(double *DrxDryD) {

	double result = 0.0;

	const int a = file->squareMesh->getCurrentZoneX();
	const int b = file->squareMesh->getCurrentZoneY();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->getCellCount(a,b); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (file->squareMesh->getCell(a,b)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)] == file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)+1]) {

				const double dx = file->squareMesh->getDx(a,b,i);
				const double dy = file->squareMesh->getDy(a,b,i);
				const double dt = file->squareMesh->getDt(a,b,i);

				const double D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;

//				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt) - ((dx-FxFyD[2]*FxFyD[0]*dt)*(dx-FxFyD[2]*FxFyD[0]*dt))/(4.0*(FxFyD[2]+D_bruit)*dt) - ((dy-FxFyD[2]*FxFyD[1]*dt)*(dy-FxFyD[2]*FxFyD[1]*dt))/(4.0*(FxFyD[2]+D_bruit)*dt);
				result += - log(4.0*PI*(DrxDryD[2]+D_bruit)*dt) - (dx-DrxDryD[0]*dt)*(dx-DrxDryD[0]*dt)/(4*(DrxDryD[2]+D_bruit)*dt) - (dy-DrxDryD[1]*dt)*(dy-DrxDryD[1]*dt)/(4*(DrxDryD[2]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->squareMesh->jeffreysPriorEnabled()) {
		result += - 2.0*log(DrxDryD[2]*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
	}
	return -result;
}

double ddrPosteriorSmoothingSquare(double *DDrxDry) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc, Drx, Dry;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(ddrGradDxSquare(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(ddrGradDySquare(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];
				Drx  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+1];
				Dry  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+2];

//				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt),2.0)+pow((dy-D_loc*gradVy*dt),2.0))/(4.0*(D_loc+D_bruit)*dt);
//				result += - log(4.0*PI*(DrxDryD[2]+D_bruit)*dt) - (dx-DrxDryD[0]*dt)*(dx-DrxDryD[0]*dt)/(4*(DrxDryD[2]+D_bruit)*dt) - (dy-DrxDryD[1]*dt)*(dy-DrxDryD[1]*dt)/(4*(DrxDryD[2]+D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingSquareSelection(double *DDrxDry) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc, Drx, Dry;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dfGradDxSquare(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dfGradDySquare(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];
				Drx  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+1];
				Dry  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+2];

//				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Drx*dt)*(dx-D_loc*Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Dry*dt)*(dy-D_loc*Dry*dt))/(4.0*(D_loc+D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrGradDxSquare(double *DDrxDry, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = DDrxDry[3*file->squareMesh->getIdentifier(x-1,y)];
		q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
		q[2] = DDrxDry[3*file->squareMesh->getIdentifier(x+1,y)];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
			q[2] = DDrxDry[3*file->squareMesh->getIdentifier(x+1,y)];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = DDrxDry[3*file->squareMesh->getIdentifier(x-1,y)];
			q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double ddrGradDySquare(double *DDrxDry, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = DDrxDry[3*file->squareMesh->getIdentifier(x,y-1)];
		q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
		q[2] = DDrxDry[3*file->squareMesh->getIdentifier(x,y+1)];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
			q[2] = DDrxDry[3*file->squareMesh->getIdentifier(x,y+1)];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = DDrxDry[3*file->squareMesh->getIdentifier(x,y-1)];
			q[1] = DDrxDry[3*file->squareMesh->getIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}


double ddrPosteriorSmoothingSquareRandomizedOptimization(double *DDrxDry) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit, Drx, Dry;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(ddrGradDxSquareRandomizedOptimization(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(ddrGradDySquareRandomizedOptimization(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];
				Drx  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+1];
				Dry  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DDrxDry) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit, Drx, Dry;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(ddrGradDxSquareRandomizedOptimization(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(ddrGradDySquareRandomizedOptimization(DDrxDry,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];
				Drx  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+1];
				Dry  = DDrxDry[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double ddrGradDxSquareRandomizedOptimization(double *DDrxDry, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
		q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
			q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double ddrGradDySquareRandomizedOptimization(double *DDrxDry, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
		q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DDrxDry[3*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
			q[1] = DDrxDry[3*file->squareMesh->getRoIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dfPosteriorSquare(double *FxFyD){

	double result = 0.0;

	const int a = file->squareMesh->getCurrentZoneX();
	const int b = file->squareMesh->getCurrentZoneY();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->getCellCount(a,b); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (file->squareMesh->getCell(a,b)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)] == file->nPointer[file->squareMesh->getCell(a,b)->getIndex(i)+1]) {

				const double dx = file->squareMesh->getDx(a,b,i);
				const double dy = file->squareMesh->getDy(a,b,i);
				const double dt = file->squareMesh->getDt(a,b,i);

				const double D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;

				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt) - ((dx-FxFyD[2]*FxFyD[0]*dt)*(dx-FxFyD[2]*FxFyD[0]*dt))/(4.0*(FxFyD[2]+D_bruit)*dt) - ((dy-FxFyD[2]*FxFyD[1]*dt)*(dy-FxFyD[2]*FxFyD[1]*dt))/(4.0*(FxFyD[2]+D_bruit)*dt);

			}
		}
	}

	// Jeffrey's Prior
	if (file->squareMesh->jeffreysPriorEnabled()) {
		result += 2.0*log(FxFyD[2]) - 2.0*log(FxFyD[2]*file->squareMesh->getDtMean() + file->squareMesh->getSigma()*file->squareMesh->getSigma());
	}
	return -result;
}

double dfPosteriorSmoothingSquare(double *DFxFy) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc, Fx, Fy;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dfGradDxSquare(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dfGradDySquare(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];
				Fx  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+1];
				Fy  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+2];

//				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt),2.0)+pow((dy-D_loc*gradVy*dt),2.0))/(4.0*(D_loc+D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingSquareSelection(double *DFxFy) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc, Fx, Fy;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dfGradDxSquare(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dfGradDySquare(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];
				Fx  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+1];
				Fy  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getIdentifier()+2];

//				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfGradDxSquare(double *DFxFy, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = DFxFy[3*file->squareMesh->getIdentifier(x-1,y)];
		q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
		q[2] = DFxFy[3*file->squareMesh->getIdentifier(x+1,y)];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
			q[2] = DFxFy[3*file->squareMesh->getIdentifier(x+1,y)];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = DFxFy[3*file->squareMesh->getIdentifier(x-1,y)];
			q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dfGradDySquare(double *DFxFy, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = DFxFy[3*file->squareMesh->getIdentifier(x,y-1)];
		q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
		q[2] = DFxFy[3*file->squareMesh->getIdentifier(x,y+1)];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
			q[2] = DFxFy[3*file->squareMesh->getIdentifier(x,y+1)];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = DFxFy[3*file->squareMesh->getIdentifier(x,y-1)];
			q[1] = DFxFy[3*file->squareMesh->getIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double polynomialFxValueSquare(double *coeffff,  double x, double y ){
  int k,l,indiceb;
  double valeur1;

  indiceb = 0;
  valeur1 = 0.0;
  for ( k=1; k <= file->squareMesh->getPolynomialOrder(); k++){
    for (l=0; l<= k; l++){
      valeur1 +=  -coeffff[indiceb]*(k-l)*pow(x-file->squareMesh->getXCentre(),FMAX((double) k-l-1.,0.0))*pow(y-file->squareMesh->getYCentre(),FMAX((double) l,0.0));
      indiceb++;
    }
  }

  return valeur1;
}

double polynomialFyValueSquare(double *coefffff, double x, double y ){
  int k,l,indicec;
  double valeur2;

  indicec = 0;
  valeur2 = 0.0;
  for (k = 1; k <= file->squareMesh->getPolynomialOrder(); k++){
    for (l=0; l <= k; l++){
      valeur2 +=  -coefffff[indicec]*(l)*pow(x-file->squareMesh->getXCentre(),FMAX((double) k-l,0.0))*pow(y-file->squareMesh->getYCentre(),FMAX((double) l-1.,0.0));
      indicec++;
    }
  }

  return valeur2;
}

double polynomialDValueSquare(double *coeffff, int x, int y) {
  double valeur3;

  if (file->squareMesh->active(x,y)) {
	  valeur3 = coeffff[file->squareMesh->getCoefficients() + file->squareMesh->getIdentifier(x,y)];
  }
  else {
	  valeur3 = 0.0;
  }

  return valeur3;
}

double polynomialVValueSquare(double *coeff, double x, double y ) {
	int k,l,indicee;
	double valeur0;

	indicee = 0;
	valeur0 =0.;
	for (k=1; k<= file->squareMesh->getPolynomialOrder();k++) {
		for (l=0; l<= k; l++){
			valeur0 +=  coeff[indicee]*pow((x-file->getXCentre())*1.0e-6,(double)(k-l))*pow((y-file->getYCentre())*1.0e-6,(double)l)*4.e-9*pow(10.,(double)(6*(k-2)));
			indicee += 1;
		}
	}

	return valeur0/4.e-21;
}

double polynomialFxValueSquareSelection(double *coeffff,  double x, double y ){
  int k,l,indiceb;
  double valeur1;

  indiceb = 0;
  valeur1 = 0.0;
  for ( k=1; k <= file->squareMesh->getPolynomialOrder(); k++){
    for (l=0; l<= k; l++){
      valeur1 +=  -coeffff[indiceb]*(k-l)*pow(x-file->squareMesh->getXCentre(),FMAX((double) k-l-1.,0.0))*pow(y-file->squareMesh->getYCentre(),FMAX((double) l,0.0));
      indiceb++;
    }
  }

  return valeur1;
}

double polynomialFyValueSquareSelection(double *coefffff, double x, double y ){
  int k,l,indicec;
  double valeur2;

  indicec = 0;
  valeur2 = 0.0;
  for (k = 1; k <= file->squareMesh->getPolynomialOrder(); k++){
    for (l=0; l <= k; l++){
      valeur2 +=  -coefffff[indicec]*(l)*pow(x-file->squareMesh->getXCentre(),FMAX((double) k-l,0.0))*pow(y-file->squareMesh->getYCentre(),FMAX((double) l-1.,0.0));
      indicec++;
    }
  }

  return valeur2;
}

double polynomialDValueSquareSelection(double *coeffff, int x, int y) {
  double valeur3;

  if (file->squareMesh->active(x,y)) {
	  valeur3 = coeffff[file->squareMesh->getCoefficients() + file->squareMesh->getIdentifier(x,y)];
  }
  else {
	  valeur3 = 0.0;
  }

  return valeur3;
}

double polynomialVValueSquareSelection(double *coeff, double x, double y ) {
	int k,l,indicee;
	double valeur0;

	indicee = 0;
	valeur0 =0.;
	for (k=1; k<= file->squareMesh->getPolynomialOrder();k++) {
		for (l=0; l<= k; l++){
			valeur0 +=  coeff[indicee]*pow((x-file->squareMesh->getXCentre())*1.0e-6,(double)(k-l))*pow((y-file->squareMesh->getYCentre())*1.0e-6,(double)l)*4.e-9*pow(10.,(double)(6*(k-2)));
			indicee += 1;
		}
	}

	return valeur0/4.e-21;
}

double dfFxValueSquare(double *V_ij, int x, int y ) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y) == 2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = V_ij[ file->squareMesh->getIdentifier(x-1,y) ];
		q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
		q[2] = V_ij[ file->squareMesh->getIdentifier(x+1,y) ];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];

	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
			q[2] = V_ij[ file->squareMesh->getIdentifier(x+1,y) ];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = V_ij[ file->squareMesh->getIdentifier(x-1,y) ];
			q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}

	return valeur1;
}

double dfFyValueSquare(double *V_ij, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = V_ij[ file->squareMesh->getIdentifier(x,y-1) ];
		q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
		q[2] = V_ij[ file->squareMesh->getIdentifier(x,y+1) ];
		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1) {

			q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
			q[2] = V_ij[ file->squareMesh->getIdentifier(x,y+1) ];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {

			q[0] = V_ij[ file->squareMesh->getIdentifier(x,y-1) ];
			q[1] = V_ij[ file->squareMesh->getIdentifier(x,y) ];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}

//	fprintf(stderr,"[%i,%i]\t%f\n",x,y,valeur2);
	return valeur2;
}

//double dfFxValueSquareSelection(double *V_ij, int x, int y);
//double dfFyValueSquareSelection(double *V_ij, int x, int y);

double squareDifferenceSquare(double *V_ij) {
	double squareDifference = 0.0;

	/* to pass to kernel
	 * squareDifference
	 * file->regularMesh->getXCells()
	 * file->regularMesh->getYCells()
	 * dfDirectFxValueRegular
	 * dfDirectFyValueRegular
	 */

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dfFxValueSquare(V_ij,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dfFyValueSquare(V_ij,a,b));
			}
		}
	}

	#pragma omp for
	for (int i = 0; i < file->squareMesh->getXCells(); i++) {
		for (int j = 0; j < file->squareMesh->getYCells(); j++) {
			if ( (file->squareMesh->active(i,j)) /*&& (file->regularMesh->getXNeighbours(i,j)!=0) && (file->regularMesh->getYNeighbours(i,j)!=0)*/ ) {
				const double fx = file->squareMesh->getCell(i,j)->getGradVx();
				const double fy = file->squareMesh->getCell(i,j)->getGradVy();
				squareDifference += pow(fabs(fx-file->squareMesh->getForceX(i,j)),2.0);
				squareDifference += file->squareMesh->getBeta()*fx*fx;
				squareDifference += pow(fabs(fy-file->squareMesh->getForceY(i,j)),2.0);
				squareDifference += file->squareMesh->getBeta()*fy*fy;
			}
		}
	}

	return squareDifference;
}
double dfPosteriorSmoothingSquareRandomizedOptimization(double *DFxFy) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit, Fx, Fy;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dfGradDxSquareRandomizedOptimization(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dfGradDySquareRandomizedOptimization(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];
				Fx  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+1];
				Fy  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+2];

//				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt),2.0)+pow((dy-D_loc*gradVy*dt),2.0))/(4.0*(D_loc+D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DFxFy) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit, Fx, Fy;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradDx(dfGradDxSquareRandomizedOptimization(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dfGradDySquareRandomizedOptimization(DFxFy,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];
				Fx  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+1];
				Fy  = DFxFy[3*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()+2];

//				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt),2.0)+pow((dy-D_loc*gradVy*dt),2.0))/(4.0*(D_loc+D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dfGradDxSquareRandomizedOptimization(double *DFxFy, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
		q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
			q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dfGradDySquareRandomizedOptimization(double *DFxFy, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
		q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DFxFy[3*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
			q[1] = DFxFy[3*file->squareMesh->getRoIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dvPosteriorSquare(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Jeffrey's Prior
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSquareSelection(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Jeffrey's Prior
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSquareRandomizedOptimization(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
					}
					file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingSquareRandomizedOptimization(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDx(dvGradDxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dvGradDySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double lambda = (double)file->squareMeshGui->lambdaSlider->value();
	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();
				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*area*(gradVx*gradVx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									 	   gradVy*gradVy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dvPosteriorSquareRandomizedOptimizationSelection(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingSquareRandomizedOptimizationSelection(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->roActive(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDx(dvGradDxSquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dvGradDySquareRandomizedOptimization(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double lambda = (double)file->squareMeshGui->lambdaSlider->value();
	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->roActive(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();
				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getRoIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*area*(gradVx*gradVx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									 	   gradVy*gradVy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingSquare(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDx(dvGradDxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dvGradDySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double lambda = (double)file->squareMeshGui->lambdaSlider->value();
	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			x_debut = file->xPointer[i];
			y_debut = file->yPointer[i];

			xPos = (int)floor( (x_debut - file->xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - x_debut;
				dy = file->yPointer[i+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();
				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt),2.0)+pow((dy-D_loc*gradVy*dt),2.0))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*area*(gradVx*gradVx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									 	   gradVy*gradVy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingSquareSelection(double *DV) {
	int xPos, yPos;
	double x_debut, y_debut, D_loc;
	double result, dt, dx, dy, D_bruit;

	result = 0.0;

	// initialize dvGrad array
	#pragma omp for
	for (int a = 0; a < file->squareMesh->getXCells(); a++) {
		for (int b = 0; b < file->squareMesh->getYCells(); b++) {
			if (file->squareMesh->active(a,b)) {
				file->squareMesh->getCell(a,b)->setGradVx(dvGradVxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradVy(dvGradVySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDx(dvGradDxSquare(DV,a,b));
				file->squareMesh->getCell(a,b)->setGradDy(dvGradDySquare(DV,a,b));
				file->squareMesh->getCell(a,b)->resetPriors();
			}
		}
	}

	const double lambda = (double)file->squareMeshGui->lambdaSlider->value();
	const double mu = (double)file->squareMeshGui->muSlider->value();
	const double area = file->squareMesh->getDx()*file->squareMesh->getDx();

	#pragma omp for
	for (int i = 0; i < file->squareMesh->selection.count - 1; i++) {
		const int ii = file->squareMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			x_debut = file->xPointer[ii];
			y_debut = file->yPointer[ii];

			xPos = (int)floor( (x_debut - file->squareMesh->selection.xMin)/file->squareMesh->getDx() );
			yPos = (int)floor( (y_debut - file->squareMesh->selection.yMin)/file->squareMesh->getDy() );

			if (file->squareMesh->active(xPos,yPos)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - x_debut;
				dy = file->yPointer[ii+1] - y_debut;

				const double gradVx = file->squareMesh->getCell(xPos,yPos)->getGradVx();
				const double gradVy = file->squareMesh->getCell(xPos,yPos)->getGradVy();
				const double gradDx = file->squareMesh->getCell(xPos,yPos)->getGradDx();
				const double gradDy = file->squareMesh->getCell(xPos,yPos)->getGradDy();

				D_bruit = file->squareMesh->getSigma()*file->squareMesh->getSigma()/dt;
				D_loc  = DV[2*file->squareMesh->getCell(xPos,yPos)->getIdentifier()];

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->squareMesh->getCell(xPos,yPos)->priorsApplied() == false) {
					// Jeffreys
					if (file->squareMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->squareMesh->getSigma()*file->squareMesh->getSigma());
						file->squareMesh->getCell(xPos,yPos)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*area*(gradVx*gradVx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									 	   gradVy*gradVy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					result -= mu*area*(gradDx*gradDx*(double)file->squareMesh->getXNeighbours(xPos,yPos) +
									   gradDy*gradDy*(double)file->squareMesh->getYNeighbours(xPos,yPos));
					file->squareMesh->getCell(xPos,yPos)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvGradVxSquare(double *DV, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)+1];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
		q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)+1];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
			q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)+1];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)+1];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dvGradVySquare(double *DV, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)+1];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
		q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)+1];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
			q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)+1];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)+1];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)+1];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dvGradDxSquare(double *DV, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
		q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dvGradDySquare(double *DV, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
		q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dvGradDxRegular(const double *DV, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
		q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)];

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			q[2] = DV[2*file->squareMesh->getIdentifier(x+1,y)];
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = DV[2*file->squareMesh->getIdentifier(x-1,y)];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dvGradDyRegular(const double *DV, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)];
		q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
		q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)];

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			q[2] = DV[2*file->squareMesh->getIdentifier(x,y+1)];
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = DV[2*file->squareMesh->getIdentifier(x,y-1)];
			q[1] = DV[2*file->squareMesh->getIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dvGradVxSquareRandomizedOptimization(double *DV, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y) == 2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x-1,y)+1] : file->squareMesh->getCell(x-1,y)->getPotential();
		q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
		q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x+1,y)+1] : file->squareMesh->getCell(x+1,y)->getPotential();

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
			q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x+1,y)+1] : file->squareMesh->getCell(x+1,y)->getPotential();
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x-1,y)+1] : file->squareMesh->getCell(x-1,y)->getPotential();
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dvGradVySquareRandomizedOptimization(double *DV, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y-1)+1] : file->squareMesh->getCell(x,y-1)->getPotential();
		q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
		q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y+1)+1] : file->squareMesh->getCell(x,y+1)->getPotential();

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
			q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y+1)+1] : file->squareMesh->getCell(x,y+1)->getPotential();
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y-1)+1] : file->squareMesh->getCell(x,y-1)->getPotential();
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)+1];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

double dvGradDxSquareRandomizedOptimization(double *DV, int x, int y) {
	double xx[3],w[3],q[3];
	double valeur1 = 0.0;

	if (file->squareMesh->getXNeighbours(x,y)==2) {
		xx[0] = file->squareMesh->getCell(x-1,y)->getXCentre();
		xx[1] = file->squareMesh->getCell(x,y)->getXCentre();
		xx[2] = file->squareMesh->getCell(x+1,y)->getXCentre();

		q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
		q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();

		vander(xx, w, q, 2);

		valeur1 = - w[1] - 2.0*w[2]*xx[1];
	} else {
		if (file->squareMesh->getXNeighbourPos(x,y) == 1) {

			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x+1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x+1,y)] : file->squareMesh->getCell(x+1,y)->getDiffusion();
			valeur1 = -(q[2] -q[1] )/(file->squareMesh->getDx());

		} else if (file->squareMesh->getXNeighbourPos(x,y) == -1) {

			q[0] = file->squareMesh->getCell(x-1,y)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x-1,y)] : file->squareMesh->getCell(x-1,y)->getDiffusion();
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
			valeur1 = -(q[1] -q[0] )/(file->squareMesh->getDx());
		}
	}
	return valeur1;
}

double dvGradDySquareRandomizedOptimization(double *DV, int x, int y) {
	double yy[3], w[3], q[3] ;
	double valeur2 = 0.0;
	if (file->squareMesh->getYNeighbours(x,y) == 2) {

		yy[0] = file->squareMesh->getCell(x,y-1)->getYCentre();
		yy[1] = file->squareMesh->getCell(x,y)->getYCentre();
		yy[2] = file->squareMesh->getCell(x,y+1)->getYCentre();

		q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
		q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
		q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();

		vander(yy, w, q, 2);

		valeur2 = - w[1] - 2*w[2]*yy[1];
	} else {
		if (file->squareMesh->getYNeighbourPos(x,y) == 1){
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
			q[2] = file->squareMesh->getCell(x,y+1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y+1)] : file->squareMesh->getCell(x,y+1)->getDiffusion();
			valeur2 = -(q[2] -q[1] )/(file->squareMesh->getDy());
		} else if (file->squareMesh->getYNeighbourPos(x,y) == -1) {
			q[0] = file->squareMesh->getCell(x,y-1)->roActive() ? DV[2*file->squareMesh->getRoIdentifier(x,y-1)] : file->squareMesh->getCell(x,y-1)->getDiffusion();
			q[1] = DV[2*file->squareMesh->getRoIdentifier(x,y)];
			valeur2 = -(q[1] -q[0] )/(file->squareMesh->getDy());
		}
	}
	return valeur2;
}

/********** VORONOI MESH **********/

double polynomialPosteriorVoronoi(double *coeff_f) {

	double D_loc, D_bruit, fx, fy;
	double result, dt, deltaa_x, deltaa_y;
	int idx;

	result = 0.0;

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {
			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				deltaa_x = file->xPointer[i+1] - file->xPointer[i];
				deltaa_y = file->yPointer[i+1] - file->yPointer[i];

				fx = polynomialFxValueVoronoi(coeff_f,
									file->voronoiMesh->getCell(idx)->getXCentre(),
									file->voronoiMesh->getCell(idx)->getYCentre() );
				fy = polynomialFyValueVoronoi(coeff_f,
									file->voronoiMesh->getCell(idx)->getXCentre(),
									file->voronoiMesh->getCell(idx)->getYCentre() );

				D_loc  = polynomialDValueVoronoi(coeff_f,idx);

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double polynomialPosteriorVoronoiSelection(double *coeff_f) {

	int i;
	double D_loc, D_bruit, fx, fy;
	double result, dt, deltaa_x, deltaa_y;
	int idx;

	result = 0.0;

	#pragma omp for
	for (i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				deltaa_x = file->xPointer[ii+1] - file->xPointer[ii];
				deltaa_y = file->yPointer[ii+1] - file->yPointer[ii];

				fx = polynomialFxValueVoronoiSelection(coeff_f,file->voronoiMesh->getCell(idx)->getXCentre(),file->voronoiMesh->getCell(idx)->getYCentre() );
				fy = polynomialFyValueVoronoiSelection(coeff_f,file->voronoiMesh->getCell(idx)->getXCentre(),file->voronoiMesh->getCell(idx)->getYCentre() );

				D_loc  = polynomialDValueVoronoiSelection(coeff_f,idx);

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double polynomialFxValueVoronoi(double *coeffff,  double x, double y ) {
  int k,l,indiceb;
  double valeur1;

  indiceb = 0;
  valeur1=0.;
  for ( k=1; k <= file->voronoiMesh->getPolynomialOrder(); k++){
    for (l=0; l<= k; l++){
      valeur1 +=  -coeffff[indiceb]*(k-l)*pow(x-file->voronoiMesh->getXCentre(),(double) FMAX((double) (k-l-1),0.0))*pow(y-file->voronoiMesh->getYCentre(),(double) FMAX((double) l,0.0));
      indiceb += 1;
    }
  }

  return valeur1;
}

double polynomialFyValueVoronoi(double *coefffff, double x, double y ) {
  int k,l,indicec;
  double valeur2;

  indicec = 0;
  valeur2 =0.;
  for (k = 1; k <= file->voronoiMesh->getPolynomialOrder(); k++){
    for (l=0; l <= k; l++){
      valeur2 +=  -coefffff[indicec]*(l)*pow(x-file->voronoiMesh->getXCentre(),(double) FMAX((double)(k-l),0.0))*pow(y-file->voronoiMesh->getXCentre(),(double) FMAX((double)(l-1),0.0));
      indicec += 1;
    }
  }

  return valeur2;
}

double polynomialDValueVoronoi(double *coeffff, int i) {
  double valeur3;

  if (file->voronoiMesh->active(i)) {
	  valeur3 = coeffff[file->voronoiMesh->getCoefficients() + file->voronoiMesh->getIdentifier(i)];
  }
  else {
	  valeur3 = 0.0;
  }

  return valeur3;
}

double polynomialVValueVoronoi(double *coeff, double x, double y ) {
	int k,l,indicee;
	double valeur0;

	indicee = 0;
	valeur0 =0.;
	for (k=1; k<= file->voronoiMesh->getPolynomialOrder();k++) {
		for (l=0; l<= k; l++){
			valeur0 +=  coeff[indicee]*pow((x-file->getXCentre())*1e-6,double(k-l))*pow((y-file->getYCentre())*1e-6,l)*4.e-9*pow(10.0,(double)(6*(k-2)));
			indicee += 1;
		}
	}

	return valeur0/4.e-21;
}

double polynomialFxValueVoronoiSelection(double *coeffff,  double x, double y ) {
  int k,l,indiceb;
  double valeur1;

  indiceb = 0;
  valeur1=0.;
  for ( k=1; k <= file->voronoiMesh->getPolynomialOrder(); k++){
    for (l=0; l<= k; l++){
      valeur1 +=  -coeffff[indiceb]*(k-l)*pow(x-(file->voronoiMesh->selection.xMax-file->voronoiMesh->selection.xMin)/2.0,(double) FMAX((double) (k-l-1),0.0))*pow(y-(file->voronoiMesh->selection.yMax-file->voronoiMesh->selection.yMin)/2.0,(double) FMAX((double) l,0.0));
      indiceb += 1;
    }
  }

  return valeur1;
}

double polynomialFyValueVoronoiSelection(double *coefffff, double x, double y ) {
  int k,l,indicec;
  double valeur2;

  indicec = 0;
  valeur2 =0.;
  for (k = 1; k <= file->voronoiMesh->getPolynomialOrder(); k++){
    for (l=0; l <= k; l++){
      valeur2 +=  -coefffff[indicec]*(l)*pow(x-(file->voronoiMesh->selection.xMax-file->voronoiMesh->selection.xMin)/2.0,(double) FMAX((double)(k-l),0.0))*pow(y-(file->voronoiMesh->selection.yMax-file->voronoiMesh->selection.yMin)/2.0,(double) FMAX((double)(l-1),0.0));
      indicec += 1;
    }
  }

  return valeur2;
}

double polynomialDValueVoronoiSelection(double *coeffff, int i) {
  double valeur3;

  if (file->voronoiMesh->active(i)) {
	  valeur3 = coeffff[file->voronoiMesh->getCoefficients() + file->voronoiMesh->getIdentifier(i)];
  }
  else {
	  valeur3 = 0.0;
  }

  return valeur3;
}

double polynomialVValueVoronoiSelection(double *coeff, double x, double y ) {
	int k,l,indicee;
	double valeur0;

	indicee = 0;
	valeur0 =0.;
	for (k=1; k<= file->voronoiMesh->getPolynomialOrder();k++) {
		for (l=0; l<= k; l++){
			valeur0 +=  coeff[indicee]*pow((x-(file->voronoiMesh->selection.xMax-file->voronoiMesh->selection.xMin)/2.0)*1e-6,double(k-l))*pow((y-(file->voronoiMesh->selection.yMax-file->voronoiMesh->selection.yMin)/2.0)*1e-6,l)*4.e-9*pow(10.0,(double)(6*(k-2)));
			indicee += 1;
		}
	}

	return valeur0/4.e-21;
}

// zonal optimization
double dfPosteriorVoronoi(double *FxFyD) {

	double D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	// current zone index
	const int a = file->voronoiMesh->getCurrentZone();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->getCellCount(a); i++) {
		if (file->voronoiMesh->getCell(a)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)] == file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
					 file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dx = file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dy = file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt ) - pow(fabs(dx-FxFyD[2]*FxFyD[0]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt) - pow(fabs(dy-FxFyD[2]*FxFyD[1]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->voronoiMesh->jeffreysPriorEnabled()) {
		result += 2.0*log(FxFyD[2]) - 2.0*log(FxFyD[2]*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
	}
	return -result;
}

double dfPosteriorSmoothingVoronoi(double *DFxFy) {

	double D_loc, D_bruit, Fx, Fy;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dfGradDxVoronoi(DFxFy,a));
			file->voronoiMesh->getCell(a)->setGradDy(dfGradDyVoronoi(DFxFy,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Fx = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()+1];
				Fy = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()+2];

//				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingVoronoiSelection(double *DFxFy) {

	double D_loc, D_bruit, Fx, Fy;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dfGradDxVoronoi(DFxFy,a));
			file->voronoiMesh->getCell(a)->setGradDy(dfGradDyVoronoi(DFxFy,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Fx  = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()+1];
				Fy  = DFxFy[3*file->voronoiMesh->getCell(idx)->getIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dfGradDxVoronoi(double *DFxFy, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DFxFy[3*file->voronoiMesh->getCell(i)->getIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += weight*DFxFy[3*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += weight*DFxFy[3*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dfGradDyVoronoi(double *DFxFy, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DFxFy[3*file->voronoiMesh->getCell(i)->getIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += weight*DFxFy[3*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += weight*DFxFy[3*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

// zonal optimization
double dPosteriorVoronoi(double *D) {

	double D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	// current zone index
	const int a = file->voronoiMesh->getCurrentZone();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->getCellCount(a); i++) {
		if (file->voronoiMesh->getCell(a)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)] == file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
					 file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dx = file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dy = file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

//				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt ) - pow(fabs(dx-FxFyD[2]*FxFyD[0]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt) - pow(fabs(dy-FxFyD[2]*FxFyD[1]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt);
				result += - log(4.0*PI*(D[0]+D_bruit)*dt) - (dx*dx)/(4.0*(D[0]+D_bruit)*dt) - (dy*dy)/(4.0*(D[0]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->voronoiMesh->jeffreysPriorEnabled()) {
		result += - (D[0]*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
	}
	return -result;
}

double dPosteriorSmoothingVoronoiSelection(double *D) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dGradDxVoronoi(D,a));
			file->voronoiMesh->getCell(a)->setGradDy(dGradDyVoronoi(D,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			const int idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				D_loc  = D[file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

//				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}
	return -result;
}

double dPosteriorSmoothingVoronoi(double *D) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dGradDxVoronoi(D,a));
			file->voronoiMesh->getCell(a)->setGradDy(dGradDyVoronoi(D,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = D[file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;

}

double dGradDxVoronoi(double *D, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = D[file->voronoiMesh->getCell(i)->getIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += weight*D[file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += weight*D[file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dGradDyVoronoi(double *D, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = D[file->voronoiMesh->getCell(i)->getIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += weight*D[file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += weight*D[file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dPosteriorSmoothingVoronoiRandomizedOptimization(double *D) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dGradDxVoronoiRandomizedOptimization(D,a));
			file->voronoiMesh->getCell(a)->setGradDy(dGradDyVoronoiRandomizedOptimization(D,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = D[file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += - (D_loc*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *D) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dGradDxVoronoiRandomizedOptimization(D,a));
			file->voronoiMesh->getCell(a)->setGradDy(dGradDyVoronoiRandomizedOptimization(D,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = D[file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += - (D_loc*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dGradDxVoronoiRandomizedOptimization(double *D, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = D[file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->roActive() ? weight*D[file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getDiffusion();
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += file->voronoiMesh->getCell(i)->getRightNeighbours(p)->roActive() ? weight*D[file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getDiffusion();
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dGradDyVoronoiRandomizedOptimization(double *D, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = D[file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += file->voronoiMesh->getCell(i)->getTopNeighbours(p)->roActive() ? weight*D[file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getDiffusion();
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->roActive() ? weight*D[file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getDiffusion();
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

// zonal optimization
double ddrPosteriorVoronoi(double *DrxDryD) {

	double D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	// current zone index
	const int a = file->voronoiMesh->getCurrentZone();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->getCellCount(a); i++) {
		if (file->voronoiMesh->getCell(a)->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)] == file->nPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
					 file->tPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dx = file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->xPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				dy = file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)+1] -
						   file->yPointer[file->voronoiMesh->getCell(a)->getIndex(i)];

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

//				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt ) - pow(fabs(dx-FxFyD[2]*FxFyD[0]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt) - pow(fabs(dy-FxFyD[2]*FxFyD[1]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt);
				result += - log(4.0*PI*(DrxDryD[2]+D_bruit)*dt) - (dx-DrxDryD[0]*dt)*(dx-DrxDryD[0]*dt)/(4*(DrxDryD[2]+D_bruit)*dt) - (dy-DrxDryD[1]*dt)*(dy-DrxDryD[1]*dt)/(4*(DrxDryD[2]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->voronoiMesh->jeffreysPriorEnabled()) {
		result += - 2.0*log(DrxDryD[2]*file->voronoiMesh->getDtMean() + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
	}
	return -result;
}

double ddrPosteriorSmoothingVoronoi(double *DDrxDry) {

	double D_loc, D_bruit, Drx, Dry;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(ddrGradDxVoronoi(DDrxDry,a));
			file->voronoiMesh->getCell(a)->setGradDy(ddrGradDyVoronoi(DDrxDry,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Drx = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()+1];
				Dry = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()+2];

//				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingVoronoiSelection(double *DDrxDry) {

	double D_loc, D_bruit, Drx, Dry;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(ddrGradDxVoronoi(DDrxDry,a));
			file->voronoiMesh->getCell(a)->setGradDy(ddrGradDyVoronoi(DDrxDry,a));
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Drx  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()+1];
				Dry  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double ddrGradDxVoronoi(double *DDrxDry, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DDrxDry[3*file->voronoiMesh->getCell(i)->getIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double ddrGradDyVoronoi(double *DDrxDry, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DDrxDry[3*file->voronoiMesh->getCell(i)->getIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}


double ddrPosteriorSmoothingVoronoiRandomizedOptimization(double *DDrxDry) {

	double D_loc, D_bruit;
	double result, dt, dx, dy, Drx, Dry;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(ddrGradDxVoronoiRandomizedOptimization(DDrxDry,a));
			file->voronoiMesh->getCell(a)->setGradDy(ddrGradDyVoronoiRandomizedOptimization(DDrxDry,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double lambda = (double)file->voronoiMeshGui->lambdaSlider->value();
	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Drx = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+1];
				Dry = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DDrxDry) {

	double D_loc, D_bruit;
	double result, dt, dx, dy, Drx, Dry;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(ddrGradDxVoronoiRandomizedOptimization(DDrxDry,a));
			file->voronoiMesh->getCell(a)->setGradDy(ddrGradDyVoronoiRandomizedOptimization(DDrxDry,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;
				Drx = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+1];
				Dry = DDrxDry[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrGradDxVoronoiRandomizedOptimization(double *DDrxDry, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DDrxDry[3*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->roActive() ? weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getDiffusion();
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += file->voronoiMesh->getCell(i)->getRightNeighbours(p)->roActive() ? weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getDiffusion();
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double ddrGradDyVoronoiRandomizedOptimization(double *DDrxDry, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DDrxDry[3*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += file->voronoiMesh->getCell(i)->getTopNeighbours(p)->roActive() ? weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getDiffusion();
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->roActive() ? weight*DDrxDry[3*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getDiffusion();
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dfFxValueVoronoi(double *V_ij, int i ) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double potential = V_ij[file->voronoiMesh->getCell(i)->getIdentifier()];
	double rCentre, lCentre;
	double rPotential, lPotential;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lPotential += weight*V_ij[file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rPotential += weight*V_ij[file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= (double)lNumber;
		lPotential /= (double)lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= (double)rNumber;
		rPotential /= (double)rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dfFyValueVoronoi(double *V_ij, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double potential = V_ij[file->voronoiMesh->getCell(i)->getIdentifier()];
	double tCentre, bCentre;
	double tPotential, bPotential;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bPotential = tPotential = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tPotential += weight*V_ij[file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bPotential += weight*V_ij[file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= (double)tNumber;
		tPotential /= (double)tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= (double)bNumber;
		bPotential /= (double)bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double squareDifferenceVoronoi(double *V_ij) {
	int i;
	double result=0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dfFxValueVoronoi(V_ij,a));
			file->voronoiMesh->getCell(a)->setGradVy(dfFyValueVoronoi(V_ij,a));
		}
	}

	#pragma omp for
	for (i = 0; i < file->voronoiMesh->getNumberOfClusters(); i++){
		if ( (file->voronoiMesh->active(i)) && (file->voronoiMesh->getNumberOfNeighbours(i)!=0) ) {
			const double fx = file->voronoiMesh->getCell(i)->getGradVx();
			const double fy = file->voronoiMesh->getCell(i)->getGradVy();
			result += pow( fx - file->voronoiMesh->getForceX(i) ,2.0);
			result += file->voronoiMesh->getBeta()*fx*fx;
			result += pow( fy - file->voronoiMesh->getForceY(i) ,2.0);
			result += file->voronoiMesh->getBeta()*fy*fy;
//			fprintf(stderr,"%f\t[%i]\t%f\t%f\t%f\t%f\n",result,i,fx,fy,file->voronoiMesh->getForceX(i),file->voronoiMesh->getForceY(i));
		}
	}

	return result;
}

double dfPosteriorSmoothingVoronoiRandomizedOptimization(double *DFxFy) {

	double D_loc, D_bruit;
	double result, dt, dx, dy, Fx, Fy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dfGradDxVoronoiRandomizedOptimization(DFxFy,a));
			file->voronoiMesh->getCell(a)->setGradDy(dfGradDyVoronoiRandomizedOptimization(DFxFy,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				Fx = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+1];
				Fy = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+2];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

//				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DFxFy) {

	double D_loc, D_bruit;
	double result, dt, dx, dy, Fx, Fy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradDx(dfGradDxVoronoiRandomizedOptimization(DFxFy,a));
			file->voronoiMesh->getCell(a)->setGradDy(dfGradDyVoronoiRandomizedOptimization(DFxFy,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				Fx = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+1];
				Fy = DFxFy[3*file->voronoiMesh->getCell(idx)->getRoIdentifier()+2];

				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

//				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfGradDxVoronoiRandomizedOptimization(double *DFxFy, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DFxFy[3*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->roActive() ? weight*DFxFy[3*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getDiffusion();
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += file->voronoiMesh->getCell(i)->getRightNeighbours(p)->roActive() ? weight*DFxFy[3*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getDiffusion();
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dfGradDyVoronoiRandomizedOptimization(double *DFxFy, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DFxFy[3*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += file->voronoiMesh->getCell(i)->getTopNeighbours(p)->roActive() ? weight*DFxFy[3*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getDiffusion();
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->roActive() ? weight*DFxFy[3*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getDiffusion();
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double quadraticTest(double *coeff) {

    double res = 0.0;
    res = 0.5*coeff[0]*coeff[0] + 0.5*(coeff[1]-0.5)*(coeff[1]-0.5);
    return res;
}

double dvPosteriorVoronoi(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			const int idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}

			}
		}
	}

	return -result;
}

double dvPosteriorVoronoiSelection(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			const int idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;

}

double dvPosteriorSmoothingVoronoi(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradDx(dvGradDxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradDy(dvGradDyVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double lambda = (double)file->voronoiMeshGui->lambdaSlider->value();
	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();
				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*file->voronoiMesh->getCell(idx)->getAreaX() +
									  gradVy*gradVy*file->voronoiMesh->getCell(idx)->getAreaY());
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingVoronoiSelection(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->active(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradDx(dvGradDxVoronoi(DV,a));
			file->voronoiMesh->getCell(a)->setGradDy(dvGradDyVoronoi(DV,a));
		}
	}

	const double lambda = (double)file->voronoiMeshGui->lambdaSlider->value();
	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->active(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();
				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*file->voronoiMesh->getCell(idx)->getAreaX() +
									  gradVy*gradVy*file->voronoiMesh->getCell(idx)->getAreaY());
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dvGradVxVoronoi(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double potential = DV[2*file->voronoiMesh->getCell(i)->getIdentifier()+1];
	double rCentre, lCentre;
	double rPotential, lPotential;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lPotential += weight*DV[2*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()+1];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rPotential += weight*DV[2*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()+1];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lPotential /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rPotential /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dvGradVyVoronoi(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double potential = DV[2*file->voronoiMesh->getCell(i)->getIdentifier()+1];
	double tCentre, bCentre;
	double tPotential, bPotential;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bPotential = tPotential = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tPotential += weight*DV[2*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()+1];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bPotential += weight*DV[2*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()+1];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tPotential /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bPotential /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dvGradDxVoronoi(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DV[2*file->voronoiMesh->getCell(i)->getIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += weight*DV[2*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getIdentifier()];
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += weight*DV[2*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getIdentifier()];
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dvGradDyVoronoi(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DV[2*file->voronoiMesh->getCell(i)->getIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += weight*DV[2*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getIdentifier()];
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += weight*DV[2*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getIdentifier()];
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dvPosteriorVoronoiRandomizedOptimization(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			const int idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - ((dx-D_loc*gradVx*dt)*(dx-D_loc*gradVx*dt) + (dy-D_loc*gradVy*dt)*(dy-D_loc*gradVy*dt))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}

			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingVoronoiRandomizedOptimization(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradDx(dvGradDxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradDy(dvGradDyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double lambda = (double)file->voronoiMeshGui->lambdaSlider->value();
	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();
				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*file->voronoiMesh->getCell(idx)->getAreaX() +
									  gradVy*gradVy*file->voronoiMesh->getCell(idx)->getAreaY());
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvPosteriorVoronoiRandomizedOptimizationSelection(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			const int idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingVoronoiRandomizedOptimizationSelection(double *DV) {

	double D_loc, D_bruit;
	double result, dt, dx, dy;

	int idx;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
		if (file->voronoiMesh->roActive(a)) {
			file->voronoiMesh->getCell(a)->setGradVx(dvGradVxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradVy(dvGradVyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradDx(dvGradDxVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->setGradDy(dvGradDyVoronoiRandomizedOptimization(DV,a));
			file->voronoiMesh->getCell(a)->resetPriors();
		}
	}

	const double lambda = (double)file->voronoiMeshGui->lambdaSlider->value();
	const double mu = (double)file->voronoiMeshGui->muSlider->value();

	#pragma omp for
	for (int i = 0; i < file->voronoiMesh->selection.count - 1; i++) {
		const int ii = file->voronoiMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			idx = file->voronoiMesh->getClusterIndex(i);

			if (file->voronoiMesh->roActive(idx)) {

				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = file->voronoiMesh->getCell(idx)->getGradVx();
				const double gradVy = file->voronoiMesh->getCell(idx)->getGradVy();
				const double gradDx = file->voronoiMesh->getCell(idx)->getGradDx();
				const double gradDy = file->voronoiMesh->getCell(idx)->getGradDy();

				D_loc  = DV[2*file->voronoiMesh->getCell(idx)->getRoIdentifier()];
				D_bruit = file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (file->voronoiMesh->getCell(idx)->priorsApplied() == false) {
					// Jeffreys
					if (file->voronoiMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->voronoiMesh->getSigma()*file->voronoiMesh->getSigma());
						file->voronoiMesh->getCell(idx)->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*file->voronoiMesh->getCell(idx)->getAreaX() +
									  gradVy*gradVy*file->voronoiMesh->getCell(idx)->getAreaY());
					result -= mu*(gradDx*gradDx*file->voronoiMesh->getCell(idx)->getAreaX() +
								  gradDy*gradDy*file->voronoiMesh->getCell(idx)->getAreaY());
					file->voronoiMesh->getCell(idx)->applySmoothingPrior();
				}

			}
		}
	}

	return -result;
}

double dvGradVxVoronoiRandomizedOptimization(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double potential = DV[2*file->voronoiMesh->getCell(i)->getRoIdentifier()+1];
	double rCentre, lCentre;
	double rPotential, lPotential;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lPotential += file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getRoIdentifier()+1] : weight*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getPotential();
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rPotential += file->voronoiMesh->getCell(i)->getRightNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getRoIdentifier()+1] : weight*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getPotential();
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lPotential /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rPotential /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dvGradVyVoronoiRandomizedOptimization(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double potential = DV[2*file->voronoiMesh->getCell(i)->getRoIdentifier()+1];
	double tCentre, bCentre;
	double tPotential, bPotential;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bPotential = tPotential = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tPotential += file->voronoiMesh->getCell(i)->getTopNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getRoIdentifier()+1] : weight*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getPotential();
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bPotential += file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getRoIdentifier()+1] : weight*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getPotential();
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tPotential /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bPotential /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dvGradDxVoronoiRandomizedOptimization(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getXCentre();
	const double diffusion = DV[2*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNLeftNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getLeftNeighbourWeight(p);
			lCentre += (file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getXCentre()-centre);
			lDiffusion += file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getLeftNeighbours(p)->getDiffusion();
			lNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNRightNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getRightNeighbourWeight(p);
			rCentre += (file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getXCentre()-centre);
			rDiffusion += file->voronoiMesh->getCell(i)->getRightNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getRightNeighbours(p)->getDiffusion();
			rNumber++;
		}
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dvGradDyVoronoiRandomizedOptimization(double *DV, int i) {
	const double centre = file->voronoiMesh->getCell(i)->getYCentre();
	const double diffusion = DV[2*file->voronoiMesh->getCell(i)->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNTopNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getTopNeighbourWeight(p);
			tCentre += (file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getYCentre()-centre);
			tDiffusion += file->voronoiMesh->getCell(i)->getTopNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getTopNeighbours(p)->getDiffusion();
			tNumber++;
		}
	}
	for (int p = 0; p < file->voronoiMesh->getCell(i)->getNBottomNeighbours(); p++) {
		if (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->active()) {
			const double weight = file->voronoiMesh->getCell(i)->getBottomNeighbourWeight(p);
			bCentre += (file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getYCentre()-centre);
			bDiffusion += file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->roActive() ? weight*DV[2*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getRoIdentifier()] : weight*file->voronoiMesh->getCell(i)->getBottomNeighbours(p)->getDiffusion();
			bNumber++;
		}
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

/********* QUAD TREE MESH **********/

double dPosteriorQuadTree(double *D) {

	double D_bruit;
	double result, dt, dx, dy, dtAverage;

//	int steps = 0;
	result = 0.0;
	dtAverage = 0.0;

	QuadTree* tree = file->treeMesh->getCurrentQuadZone();

	for (int i = 0; i < tree->getCount(); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (tree->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[tree->getIndex(i)+1] - file->tPointer[tree->getIndex(i)];

				dtAverage += dt;
				dx = file->xPointer[tree->getIndex(i)+1] - file->xPointer[tree->getIndex(i)];
				dy = file->yPointer[tree->getIndex(i)+1] - file->yPointer[tree->getIndex(i)];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D[0]+D_bruit)*dt) - (dx*dx)/(4.0*(D[0]+D_bruit)*dt) - (dy*dy)/(4.0*(D[0]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->treeMesh->jeffreysPriorEnabled()) {
		result += -log(D[0]*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
	}
	return -result;

}

void dPosteriorSmoothingSetupQuadTree(double *D, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(dGradDxQuadTree(D,tree));
			tree->setGradDy(dGradDyQuadTree(D,tree));
			tree->resetPriors();
		}
		return;
	}
	dPosteriorSmoothingSetupQuadTree(D,tree->nw);
	dPosteriorSmoothingSetupQuadTree(D,tree->ne);
	dPosteriorSmoothingSetupQuadTree(D,tree->sw);
	dPosteriorSmoothingSetupQuadTree(D,tree->se);
	return;
}

double dPosteriorSmoothingQuadTree(double *D) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dPosteriorSmoothingSetupQuadTree(D,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = D[currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

//				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dPosteriorSmoothingQuadTreeSelection(double *D) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dPosteriorSmoothingSetupQuadTree(D,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = D[currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

//				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dGradDxQuadTree(double *D, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = D[tree->identifier];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += D[tree->leftNeighbours[l]->identifier];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += D[tree->rightNeighbours[r]->identifier];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dGradDyQuadTree(double *D, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = D[tree->identifier];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += D[tree->topNeighbours[t]->identifier];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += D[tree->bottomNeighbours[b]->identifier];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}

void dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *D, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roActive()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(dGradDxQuadTreeRandomizedOptimization(D,tree));
			tree->setGradDy(dGradDyQuadTreeRandomizedOptimization(D,tree));
			tree->resetPriors();
		}
		return;
	}
	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,tree->nw);
	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,tree->ne);
	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,tree->sw);
	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,tree->se);
	return;
}

double dPosteriorSmoothingQuadTreeRandomizedOptimization(double *D) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = D[currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *D) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dPosteriorSmoothingSetupQuadTreeRandomizedOptimization(D,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = D[currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - (dx*dx)/(4.0*(D_loc+D_bruit)*dt) - (dy*dy)/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += -log(D_loc*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dGradDxQuadTreeRandomizedOptimization(double *D, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = D[tree->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += tree->leftNeighbours[l]->roActive() ? D[tree->leftNeighbours[l]->getRoIdentifier()] : tree->leftNeighbours[l]->getDiffusion();
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += tree->rightNeighbours[r]->roActive() ? D[tree->rightNeighbours[r]->getRoIdentifier()] : tree->rightNeighbours[r]->getDiffusion();
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dGradDyQuadTreeRandomizedOptimization(double *D, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = D[tree->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += tree->topNeighbours[t]->roActive() ? D[tree->topNeighbours[t]->getRoIdentifier()] : tree->topNeighbours[t]->getDiffusion();
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += tree->bottomNeighbours[b]->roActive() ? D[tree->bottomNeighbours[b]->getRoIdentifier()] : tree->bottomNeighbours[b]->getDiffusion();
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double ddrPosteriorQuadTree(double *DrxDryD) {

	double D_bruit;
	double result, dt, dx, dy, dtAverage;

//	int steps = 0;
	result = 0.0;
	dtAverage = 0.0;

	QuadTree* tree = file->treeMesh->getCurrentQuadZone();

	for (int i = 0; i < tree->getCount(); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (tree->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[tree->getIndex(i)+1] - file->tPointer[tree->getIndex(i)];

				dtAverage += dt;
				dx = file->xPointer[tree->getIndex(i)+1] - file->xPointer[tree->getIndex(i)];
				dy = file->yPointer[tree->getIndex(i)+1] - file->yPointer[tree->getIndex(i)];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(DrxDryD[2]+D_bruit)*dt) - (dx-DrxDryD[0]*dt)*(dx-DrxDryD[0]*dt)/(4*(DrxDryD[2]+D_bruit)*dt) - (dy-DrxDryD[1]*dt)*(dy-DrxDryD[1]*dt)/(4*(DrxDryD[2]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->treeMesh->jeffreysPriorEnabled()) {
		result += - 2.0*log(DrxDryD[2]*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
	}
	return -result;

}

double dfPosteriorQuadTree(double *FxFyD) {

	double D_bruit;
	double result, dt, dx, dy, dtAverage;

//	int steps = 0;
	result = 0.0;
	dtAverage = 0.0;

	QuadTree* tree = file->treeMesh->getCurrentQuadZone();

	for (int i = 0; i < tree->getCount(); i++) {

		// check that index is not at end of list
		// check that deltas between data points in same file
		if (tree->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[tree->getIndex(i)] == file->nPointer[tree->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[tree->getIndex(i)+1] - file->tPointer[tree->getIndex(i)];

				dtAverage += dt;

				dx = file->xPointer[tree->getIndex(i)+1] - file->xPointer[tree->getIndex(i)];

				dy = file->yPointer[tree->getIndex(i)+1] - file->yPointer[tree->getIndex(i)];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt ) - pow(fabs(dx-FxFyD[2]*FxFyD[0]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt) - pow(fabs(dy-FxFyD[2]*FxFyD[1]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt);
			}
		}
	}

	// Jeffrey's Prior
	if (file->treeMesh->jeffreysPriorEnabled()) {
		result += 2.0*log(FxFyD[2]) - 2.0*log(FxFyD[2]*file->treeMesh->getDtMean() + file->treeMesh->getSigma()*file->treeMesh->getSigma());
	}
	return -result;

}

void dfPosteriorSmoothingSetupQuadTree(double *DFxFy, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(dfGradDxQuadTree(DFxFy,tree));
			tree->setGradDy(dfGradDyQuadTree(DFxFy,tree));
			tree->resetPriors();
		}
		return;
	}
	dfPosteriorSmoothingSetupQuadTree(DFxFy,tree->nw);
	dfPosteriorSmoothingSetupQuadTree(DFxFy,tree->ne);
	dfPosteriorSmoothingSetupQuadTree(DFxFy,tree->sw);
	dfPosteriorSmoothingSetupQuadTree(DFxFy,tree->se);
	return;
}

void ddrPosteriorSmoothingSetupQuadTree(double *DDrxDry, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(ddrGradDxQuadTree(DDrxDry,tree));
			tree->setGradDy(ddrGradDyQuadTree(DDrxDry,tree));
			tree->resetPriors();
		}
		return;
	}
	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,tree->nw);
	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,tree->ne);
	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,tree->sw);
	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,tree->se);
	return;
}

double dfPosteriorSmoothingQuadTree(double *DFxFy) {
	int i;
	double D_loc, D_bruit, Fx, Fy;
	double result, dt, dx, dy;

	result = 0.0;

	dfPosteriorSmoothingSetupQuadTree(DFxFy,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DFxFy[3*currentLeaf->identifier];
				Fx = DFxFy[3*currentLeaf->identifier+1];
				Fy = DFxFy[3*currentLeaf->identifier+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingQuadTree(double *DDrxDry) {
	int i;
	double D_loc, D_bruit, Drx, Dry;
	double result, dt, dx, dy;

	result = 0.0;

	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DDrxDry[3*currentLeaf->identifier];
				Drx = DDrxDry[3*currentLeaf->identifier+1];
				Dry = DDrxDry[3*currentLeaf->identifier+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingQuadTreeSelection(double *DFxFy) {
	int i;
	double D_loc, D_bruit, Fx, Fy;
	double result, dt, dx, dy;

	result = 0.0;

	dfPosteriorSmoothingSetupQuadTree(DFxFy,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc  = DFxFy[3*currentLeaf->identifier];
				Fx  = DFxFy[3*currentLeaf->identifier+1];
				Fy  = DFxFy[3*currentLeaf->identifier+2];

//				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);
				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingQuadTreeSelection(double *DDrxDry) {
	int i;
	double D_loc, D_bruit, Drx, Dry;
	double result, dt, dx, dy;

	result = 0.0;

	ddrPosteriorSmoothingSetupQuadTree(DDrxDry,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DDrxDry[3*currentLeaf->identifier];
				Drx = DDrxDry[3*currentLeaf->identifier+1];
				Dry = DDrxDry[3*currentLeaf->identifier+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfGradDxQuadTree(double *DFxFy, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DFxFy[3*tree->identifier];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += DFxFy[3*tree->leftNeighbours[l]->identifier];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += DFxFy[3*tree->rightNeighbours[r]->identifier];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double ddrGradDxQuadTree(double *DDrxDry, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DDrxDry[3*tree->identifier];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += DDrxDry[3*tree->leftNeighbours[l]->identifier];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += DDrxDry[3*tree->rightNeighbours[r]->identifier];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dfGradDyQuadTree(double *DFxFy, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DFxFy[3*tree->identifier];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += DFxFy[3*tree->topNeighbours[t]->identifier];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += DFxFy[3*tree->bottomNeighbours[b]->identifier];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}

double ddrGradDyQuadTree(double *DDrxDry, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DDrxDry[3*tree->identifier];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += DDrxDry[3*tree->topNeighbours[t]->identifier];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += DDrxDry[3*tree->bottomNeighbours[b]->identifier];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}


void ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roActive()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(ddrGradDxQuadTreeRandomizedOptimization(DDrxDry,tree));
			tree->setGradDy(ddrGradDyQuadTreeRandomizedOptimization(DDrxDry,tree));
			tree->resetPriors();
		}
		return;
	}
	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,tree->nw);
	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,tree->ne);
	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,tree->sw);
	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,tree->se);
	return;
}

double ddrPosteriorSmoothingQuadTreeRandomizedOptimization(double *DDrxDry) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy, Drx, Dry;

	result = 0.0;

	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DDrxDry[3*currentLeaf->getRoIdentifier()];
				Drx = DDrxDry[3*currentLeaf->getRoIdentifier()+1];
				Dry = DDrxDry[3*currentLeaf->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DDrxDry) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy, Drx, Dry;

	result = 0.0;

	ddrPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DDrxDry,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DDrxDry[3*currentLeaf->getRoIdentifier()];
				Drx = DDrxDry[3*currentLeaf->getRoIdentifier()+1];
				Dry = DDrxDry[3*currentLeaf->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-Drx*dt)*(dx-Drx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-Dry*dt)*(dy-Dry*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double ddrGradDxQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DDrxDry[3*tree->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += tree->leftNeighbours[l]->roActive() ? DDrxDry[3*tree->leftNeighbours[l]->getRoIdentifier()] : tree->leftNeighbours[l]->getDiffusion();
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += tree->rightNeighbours[r]->roActive() ? DDrxDry[3*tree->rightNeighbours[r]->getRoIdentifier()] : tree->rightNeighbours[r]->getDiffusion();
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double ddrGradDyQuadTreeRandomizedOptimization(double *DDrxDry, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DDrxDry[3*tree->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += tree->topNeighbours[t]->roActive() ? DDrxDry[3*tree->topNeighbours[t]->getRoIdentifier()] : tree->topNeighbours[t]->getDiffusion();
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += tree->bottomNeighbours[b]->roActive() ? DDrxDry[3*tree->bottomNeighbours[b]->getRoIdentifier()] : tree->bottomNeighbours[b]->getDiffusion();
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double squareDifferenceQuadTree(double *V_ij) {
	double squareDifference = 0.0;

	squareDifferenceSetupQuadTree(V_ij,file->treeMesh->quadTree);
	squareDifferenceQuadTree(V_ij,file->treeMesh->quadTree,&squareDifference);

	return squareDifference;
}

void squareDifferenceSetupQuadTree(double *V_ij,QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradVx(dfFxValueQuadTree(V_ij,tree));
			tree->setGradVy(dfFyValueQuadTree(V_ij,tree));
		}
		return;
	}

	squareDifferenceSetupQuadTree(V_ij,tree->nw);
	squareDifferenceSetupQuadTree(V_ij,tree->ne);
	squareDifferenceSetupQuadTree(V_ij,tree->sw);
	squareDifferenceSetupQuadTree(V_ij,tree->se);

	return;
}

void squareDifferenceQuadTree(double *V_ij, QuadTree *tree, double *squareDifference) {

	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			const double fx = tree->getGradVx();
			const double fy = tree->getGradVy();

			*squareDifference += pow( fabs(fx-tree->getFx()),2.0);
			*squareDifference += file->treeMesh->getBeta()*pow(fx,2.0);
			*squareDifference += pow( fabs(fy-tree->getFy()),2.0);
			*squareDifference += file->treeMesh->getBeta()*pow(fy,2.0);
		}
		return;
	}

	squareDifferenceQuadTree(V_ij,tree->nw,squareDifference);
	squareDifferenceQuadTree(V_ij,tree->ne,squareDifference);
	squareDifferenceQuadTree(V_ij,tree->sw,squareDifference);
	squareDifferenceQuadTree(V_ij,tree->se,squareDifference);

	return;
}

void dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roActive()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradDx(dfGradDxQuadTreeRandomizedOptimization(DFxFy,tree));
			tree->setGradDy(dfGradDyQuadTreeRandomizedOptimization(DFxFy,tree));
			tree->resetPriors();
		}
		return;
	}
	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,tree->nw);
	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,tree->ne);
	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,tree->sw);
	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,tree->se);
	return;
}

double dfPosteriorSmoothingQuadTreeRandomizedOptimization(double *DFxFy) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy, Fx, Fy;

	result = 0.0;

	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,file->treeMesh->quadTree);

	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DFxFy[3*currentLeaf->getRoIdentifier()];
				Fx = DFxFy[3*currentLeaf->getRoIdentifier()+1];
				Fy = DFxFy[3*currentLeaf->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DFxFy) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy, Fx, Fy;

	result = 0.0;

	dfPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DFxFy,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				D_loc = DFxFy[3*currentLeaf->getRoIdentifier()];
				Fx = DFxFy[3*currentLeaf->getRoIdentifier()+1];
				Fy = DFxFy[3*currentLeaf->getRoIdentifier()+2];

				result += - log(4.0*PI*(D_loc+D_bruit)*dt) - ((dx-D_loc*Fx*dt)*(dx-D_loc*Fx*dt))/(4.0*(D_loc+D_bruit)*dt) - ((dy-D_loc*Fy*dt)*(dy-D_loc*Fy*dt))/(4.0*(D_loc+D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= mu*(gradDx*gradDx*currentLeaf->xArea + gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dfGradDxQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DFxFy[3*tree->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += tree->leftNeighbours[l]->roActive() ? DFxFy[3*tree->leftNeighbours[l]->getRoIdentifier()] : tree->leftNeighbours[l]->getDiffusion();
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += tree->rightNeighbours[r]->roActive() ? DFxFy[3*tree->rightNeighbours[r]->getRoIdentifier()] : tree->rightNeighbours[r]->getDiffusion();
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dfGradDyQuadTreeRandomizedOptimization(double *DFxFy, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DFxFy[3*tree->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += tree->topNeighbours[t]->roActive() ? DFxFy[3*tree->topNeighbours[t]->getRoIdentifier()] : tree->topNeighbours[t]->getDiffusion();
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += tree->bottomNeighbours[b]->roActive() ? DFxFy[3*tree->bottomNeighbours[b]->getRoIdentifier()] : tree->bottomNeighbours[b]->getDiffusion();
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dfFxValueQuadTree(double *V_ij, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double potential = V_ij[tree->identifier];
	double rCentre, lCentre;
	double rPotential, lPotential;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lPotential += V_ij[tree->leftNeighbours[l]->identifier];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rPotential += V_ij[tree->rightNeighbours[r]->identifier];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lPotential /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rPotential /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dfFyValueQuadTree(double *V_ij, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double potential = V_ij[tree->identifier];
	double tCentre, bCentre;
	double tPotential, bPotential;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tPotential = bPotential = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tPotential += V_ij[tree->topNeighbours[t]->identifier];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bPotential += V_ij[tree->bottomNeighbours[b]->identifier];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tPotential /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bPotential /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}

double polynomialPosteriorQuadTree(double *coeff_f) {

	int i;
	double D_loc, D_bruit, fx, fy;
	double result, dt, deltaa_x, deltaa_y;

	QuadTree* currentLeaf;
	result = 0.0;

	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				deltaa_x = file->xPointer[i+1] - file->xPointer[i];
				deltaa_y = file->yPointer[i+1] - file->yPointer[i];

				fx = polynomialFxValueQuadTree(coeff_f,currentLeaf->getXCentre(),currentLeaf->getYCentre());
				fy = polynomialFyValueQuadTree(coeff_f,currentLeaf->getXCentre(),currentLeaf->getYCentre());

				D_loc  = polynomialDValueQuadTree(coeff_f,currentLeaf);

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}
	return -result;

}

double polynomialPosteriorQuadTreeSelection(double *coeff_f) {

	int i;
	double D_loc, D_bruit, fx, fy;
	double result, dt, deltaa_x, deltaa_y;

	QuadTree* currentLeaf;
	result = 0.0;

	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {
			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				deltaa_x = file->xPointer[ii+1] - file->xPointer[ii];
				deltaa_y = file->yPointer[ii+1] - file->yPointer[ii];

				fx = polynomialFxValueQuadTree(coeff_f,currentLeaf->getXCentre(),currentLeaf->getYCentre());
				fy = polynomialFyValueQuadTree(coeff_f,currentLeaf->getXCentre(),currentLeaf->getYCentre());

				D_loc  = polynomialDValueQuadTree(coeff_f,currentLeaf);

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;
				result += -log(4.*PI*(D_loc + D_bruit)*dt) - pow((deltaa_x-D_loc*fx*dt), 2.0)/(4.*(D_loc + D_bruit)*dt) - pow((deltaa_y-D_loc*fy*dt), 2.0)/(4.*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}
	return -result;

}

double polynomialFxValueQuadTree(double *coeffff, double x, double y) {
	  int k,l,indiceb;
	  double valeur1;

	  indiceb = 0;
	  valeur1=0.;
	  for ( k=1; k <= file->treeMesh->getPolynomialOrder(); k++){
	    for (l=0; l<= k; l++){
	      valeur1 +=  -coeffff[indiceb]*(k-l)*pow(x-file->treeMesh->getXCentre(),(double) FMAX((double) (k-l-1),0.0))*pow(y-file->treeMesh->getYCentre(),(double) FMAX((double) l,0.0));
	      indiceb += 1;
	    }
	  }

	  return valeur1;

}

double polynomialFyValueQuadTree(double *coefffff, double x, double y) {
	  int k,l,indicec;
	  double valeur2;

	  indicec = 0;
	  valeur2 =0.;
	  for (k = 1; k <= file->treeMesh->getPolynomialOrder(); k++){
	    for (l=0; l <= k; l++){
	      valeur2 +=  -coefffff[indicec]*(l)*pow(x-file->treeMesh->getXCentre(),(double) FMAX((double)(k-l),0.0))*pow(y-file->treeMesh->getXCentre(),(double) FMAX((double)(l-1),0.0));
	      indicec += 1;
	    }
	  }

	  return valeur2;
}

double polynomialDValueQuadTree(double *coeffff, QuadTree *tree) {
	  double valeur3;

	  if (tree->active()) {
		  valeur3 = coeffff[file->treeMesh->getCoefficients() + tree->identifier];
	  }
	  else {
		  valeur3 = 0.0;
	  }

	  return valeur3;

}

double polynomialEpValueQuadTree(double *coeff, double x, double y ) {
	int k,l,indicee;
	double valeur0;

	indicee = 0;
	valeur0 =0.;
	for (k=1; k<= file->treeMesh->getPolynomialOrder();k++) {
		for (l=0; l<= k; l++){
			valeur0 +=  coeff[indicee]*pow((x-file->getXCentre())*1e-6,double(k-l))*pow((y-file->getYCentre())*1e-6,l)*4.e-9*pow(10.0,(double)(6*(k-2)));
			indicee += 1;
		}
	}

	return valeur0/4.e-21;

}

void dvPosteriorSetupQuadTree(double *DV, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradVx(dvGradVxQuadTree(DV,tree));
			tree->setGradVy(dvGradVyQuadTree(DV,tree));
			tree->resetPriors();
		}
		return;
	}

	dvPosteriorSetupQuadTree(DV,tree->nw);
	dvPosteriorSetupQuadTree(DV,tree->ne);
	dvPosteriorSetupQuadTree(DV,tree->sw);
	dvPosteriorSetupQuadTree(DV,tree->se);
	return;
}

void dvPosteriorSmoothingSetupQuadTree(double *DV, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->active()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradVx(dvGradVxQuadTree(DV,tree));
			tree->setGradVy(dvGradVyQuadTree(DV,tree));
			tree->setGradDx(dvGradDxQuadTree(DV,tree));
			tree->setGradDy(dvGradDyQuadTree(DV,tree));
			tree->resetPriors();
		}
		return;
	}
	dvPosteriorSmoothingSetupQuadTree(DV,tree->nw);
	dvPosteriorSmoothingSetupQuadTree(DV,tree->ne);
	dvPosteriorSmoothingSetupQuadTree(DV,tree->sw);
	dvPosteriorSmoothingSetupQuadTree(DV,tree->se);
	return;
}

double dvPosteriorQuadTree(double *DV) {
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSetupQuadTree(DV,file->treeMesh->quadTree);
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();

				D_loc  = DV[2*currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorQuadTreeSelection(double *DV) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSetupQuadTree(DV,file->treeMesh->quadTree);
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();

				D_loc  = DV[2*currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingQuadTree(double *DV) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSmoothingSetupQuadTree(DV,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();
				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = DV[2*currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*currentLeaf->xArea +
									  gradVy*gradVy*currentLeaf->yArea);
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingQuadTreeSelection(double *DV) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSmoothingSetupQuadTree(DV,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->active()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();
				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = DV[2*currentLeaf->identifier];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*currentLeaf->xArea +
									  gradVy*gradVy*currentLeaf->yArea);
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvGradVxQuadTree(double *DV, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double potential = DV[2*tree->identifier+1];
	double rCentre, lCentre;
	double rPotential, lPotential;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lPotential += DV[2*tree->leftNeighbours[l]->identifier+1];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rPotential += DV[2*tree->rightNeighbours[r]->identifier+1];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lPotential /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rPotential /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dvGradVyQuadTree(double *DV, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double potential = DV[2*tree->identifier+1];
	double tCentre, bCentre;
	double tPotential, bPotential;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tPotential = bPotential = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tPotential += DV[2*tree->topNeighbours[t]->identifier+1];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bPotential += DV[2*tree->bottomNeighbours[b]->identifier+1];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tPotential /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bPotential /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}

double dvGradDxQuadTree(double *DV, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DV[2*tree->identifier];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += DV[2*tree->leftNeighbours[l]->identifier];
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += DV[2*tree->rightNeighbours[r]->identifier];
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dvGradDyQuadTree(double *DV, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DV[2*tree->identifier];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += DV[2*tree->topNeighbours[t]->identifier];
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += DV[2*tree->bottomNeighbours[b]->identifier];
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
//		printf("zero y force : %i\n",tree->count);
	}

	return yForce;
}

void dvPosteriorSetupQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roActive()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradVx(dvGradVxQuadTreeRandomizedOptimization(DV,tree));
			tree->setGradVy(dvGradVyQuadTreeRandomizedOptimization(DV,tree));
			tree->resetPriors();
		}
		return;
	}
	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,tree->nw);
	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,tree->ne);
	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,tree->sw);
	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,tree->se);
	return;
}

void dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
		if (tree->roActive()/* && (tree->nLeftNeighbours + tree->nRightNeighbours + tree->nTopNeighbours + tree->nBottomNeighbours > 0)*/) {
			tree->setGradVx(dvGradVxQuadTreeRandomizedOptimization(DV,tree));
			tree->setGradVy(dvGradVyQuadTreeRandomizedOptimization(DV,tree));
			tree->setGradDx(dvGradDxQuadTreeRandomizedOptimization(DV,tree));
			tree->setGradDy(dvGradDyQuadTreeRandomizedOptimization(DV,tree));
			tree->resetPriors();
		}
		return;
	}
	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,tree->nw);
	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,tree->ne);
	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,tree->sw);
	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,tree->se);
	return;
}

double dvPosteriorQuadTreeRandomizedOptimization(double *DV) {
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,file->treeMesh->quadTree);
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (int i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {

				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();

				D_loc  = DV[2*currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingQuadTreeRandomizedOptimization(double *DV) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->localizationCount - 1; i++) {
		if (file->nPointer[i] == file->nPointer[i+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[i+1] - file->tPointer[i];
				dx = file->xPointer[i+1] - file->xPointer[i];
				dy = file->yPointer[i+1] - file->yPointer[i];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();
				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = DV[2*currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*currentLeaf->xArea +
									  gradVy*gradVy*currentLeaf->yArea);
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvPosteriorQuadTreeRandomizedOptimizationSelection(double *DV) {
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSetupQuadTreeRandomizedOptimization(DV,file->treeMesh->quadTree);
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (int i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();

				D_loc  = DV[2*currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += -log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
				}
			}
		}
	}

	return -result;
}

double dvPosteriorSmoothingQuadTreeRandomizedOptimizationSelection(double *DV) {
	int i;
	double D_loc, D_bruit;
	double result, dt, dx, dy;

	result = 0.0;

	dvPosteriorSmoothingSetupQuadTreeRandomizedOptimization(DV,file->treeMesh->quadTree);

	const double lambda = (double)file->treeMeshGui->lambdaSlider->value();
	const double mu = (double)file->treeMeshGui->muSlider->value();
	QuadTree* currentLeaf = NULL;

	#pragma omp for
	for (i = 0; i < file->treeMesh->selection.count - 1; i++) {
		const int ii = file->treeMesh->selection.indexArray[i];
		if (file->nPointer[ii] == file->nPointer[ii+1]) {

			currentLeaf = file->treeMesh->quadTreeLocalizationPointer[i];

			if (currentLeaf->roActive()) {
				dt = file->tPointer[ii+1] - file->tPointer[ii];
				dx = file->xPointer[ii+1] - file->xPointer[ii];
				dy = file->yPointer[ii+1] - file->yPointer[ii];

				const double gradVx = currentLeaf->getGradVx();
				const double gradVy = currentLeaf->getGradVy();
				const double gradDx = currentLeaf->getGradDx();
				const double gradDy = currentLeaf->getGradDy();

				D_loc  = DV[2*currentLeaf->getRoIdentifier()];

				D_bruit = file->treeMesh->getSigma()*file->treeMesh->getSigma()/dt;

				result += - log(4.0*PI*(D_loc + D_bruit)*dt) - (pow((dx-D_loc*gradVx*dt), 2.0) + pow((dy-D_loc*gradVy*dt), 2.0))/(4.0*(D_loc + D_bruit)*dt);

				// Priors
				if (currentLeaf->priorsApplied() == false) {
					// Jeffreys
					if (file->treeMesh->jeffreysPriorEnabled()) {
						result += 2.0*log(D_loc) - 2.0*log(D_loc*dt + file->treeMesh->getSigma()*file->treeMesh->getSigma());
						currentLeaf->applyJeffreysPrior();
					}
					// Smoothing
					result -= lambda*(gradVx*gradVx*currentLeaf->xArea +
									  gradVy*gradVy*currentLeaf->yArea);
					result -= mu*(gradDx*gradDx*currentLeaf->xArea +
								  gradDy*gradDy*currentLeaf->yArea);
					currentLeaf->applySmoothingPrior();
				}
			}
		}
	}

	return -result;
}

double dvGradVxQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double potential = DV[2*tree->getRoIdentifier()+1];
	double rCentre, lCentre;
	double rPotential, lPotential;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lPotential += tree->leftNeighbours[l]->roActive() ? DV[2*tree->leftNeighbours[l]->getRoIdentifier()+1] : tree->leftNeighbours[l]->getPotential();
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rPotential += tree->rightNeighbours[r]->roActive() ? DV[2*tree->rightNeighbours[r]->getRoIdentifier()+1] : tree->rightNeighbours[r]->getPotential();
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lPotential /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rPotential /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}

double dvGradVyQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double potential = DV[2*tree->getRoIdentifier()+1];
	double tCentre, bCentre;
	double tPotential, bPotential;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tPotential = bPotential = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tPotential += tree->topNeighbours[t]->roActive() ? DV[2*tree->topNeighbours[t]->getRoIdentifier()+1] : tree->topNeighbours[t]->getPotential();
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bPotential += tree->bottomNeighbours[b]->roActive() ? DV[2*tree->bottomNeighbours[b]->getRoIdentifier()+1] : tree->bottomNeighbours[b]->getPotential();
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tPotential /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bPotential /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double dvGradDxQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	const double centre = tree->getXCentre();
	const double diffusion = DV[2*tree->getRoIdentifier()];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;

	// left xCentres
	for (int l = 0; l < tree->nLeftNeighbours; l++) {
		lCentre += tree->leftNeighbours[l]->getXCentre();
		lDiffusion += tree->leftNeighbours[l]->roActive() ? DV[2*tree->leftNeighbours[l]->getRoIdentifier()] : tree->leftNeighbours[l]->getDiffusion();
	}
	// right xCentres
	for (int r = 0; r < tree->nRightNeighbours; r++) {
		rCentre += tree->rightNeighbours[r]->getXCentre();
		rDiffusion += tree->rightNeighbours[r]->roActive() ? DV[2*tree->rightNeighbours[r]->getRoIdentifier()] : tree->rightNeighbours[r]->getDiffusion();
	}

	// average left xCentres
	if (tree->nLeftNeighbours > 0) {
		lCentre /= (double)tree->nLeftNeighbours;
		lDiffusion /= (double)tree->nLeftNeighbours;
	}

	// average right xCentres
	if (tree->nRightNeighbours > 0) {
		rCentre /= (double)tree->nRightNeighbours;
		rDiffusion /= (double)tree->nRightNeighbours;
	}

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (tree->nLeftNeighbours != 0 && tree->nRightNeighbours != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (tree->nLeftNeighbours != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nRightNeighbours != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
//		printf("zero x force : %i\n",tree->count);
	}

	return xForce;
}

double dvGradDyQuadTreeRandomizedOptimization(double *DV, QuadTree *tree) {
	const double centre = tree->getYCentre();
	const double diffusion = DV[2*tree->getRoIdentifier()];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	tCentre = bCentre = 0.0;
	tDiffusion = bDiffusion = 0.0;

	// top yCentres
	for (int t = 0; t < tree->nTopNeighbours; t++) {
		tCentre += tree->topNeighbours[t]->getYCentre();
		tDiffusion += tree->topNeighbours[t]->roActive() ? DV[2*tree->topNeighbours[t]->getRoIdentifier()] : tree->topNeighbours[t]->getDiffusion();
	}
	// bottom yCentres
	for (int b = 0; b < tree->nBottomNeighbours; b++) {
		bCentre += tree->bottomNeighbours[b]->getYCentre();
		bDiffusion += tree->bottomNeighbours[b]->roActive() ? DV[2*tree->bottomNeighbours[b]->getRoIdentifier()] : tree->bottomNeighbours[b]->getDiffusion();
	}

	// average top yCentres
	if (tree->nTopNeighbours > 0) {
		tCentre /= (double)tree->nTopNeighbours;
		tDiffusion /= (double)tree->nTopNeighbours;
	}

	// average bottom yCentres
	if (tree->nBottomNeighbours > 0) {
		bCentre /= (double)tree->nBottomNeighbours;
		bDiffusion /= (double)tree->nBottomNeighbours;
	}

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tree->nTopNeighbours != 0 && tree->nBottomNeighbours != 0) {
		vander(y,w,q,2);
		// return derivative of polynomial
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on left side
	else if (tree->nBottomNeighbours != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on right side
	else if (tree->nTopNeighbours != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}

double singleTrajectoryPosterior(double *coeff) {

	double result = 0.0;
	double dt;

	for (int i = iMAP->singleTrajectoryInferenceGui->startIndex; i < iMAP->singleTrajectoryInferenceGui->endIndex-1; i++) {

		dt = file->tPointer[i+1] - file->tPointer[i];
		const double dx = file->xPointer[i+1] - file->xPointer[i];
		const double dy = file->yPointer[i+1] - file->yPointer[i];

		const double fx = coeff[1];
		const double fy = coeff[2];

		const double D_bruit = (double)(iMAP->singleTrajectoryInferenceGui->noiseSigmaSlider->value()/1000.0)*(iMAP->singleTrajectoryInferenceGui->noiseSigmaSlider->value()/1000.0)/dt;

		result += -log(4.0*PI*(coeff[0] + D_bruit)*dt) - pow((dx-coeff[0]*fx*dt), 2.0)/(4.0*(coeff[0] + D_bruit)*dt) - pow((dy-coeff[0]*fy*dt), 2.0)/(4.0*(coeff[0] + D_bruit)*dt);

	}
	result += 2.0*log(coeff[0]) - 2.0*log(coeff[0]*dt + (double)(iMAP->singleTrajectoryInferenceGui->noiseSigmaSlider->value()*iMAP->singleTrajectoryInferenceGui->noiseSigmaSlider->value()/1000.0/1000.0));
	
	return -result;
}

double customSelectionDFPosterior(double *DFxFy) {

	double res, dt;

	res = 0.0;

	for (int i = 0; i < file->selection->cell.count; i++) {
		if (file->selection->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->selection->getIndex(i)] == file->nPointer[file->selection->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[file->selection->getIndex(i)+1] - file->tPointer[file->selection->getIndex(i)];

				const double dx = file->xPointer[file->selection->getIndex(i)+1] - file->xPointer[file->selection->getIndex(i)];
				const double dy = file->yPointer[file->selection->getIndex(i)+1] - file->yPointer[file->selection->getIndex(i)];
				const double D_bruit = iMAP->customSelectionInferenceGui->getSigma()*iMAP->customSelectionInferenceGui->getSigma()/dt;

				res += -log(4.0*PI*(DFxFy[0] + D_bruit)*dt) - pow( fabs(dx-DFxFy[0]*DFxFy[1]*dt), 2.0)/(4.0*(DFxFy[0] + D_bruit)*dt) - pow( fabs(dy-DFxFy[0]*DFxFy[2]*dt), 2.0)/(4.0*(DFxFy[0] + D_bruit)*dt);
			}
		}
	}
	res += 2.0*log( DFxFy[0] ) - 2.0*log( DFxFy[0]*dt + iMAP->customSelectionInferenceGui->getSigma()*iMAP->customSelectionInferenceGui->getSigma() );

	return -res;
}

double customSelectionDDrPosterior(double *DFxFy) {

	double result, dt;

	result = 0.0;

	for (int i = 0; i < file->selection->cell.count; i++) {
		if (file->selection->getIndex(i) < file->localizationCount-1) {
			if (file->nPointer[file->selection->getIndex(i)] == file->nPointer[file->selection->getIndex(i)+1]) {
				// "next" index does not necessarily lie in the "current zone"
				dt = file->tPointer[file->selection->getIndex(i)+1] - file->tPointer[file->selection->getIndex(i)];

				const double dx = file->xPointer[file->selection->getIndex(i)+1] - file->xPointer[file->selection->getIndex(i)];
				const double dy = file->yPointer[file->selection->getIndex(i)+1] - file->yPointer[file->selection->getIndex(i)];
				const double D_bruit = iMAP->customSelectionInferenceGui->getSigma()*iMAP->customSelectionInferenceGui->getSigma()/dt;

				result += - log(4.0*PI*(DFxFy[0]+D_bruit)*dt) - (dx-DFxFy[1]*dt)*(dx-DFxFy[1]*dt)/(4*(DFxFy[0]+D_bruit)*dt) - (dy-DFxFy[2]*dt)*(dy-DFxFy[2]*dt)/(4*(DFxFy[0]+D_bruit)*dt);
//				res += -log(4.0*PI*(DFxFy[0] + D_bruit)*dt) - pow( fabs(dx-DFxFy[0]*DFxFy[1]*dt), 2.0)/(4.0*(DFxFy[0] + D_bruit)*dt) - pow( fabs(dy-DFxFy[0]*DFxFy[2]*dt), 2.0)/(4.0*(DFxFy[0] + D_bruit)*dt);
			}
		}
	}
	result += - 2.0*log( DFxFy[0]*dt + iMAP->customSelectionInferenceGui->getSigma()*iMAP->customSelectionInferenceGui->getSigma() );

	return -result;
}

void textDisplayUpdate(char *update) {
	Fl::lock();

	sprintf(iMAP->outputStream,"%i   %s",iMAP->textBuffer->count_lines(0,iMAP->textBuffer->length()),update);
	iMAP->textBuffer->append(iMAP->outputStream);
	iMAP->textDisplay->insert_position(iMAP->textBuffer->length());
	iMAP->textDisplay->show_insert_position();

	Fl::unlock();
	Fl::check();
}
