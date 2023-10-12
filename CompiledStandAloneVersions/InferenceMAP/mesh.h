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

#ifndef MESH_H_
#define MESH_H_

#include <time.h>

// FLTK Libraries
#include <FL/Fl_Group.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl.H>
#include <FL/Fl_Native_File_Chooser.h>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Group.H>

#include "tessellation.h"
#include "trees.h"
#include "draw.h"
#include "structs.h"
#include "selection.h"

class SquareCell;
class VoronoiCell;
class GaussianMixtureCell;

class TreeMesh;
class GaussianMixtureMesh;
class RandomizedOptimizationGui;

void randomizedOptimizationSaveCallback(Fl_Button*w,int*v);

class SquareMesh {

	double dx;
	double dy;
	int xCells;
	int yCells;
	double xMin,xMax;
	double yMin,yMax;
	double tMin,tMax;
	double xRange,yRange;
	int dimensions;
	int totalVariables;
	bool selectionEnable;
	int iterations[1];
	double fret[1];
	double sigma;
	double xMean,yMean,dtMean;
	long seed;

	double *potentialArray;

	SquareCell *maxCell;
	SquareCell *minCell;

	// stores cells in mesh
	SquareCell **cells;

	// inference
	int polynomialOrder;
	int numberOfCoefficients;
	int activeDetections;
	int activeCells;
	double forceMax;
	double forceMin;
	double diffusionMax;
	double diffusionMin;
	double diffusionLogMax;
	double diffusionLogMin;
	double potentialMax;
	double potentialMin;
	double speedMin;
	double speedMax;
	int optimizationMode;
	int currentZoneX,currentZoneY;
	int currentZone2X,currentZone2Y;
	double beta;
	bool zonalPotentialsCalculated;
	float maximumNeighbourDistance;
	bool smoothingPriorEnable;
	bool jeffreysPriorEnable;
	int minPointsPerCell;

	public:
		SquareMesh(double dx);
		SquareMesh(double dx, SelectionCell selection);

		SelectionCell selection;
		double **hessian;
		double cost;

		// landscape view variables
		float *landscapeTriangles;
		float *landscapeNormals;
		float *landscapeColors;

		// for D,gradV calculations
		int **localizationAssignment;

		RandomizedOptimizationGui *randomizedOptimizationGui;

		// stochastic gradient descent
		int roZones;
		int *roXCoords;
		int *roYCoords;
		bool *previousRoSelection;
		bool *currentRoSelection;
		bool roEnable;

		bool consecutiveRoSelection(int id) {
			for (int q = 0; q < roZones; q++) {
				if (previousRoSelection[id]) {
					return true;
				}
			}
			return false;
		}

		float *dvGradxArray;
		float *dvGradyArray;
		int *identityArray;
		float *optimizationArrayFloat;
		double *optimizationArray;
		double *intermediateOptimizationArray;
		float *costArray;

		// accessors
		double getDx() { return dx; }
		double getDy() { return dy; }
		int getXCells() { return xCells; }
		int getYCells() { return yCells; }
		double getXCentre() { return (xMin+xMax)/2.0; }
		double getYCentre() { return (yMin+yMax)/2.0; }
		SquareCell* getMaxCell() { return maxCell; }
		SquareCell* getMinCell() { return minCell; }
		int getPolynomialOrder() { return polynomialOrder; }
		int getCoefficients() { return numberOfCoefficients; }
		SquareCell* getCell(int x, int y);
		SquareCell* getCell(int id);
		SquareCell* getCellRo(int id);
		bool roActive(int x, int y);
		bool active(int x, int y);
		int getCellCount(int x, int y);
		double getTMin() { return tMin; }
		double getTMax() { return tMax; }
		int getRoIdentifier(int x, int y);
		int getIdentifier(int x, int y);
		int getDimensions() { return dimensions; }
		double getSigma() { return sigma; }
		double getOptimizationArray(int index) {
			if (index < dimensions) { return optimizationArray[index]; }
			else { return 0.0; }
		}
		double getForceMax() { return forceMax; }
		double getForceMin() { return forceMin; }
		double getDiffusionMax() { return diffusionMax; }
		double getDiffusionMin() { return diffusionMin; }
		double getDiffusionLogMax() { return diffusionLogMax; }
		double getDiffusionLogMin() { return diffusionLogMin; }
		double getPotentialMax() { return potentialMax; }
		double getPotentialMin() { return potentialMin; }
		double getSpeedMin() { return speedMin; }
		double getSpeedMax() { return speedMax; }
		int getVariables() { return totalVariables; }
		int getCount(int x, int y);
		double getForceX(int x, int y);
		double getForceY(int x, int y);
		double getForce(int x, int y);
		double getDiffusion(int x, int y);
		double getDiffusionLog(int x, int y);
		double getPotential(int x, int y);
		double getSpeed(int x, int y);
		int getActiveDetections() { return activeDetections; }
		int getActiveCells() { return activeCells; }
		int getCurrentZoneX() { return currentZoneX; }
		int getCurrentZoneY() { return currentZoneY; }
		int getCurrentZone2X() { return currentZone2X; }
		int getCurrentZone2Y() { return currentZone2Y; }
		int getOptimizationMode() { return optimizationMode; }
		int getXNeighbours(int x, int y);
		int getYNeighbours(int x, int y);
		double getXMean(int x, int y);
		double getYMean(int x, int y);
		int getXNeighbourPos(int x, int y);
		int getYNeighbourPos(int x, int y);
		double getBeta() { return beta; }
		double getDx(int x, int y, int i);
		double getDy(int x, int y, int i);
		double getDt(int x, int y, int i);
		double getDtMean() { return dtMean; }
		bool selectionMode() { return selectionEnable; }
		float getMaximumNeighbourDistance() { return maximumNeighbourDistance; }
		bool smoothingPriorEnabled() { return smoothingPriorEnable; }
		bool jeffreysPriorEnabled() { return jeffreysPriorEnable; }
		int getMinPointsPerCell() { return minPointsPerCell; }
		void enableRo() { roEnable = true; }
		void disableRo() { roEnable = false; }
		double getCentroidX(int x, int y);
		double getCentroidY(int x, int y);
		int getTotalVariables() { return totalVariables; }

		// mutators
		void setXMin(float value) { xMin = value; }
		void setXMax(float value) { xMax = value; }
		void setYMin(float value) { yMin = value; }
		void setYMax(float value) { yMax = value; }
		void setTMin(float value) { tMin = value; }
		void setTMax(float value) { tMax = value; }
		void activateCell(int x, int y);
		void deactivateCell(int x, int y);
		void setCurrentZone(int x, int y) { currentZoneX = x; currentZoneY = y; }
		void setCurrentZone2(int x, int y) { currentZone2X = x; currentZone2Y = y; }
		void clearArrays();
		void setMaximumNeighbourDistance(float mnd) { maximumNeighbourDistance = mnd; }
		void enableSmoothingPrior() { smoothingPriorEnable = true; }
		void disableSmoothingPrior() { smoothingPriorEnable = false; }
		void enableJeffreysPrior() { jeffreysPriorEnable = true; }
		void disableJeffreysPrior() { jeffreysPriorEnable = false; }
		void setOptimizationMode(int m) { optimizationMode = m; }

		void updateNeighbours();

		void saveWork();
		void savePosterior();

		void offsetPotentials();

		void findExtremesD();
		void findExtremesDF();
		void findExtremesDDr();
		void findExtremesDV();

		// inference methods
		void infer();

		// D Inference
		void preInferD();
		void inferD();
		int postInferD();
		void preInferDSmoothing();
		void inferDSmoothing();
		int postInferDSmoothing();

		// DF Inference
		void preInferDF();
		void inferDF();
		int postInferDF();
		void inferPotentialsDF();

		void preInferDFSmoothing();
		void inferDFSmoothing();
		int postInferDFSmoothing();

		// DDr Inference
		void preInferDDr();
		void inferDDr();
		int postInferDDr();

		void preInferDDrSmoothing();
		void inferDDrSmoothing();
		int postInferDDrSmoothing();

		// DV Inference
		void preInferDV();
		void inferDV();
		int postInferDV();

		// global polynomial potential optimization
		void preInferPolynomial();
		void inferPolynomial();
		void postInferPolynomial();

		void inferRandomizedOptimizationD();
		int inferRandomizedOptimizationDF();
		void inferRandomizedOptimizationDDr();
		void inferRandomizedOptimizationDV();

		// Import/Export
		void exportMesh();

		virtual ~SquareMesh();

		friend class SquareMeshGui;
		friend class Simulation;
		friend class RandomizedOptimizationGui;
};


class VoronoiMesh {

	double xMin,xMax;
	double yMin,yMax;
	double tMin,tMax;
	int dimensions;
	int totalVariables;
	bool slidingWindowEnable;
	int iterations[1];
	double fret[1];
	double sigma;
	double xMean,yMean,dtMean;

	double *intermediateOptimizationArray;
	double *optimizationArray;
	double *potentialArray;

	VoronoiCell *maxCell;
	VoronoiCell *minCell;

	// stores cells in mesh
	VoronoiCell *cells;
	int cornerCellIndices[4];
	float charDistance;

	// inference
	int polynomialOrder;
	int numberOfCoefficients;
	int activeDetections;
	int activeCells;
	double forceMax;
	double forceMin;
	double diffusionMax;
	double diffusionMin;
	double diffusionLogMax;
	double diffusionLogMin;
	double potentialMax;
	double potentialMin;
	float areaMax;
	float areaMin;
	float perimeterMax;
	float perimeterMin;
	int optimizationMode;
	int currentZone;
	int currentZone2;
	double beta;
	bool zonalPotentialsCalculated;
	bool selectionEnable;
	double maximumNeighbourDistance;
	bool smoothingPriorEnable;
	bool jeffreysPriorEnable;
	int minPointsPerCell;

	// kmeans clustering
	// inputs
	int spatialDimensions;
	int maxIterations;
	int numberOfClusters;
	double *point;
	int distanceType;
	// outputs
	double *clusterCentres;
	int *clusterIndices;
	double *clusterEnergies;
	int *clusterPopulations;
	double *clusterVariances;
	double minVariance,maxVariance;
	int numberOfIterations;
	long seed;
	int clusteringMethod;

	public:
		VoronoiMesh();
		VoronoiMesh(SelectionCell selection);

		// rendering functions
		float *colorArray;
		double **hessian;
		float *linesOverlay;
		float *polygonsOverlay;

		// landscape view variables
		int landscapeVertices;
		float *landscapeTriangles;
		float *landscapeNormals;
		float *landscapeColors;

		SelectionCell selection;

		bool roEnable;
		int roZones;
		int *roClusters;

		void saveWork();
		void savePosterior();

		// accessors
		int getActiveCellId(int counter);
		void enableRo() { roEnable = true; }
		void disableRo() { roEnable = false; }
		bool roActive(int c);
		bool selectionMode() { return selectionEnable; }
		double getXCentre() { return (xMin+xMax)/2.0; }
		double getYCentre() { return (yMin+yMax)/2.0; }
		VoronoiCell* getMaxCell() { return maxCell; }
		VoronoiCell* getMinCell() { return minCell; }
		int getPolynomialOrder() { return polynomialOrder; }
		int getCoefficients() { return numberOfCoefficients; }
		VoronoiCell* getCell(int i);
		bool active(int i);
		int getCellCount(int xi);
		double getTMin() { return tMin; }
		double getTMax() { return tMax; }
		int getIdentifier(int i);
		bool slidingWindow() { return slidingWindowEnable; }
		int getDimensions() { return dimensions; }
		double getSigma() { return sigma; }
		double getOptimizationArray(int index) {
			if (index < dimensions) { return optimizationArray[index]; }
			else { return 0.0; }
		}
		double getForceMax() { return forceMax; }
		double getForceMin() { return forceMin; }
		double getDiffusionMax() { return diffusionMax; }
		double getDiffusionMin() { return diffusionMin; }
		double getDiffusionLogMax() { return diffusionLogMax; }
		double getDiffusionLogMin() { return diffusionLogMin; }
		double getPotentialMax() { return potentialMax; }
		double getPotentialMin() { return potentialMin; }
		double getDiffusion(int i);
		double getDiffusionLog(int i);
		double getPotential(int i);
		double getForceX(int i);
		double getForceY(int i);
		double getForce(int i);
		int getCount(int i);
		int getActiveDetections() { return activeDetections; }
		int getActiveCells() { return activeCells; }
		int getCurrentZone() { return currentZone; }
		int getCurrentZone2() { return currentZone2; }
		double getCentroidX(int i);
		double getCentroidY(int i);
		int getOptimizationMode() { return optimizationMode; }
		double getMinVariance() { return minVariance; }
		double getMaxVariance() { return maxVariance; }
		double getBeta() { return beta; }
		int getDistanceType() { return distanceType; }
		int getNumberOfClusters() { return numberOfClusters; }
		int getClusterIndex(int i) { return clusterIndices[i]; }
		int getClusteringMethod() { return clusteringMethod; }
		float getCharDistance() { return charDistance; }
		int getCellSelection(double x, double y);
		int getNumberOfNeighbours(int i);
		int getVariables() { return totalVariables; }
		float getMaxArea() { return areaMax; }
		float getMinArea() { return areaMin; }
		float getPerimeterMax() { return perimeterMax; }
		float getPerimeterMin() { return perimeterMin; }
		float getAreaMax() { return areaMax; }
		float getAreaMin() { return areaMin; }
		double getMaximumNeighbourDistance() { return maximumNeighbourDistance/1000.0; }
		bool smoothingPriorEnabled() { return smoothingPriorEnable; }
		bool jeffreysPriorEnabled() { return jeffreysPriorEnable; }
		int getMinPointsPerCell() { return minPointsPerCell; }
		float getMaxPerimeter() { return perimeterMax; }
		float getMinPerimeter() { return perimeterMin; }
		double getDtMean() { return dtMean; }
		int getTotalVariables() { return totalVariables; }

		// mutators
		void setXMin(float value) { xMin = value; }
		void setXMax(float value) { xMax = value; }
		void setYMin(float value) { yMin = value; }
		void setYMax(float value) { yMax = value; }
		void setTMin(float value) { tMin = value; }
		void setTMax(float value) { tMax = value; }
		void activateCellClick(int i);
		void deactivateCellClick(int i);
		void setCurrentZone(int z) { currentZone = z; }
		void setCurrentZone2(int z) { currentZone2 = z; }
		void clearArrays();
		void setDistanceType(int d) { distanceType = d; }
		void setClusteringMethod(int i) { clusteringMethod = i; }
		void setMaximumNeighbourDistance(double mnd) { maximumNeighbourDistance = mnd; }
		void enableSmoothingPrior() { smoothingPriorEnable = true; }
		void disableSmoothingPrior() { smoothingPriorEnable = false; }
		void enableJeffreysPrior() { jeffreysPriorEnable = true; }
		void disableJeffreysPrior() { jeffreysPriorEnable = false; }
		void setOptimizationMode(int m) { optimizationMode = m; }

		// misc
		void updateNeighbours();

		void findExtremesD();
		void findExtremesDF();
		void findExtremesDDr();
		void findExtremesDV();

		void offsetPotentials();

		// inference methods
		void infer();

		// kmeans clustering
		void createClusters();
		void createClustersSelection();
		double* constantPointsInitialization();
		double* constantPointsInitializationSelection();

		// global polynomial potential optimization
		void preInferPolynomial();
		void inferPolynomial();
		void postInferPolynomial();

		// D Inference
		void preInferD();
		void inferD();
		int postInferD();
		void preInferDSmoothing();
		void inferDSmoothing();
		int postInferDSmoothing();

		// DF Inference
		void preInferDF();
		void inferDF();
		int postInferDF();
		void inferPotentialsDF();

		void preInferDFSmoothing();
		void inferDFSmoothing();
		int postInferDFSmoothing();

		// DDr Inference
		void preInferDDr();
		void inferDDr();
		int postInferDDr();

		void preInferDDrSmoothing();
		void inferDDrSmoothing();
		int postInferDDrSmoothing();

		// DV Inference
		void preInferDV();
		void inferDV();
		int postInferDV();
		void inferRandomizedOptimizationDV();
		void inferRandomizedOptimizationD();
		void inferRandomizedOptimizationDDr();
		int inferRandomizedOptimizationDF();
		RandomizedOptimizationGui *randomizedOptimizationGui;

		void findNeighbours(int c);

		// voronoi diagram
		float *xCentreList;
		float *yCentreList;
		VoronoiDiagramGenerator voronoiDiagram;

		void exportMesh();
		void exportMeshBatch();

		virtual ~VoronoiMesh();

		friend class VoronoiCell;
		friend class Simulation;
		friend class RandomizedOptimizationGui;
};

class TreeMesh {

	double xMin,xMax;
	double yMin,yMax;
	double zMin,zMax;
	double tMin,tMax;
	int dimensions;
	bool slidingWindowEnable;
	bool selectionEnable;
	int iterations[1];
	double fret[1];
	double sigma;
	double xMean,yMean,zMean,dtMean;

	double *intermediateOptimizationArray;
	double *optimizationArray;
	double *potentialArray;

	// inference
	int polynomialOrder;
	int numberOfCoefficients;
	int activeDetections;
	int activeCells;
	double forceMax;
	double forceMin;
	double diffusionMax;
	double diffusionMin;
	double diffusionLogMax;
	double diffusionLogMin;
	double potentialMax;
	double potentialMin;
	int optimizationMode;
	QuadTree* currentQuadZone;
	double beta;
	bool zonalPotentialsCalculated;
	double maximumNeighbourDistance;
	bool smoothingPriorEnable;
	bool jeffreysPriorEnable;
	int minPointsPerCell;
	int minLeafPower;
	int minCapacity;


	bool roEnable;
	int roZones;

	QuadTree* maxQuadTree;
	QuadTree* minQuadTree;

	public:
		TreeMesh();
		TreeMesh(SelectionCell selection);
		int totalVariables;
		int allVariables;

		SelectionCell selection;
		RandomizedOptimizationGui *randomizedOptimizationGui;

		// rendering functions
		float *linesOverlay;
		float *quadsOverlay;
		float *quadsColor;

		// landscape view variables
		int landscapeVertices;
		float *landscapeTriangles;
		float *landscapeNormals;
		float *landscapeColors;
		int landscapeBorderVertices;
		float *landscapeBorderTriangles;
		float *landscapeBorderNormals;
		float *landscapeBorderColors;

		// misc
		void updateNeighbours(QuadTree* tree);

		QuadTree *quadTree;
		QuadTree **quadTreeLocalizationPointer;
		void initializeTreeMesh(QuadTree* tree);
		void applyTreeMesh(QuadTree* tree);
		void updateTreeMesh(QuadTree* tree);
		void initializePotentials(QuadTree* tree);
		void resetNeighbourCount(QuadTree* tree);
		void findNeighbours(QuadTree* tree);
		void findNeighbours(QuadTree* tree, QuadTree* root);
		void assignNeighbours(QuadTree* tree);
		void assignNeighbours(QuadTree* tree, QuadTree* root);
		void assignPotentials(QuadTree* tree);
		void assignAreas(QuadTree* tree);
		void roDeactivateZones(QuadTree* tree);
		void enableRo() { roEnable = true; }
		void disableRo() { roEnable = false; }
		void getTreeRo(QuadTree *tree,int identifier);
		void adjustPotentials(QuadTree* tree);
		void initializeTreeMesh(QuadTree* tree,SelectionCell selection);
		void initializePosterioriSelectionZones(QuadTree* tree, int *sel);

		double **hessian;

		int minCount;
		int maxCount;

		int maxQuadTreePower;
		int maxOctTreePower;
		QuadTree* selectedQuadLeaf;
		QuadTree* selectedQuadLeaf2;
		bool leafSelected;
		bool leaf2Selected;
		int identifierCounter;

		QuadTree* identifiedQuadLeaf;

		void saveWork();
		void savePosterior();

		// accessors
		double getXCentre() { return (xMin+xMax)/2.0; }
		double getYCentre() { return (yMin+yMax)/2.0; }
		double getZCentre() { return (zMin+zMax)/2.0; }
		double getXMax() { return xMax; }
		double getXMin() { return xMin; }
		double getYMax() { return yMax; }
		double getYMin() { return yMin; }
		QuadTree* getMaxQuadTree() { return maxQuadTree; }
		QuadTree* getMinQuadTree() { return minQuadTree; }

		int getPolynomialOrder() { return polynomialOrder; }
		int getCoefficients() { return numberOfCoefficients; }

		bool activated(int x, int y);
		int getCellCount(int x, int y);
		double getTMin() { return tMin; }
		double getTMax() { return tMax; }
		int getIdentifier(int x, int y);
		bool slidingWindow() { return slidingWindowEnable; }
		int getDimensions() { return dimensions; }
		double getSigma() { return sigma; }
		double getOptimizationArray(int index) {
			if (index < dimensions) { return optimizationArray[index]; }
			else { return 0.0; }
		}
		double getForceMax() { return forceMax; }
		double getForceMin() { return forceMin; }
		double getDiffusionMax() { return diffusionMax; }
		double getDiffusionMin() { return diffusionMin; }
		double getDiffusionLogMax() { return diffusionLogMax; }
		double getDiffusionLogMin() { return diffusionLogMin; }
		double getPotentialMax() { return potentialMax; }
		double getPotentialMin() { return potentialMin; }
		int getVariables() { return totalVariables; }

		double getForceX(QuadTree *tree) { return tree->getFx(); }
		double getForceY(QuadTree *tree) { return tree->getFy(); }
		double getDiffusion(QuadTree *tree) { return tree->getD(); }
		double getDiffusionLog(QuadTree *tree) { return tree->getDlog(); }
		double getPotential(QuadTree *tree) { return tree->getPotential(); }
		int getActiveDetections() { return activeDetections; }
		int getActiveCells() { return activeCells; }
		QuadTree* getCurrentQuadZone() { return currentQuadZone; }
		void getTree(QuadTree *tree, int identifier);
		int getOptimizationMode() { return optimizationMode; }
		int getXNeighbours(int x, int y);
		int getYNeighbours(int x, int y);
		double getXMean(int x, int y);
		double getYMean(int x, int y);
		int getXNeighbourPos(int x, int y);
		int getYNeighbourPos(int x, int y);
		double getBeta() { return beta; }
		void getLeafSelection(double x, double y, QuadTree *tree, bool *found);
		void getLeafSelection2(double x, double y, QuadTree *tree, bool *found);
		bool selectionMode() { return selectionEnable; }
		double getMaximumNeighbourDistance() { return maximumNeighbourDistance; }
		bool smoothingPriorEnabled() { return smoothingPriorEnable; }
		bool jeffreysPriorEnabled() { return jeffreysPriorEnable; }
		int getMinNumberOfPointsPerCell() { return minPointsPerCell; }
		int getMinLeafPower() { return minLeafPower; }
		int getMinCapacity() { return minCapacity; }
		double getDtMean() { return dtMean; }
		int getTotalVariables() { return totalVariables; }

		// mutators
		void setXMin(float value) { xMin = value; }
		void setXMax(float value) { xMax = value; }
		void setYMin(float value) { yMin = value; }
		void setYMax(float value) { yMax = value; }
		void setZMin(float value) { zMin = value; }
		void setZMax(float value) { zMax = value; }
		void setTMin(float value) { tMin = value; }
		void setTMax(float value) { tMax = value; }
		void activateCell(int x, int y);
		void deactivateCell(int x, int y);
		void activateCell(int x, int y, int z);
		void deactivateCell(int x, int y, int z);
		void setCurrentZone(QuadTree* tree) { currentQuadZone = tree; }
		void clearArrays();
		void setMinTree(QuadTree* tree) { minQuadTree = tree; }
		void setMaxTree(QuadTree* tree) { maxQuadTree = tree; }
		void setMaximumNeighbourDistance(double mnd) { maximumNeighbourDistance = mnd; }
		void enableSmoothingPrior() { smoothingPriorEnable = true; }
		void disableSmoothingPrior() { smoothingPriorEnable = false; }
		void enableJeffreysPrior() { jeffreysPriorEnable = true; }
		void disableJeffreysPrior() { jeffreysPriorEnable = false; }
		void setMinLeafPower(int m) { minLeafPower = m; }
		void setMinNumberOfPointsPerCell(int m) { minPointsPerCell = m; }
		void setMinCapacity(int c) { minCapacity = c; }
		void setOptimizationMode(int m) { optimizationMode = m; }
		void offsetPotentials();

		// inference methods
		void infer();

		void findExtremesD();
		void findExtremesDF();
		void findExtremesDDr();
		void findExtremesDV();

		// global polynomial potential optimization
		void preInferPolynomial();
		void inferPolynomial();
		void postInferPolynomial();
		void findExtremesPolynomial(QuadTree* tree);

		// D Inference
		void preInferD();
		void inferD(QuadTree* tree,int *clusterNumber);
		int postInferD();
		void preInferDSmoothing();
		void inferDSmoothing(QuadTree* tree);
		int postInferDSmoothing();

		void findExtremesD(QuadTree* tree);
		void findExtremesDSmoothing(QuadTree* tree);

		// DF Inference
		void preInferDF();
		void inferDF(QuadTree* tree,int *clusterNumber);
		int postInferDF();
		void inferPotentialsDF();
		void inferPotentialsDF(QuadTree* tree);

		void preInferDFSmoothing();
		void inferDFSmoothing(QuadTree* tree);
		int postInferDFSmoothing();

		void findExtremesDF(QuadTree* tree);
		void findExtremesDFSmoothing(QuadTree* tree);

		// DDr Inference
		void preInferDDr();
		void inferDDr(QuadTree* tree,int *clusterNumber);
		int postInferDDr();

		void preInferDDrSmoothing();
		void inferDDrSmoothing(QuadTree* tree);
		int postInferDDrSmoothing();

		void findExtremesDDr(QuadTree* tree);
		void findExtremesDDrSmoothing(QuadTree* tree);

		// DV Inference
		void preInferDV();
		void inferDV(QuadTree* tree);
		int postInferDV();

		void findExtremesDV(QuadTree* tree);

		void inferRandomizedOptimizationDV();
		void inferRandomizedOptimizationD();
		void inferRandomizedOptimizationDDr();
		int inferRandomizedOptimizationDF();

		void countWithinCircle(QuadTree *tree, float xCentre, float yCentre, float radius, int *count);
		void assignWithinCircleD(QuadTree *tree, float xCentre, float yCentre, float radius, int *count);
		void assignWithinCircleDF(QuadTree *tree, float xCentre, float yCentre, float radius, int *count);
		void assignWithinCircleDDr(QuadTree *tree, float xCentre, float yCentre, float radius, int *count);
		void assignWithinCircleDV(QuadTree *tree, float xCentre, float yCentre, float radius, int *count);

		int getRoZones() { return roZones; }
		QuadTree* getROQuadTree(int q) { return roQuadTrees[q]; }
		QuadTree **roQuadTrees;
		int *roIdentifiers;

		void exportMesh();
		void exportMesh(QuadTree* tree,FILE * writeFile);

		virtual ~TreeMesh();

		friend class TreeMeshGui;
		friend class Simulation;
};

class RandomizedOptimizationGui {
	private:

		Fl_Button *saveButton;
		DataNode *head;
		DataNode *tail;
		double maxCost;
		double minCost;
		int n;
	public:
		Fl_Window *window;

		RandomizedOptimizationGui();

		void finish() {
			saveButton->activate();
			saveButton->show();
			window->hide();
		}

		void save() {
			Fl_Native_File_Chooser *native = new Fl_Native_File_Chooser();
			char nativeFilename[FILENAME_MAX];

			native->title("Save Randomized Optimization Data");
			native->type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
			native->filter("Text File\t*.txt\n");
		    native->options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);

			if (native->show() == 0) {
				// store filename
				sprintf(nativeFilename,"%s%s",native->filename(),".txt");
				FILE * writeFile;

				writeFile = fopen (nativeFilename,"wb");

				// write file
				if (n > 1) {
					minCost = 100000000.0;
					maxCost = -100000000.0;
					DataNode *temp = head;
					while (temp != NULL) {
						if (minCost > temp->data) { minCost = temp->data; }
						if (maxCost < temp->data) { maxCost = temp->data; }
						temp = temp->next;
					}
				}

				int its = 1;
				if (n > 1) {
					DataNode *temp = head;
					while (temp != NULL) {
						fprintf(writeFile,"%i\t%f\n",its,temp->data);
						temp = temp->next;
						its++;
					}
				}

				fclose(writeFile);
			}

		}

		void update(double cost, int zones);

		void hide() { window->hide(); }

		void addNode(double data) {
			if (n == 0) {
				head->data = data;
				n++;
			} else {
				DataNode *newNode = new DataNode;
				newNode->next = NULL;
				newNode->data = data;

				DataNode *currentNode = head;
				while (currentNode) {
					if (currentNode->next == NULL) {
						currentNode->next = newNode;
						tail = currentNode->next;
						n++;
						return;
					}
					currentNode = currentNode->next;
				}
			}
		}

		void findExtremes() {
			if (n > 1) {
				minCost = 100000000.0;
				maxCost = -100000000.0;
				DataNode *temp = head;
				while (temp != NULL) {
					if (minCost > temp->data) { minCost = temp->data; }
					if (maxCost < temp->data) { maxCost = temp->data; }
					temp = temp->next;
				}
			}
		}
};

#endif /* MESH_H_ */
