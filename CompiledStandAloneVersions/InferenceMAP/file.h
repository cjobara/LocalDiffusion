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

#ifndef FILE_H_
#define FILE_H_

class File;

#include "simulation.h"
#include "structs.h"
#include "tabs.h"
#include "mesh.h"
#include "cell.h"
#include "trees.h"
#include "selection.h"
#include "gui.h"

class File {

	public:

	// constructor
	File(int files, char **list, int fileType);
	~File();

	// Density Calculation
	Partition **partitionMatrix;
	int xParts;
	int yParts;
	double unitLength;
	bool densityCalculated;

	// Square Meshing GUI
	SquareMeshGui *squareMeshGui;
	bool squareMeshOverlay;
	bool squareMeshApply;

	// Voronoi Meshing GUI
	VoronoiMeshGui *voronoiMeshGui;
	bool voronoiMeshOverlay;
	bool voronoiMeshApply;

	// Quad-Tree Meshing GUI
	TreeMeshGui *treeMeshGui;
	bool treeMeshOverlay;
	bool treeMeshApply;

	// menu variables
	Tabs *tabButton;

	// file storage variables
	int fileType;
	char **fileNameList;
	int numberOfFiles;
	char fileName [260];
	char filePath [260];
	int localizationCount;
	int localizationCountDisplay;
	int localizationCountVoronoi;
	int *localizationCountFile;
	double *xPointer;
	double *yPointer;
	double *iPointer;
	double *tPointer;
	double *xPointerHolder;
	double *yPointerHolder;
	double *iPointerHolder;
	double *tPointerHolder;
	int *nPointerHolder;
	bool *insideInterval;
	int *nPointer;
	double exposureTime;

	// file-specific data
	double xMin,xMax;
	double yMin,yMax;
	double iMin,iMax;
	double tMin,tMax;
	double tMinInitial,tMaxInitial;
	double fRange;
	double xFileOffset,yFileOffset;
	double xRange,yRange,zRange;
	double maxRange;
	int nMin,nMax;

	// 3d variables
	bool file3D;
	double zMin;
	double zMax;
	double *zPointer;
	double averageDz;
	double zFileOffset;
	int *indexArray; // for depth sorting
	float *coordsTemp;
	int *linesIndexArray;

	// render arrays
	float *pointsArray;
	float *linesArray;
	float *linesColorArray;
	float *colorArray;

	// trajectory storage variables
	int trackNumber;
	float *trackLengths;
	int *trackFrames;
	Trajectory *tracks;

	// display options
	int plotType;
	bool fileLoaded;
	float detectionSize;
	bool plotChanged;
	float orthoLimit;
	float xTranslate,yTranslate,zTranslate;
	float xRotate,yRotate,zRotate;
	bool drawTrajectories;
	bool drawLocalizations;
	int animationIndex;
	int animationProgress;
	bool animateTrajectories;
	float gridSpacing;
	bool flipColormap;
	bool landscapePlot;

	// intervals
	double startInterval;
	double endInterval;

	// annotations
	bool boundingBoxEnable;
	bool dimensionsEnable;
	bool gridEnable;
	bool ticksEnable;

	// thresholding parameters
	void setMinFrame(float minimum);
	void setMaxFrame(float maximum);
	int getVisibility(int index);
	int *tMinVisibility;
	int *tMaxVisibility;

	// colormap options
	int colormapType;
	bool colormapFlip;
	float cMin,cMax;
	float iOffset;

	// meshing pointers
	SquareMesh *squareMesh;
	VoronoiMesh *voronoiMesh;
	TreeMesh *treeMesh;
	int meshType;
	double averageDx,averageDy;
	double averageDt;

	// inference-related
	bool inferred;
	int optimizationMode;
	int optimizationFunction;

	// simulated trajectories
	Simulation *simulation;

	// selection
	Selection *selection;

	// random optimization
	float randomizedOptimizationTolerance;
	bool drawRandomizedOptimizationZones;

	// accessors

	void applySquareMesh() { squareMeshApply = true; }
	void unapplySquareMesh() { squareMeshApply = false; }
	bool squareMeshApplied() { return squareMeshApply; }
	void applyVoronoiMesh() { voronoiMeshApply = true; }
	void unapplyVoronoiMesh() { voronoiMeshApply = false; }
	bool voronoiMeshApplied() { return voronoiMeshApply; }
	void applyTreeMesh() { treeMeshApply = true; }
	void unapplyTreeMesh() { treeMeshApply = false; }
	bool treeMeshApplied() { return treeMeshApply; }

	double getXCentre() { return (xMin+xMax)/2.0; }
	double getYCentre() { return (yMin+yMax)/2.0; }

	int getMaxPoints() {
		int max = -1;
		for (int h = 0; h < numberOfFiles; h++) {
			if (max < localizationCountFile[h]) { max = localizationCountFile[h]; }
		}
		return max;
	}
	int getMinPoints() {
		int min = 10000000;
		for (int h = 0; h < numberOfFiles; h++) {
			if (min > localizationCountFile[h]) { min = localizationCountFile[h]; }
		}
		return min;
	}
	void updateAlpha(float alpha);
	void updateIntervalDisplay();

	int updateIntervalData();
	void reloadOriginalData();

	// file export options
	void saveFile(char* saveFileName, int fileType);
};

#endif /* FILE_H_ */
