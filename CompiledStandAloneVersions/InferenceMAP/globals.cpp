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

#include "globals.h"

Globals::Globals() {

	smoothingPriorActive = false;

	whiteBackground = false;
	file3D = false;

	stereoVision = false;
	stereoVisionGui = NULL;

	saveTiffBatch = false;

	showLogo = false;
	bgColor = fl_darker(FL_DARK3);

	normalFont = FL_HELVETICA;
	boldFont = FL_HELVETICA_BOLD;

	visualizationBox = NULL;
	trajectoriesBox = NULL;
	textDisplayBox = NULL;

	landscapeAxesEnable = true;
	updateDisplayButton = NULL;
	updateDisplay = true;
	updatePotentialCalculation = false;

	startTime = 0.0;
	endTime = 0.0;
	pauseCalculation = false;
	stopCalculation = false;

	sprintf(outputStream,"InferenceMAP v1.0\n\nMohamed El Beheiry\nJean-Baptiste Masson\n\nInstitut Curie\nInstitut Pasteur\n\nCopyright (c) 2015\n\nFor academic use only.");

	movieMakerGui = NULL;
	movieMakerPlayEnable = false;
	movieMakerRecordEnable = false;
	movieMakerOpened = false;
	currentRow = 0;
	landscapeAlpha = 0.0;
	landscapeScale = 0.0;

	textDisplay = NULL;
	textBuffer = NULL;

	defaultZoom = 0.0;
	defaultRotate = 0.0;
	gaussianBell = false;
	bellDimensions = 0;
	bellQuads = NULL;
	bellColors = NULL;
	bellNormals = NULL;

	movieMakerOpened = false;

	backgroundRGB[0] = 0.0;
	backgroundRGB[1] = 0.0;
	backgroundRGB[2] = 0.0;

	linesOverlay = NULL;
	quadsOverlay = NULL;
	quadsColor = NULL;

	overlayAdjustment = true;
	fovy = 0.0;
	initialLookAt = 23.5;

	fileParametersButton = NULL;
	exposureTimeInput = NULL;
	fileParametersWindow = NULL;
	tiffSequence = NULL;
	pixelSizeInput = NULL;
	colorbarCanvas = NULL;

	pixelSize = 160.0;
	exposureTime = 0.025;

	activeZones = 0;
	colorbar = NULL;
	macOffset = 0;
	fileInfo = NULL;
	mainGroup = NULL;
	mainWindow = NULL;
	menubar = NULL;
	glWindow = NULL;
	tabBarWidth = 680;
	fullScreenButton = NULL;
	screenDimensions[2] = 1110;
	screenDimensions[3] = 730;
	units = 0; // um by default

	iOffsetSlider = NULL;
	cMinSlider = NULL;
	cMaxSlider = NULL;

	localizationsButton = NULL;
	drawTrajectoriesButton = NULL;
	animateTrajectoriesButton = NULL;

	annotationsGui = NULL;

//	batchTrajectoryInferenceGui = NULL;
//	batchTrajectoryInferenceMode = false;
	singleTrajectoryInferenceGui = NULL;
	singleTrajectoryInferenceMode = false;
	customSelectionInferenceGui = NULL;
	customSelectionInferenceMode = false;
	selectionMade = false;
	selectionButtonPressed = false;

	zoomBox[0] = 0.0;
	zoomBox[1] = 0.0;
	zoomBox[2] = 0.0;
	zoomBox[3] = 0.0;

	xMinRange = -200000.0;
	xMaxRange = 200000.0;
	yMinRange = -200000.0;
	yMaxRange = 200000.0;

	fileNumber = 0;
	currentFile = 0;

	tiffCount = 0;

	seed = -10678;

	overlayTiffGui = NULL;
	tiffOverlay = false;
	synchronizeTiff = false;

	gaussianTex = 0;
	fontTex = 0;

	#ifdef __APPLE__
		fontTex = txfLoadFont("helvetica48.txf");
	#endif

	zoomBoxEnable = false;

	animateAccumulationButton = NULL;
	animateDelayedButton = NULL;
	fpsSlider = NULL;
	animateDelayedStepsSlider = NULL;

	intervalBox = NULL;
	startIntervalSlider = NULL;
	endIntervalSlider = NULL;
}
