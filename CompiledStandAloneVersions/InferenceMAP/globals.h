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

#ifndef GLOBALS_H_
#define GLOBALS_H_

// Standard Libraries
#include <stdio.h>
#include <iostream>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>

// OpenGL Libraries
#ifdef __APPLE__
	#include <GL/glew.h>
	#include <OpenGL/OpenGL.h>
//	#include <OpenGL/gl.h>
//	#include <OpenGL/glu.h>
	#include <GLUT/glut.H>
	#include <FL/gl.h>
#elif _WIN32
	#include <GL/glew.h>
	#include <GL/glext.h>
	#include <glut.h>
#endif

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

// Image Libraries
#include <tiffio.h>
#include <omp.h>

#define MAX_FILES 10

class Globals;

// User Header Files
#include "graphics.h"
#include "file.h"
#include "tabs.h"
#include "gui.h"
#include "text.h"
#include "annotations.h"
#include "movie.h"
#include "draw.h"

class Globals {

	public:

		// 3D Trajectory File
		bool file3D;

		// Stereo Vision
		StereoVisionGui *stereoVisionGui;
		bool stereoVision;

		// Inference Calculation
		bool pauseCalculation;
		bool stopCalculation;
		bool smoothingPriorActive;

		// Movie Maker
		MovieMakerGui *movieMakerGui;
		int operationArray[MAX_ROWS];
		int axisArray[MAX_ROWS];
		float magnitudeArray[MAX_ROWS];
		float incrementArray[MAX_ROWS];
		int framesArray[MAX_ROWS];
		float overlapArray[MAX_ROWS];
		int framesPlayedArray[MAX_ROWS];
		bool movieMakerPlayEnable;
		bool movieMakerRecordEnable;
		bool movieMakerOpened;
		int currentRow;
		float landscapeAlpha;
		float landscapeScale;

		// Intro Window
		bool gaussianBell;
		int bellDimensions;
		float *bellQuads;
		float *bellColors;
		float *bellNormals;
		float defaultRotate;
		float defaultZoom;
		bool showLogo;

		// Main Window
		Fl_Group *mainGroup;
		Fl_Window *mainWindow;
		Fl_Menu_Bar *menubar;
		iMAPGlWindow *glWindow;
		FileInfo *fileInfo;
		int screenDimensions[4];
		int tabBarWidth;
		Fl_Button *fullScreenButton;
		int macOffset;
		ColorBar *colorbar;
		Fl_Box *colorbarCanvas;
		clock_t startTime,endTime;
		Fl_Box *textDisplayBox;
		Fl_Text_Display *textDisplay;
		Fl_Text_Buffer *textBuffer;
		char outputStream[1024];
		Fl_Check_Button *updateDisplayButton;
		bool updateDisplay;
		bool updatePotentialCalculation;
		Fl_Box *visualizationBox;
		Fl_Box *trajectoriesBox;
		Fl_Font normalFont;
		Fl_Font boldFont;
		Fl_Color bgColor;

		// Trajectory Inference
//		BatchTrajectoryInferenceGui *batchTrajectoryInferenceGui;
//		bool batchTrajectoryInferenceMode;
		SingleTrajectoryInferenceGui *singleTrajectoryInferenceGui;
		bool singleTrajectoryInferenceMode;
		CustomSelectionInferenceGui *customSelectionInferenceGui;
		bool customSelectionInferenceMode;
		bool selectionMade;
		bool selectionButtonPressed;
		bool landscapeAxesEnable;
		DensityGui *densityGui;

		// Annotations
		int activeZones;
		bool overlayAdjustment;
		float *linesOverlay;
		float *quadsOverlay;
		float *quadsColor;
		bool zoomBoxEnable;
		bool whiteBackground;

		// Visualization widgets
		Slider *iOffsetSlider;
		Slider *cMinSlider;
		Slider *cMaxSlider;
		Fl_Check_Button *localizationsButton;
		Fl_Check_Button *drawTrajectoriesButton;
		Fl_Check_Button *animateTrajectoriesButton;
		Fl_Check_Button *animateAccumulationButton;
		Fl_Check_Button *animateDelayedButton;
		Slider *fpsSlider;
		Slider *animateDelayedStepsSlider;
		AnnotationsGui *annotationsGui;
		Slider *startIntervalSlider;
		Slider *endIntervalSlider;
		Fl_Box *intervalBox;
		double backgroundRGB[3];

		// Global file parameters
		int fileNumber;
		int currentFile;
		char fileNameStorage[MAX_FILES][FILENAME_MAX];
		Tabs *tabButtonArray[MAX_FILES];
		File *fileObjectArray[MAX_FILES];
		float xMinRange,xMaxRange;
		float yMinRange,yMaxRange;

		// For file formats where exposure time is not specified
		Fl_Window *fileParametersWindow;
		Fl_Float_Input *exposureTimeInput;
		Fl_Button *fileParametersButton;
		Fl_Float_Input *pixelSizeInput;
		double pixelSize;
		double exposureTime;

		// Texture ID
		GLuint gaussianTex;
		TexFont *fontTex;
		char cwd[FILENAME_MAX];

		// OpenGL Display Parameters
		double zoomBox[4];
		double mv[16];
		short units;

		// 3D file format
		double fovy;
		double initialLookAt;

		// Monte Carlo Sampling
		long seed;

		// Recording trajectory
		char captureScreenName[FILENAME_MAX];
		int tiffCount;

		// Overlay TIFF
		OverlayTiffGui *overlayTiffGui;
		TiffSequence *tiffSequence;
		bool tiffOverlay;
		bool synchronizeTiff;

		bool saveTiffBatch;

		// Constructor
		Globals();

};

#endif /* GLOBALS_H_ */
