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

#ifndef GRAPHICS_H_
#define GRAPHICS_H_

#include <FL/Fl_Gl_Window.H>
#include <FL/Fl.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Menu.H>
#include <FL/FL_Menu_Window.H>

#include "file.h"
#include "tiff.h"

class iMAPGlWindow;

void setupGLExtensions();
float* detectionImage(float sigma, int dimension);
void colormap(File *fileObject, int index, float value);
float* colormap(int colormapType, float value, float cMin, float cMax, bool flip);
void colormap(float *rgb, int colormapType, float value, float cMin, float cMax, bool flip);
void defaultScreen();
void drawGaussianBell();
void drawLandscapeAxes(int screen);
void drawAxes(int screen);
void drawZoomBox(float xRange, float yRange);
void Normalise(xyz *p);

class iMAPGlWindow : public Fl_Gl_Window {

	// handle mouse events
	int handle(int e);
	void GlInit();
	void draw();
	void drawMaster();
	void drawSlave();

	static void TimerCallback(void *userdata);

	// for handling translations
	double modelMatrix[16];
	double projMatrix[16];
	int viewport[4];

    // Constructor
    public:
	iMAPGlWindow(int x,int y,int w,int h,const char*l=0) : Fl_Gl_Window(x,y,w,h,l) {
		boxChanged = false;
		Fl::add_timeout(1.0f/40.0f,TimerCallback,(void*)this);
		// to avoid weird lines upon loading
		dragging = false;
		pushPosition[0] = pushPosition[1] = pushPosition[2] = pushPosition[3] = 0.0f;
		releasePosition[0] = releasePosition[1] = releasePosition[2] = releasePosition[3] = 0.0f;
		end();
	}

	// Class variables
	int currentX,currentY;
	int lastX,lastY;
	float lastZ,currentZ;
	float x1,y1,x2,y2;

	// position variables
	double pushPosition[4];
	double releasePosition[4];
	double cellSelect[4];

	// update variables
	bool doubleClick;
	bool boxChanged;
	bool dragging;

};

#endif /* GRAPHICS_H_ */
