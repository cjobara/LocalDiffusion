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

#ifndef TIFF_H_
#define TIFF_H_

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

#include "gui.h"
#include "draw.h"

class OverlayTiffGui;
class TiffSequence;

int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression);
void saveScreenshot();
void captureScreen(char* captureScreenName);
void writetiffOffscreen(char *offscreenFileName, Fl_Widget* widget, char* extension);
void saveScreenSequence();
void readTiffImage(char *name, int width, int height, unsigned short bitsPerSample, unsigned int samplesPerPixel, int slice, GLuint *texture);
void getTiffAttributes(char *name, int *width, int *height, unsigned short *bitsPerSample, unsigned short *samplesPerPixel, int *n);

// Overlay TIFF GUI related
void overlayTiffWindowCallback(Fl_Window*w,void*);
void adjustImageCallback(Fl_Widget*,void*);
void yFlipButtonCallback(Fl_Button*w, void*v);
void xFlipButtonCallback(Fl_Button*w, void*v);
void clampButtonCallback(Fl_Check_Button*w,void*v);
void sequenceSliderCallback(Slider*w,void*v);
void syncButtonCallback(Fl_Check_Button*w,void*v);

class OverlayTiffGui {
	public:
		OverlayTiffGui();

		Fl_Window *window;
		Fl_Tabs *tabs;

		// geometry
		Fl_Group *geometryGroup;
		Slider *resolutionSlider;
		Slider *xAlignSlider;
		Slider *yAlignSlider;
//		Slider *zAlignSlider;
		Fl_Button *xFlipButton;
		Fl_Button *yFlipButton;
//		Fl_Button *zFlipButton;
		Fl_Button *transposeButton;

		// image
		Fl_Group *imageGroup;
		Slider *alphaSlider;
		Slider *contrastSlider;
		Slider *gammaSlider;
		Slider *brightnessSlider;
		Slider *maxSlider;
		Slider *minSlider;
		Fl_Button *normalizeButton;
		Fl_Check_Button *clampButton;

		// sequence
		Fl_Group *sequenceGroup;
		Fl_Check_Button *synchronizeButton;
		Slider* sequenceSlider;
		Fl_Button *playButton;
		Fl_Button *stopButton;
		bool playEnable;

		int count;

		void show() { window->show(); }
		void hide() { window->hide(); }
		int visible() { return window->visible(); }
};

class TiffSequence {
	public:
		TiffSequence(char *loadedFilename);
		char filename[FILENAME_MAX];
		int stackSize;
		int width, height;
		unsigned short bitsPerSample;
		unsigned short samplesPerPixel;
		GLuint texture;
		int currentImage;
		int slice;
		float resolution;
		int saturation;
		float xMin,xMax,yMin,yMax;
		float xFlip,yFlip;

		void adjust();
		void draw(int slice);
		void setSlice(int s) { slice = s; }
		void clip() const;
		void unclip() const;
};

#endif /* TIFF_H_ */
