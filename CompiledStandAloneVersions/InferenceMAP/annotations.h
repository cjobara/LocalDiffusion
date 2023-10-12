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

#ifndef ANNOTATIONS_H_
#define ANNOTATIONS_H_

#include <FL/Fl_Int_Input.H>

class AnnotationsGui;

// 2D Annotations
void drawBoundingBox();
void drawDimensions();
void drawGrid();
void drawTicks();

// 3D Annotations
void drawBoundingBox3D();

void annotationsWindowCallback(Fl_Widget*w,void*);
void boundingBoxCallback(Fl_Widget*,void*);
void dimensionsCallback(Fl_Widget*,void*);
void gridCallback(Fl_Widget*,void*);
void ticksCallback(Fl_Widget*,void*);
void displaySizeCallback(Fl_Widget*,void*);
void updateDisplaySizeCallback(Fl_Button*b,void*);
void umUnitsButtonCallback(Fl_Check_Button*w,int*v);
void nmUnitsButtonCallback(Fl_Check_Button*w,int*v);
void backgroundColorCallback(Fl_Button*w,void*);
int countDigits(int number);

class AnnotationsGui {

	private:
		Fl_Window *window;
		Fl_Tabs *tabs;

		Fl_Group *overlayGroup;
		Fl_Check_Button *ticksButton;
		Fl_Check_Button *dimensionsButton;
		Fl_Check_Button *boundingBoxButton;
		Fl_Check_Button *gridButton;
		Fl_Check_Button *umUnitsButton;
		Fl_Check_Button *nmUnitsButton;
		Slider *gridSpacingSlider;
		Fl_Check_Button *relativeSpacingButton;

	public:
		AnnotationsGui();

		int relativeSpacing() { return relativeSpacingButton->value(); }
		float getSpacing() { return gridSpacingSlider->value()/1000.0; }
		void setNanometres() {
			nmUnitsButton->value(1);
			umUnitsButton->value(0);
		}
		void setMicrometres() {
			umUnitsButton->value(1);
			nmUnitsButton->value(0);
		}
		void show() { window->show(); }
		void hide() { window->hide(); }

};

#endif /* ANNOTATIONS_H_ */
