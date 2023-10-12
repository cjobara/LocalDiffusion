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

#ifndef MOVIE_H_
#define MOVIE_H_

#define MAX_COLS 7
#define MAX_ROWS 21

#ifdef _WIN32
	#include <FL\Fl_Window.H>
	#include <FL\Fl_Box.H>
	#include <FL\Fl_Button.H>
	#include <FL\Fl_Float_Input.H>
	#include <FL\Fl_Choice.H>
	#include <FL\Fl_Scroll.H>
#elif __APPLE__
	#include <FL/Fl_Window.H>
	#include <FL/Fl_Box.H>
	#include <FL/Fl_Button.H>
	#include <FL/Fl_Float_Input.H>
	#include <FL/Fl_Choice.H>
	#include <FL/Fl_Scroll.H>
#endif

#include <FL/fl_message.H>
#include "gui.h"

class MovieMakerTable;
class MovieMakerGui;

void movieMakerAddRowCallback(Fl_Widget*, void*);
void movieMakerRemoveRowCallback(Fl_Widget*, void*);
void movieMakerRecordCallback(Fl_Widget*, void*);
void movieMakerPlayCallback(Fl_Widget*, void*);
void movieMakerStopCallback(Fl_Widget*, void*);
void movieMakerOperationOverlapCallback(Fl_Widget* w, void*);
void movieMakerOperation( int operation, int axis, float increment );
void movieMakerUpdateUnitsCallback(Fl_Widget* w, void*);
void movieMakerWindowCallback(Fl_Widget*, void*);
void movieMakerSaveScriptCallback(Fl_Widget*,void*);
void movieMakerLoadScriptCallback(Fl_Widget*,void*);

const Fl_Menu_Item movieMakerOperationItems[] = {
    {"Rotation"},
    {"Translation"},
    {"Zoom"},
    {"Alpha"},
    {"Scale"},
    {0},
};

const Fl_Menu_Item movieMakerAxisItems[] = {
    {"x"},
    {"y"},
    {"z"},
    {0},
};

class MovieMakerGui {

	public:
		MovieMakerGui();
		MovieMakerTable *table;

		void show() { movieMakerWindow->show(); }
		void hide() { movieMakerWindow->hide(); }

		Fl_Button* stopButton;

	private:

		Fl_Window* movieMakerWindow;
		Fl_Button* playButton;
		Fl_Button* recordButton;
		Fl_Button* addRowButton;
		Fl_Button* removeButton;
		Fl_Scroll* movieMakerScroll;

		Fl_Button* saveScriptButton;
		Fl_Button* loadScriptButton;

};

class MovieMakerTable : public Fl_Scroll {
    void *w[MAX_ROWS][MAX_COLS];        // widget pointers
    int xPassed;
    int yPassed;
    int rows;
    int cellh;
    int yy;

public:
    MovieMakerTable(int X, int Y, int W, int H, const char*L=0) : Fl_Scroll(X,Y,W,H,L) {

        rows = 1;
        xPassed = X;
        yPassed = Y;
        cellh = 20;
        yy = yPassed;

        int xx = xPassed;
        yy = Y;

        int cellWidth[MAX_COLS] = { 35,100,45,70,70,45,75 };

        xx = xPassed;
        yy = Y;

        begin();

        // Create widgets
        for ( int r=1; r<=rows; r++ ) {
            for ( int c=0; c<MAX_COLS; c++ ) {
            	if (c==0) {
                    Fl_Button *box = new Fl_Button(xx,yy,cellWidth[c],cellh);
                    box->label("1");
					box->labelsize(12);
            		box->align(FL_ALIGN_CENTER);
                    w[r][c] = (void*)box;
                } else {
                	if (c==1) {
                		Fl_Choice *operation = new Fl_Choice(xx,yy,cellWidth[c],cellh);
                		operation->menu(movieMakerOperationItems);
                		operation->box(FL_THIN_DOWN_BOX);
                		operation->value(0); // initially
                		operation->when(FL_WHEN_CHANGED);
						operation->labelsize(12);
						operation->textsize(12);
                		operation->callback((Fl_Callback*)movieMakerUpdateUnitsCallback,(void*)(0));
                        w[r][c] = (void*)operation;
                	}
                	else if (c==2) {
                		Fl_Choice *axis = new Fl_Choice(xx,yy,cellWidth[c],cellh);
                		axis->menu(movieMakerAxisItems);
                		axis->box(FL_THIN_DOWN_BOX);
                		axis->value(0); // initially
						axis->labelsize(12);
						axis->textsize(12);
                		axis->when(FL_WHEN_CHANGED);
                        w[r][c] = (void*)axis;
                	}
                	else if (c==3) {
                		Fl_Float_Input *magnitude = new Fl_Float_Input(xx,yy,cellWidth[c],cellh);
                		magnitude->box(FL_THIN_DOWN_BOX);
                		magnitude->value("0"); // initially
						magnitude->textsize(12);
						magnitude->labelsize(12);
                		magnitude->align(FL_ALIGN_CENTER);
                        w[r][c] = (void*)magnitude;
                	}
                	else if (c==4) {
                		Fl_Float_Input *increment = new Fl_Float_Input(xx,yy,cellWidth[c],cellh);
                		increment->box(FL_THIN_DOWN_BOX);
                		increment->value("0"); // initially
						increment->labelsize(12);
						increment->textsize(12);
                		increment->align(FL_ALIGN_CENTER);
                        w[r][c] = (void*)increment;
                	}
                	else if (c==5) {
                		Fl_Box *unit = new Fl_Box(xx,yy,cellWidth[c],cellh,"");
                		unit->box(FL_THIN_DOWN_BOX);
                		unit->align(FL_ALIGN_CENTER);
						unit->labelsize(12);
                		unit->label("deg"); // initially
                        w[r][c] = (void*)unit;
                	}
                	else if (c==6) {
                		Fl_Float_Input *overlap = new Fl_Float_Input(xx,yy,cellWidth[c],cellh);
                		overlap->box(FL_THIN_DOWN_BOX);
                		overlap->value(""); // initially
						overlap->labelsize(12);
						overlap->textsize(12);
						overlap->deactivate();
                		overlap->align(FL_ALIGN_CENTER);
                        w[r][c] = (void*)overlap;
                	}
                }
                xx += cellWidth[c];
            }
            xx = xPassed;
            yy += cellh;
        }
        end();
    }

    int getOperation(int row);
	void setOperation(int row, int value);
    int getAxis(int row);
	void setAxis(int row, int value);
    float getMagnitude(int row);
	void setMagnitude(int row, float value);
    float getIncrement(int row);
	void setIncrement(int row, float value);
	float getOverlap(int row);
	void setOverlap(int row, float value);
    void addRow();
    void removeRow();
    int getRows();
    void assignValues();
    Fl_Box* getUnitsPointer(Fl_Choice* op);
    Fl_Choice* getAxisPointer(Fl_Choice* op);
	Fl_Float_Input* getMagnitudePointer(Fl_Choice* op);
    void highlightRow(int row);
};


#endif /* MOVIE_H_ */
