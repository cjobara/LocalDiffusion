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

#ifdef _WIN32
	#include <FL\Fl_Double_Window.H>
	#include <FL\Fl_Native_File_Chooser.H>
	#include <FL\Fl_Check_Button.H>
	#include <FL\Fl_Menu_Item.H>
#elif __APPLE__
	#include <FL/Fl_Double_Window.H>
	#include <FL/Fl_Native_File_Chooser.H>
	#include <FL/Fl_Check_Button.H>
#endif

#include "movie.h"
#include "globals.h"

extern Globals *iMAP;
extern File *file;

MovieMakerGui::MovieMakerGui() {

    const int cellWidth[MAX_COLS] = { 35,100,45,70,70,45,75 };
    const char *header[MAX_COLS] = { "#","Operation", "Axis", "Magnitude", "Increment", "Units","% Overlap" };
    int xx = 10;

	movieMakerWindow = new Fl_Window(485,380, "Movie Maker");
	movieMakerWindow->begin();
	movieMakerWindow->callback(movieMakerWindowCallback);

	// draw header
	for ( int c=0; c<MAX_COLS; c++ ) {
		Fl_Box *box = new Fl_Box(xx,10,cellWidth[c],20,header[c]);
		box->box(FL_THIN_UP_BOX);
		box->labelsize(12);
		box->show();
		xx += cellWidth[c];
	}

	movieMakerScroll = new Fl_Scroll(10,30,465,280);
		table = new MovieMakerTable(10,30,440,280);
	movieMakerScroll->end();

	playButton = new Fl_Button(110,345,50,25,"Play");
	playButton->labelcolor(iMAP->bgColor);
	playButton->label("@>");
	playButton->labelsize(12);
	playButton->callback(movieMakerPlayCallback, (void*)this);

	stopButton = new Fl_Button(165,345,50,25,"Stop");
	stopButton->labelcolor(iMAP->bgColor);
	stopButton->label("@square");
	stopButton->labelsize(12);
	stopButton->callback(movieMakerStopCallback, (void*)this);

	recordButton = new Fl_Button(220,345,50,25,"Record");
	recordButton->labelcolor(FL_RED);
	recordButton->labelsize(12);
	recordButton->label("@circle");
	recordButton->callback(movieMakerRecordCallback, (void*)this);

	addRowButton = new Fl_Button(87,315,100,25,"Add Row");
	addRowButton->labelsize(12);
	addRowButton->callback(movieMakerAddRowCallback, (void*)this);

	removeButton = new Fl_Button(193,315,100,25,"Remove Row");
	removeButton->labelsize(12);
	removeButton->callback(movieMakerRemoveRowCallback, (void*)this);

	saveScriptButton = new Fl_Button(310,315,100,25,"Save Script");
	saveScriptButton->labelsize(12);
	saveScriptButton->callback(movieMakerSaveScriptCallback,(void*)this);

	loadScriptButton = new Fl_Button(310,345,100,25,"Load Script");
	loadScriptButton->labelsize(12);
	loadScriptButton->callback(movieMakerLoadScriptCallback,(void*)this);

	movieMakerWindow->end();
	movieMakerWindow->show();
}

int MovieMakerTable::getOperation(int row) {
	Fl_Choice *operationPointer = (Fl_Choice*)w[row][1];
	return operationPointer->value();
}

int MovieMakerTable::getAxis(int row) {
	Fl_Choice *axisPointer = (Fl_Choice*)w[row][2];
	return axisPointer->value();
}

float MovieMakerTable::getMagnitude(int row) {
	Fl_Float_Input *magnitudePointer = (Fl_Float_Input*)w[row][3];
	return atof(magnitudePointer->value());
}

float MovieMakerTable::getIncrement(int row) {
	Fl_Float_Input *incrementPointer = (Fl_Float_Input*)w[row][4];
	return atof(incrementPointer->value());
}

float MovieMakerTable::getOverlap(int row) {
	Fl_Float_Input *overlapPointer = (Fl_Float_Input*)w[row][6];
	return atof(overlapPointer->value());
}

void MovieMakerTable::addRow() {

    const int cellWidth[MAX_COLS] = { 35,100,45,70,70,45,75 };
    int xx = xPassed;
    rows++;

    char rowNumber[20];
    sprintf(rowNumber," %i",rows);

    for ( int r=rows; r<=rows; r++ ) {
        for ( int c=0; c<MAX_COLS; c++ ) {
        	if (c==0) {
                Fl_Button *box = new Fl_Button(xx,yy,cellWidth[c],cellh);
                box->copy_label(rowNumber);
				box->labelsize(12);
                w[r][c] = (void*)box;
            } else {
            	if (c==1) {
            		Fl_Choice *operation = new Fl_Choice(xx,yy,cellWidth[c],cellh);
            		operation->menu(movieMakerOperationItems);
            		operation->box(FL_THIN_DOWN_BOX);
            		operation->value(0); // initially
					operation->labelsize(12);
					operation->textsize(12);
            		operation->when(FL_WHEN_CHANGED);
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
					magnitude->labelsize(12);
					magnitude->textsize(12);
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
            		overlap->value("0"); // initially
					overlap->labelsize(12);
					overlap->textsize(12);
            		overlap->align(FL_ALIGN_CENTER);
                    w[r][c] = (void*)overlap;
            	}
            }
            xx += cellWidth[c];
        }
        xx = xPassed;
        yy += cellh;
        MovieMakerTable::redraw();
    }
    end();
}

void MovieMakerTable::removeRow() {
	if (rows > 1) {
		if (yy > 60) { yy-=cellh; }

		Fl_Button *rowNumber = (Fl_Button*)w[rows][0];
		Fl_Choice *operationPointer = (Fl_Choice*)w[rows][1];
		Fl_Choice *axisPointer = (Fl_Choice*)w[rows][2];
		Fl_Float_Input *magnitudePointer = (Fl_Float_Input*)w[rows][3];
		Fl_Float_Input *incrementPointer = (Fl_Float_Input*)w[rows][4];
		Fl_Box *unitPointer = (Fl_Box*)w[rows][5];
		Fl_Float_Input *overlapPointer = (Fl_Float_Input*)w[rows][6];

		rowNumber->hide();
		operationPointer->hide();
		axisPointer->hide();
		magnitudePointer->hide();
		unitPointer->hide();
		incrementPointer->hide();
		overlapPointer->hide();

		rows--;
	}
}

int MovieMakerTable::getRows() {
	return rows;
}

// must be run before playing/recording movie
void MovieMakerTable::assignValues() {
	for (int r = 1; r <= rows; r++) {
		Fl_Choice *operationPointer = (Fl_Choice*)w[r][1];
		Fl_Choice *axisPointer = (Fl_Choice*)w[r][2];
		Fl_Float_Input *magnitudePointer = (Fl_Float_Input*)w[r][3];
		Fl_Float_Input *incrementPointer = (Fl_Float_Input*)w[r][4];
		Fl_Float_Input *overlapPointer = (Fl_Float_Input*)w[r][6];

		iMAP->operationArray[r] = operationPointer->value();
		iMAP->axisArray[r] = axisPointer->value();
		iMAP->magnitudeArray[r] = atof(magnitudePointer->value());
		iMAP->incrementArray[r] = atof(incrementPointer->value());
		iMAP->overlapArray[r] = atof(overlapPointer->value())/100;

		// error proofing
		if (atof(magnitudePointer->value()) == 0) {
			magnitudePointer->value("1");
			iMAP->magnitudeArray[r] = atof(magnitudePointer->value());
		}
		else if (atof(magnitudePointer->value()) < 0) {
			char mag [30];
			sprintf(mag,"%f",fabsf(atof(magnitudePointer->value())));
			magnitudePointer->value(mag);
			iMAP->magnitudeArray[r] = fabsf(atof(magnitudePointer->value()));
		}
		if (atof(incrementPointer->value()) == 0) {
			incrementPointer->value("1");
			iMAP->incrementArray[r] = atof(incrementPointer->value());
		}
		if (atof(overlapPointer->value()) > 100) {
			overlapPointer->value("100");
			iMAP->overlapArray[r] = atof(overlapPointer->value());
		}

		// calculate the number of frames for each step
		switch (iMAP->operationArray[r]) {
			case 0: // case of rotation
				iMAP->framesArray[r] = (int)(iMAP->magnitudeArray[r]/fabsf(iMAP->incrementArray[r]));
				break;
			case 1: // case for translation
				iMAP->framesArray[r] = (int)(iMAP->magnitudeArray[r]/fabsf(iMAP->incrementArray[r]));
				break;
			case 2: // case for zoom
				iMAP->framesArray[r] = (int)(iMAP->magnitudeArray[r]/fabsf(iMAP->incrementArray[r]));
				break;
			case 3: // case for alpha
				iMAP->framesArray[r] = (int)(iMAP->magnitudeArray[r]/fabsf(iMAP->incrementArray[r]));
				break;
			case 4: // case for scale
				iMAP->framesArray[r] = (int)(iMAP->magnitudeArray[r]/fabsf(iMAP->incrementArray[r]));
				break;
		}
	}
}

void movieMakerOperation( int operation, int axis, float increment ) {
	switch (operation) {
	case 0: // rotation
		switch (axis) {
		case 0: // x axis
			file->xRotate += increment;
			break;
		case 1: // y axis
			file->yRotate += increment;
			break;
		case 2:
			file->zRotate += increment;
			break;
		}
		break;
	case 1: // translation
		switch (axis) {
		case 0: // x axis
			file->xTranslate -= increment;
			break;
		case 1: // y axis
			file->yTranslate -= increment;
			break;
		case 2:
			file->zTranslate += increment;
			break;
		}
		break;
	case 2: // zoom
		iMAP->fovy += increment;
		break;
	case 3: // alpha
		if (file->squareMeshOverlay == false && file->voronoiMeshOverlay == false && file->treeMeshOverlay == false) {} else {
			iMAP->landscapeAlpha += increment;
			switch(file->meshType) {
				case 0: // regular mesh
					file->squareMeshGui->landscapeAlphaSlider->value(iMAP->landscapeAlpha);
					break;
				case 1: // voronoi mesh
					file->voronoiMeshGui->landscapeAlphaSlider->value(iMAP->landscapeAlpha);
					break;
				case 2: // quad-tree
					file->treeMeshGui->landscapeAlphaSlider->value(iMAP->landscapeAlpha);
					break;
			}
		}
		break;
	case 4: // scale
		if (file->squareMeshOverlay == false && file->voronoiMeshOverlay == false && file->treeMeshOverlay == false) {} else {
			iMAP->landscapeScale += increment;
			switch(file->meshType) {
				case 0: // regular mesh
					file->squareMeshGui->scaleLandscapeSlider->value(iMAP->landscapeScale);
					break;
				case 1: // voronoi mesh
					file->voronoiMeshGui->scaleLandscapeSlider->value(iMAP->landscapeScale);
					break;
				case 2: // quad-tree
					file->treeMeshGui->scaleLandscapeSlider->value(iMAP->landscapeScale);
					break;
			}
		}
		break;
	}
}

// return units pointer for row of the passed operation (op)
Fl_Box* MovieMakerTable::getUnitsPointer(Fl_Choice* op) {
	for (int q = 1; q<=rows; q++) {
		Fl_Choice* opCompare = (Fl_Choice*)w[q][1];
		if (op == opCompare) {
			Fl_Box* unitsReturn = (Fl_Box*)w[q][5];
			return unitsReturn;
		}
	}
	return 0;
}

// return axis pointer for row of the passed operation (op)
Fl_Choice* MovieMakerTable::getAxisPointer(Fl_Choice* op) {
	for (int q = 1; q<=rows; q++) {
		Fl_Choice* opCompare = (Fl_Choice*)w[q][1];
		if (op == opCompare) {
			Fl_Choice* axisReturn = (Fl_Choice*)w[q][2];
			return axisReturn;
		}
	}
	return 0;
}

// retun magnitude pointer for row of the passed operation (op)
Fl_Float_Input* MovieMakerTable::getMagnitudePointer(Fl_Choice* op) {
	for (int q = 1; q<=rows; q++) {
		Fl_Choice* opCompare = (Fl_Choice*)w[q][1];
		if (op == opCompare) {
			Fl_Float_Input* axisReturn = (Fl_Float_Input*)w[q][3];
			return axisReturn;
		}
	}
	return 0;
}

// highlight row in movieMakerTable
void MovieMakerTable::highlightRow(int row) {
	// get current row
	Fl_Input* current = (Fl_Input*)w[row][0];
	if (row != 0) {
		current->color(FL_CYAN);
	}
	for (int q = 1; q <= rows; q++) {
		if (q != row) {
			Fl_Input* inactiveRow = (Fl_Input*)w[q][0];
			inactiveRow->color(FL_BACKGROUND_COLOR);
		}
	}
}

// add row to movie maker script
void movieMakerAddRowCallback(Fl_Widget*, void*) {
	if (iMAP->movieMakerGui->table->getRows() < MAX_ROWS-1) {
		// running addRow() once doesn't update display for whatever reason
		iMAP->movieMakerGui->table->addRow();
		iMAP->movieMakerGui->table->removeRow();
		iMAP->movieMakerGui->table->addRow();
	}
}

// remove row from movie maker script
void movieMakerRemoveRowCallback(Fl_Widget*, void*) {
	iMAP->movieMakerGui->table->removeRow();
}

void movieMakerPlayCallback(Fl_Widget*, void*) {
	iMAP->movieMakerGui->table->assignValues();
	iMAP->movieMakerPlayEnable = true;
	float alpha, scale;
	if (file->squareMeshOverlay == false && file->voronoiMeshOverlay == false && file->treeMeshOverlay == false) {
		alpha = 1.0;
		scale = 1.0;
	} else {
		switch(file->meshType) {
			case 0:
				alpha = file->squareMeshGui->landscapeAlphaSlider->value();
				scale = file->squareMeshGui->scaleLandscapeSlider->value();
				break;
			case 1:
				alpha = file->voronoiMeshGui->landscapeAlphaSlider->value();
				scale = file->voronoiMeshGui->scaleLandscapeSlider->value();
				break;
			case 2:
				alpha = file->treeMeshGui->landscapeAlphaSlider->value();
				scale = file->treeMeshGui->scaleLandscapeSlider->value();
				break;
			default:
				alpha = 1.0;
				scale = 1.0;
				break;
			}
	}
	iMAP->landscapeAlpha = alpha;
	iMAP->landscapeScale = scale;
}

void movieMakerStopCallback(Fl_Widget*, void*) {
	iMAP->movieMakerPlayEnable = false;
	iMAP->movieMakerRecordEnable = false;
	for (int p = 0; p < MAX_COLS; p++) { iMAP->framesPlayedArray[iMAP->currentRow] = 0; }
	iMAP->currentRow = 1;
}

void movieMakerRecordCallback(Fl_Widget*, void*) {
	iMAP->movieMakerGui->table->assignValues();
	iMAP->tiffCount = 1;

	Fl_Native_File_Chooser native;
	native.title("Save sequence as");
	native.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native.filter("TIFF Files\t*.tif\n");
	native.show();

	char nativeFilename[FILENAME_MAX];
	strcpy(nativeFilename,native.filename());

	// in case save screen is cancelled
	if (nativeFilename[0] != (char)NULL) {
		strcpy(iMAP->captureScreenName, native.filename());

		float alpha, scale;
		if (file->squareMeshOverlay == false && file->voronoiMeshOverlay == false && file->treeMeshOverlay == false) {
			alpha = 1.0;
			scale = 1.0;
		} else {
		switch(file->meshType) {
			case 0:
				alpha = file->squareMeshGui->landscapeAlphaSlider->value();
				scale = file->squareMeshGui->scaleLandscapeSlider->value();
				break;
			case 1:
				alpha = file->voronoiMeshGui->landscapeAlphaSlider->value();
				scale = file->voronoiMeshGui->scaleLandscapeSlider->value();
				break;
			case 2:
				alpha = file->treeMeshGui->landscapeAlphaSlider->value();
				scale = file->treeMeshGui->scaleLandscapeSlider->value();
				break;
			default:
				alpha = 1.0;
				scale = 1.0;
				break;
			}
		}
		iMAP->landscapeAlpha = alpha;
		iMAP->landscapeScale = scale;

		iMAP->movieMakerPlayEnable = true;
		iMAP->movieMakerRecordEnable = true;
	}
}

// update units cell callback
void movieMakerUpdateUnitsCallback(Fl_Widget* w, void*) {
	Fl_Choice* operationPointer = (Fl_Choice*)w;
	Fl_Box* units = iMAP->movieMakerGui->table->getUnitsPointer(operationPointer);
	Fl_Choice* axis = iMAP->movieMakerGui->table->getAxisPointer(operationPointer);
	Fl_Float_Input* magnitude = iMAP->movieMakerGui->table->getMagnitudePointer(operationPointer);

	char mag [30];
	sprintf(mag,"%f",200);

	switch (operationPointer->value()) {
	case 0: // rotation
		units->label("deg");
		axis->activate();
		magnitude->activate();
		break;
	case 1: // translation
		units->label("um");
		axis->activate();
		magnitude->activate();
		break;
	case 2: // zoom
		units->label("deg");
		axis->deactivate();
		magnitude->activate();
		break;
	case 3: // alpha
		units->label("%");
		axis->deactivate();
		magnitude->activate();
		break;
	case 4: // scale
		units->label("%");
		axis->deactivate();
		magnitude->activate();
		break;
	}
}

void MovieMakerTable::setOperation(int row, int value) {
	Fl_Choice *operationPointer = (Fl_Choice*)w[row][1];
	operationPointer->value(value);

	// for zoom
	if (value == 2 || value == 3 || value == 4) {
		Fl_Choice *axisPointer = (Fl_Choice*)w[row][2];
		axisPointer->deactivate();
	}
}

void MovieMakerTable::setAxis(int row, int value) {
	Fl_Choice *axisPointer = (Fl_Choice*)w[row][2];
	axisPointer->value(value);
}

void MovieMakerTable::setMagnitude(int row, float value) {
	Fl_Float_Input *magnitudePointer = (Fl_Float_Input*)w[row][3];
	char s [30];
	sprintf(s,"%.2f",value);
	magnitudePointer->value(s);
}

void MovieMakerTable::setIncrement(int row, float value) {
	Fl_Float_Input *incrementPointer = (Fl_Float_Input*)w[row][4];
	char s [30];
	sprintf(s,"%.2f",value);
	incrementPointer->value(s);
}

void MovieMakerTable::setOverlap(int row, float value) {
	Fl_Float_Input *overlapPointer = (Fl_Float_Input*)w[row][6];
	char s [30];
	sprintf(s,"%.0f",value);
	overlapPointer->value(s);
}


void movieMakerWindowCallback(Fl_Widget*, void*) {
	iMAP->movieMakerGui->hide();
	iMAP->movieMakerOpened = false;
}

void movieMakerSaveScriptCallback(Fl_Widget*,void*) {
	Fl_Native_File_Chooser native;
	char nativeFilename[FILENAME_MAX];

	native.title("Save Script As");
	native.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
	native.filter("Script Files\t*.script\n");
	native.show();
	strcpy(nativeFilename,native.filename());
	strcat(nativeFilename,".script");

	// load file conditions (to avoid closing file if open screen is cancelled)
	if (nativeFilename[0]!=(char)NULL) {
		FILE * writeFile;

		writeFile = fopen (nativeFilename,"wb");

		for (int b = 1; b <= iMAP->movieMakerGui->table->getRows(); b++) {
			fprintf(writeFile, "%i\t%i\t%f\t%f\t%f\n",
					iMAP->movieMakerGui->table->getOperation(b),
					iMAP->movieMakerGui->table->getAxis(b),
					iMAP->movieMakerGui->table->getMagnitude(b),
					iMAP->movieMakerGui->table->getIncrement(b),
					iMAP->movieMakerGui->table->getOverlap(b)
					);
		}

		fclose (writeFile);
	}
}

void movieMakerLoadScriptCallback(Fl_Widget*,void*) {
	Fl_Native_File_Chooser native;
	char nativeFilename[FILENAME_MAX];
	native.title("Load Script File");
	native.type(Fl_Native_File_Chooser::BROWSE_FILE);
	native.filter("Script Files\t*.script\n");
	native.show();
	strcpy(nativeFilename,native.filename());

	// load file conditions (to avoid closing file if open screen is cancelled)
	if ((nativeFilename[0]!=(char)NULL) || (nativeFilename[0]!=(char)NULL && iMAP->fileNameStorage[iMAP->fileNumber][0]!=(char)NULL)) {
		// clear any existing rows
		int r = iMAP->movieMakerGui->table->getRows();
		for (int l = 0; l < r; l++) { iMAP->movieMakerGui->table->removeRow(); }

		// read file
		FILE *fileHandle = fopen(nativeFilename,"r");

		int scriptRows = 0;
		// determine number of lines in file
		int ch;
		do {
			ch = fgetc(fileHandle);
			if (ch=='\n') {
				scriptRows++;
			}
		} while( ch != EOF);

		rewind(fileHandle);

		int operation = 0;
		int axis = 0;
		float magnitude = 0;
		float increment = 0;
		float overlap = 0;

		// add row to table
		for (int e = 1; e < scriptRows; e++) {
			iMAP->movieMakerGui->table->addRow();
			iMAP->movieMakerGui->table->removeRow();
			iMAP->movieMakerGui->table->addRow();
		}

		for(int j = 1; j <= scriptRows; j++) {
			operation = 0;
			axis = 0;
			magnitude = 0;
			increment = 0;
			overlap = 0;
			// assure no more than MAX_ROWS
			if (j < MAX_ROWS-1) {

				fscanf(fileHandle,"%i\t%i\t%f\t%f\t%f\t",
						&operation,
						&axis,
						&magnitude,
						&increment,
						&overlap
						);

				// assign read values
				iMAP->movieMakerGui->table->setOperation(j,operation);
				iMAP->movieMakerGui->table->setAxis(j,axis);
				iMAP->movieMakerGui->table->setMagnitude(j,magnitude);
				iMAP->movieMakerGui->table->setIncrement(j,increment);
				iMAP->movieMakerGui->table->setOverlap(j,overlap);
			}
		}
		fclose(fileHandle);
	}
}
