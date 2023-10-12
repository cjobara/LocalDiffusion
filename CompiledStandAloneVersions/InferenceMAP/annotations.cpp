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

#include "file.h"
#include "globals.h"
#include "annotations.h"
#include "draw.h"
#include "text.h"

extern File *file;
extern Globals *iMAP;

AnnotationsGui::AnnotationsGui() {

	if (iMAP->fileNumber > 0) {

		const int w = 250;
		const int h = 215;

		window = new Fl_Window(w,h,"Annotations");
		window->callback(annotationsWindowCallback);
		window->color(iMAP->bgColor);
		window->begin();
		window->set_non_modal();

		tabs = new Fl_Tabs(0,0,w,h);
		tabs->labelsize(12);
//		tabs->color(FL_DARK1,FL_DARK1);
		tabs->box(FL_BORDER_FRAME);
		tabs->labelcolor(FL_WHITE);
		tabs->labelfont(iMAP->normalFont);
		tabs->begin();

        char length [10];
        sprintf(length,"%i",iMAP->glWindow->w());

		//displayGroup = new Fl_Group(0,25,w,h,"Display");
		//displayGroup->labelsize(12);
		//displayGroup->labelfont(1);
		//displayGroup->labelcolor(FL_DARK1);
		//displayGroup->labelfont(iMAP->normalFont);
		//displayGroup->box(FL_BORDER_FRAME);
		//displayGroup->begin();
		//{
		//	screenLengthInput = new Fl_Int_Input(75,35,80,20,"Size [px]");
		//	screenLengthInput->labelcolor(FL_WHITE);
		//	screenLengthInput->color(iMAP->bgColor);
		//	screenLengthInput->labelfont(iMAP->boldFont);
		//	screenLengthInput->textfont(iMAP->boldFont);
		//	screenLengthInput->textcolor(FL_WHITE);
		//	screenLengthInput->labelsize(12);
		//	screenLengthInput->textsize(12);
		//	screenLengthInput->value(length);
		//	screenLengthInput->show();

		//	updateDisplayButton = new Fl_Button(150+10,35,w-150-20,20,"Update");
		//	updateDisplayButton->labelsize(12);
		//	updateDisplayButton->labelfont(iMAP->boldFont);
		//	updateDisplayButton->callback((Fl_Callback*)updateDisplaySizeCallback,(int*)atoi(screenLengthInput->value()));
		//	updateDisplayButton->show();

		//	backgroundColorButton = new Fl_Button(10,65,w-20,25,"Background Colour");
		//	backgroundColorButton->labelsize(12);
		//	backgroundColorButton->labelfont(iMAP->boldFont);
		//	backgroundColorButton->callback((Fl_Callback*)backgroundColorCallback);
		//	backgroundColorButton->show();

		//}
		//displayGroup->end();

		overlayGroup = new Fl_Group(0,25,w,h,"Overlay");
		overlayGroup->labelsize(12);
		overlayGroup->labelcolor(FL_DARK1);
		overlayGroup->box(FL_BORDER_FRAME);
		overlayGroup->labelfont(iMAP->normalFont);
		overlayGroup->begin();
		{
			boundingBoxButton = new Fl_Check_Button(10,35,w-20,20,"Bounding Box");
			boundingBoxButton->labelsize(12);
			boundingBoxButton->labelcolor(FL_WHITE);
			boundingBoxButton->labelfont(iMAP->boldFont);
			boundingBoxButton->value(1);
			boundingBoxButton->callback((Fl_Callback*)boundingBoxCallback);
			boundingBoxButton->show();

			gridButton = new Fl_Check_Button(10,60,w-20,20,"Grid");
			gridButton->labelsize(12);
			gridButton->labelcolor(FL_WHITE);
			gridButton->labelfont(iMAP->boldFont);
			gridButton->value(0);
			gridButton->callback((Fl_Callback*)gridCallback);
			gridButton->show();

			dimensionsButton = new Fl_Check_Button(10,85,w-20,20,"Dimensions");
			dimensionsButton->labelsize(12);
			dimensionsButton->labelcolor(FL_WHITE);
			dimensionsButton->labelfont(iMAP->boldFont);
			dimensionsButton->value(1);
			dimensionsButton->callback((Fl_Callback*)dimensionsCallback);
			dimensionsButton->show();

			ticksButton = new Fl_Check_Button(10,110,w-20,20,"Ticks");
			ticksButton->labelsize(12);
			ticksButton->labelcolor(FL_WHITE);
			ticksButton->labelfont(iMAP->boldFont);
			ticksButton->value(1);
			ticksButton->callback((Fl_Callback*)ticksCallback);
			ticksButton->show();

			umUnitsButton = new Fl_Check_Button(10,135,w/2-40,20,"um Units");
			umUnitsButton->labelsize(12);
			umUnitsButton->labelcolor(FL_WHITE);
			umUnitsButton->labelfont(iMAP->boldFont);
			umUnitsButton->value(1);
			umUnitsButton->callback((Fl_Callback*)umUnitsButtonCallback);
			umUnitsButton->show();

			nmUnitsButton = new Fl_Check_Button(w/2-10,135,w/2-40,20,"nm Units");
			nmUnitsButton->labelsize(12);
			nmUnitsButton->labelcolor(FL_WHITE);
			nmUnitsButton->labelfont(iMAP->boldFont);
			nmUnitsButton->value(0);
			nmUnitsButton->callback((Fl_Callback*)nmUnitsButtonCallback);
			nmUnitsButton->show();

			relativeSpacingButton = new Fl_Check_Button(10,160,w-20,20,"Grid Spacing [nm]");
			relativeSpacingButton->labelsize(12);
			relativeSpacingButton->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
			relativeSpacingButton->labelcolor(FL_WHITE);
			relativeSpacingButton->labelfont(iMAP->boldFont);
			relativeSpacingButton->value(0);
			relativeSpacingButton->show();

			gridSpacingSlider = new Slider(10,185,w-20,20);
			gridSpacingSlider->labelsize(12);
			gridSpacingSlider->labelcolor(FL_WHITE);
			gridSpacingSlider->labelfont(iMAP->boldFont);
			gridSpacingSlider->value(250);
			gridSpacingSlider->bounds(100,1000);
			gridSpacingSlider->show();

		}
		overlayGroup->end();
		tabs->end();

		window->end();
		window->show();
	}
}

void backgroundColorCallback(Fl_Button*w,void*) {
	fl_color_chooser("Background Colour",iMAP->backgroundRGB[0],iMAP->backgroundRGB[1],iMAP->backgroundRGB[2],0);
}

void annotationsWindowCallback(Fl_Widget*w,void*) {
	iMAP->annotationsGui->hide();
}

void boundingBoxCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (file->boundingBoxEnable) { file->boundingBoxEnable = false;	}
		else { file->boundingBoxEnable = true; }
	}
}

void dimensionsCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (file->dimensionsEnable) { file->dimensionsEnable = false;	}
		else { file->dimensionsEnable = true; }
	}
}

void gridCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (file->gridEnable) { file->gridEnable = false;	}
		else { file->gridEnable = true; }
	}
}

void ticksCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (file->ticksEnable) { file->ticksEnable = false;	}
		else { file->ticksEnable = true; }
	}
}

void updateDisplaySizeCallback(Fl_Button*b,void*v) {

	char lengthString [10];
	sprintf(lengthString,"%i",v);
	int length = atoi(lengthString);

	int w = iMAP->screenDimensions[2]-90;
	int h = iMAP->screenDimensions[3]-90;

	const float ar = 730.0f/1110.0f;

	if (w > 1110 && h > 730) {

//		if (w-420 > h) { w = h-50+420; }
//		else { h = w-420; }

		if ((float)w*ar > (float)h) {
			w = h+380;
		} else {
			h = w-380;
		}

		iMAP->mainGroup = new Fl_Group(0,0,w,h-iMAP->macOffset);
		iMAP->mainWindow->resize(0,30-iMAP->macOffset,w,h-iMAP->macOffset);

		#ifdef _WIN32
			// create menubar
			iMAP->menubar->resize(0,0,w,30);
		#endif

		// create gl window
		iMAP->glWindow->resize(0,50-iMAP->macOffset,h-50,h-50);
		iMAP->fileInfo->resize(h-40,40-iMAP->macOffset,410,65);

		// visualization widgets
		iMAP->iOffsetSlider->resize(h-40, 130-iMAP->macOffset, 200, 20);
		iMAP->cMinSlider->resize(h-40, 170-iMAP->macOffset, 200, 20);
		iMAP->cMaxSlider->resize(h-40, 210-iMAP->macOffset, 200, 20);
		iMAP->drawTrajectoriesButton->resize(h+170,115-iMAP->macOffset,200,20);
		iMAP->localizationsButton->resize(h+170,135-iMAP->macOffset,200,20);
		iMAP->animateTrajectoriesButton->resize(h+170,155-iMAP->macOffset,200,20);
		iMAP->fullScreenButton->resize(w-35,h-35-iMAP->macOffset,25,25);
		iMAP->colorbar->resize(h-40,h-30-iMAP->macOffset,50,20);

		iMAP->tabBarWidth = h-50;
		// resize widths of any existing tabs
		for (int d = 0; d < iMAP->fileNumber; d++) {
			iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber),20);
		}
	}
	else {
		(b->value(0));
	}
	if (!b->value()) {
		iMAP->mainGroup->resize(0,0,1110,730-iMAP->macOffset);
		iMAP->mainWindow->resize(0,0,1110,730-iMAP->macOffset);

		#ifdef _WIN32
			// create menubar
			iMAP->menubar->resize(0,0,w,30);
		#endif
		// create gl window
		iMAP->glWindow->resize(0,50-iMAP->macOffset,680,680);
		iMAP->fileInfo->resize(690,40-iMAP->macOffset,410,65);

		// visualization widgets
		iMAP->iOffsetSlider->resize(690, 130-iMAP->macOffset, 200, 20);
		iMAP->cMinSlider->resize(690, 170-iMAP->macOffset, 200, 20);
		iMAP->cMaxSlider->resize(690, 210-iMAP->macOffset, 200, 20);
		iMAP->drawTrajectoriesButton->resize(900,115-iMAP->macOffset,200,20);
		iMAP->localizationsButton->resize(900,135-iMAP->macOffset,200,20);
		iMAP->animateTrajectoriesButton->resize(900,155-iMAP->macOffset,200,20);
		iMAP->fullScreenButton->resize(1075,695-iMAP->macOffset,25,25);
		iMAP->colorbar->resize(690,700-iMAP->macOffset,50,20);
		iMAP->screenDimensions[2] = 1110;
		iMAP->screenDimensions[3] = 730;

		iMAP->tabBarWidth = 680;
		// resize widths of any existing tabs
		for (int d = 0; d < iMAP->fileNumber; d++) {
			iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber),20);
		}
	}

	Fl::redraw();
}

/********** 2D ANNOTATIONS ***********/
void drawBoundingBox() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);

	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
		updateView();
		glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f-iMAP->backgroundRGB[0]);
		glLineWidth(3.0f);
		glBegin(GL_LINE_LOOP);
			glVertex2f(file->xMin,file->yMax);
			glVertex2f(file->xMax,file->yMax);
			glVertex2f(file->xMax,file->yMin);
			glVertex2f(file->xMin,file->yMin);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex2f(file->xMin,file->yMax);
			glVertex2f(file->xMax,file->yMax);
			glVertex2f(file->xMax,file->yMin);
			glVertex2f(file->xMin,file->yMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex2f(file->xMin,file->yMax);
			glVertex2f(file->xMin,file->yMax);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex2f(file->xMin,file->yMin);
			glVertex2f(file->xMin,file->yMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex2f(file->xMax,file->yMax);
			glVertex2f(file->xMax,file->yMax);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex2f(file->xMax,file->yMin);
			glVertex2f(file->xMax,file->yMin);
		glEnd();
	glPopMatrix();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
}

void drawDimensions() {

	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_TEXTURE_2D);

	const float xGridRange = file->xRange;
	const float yGridRange = file->yRange;

	float gridSpacing = 0.1;

	if (iMAP->annotationsGui != NULL) {
		if (iMAP->annotationsGui->relativeSpacing()) {
			gridSpacing = iMAP->annotationsGui->getSpacing();
		}
		else {
			gridSpacing = (file->maxRange)/15.0f;

			gridSpacing *= 1000.0;

			const int digits = countDigits((int)(gridSpacing))-1;
			const float pwr = pow((float)10.0,(float)digits);

			gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
			gridSpacing /= 1000.0;
		}
	}
	else {
		gridSpacing = (file->maxRange)/15.0f;

		gridSpacing *= 1000.0;

		const int digits = countDigits((int)(gridSpacing))-1;
		const float pwr = pow((float)10.0,(float)digits);

		gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
		gridSpacing /= 1000.0;
	}

	const float scale = gridSpacing/100.0f;

	float offset;
	if (iMAP->units == 0) {	offset = gridSpacing*1.4f; }
	else { offset = gridSpacing*1.6; }
	char label[20];

	txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

	const bool relative = true;

	glPushMatrix();
		updateView();

		glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);

		// x
		for (int j = 0; j <= (int)(xGridRange/gridSpacing); j++) {
			if (j%2 == 0) {
				if (iMAP->units) {
					if (relative) {	sprintf(label,"%6.0f",(float)j*gridSpacing*1000.0f); }
					else { sprintf(label,"%6.0f",(file->xMin+(float)j*gridSpacing)*1000.0f); }
				}
				else {
					if (relative) {	sprintf(label,"%6.1f",(float)j*gridSpacing); }
					else { sprintf(label,"%6.1f",file->xMin+(float)j*gridSpacing); }
				}
				glPushMatrix();
					glTranslatef(file->xMin+(float)j*gridSpacing,file->yMin-offset,0.0);
					glRotatef(90.0f,0.0f,0.0f,1.0f);
					glScalef(scale, scale, scale);
					txfRenderString(iMAP->fontTex, label, strlen(label));
				glPopMatrix();
			}
		}
		// y
		for (int j = 0; j <= (int)(yGridRange/gridSpacing); j++) {
			if (j%2 == 0) {
				if (iMAP->units) {
					if (relative) {	sprintf(label,"%6.0f",(float)j*gridSpacing*1000.0f); }
					else { sprintf(label,"%6.0f",(file->yMin+(float)j*gridSpacing)*1000.0f); }
				}
				else {
					if (relative) {	sprintf(label,"%6.1f",(float)j*gridSpacing); }
					else { sprintf(label,"%6.1f",file->yMin+(float)j*gridSpacing); }
				}
				glPushMatrix();
					glTranslatef(file->xMin-offset,file->yMin+(float)j*gridSpacing,0.0);
					glScalef(scale, scale, scale);
					txfRenderString(iMAP->fontTex, label, strlen(label));
				glPopMatrix();
			}
		}

	glPopMatrix();

	glBindTexture(GL_TEXTURE_2D,0);
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);

	glDisable(GL_DEPTH_TEST);

}

void drawTicks() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float xGridRange = file->xMax - file->xMin;
	float yGridRange = file->yMax - file->yMin;

	float gridSpacing = 0.1;

	if (iMAP->annotationsGui != NULL) {
		if (iMAP->annotationsGui->relativeSpacing()) {
			gridSpacing = iMAP->annotationsGui->getSpacing();
		}
		else {
			gridSpacing = (file->maxRange)/15.0f;

			gridSpacing *= 1000.0;

			const int digits = countDigits((int)(gridSpacing))-1;
			const float pwr = pow((float)10.0,(float)digits);

			gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
			gridSpacing /= 1000.0;
		}
	}
	else {
		gridSpacing = (file->maxRange)/15.0f;

		gridSpacing *= 1000.0;

		const int digits = countDigits((int)(gridSpacing))-1;
		const float pwr = pow((float)10.0,(float)digits);

		gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
		gridSpacing /= 1000.0;
	}


	const float tickSize = gridSpacing/8.0;

	glPushMatrix();
		updateView();

		glLineWidth(1);

		glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);

		// vertical segments
		for (int j = 1; j <= (int)(xGridRange/gridSpacing); j++) {
			glBegin(GL_TRIANGLES);
					// TICKS
					// bottom back
				glVertex2f(file->xMin+(float)j*gridSpacing-tickSize/2.0,file->yMin);
				glVertex2f(file->xMin+(float)j*gridSpacing,file->yMin+tickSize);
				glVertex2f(file->xMin+(float)j*gridSpacing+tickSize/2.0,file->yMin);
//					// top back
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax,file->zMin);
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax-tickSize,file->zMin);
//					// top z front
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax,file->zMin);
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax,file->zMin+tickSize);
//					// top z back
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax,file->zMax);
//					glVertex3f(file->xMin+(float)j*gridSpacing,file->yMax,file->zMax-tickSize);
			glEnd();
		}
		// horizontal segments
		for (int j = 1; j <= (int)(yGridRange/gridSpacing); j++) {
			glBegin(GL_TRIANGLES);
				// TICKS
				// bottom back
				glVertex2f(file->xMin,file->yMin+(float)j*gridSpacing-tickSize/2.0);
				glVertex2f(file->xMin+tickSize,	   file->yMin+(float)j*gridSpacing);
				glVertex2f(file->xMin,file->yMin+(float)j*gridSpacing+tickSize/2.0);
//				// top back
//				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,file->zMin);
//				glVertex3f(file->xMax-tickSize,file->yMin+(float)j*gridSpacing,file->zMin);
//				// right front
//				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,file->zMin);
//				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,file->zMin+tickSize);
//				// right back
//				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,file->zMax);
//				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,file->zMax-tickSize);
			glEnd();
		}
//		// depth segments
//		for (int j = 1; j <= (int)(zGridRange/gridSpacing); j++) {
//			glBegin(GL_LINES);
//				// TICKS
//				// top left x
//				glVertex3d(file->xMin,file->yMax,file->zMin+(float)j*gridSpacing);
//				glVertex3d(file->xMin+tickSize,file->yMax,file->zMin+(float)j*gridSpacing);
//				// bottom right y
//				glVertex3d(file->xMax,file->yMin,file->zMin+(float)j*gridSpacing);
//				glVertex3d(file->xMax,file->yMin+tickSize,file->zMin+(float)j*gridSpacing);
//				// top right y
//				glVertex3d(file->xMax,file->yMax,file->zMin+(float)j*gridSpacing);
//				glVertex3d(file->xMax,file->yMax-tickSize,file->zMin+(float)j*gridSpacing);
//				// top right x
//				glVertex3d(file->xMax,file->yMax,file->zMin+(float)j*gridSpacing);
//				glVertex3d(file->xMax-tickSize,file->yMax,file->zMin+(float)j*gridSpacing);
//			glEnd();
//		}
	glPopMatrix();

	glDisable(GL_BLEND);
}

void drawGrid() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float xGridRange = file->xMax - file->xMin;
	float yGridRange = file->yMax - file->yMin;

	float gridSpacing = 0.1;

	if (iMAP->annotationsGui != NULL) {
		if (iMAP->annotationsGui->relativeSpacing()) {
			gridSpacing = iMAP->annotationsGui->getSpacing();
		}
		else {
			gridSpacing = (file->maxRange)/15.0f;

			gridSpacing *= 1000.0;

			const int digits = countDigits((int)(gridSpacing))-1;
			const float pwr = pow((float)10.0,(float)digits);

			gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
			gridSpacing /= 1000.0;
		}
	}
	else {
		gridSpacing = (file->maxRange)/15.0f;

		gridSpacing *= 1000.0;

		const int digits = countDigits((int)(gridSpacing))-1;
		const float pwr = pow((float)10.0,(float)digits);

		gridSpacing = floor(gridSpacing/pwr+0.5)*pwr;
		gridSpacing /= 1000.0;
	}

	glPushMatrix();
		updateView();

		glLineWidth(1);

		glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f);

		// vertical segments
		for (int j = 1; j <= (int)(xGridRange/gridSpacing); j++) {
			glBegin(GL_LINES);
					// GRIDS
					// bottom back
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMin);
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
					// top back
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMin);
					// top z front
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
					// top z back
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
					glVertex2f(file->xMin+(float)j*gridSpacing,file->yMax);
			glEnd();
		}
		// horizontal segments
		for (int j = 1; j <= (int)(yGridRange/gridSpacing); j++) {
			glBegin(GL_LINES);
				// GRIDS
				// bottom back
				glVertex3f(file->xMin,file->yMin+(float)j*gridSpacing,0.0);
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
				// top back
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
				glVertex3f(file->xMin,file->yMin+(float)j*gridSpacing,0.0);
				// right front
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
				// right back
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
				glVertex3f(file->xMax,file->yMin+(float)j*gridSpacing,0.0);
			glEnd();
		}
	glPopMatrix();

	glDisable(GL_BLEND);
}

void umUnitsButtonCallback(Fl_Check_Button*w,int*v) {
	if (w->value()) {
		iMAP->annotationsGui->setMicrometres();
		iMAP->units = 0;
	}
}

void nmUnitsButtonCallback(Fl_Check_Button*w,int*v) {
	if (w->value()) {
		iMAP->annotationsGui->setNanometres();
		iMAP->units = 1;
	}
}

int countDigits(int number) {
    if (number < 10) {
        return 1;
    }
    int count = 0;
    while (number > 0) {
        number /= 10;
        count++;
    }
    return count;
}

/********** 3D ANNOTATIONS ***********/
void drawBoundingBox3D() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glPushMatrix();
		updateView();
		glLineWidth(2);
		glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],1.0f-iMAP->backgroundRGB[0]);
		glBegin(GL_LINE_LOOP);
			glVertex3d(file->xMin,file->yMax,file->zMax);
			glVertex3d(file->xMax,file->yMax,file->zMax);
			glVertex3d(file->xMax,file->yMin,file->zMax);
			glVertex3d(file->xMin,file->yMin,file->zMax);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3d(file->xMin,file->yMax,file->zMin);
			glVertex3d(file->xMax,file->yMax,file->zMin);
			glVertex3d(file->xMax,file->yMin,file->zMin);
			glVertex3d(file->xMin,file->yMin,file->zMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex3d(file->xMin,file->yMax,file->zMax);
			glVertex3d(file->xMin,file->yMax,file->zMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex3d(file->xMin,file->yMin,file->zMax);
			glVertex3d(file->xMin,file->yMin,file->zMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex3d(file->xMax,file->yMax,file->zMax);
			glVertex3d(file->xMax,file->yMax,file->zMin);
		glEnd();
		glBegin(GL_LINE_STRIP);
			glVertex3d(file->xMax,file->yMin,file->zMax);
			glVertex3d(file->xMax,file->yMin,file->zMin);
		glEnd();
	glPopMatrix();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
}

