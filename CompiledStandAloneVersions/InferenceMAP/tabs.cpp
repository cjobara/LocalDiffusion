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

#include "tabs.h"

#include "file.h"
#include "globals.h"

extern Globals *iMAP;
extern File *file;

void closeFile(Tabs*t) {

	// find out which is the corresponding file number of the closed tab
	int filePos = 0;
	while (t != iMAP->tabButtonArray[filePos]) {
		filePos++;
	}

	if (iMAP->fileNumber > 1) {
		t->hide();

		if (iMAP->currentFile == filePos) {
			int d = filePos;
			while (d < iMAP->fileNumber) {

				if (d != 0) { d--; }
				else { d++; }

				if (d != filePos) {
					break;
				}

			}
			iMAP->currentFile = d;
			file = iMAP->fileObjectArray[iMAP->currentFile];
			file->tabButton->do_callback();
		}

//		// close any file-specific guis
//		switch (iMAP->fileObjectArray[filePos]->fileType) {
//		case 0: // ROI file
//			break;
//		default:
//			printf("error\n");
//			break;
//		}

		// delete file corresponding to closed tab
		delete iMAP->tabButtonArray[filePos];
		delete iMAP->fileObjectArray[filePos];

		// adjust fileNumber correspondence
		for (int g = 0; g < iMAP->fileNumber; g++) {
			if (g > filePos) {
				// shift position over one
				iMAP->fileObjectArray[g-1] = iMAP->fileObjectArray[g];
				iMAP->tabButtonArray[g-1] = iMAP->fileObjectArray[g]->tabButton;
			}
		}
		iMAP->fileNumber--;

		// resize widths of any existing tabs
		for (int d = 0; d < iMAP->fileNumber; d++) {
			iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber),20);
		}
	}

}

void selectTab(Tabs*t) {
	// find out which is the corresponding file number of the selected tab
	int filePos = 0;
	while (t != iMAP->tabButtonArray[filePos]) {
		filePos++;
	}
	iMAP->currentFile = filePos;
	file = iMAP->fileObjectArray[iMAP->currentFile];
	file->tabButton->do_callback();
}

Tabs::Tabs(int x,int y,int w,int h,const char*l) : Fl_Button(x,y,w,h,l) {
	this->labelsize(12);
	this->align(FL_ALIGN_INSIDE|FL_ALIGN_CLIP|FL_ALIGN_LEFT);
	this->labelfont(iMAP->normalFont);
	this->color(FL_WHITE);
	this->down_color(FL_WHITE);
	this->box(FL_BORDER_FRAME);
	this->labelcolor(FL_WHITE);
}
