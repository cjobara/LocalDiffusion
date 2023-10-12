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

#include <fstream>

#include <sys/types.h>
//#include <sys/wait.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

// OpenGL Libraries
#ifdef __APPLE__
	#include <unistd.h>
#endif

#include <FL/fl_message.H>
#include <FL/Fl_Sys_Menu_Bar.H>

#include "globals.h"
#include "graphics.h"
#include "file.h"
#include "tabs.h"
#include "gui.h"
#include "tiff.h"
#include "annotations.h"
#include "density.h"

#define PI 3.1415926535897932384626433832795028841971693993751

Globals *iMAP;
File *file;

void mainWindowCallback(Fl_Widget*, void*);
void fileOpenCallback(Fl_Widget*, void*v);
void fileOpenExampleCallback(Fl_Widget*, int*v);
void fileTabsCallback(Tabs*w,void*);
void plotTypeCallback(Fl_Widget*,void*v);
void iOffsetSliderCallback(Slider*w,void*v);
void cMinSliderCallback(Slider*w,void*v);
void cMaxSliderCallback(Slider*w,void*v);
void trajectoriesButtonCallback(Fl_Check_Button*w,void*v);
void animateButtonCallback(Fl_Check_Button*w,void*v);
void aboutCallback(Fl_Widget*,void*);
void squareMeshingCallback(Fl_Widget*,void*);
void voronoiMeshingCallback(Fl_Widget*,void*);
void quadTreeMeshingCallback(Fl_Widget*,void*);
void resizeButtonCallback(Fl_Button*b,void*);
void saveScreenCallback(Fl_Widget*,void*);
void overlayTiffCallback(Fl_Widget*,void*);
void exposureTimeButtonCallback(Fl_Widget*,void*);
void movieMakerCallback(Fl_Widget*,void*);
void updateDisplayButtonCallback(Fl_Check_Button*w,void*);
void userManualCallback(Fl_Widget*,void*);
void licenseCallback(Fl_Widget*,void*);
void animateAccumulationCallback(Fl_Check_Button*w,void*v);
void animateDelayedCallback(Fl_Check_Button*w,void*v);
void startIntervalCallback(Slider*w,void*v);
void endIntervalCallback(Slider*w,void*v);
void batchVoronoiCallback(Fl_Widget*w, int*v);
void batchSingleTrajectoryCallback(Fl_Widget*w, int*v);
void stereoVisionCallback(Fl_Widget*w,void*v);
void whiteBackgroundCallback(Fl_Widget*w,void*v);

void mainWindowCallback(Fl_Widget*, void*) {
	if (Fl::event()==FL_SHORTCUT && Fl::event_key()==FL_Escape) { return; }

	if (fl_ask("Are you sure you want to quit?")) {
		iMAP->glWindow->hide();
		iMAP->mainWindow->hide();
		exit(0);
	}

	return;
}

void userManualCallback(Fl_Widget*,void*) {

	char manualLocation[FILENAME_MAX];
	strcpy(manualLocation,"file://");
	strcat(manualLocation,iMAP->cwd);

	#ifdef __APPLE__
		strcat(manualLocation,"InferenceMAP_Manual.pdf");
	#elif _WIN32
		strcat(manualLocation,"Resources/InferenceMAP_Manual.pdf");
	#endif

	fl_open_uri(manualLocation);
	
}

void licenseCallback(Fl_Widget*,void*) {

	char licenseLocation[FILENAME_MAX];
	strcpy(licenseLocation,"file://");
	strcat(licenseLocation,iMAP->cwd);

	#ifdef __APPLE__
		strcat(licenseLocation,"InferenceMAP_License.pdf");
	#elif _WIN32
		strcat(licenseLocation,"Resources/InferenceMAP_License.pdf");
	#endif

	fl_open_uri(licenseLocation);
	
}

void fileOpenCallback(Fl_Widget*, void*v) {

	/* file type legend
	 *  v == 0 -> xyif file
	 *  v == 1 -> multi csv file
	 *  v == 2 -> SLIMfast file
	 *  v == 3 -> SPTrack file
	 *  v == 4 -> tr x y t file
	 *  v == 5 -> x y t file
	 *  v == 6 -> tr x y z t file
	*/

	if (iMAP->fileNumber < 10) {

		char **list;

		char fileTypeString [1];
		sprintf(fileTypeString,"%i",v);
		const int fileType = atoi(fileTypeString);

		Fl_Native_File_Chooser native;
		char nativeFilename[FILENAME_MAX];

		bool fileOpened = false;

		switch(fileType) {
			case 0: // xyit file
			{
				native.title("Select xyif File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
				native.filter("Text Files\t*.txt");
				native.show();
				strcpy(nativeFilename,native.filename());

				const int fileCount = native.count();

				// load file conditions (to avoid closing file if open screen is cancelled)
				if ((nativeFilename[0]!=(char)NULL) || (nativeFilename[0]!=(char)NULL && iMAP->fileNameStorage[iMAP->fileNumber][0]!=(char)NULL)) {

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}
			case 1: // Multi-Trajectory CSV
			{
				native.title("Select Multi-Trajectory CSV File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_FILE);
				native.filter("Text Files\t*.txt");

				// in case save screen is canceled
				if (native.show() == 0) {

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;

				}
				break;
			}
			case 2: // SLIMfast
			{

				native.title("Select SLIMfast File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_FILE);
				native.filter("Text Files\t*.txt");

				// in case save screen is canceled
				if (native.show() == 0) {

					// user must enter exposure time

					char et [10];
					sprintf(et,"%.1f",30.0);

					iMAP->fileParametersWindow = new Fl_Window(200,105,"File Parameters");
					iMAP->fileParametersWindow->callback((Fl_Callback*)nullCallback);
					iMAP->fileParametersWindow->begin();
					iMAP->fileParametersWindow->set_non_modal();

					iMAP->exposureTimeInput = new Fl_Float_Input(10,10,70,25,"Exposure Time [ms]");
					iMAP->exposureTimeInput->align(FL_ALIGN_RIGHT);
					iMAP->exposureTimeInput->labelsize(12);
					iMAP->exposureTimeInput->textsize(12);
					iMAP->exposureTimeInput->value(et);
					iMAP->exposureTimeInput->show();

					sprintf(et,"%.1f",160.0);
					iMAP->pixelSizeInput = new Fl_Float_Input(10,40,70,25,"Pixel Size [nm/px]");
					iMAP->pixelSizeInput->align(FL_ALIGN_RIGHT);
					iMAP->pixelSizeInput->labelsize(12);
					iMAP->pixelSizeInput->textsize(12);
					iMAP->pixelSizeInput->value(et);
					iMAP->pixelSizeInput->show();

					iMAP->fileParametersButton = new Fl_Button(10,70,180,25,"Apply");
					iMAP->fileParametersButton->labelfont(1);
					iMAP->fileParametersButton->callback((Fl_Callback*)exposureTimeButtonCallback);
					iMAP->fileParametersButton->show();

					iMAP->fileParametersWindow->end();
					iMAP->fileParametersWindow->show();

					while(iMAP->fileParametersWindow->shown()) { Fl::wait(); }

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}
			case 3: // SPTrack
			{
				// user must enter exposure time

				native.title("Select SPTrack File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_FILE);
				native.filter("Text Files\t*.trc");

				// in case save screen is canceled
				if (native.show() == 0) {

					char et [10];
					sprintf(et,"%.1f",25.0);

					iMAP->fileParametersWindow = new Fl_Window(200,105,"File Parameters");
					iMAP->fileParametersWindow->callback((Fl_Callback*)nullCallback);
					iMAP->fileParametersWindow->begin();
					iMAP->fileParametersWindow->set_non_modal();

					iMAP->exposureTimeInput = new Fl_Float_Input(10,10,70,25,"Exposure Time [ms]");
					iMAP->exposureTimeInput->align(FL_ALIGN_RIGHT);
					iMAP->exposureTimeInput->labelsize(12);
					iMAP->exposureTimeInput->textsize(12);
					iMAP->exposureTimeInput->value(et);
					iMAP->exposureTimeInput->show();

					sprintf(et,"%.1f",160.0);
					iMAP->pixelSizeInput = new Fl_Float_Input(10,40,70,25,"Pixel Size [nm/px]");
					iMAP->pixelSizeInput->align(FL_ALIGN_RIGHT);
					iMAP->pixelSizeInput->labelsize(12);
					iMAP->pixelSizeInput->textsize(12);
					iMAP->pixelSizeInput->value(et);
					iMAP->pixelSizeInput->show();

					iMAP->fileParametersButton = new Fl_Button(10,70,180,25,"Apply");
					iMAP->fileParametersButton->labelfont(1);
					iMAP->fileParametersButton->callback((Fl_Callback*)exposureTimeButtonCallback);
					iMAP->fileParametersButton->show();

					iMAP->fileParametersWindow->end();
					iMAP->fileParametersWindow->show();

					while(iMAP->fileParametersWindow->shown()) { Fl::wait(); }

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}
			case 4: // tr x y t
			{
				native.title("Select tr x y t File");
				native.type(Fl_Native_File_Chooser::BROWSE_FILE);
				native.filter("Text Files\t*.trxyt");

				// in case save screen is canceled
				if (native.show() == 0) {

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}
			case 5: // x y t file
			{
				native.title("Select xyt File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
				native.filter("Text Files\t*.xyt");


				// load file conditions (to avoid closing file if open screen is cancelled)
				if (native.show() == 0) {
					
					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());
					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));

					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}

			case 6: // tr x y z t
			{
				native.title("Select tr x y z t File");
				native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
				native.filter("Text Files\t*.trxyzt");

				// in case save screen is canceled
				if (native.show() == 0) {

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					iMAP->file3D = true;
					iMAP->fovy = 2.0*atan(iMAP->fileObjectArray[iMAP->fileNumber]->maxRange/2.0/15.0)*180.0/PI;

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}
			default:
				break;
		}

		if (fileOpened) {

			// hide any existing gui windows for currentFile
			if (iMAP->fileNumber != 0) {
				// hide meshing GUIs if open
				if (file->squareMeshGui != NULL) {
					file->squareMeshGui->close();
					file->squareMeshGui = NULL;
				}
				if (file->voronoiMeshGui != NULL) {
					file->voronoiMeshGui->close();
					file->voronoiMeshGui = NULL;
				}
				if (file->treeMeshGui != NULL) {
					file->treeMeshGui->close();
					file->treeMeshGui = NULL;
				}
				if (iMAP->singleTrajectoryInferenceGui != NULL) {
					iMAP->singleTrajectoryInferenceGui->close();
					iMAP->singleTrajectoryInferenceGui = NULL;
				}
				if (iMAP->overlayTiffGui != NULL) {
					iMAP->overlayTiffGui->hide();
					iMAP->overlayTiffGui->window->do_callback();
				}
			}

			file = iMAP->fileObjectArray[iMAP->fileNumber];
			file->tabButton = new Tabs(iMAP->tabBarWidth*(iMAP->fileNumber)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			file->tabButton->labelfont(iMAP->normalFont);
			file->tabButton->labelcolor(FL_WHITE);
			file->tabButton->copy_label(file->fileName);
			file->tabButton->show();

			iMAP->mainWindow->add(file->tabButton);

			// resize widths of any existing tabs
			for (int d = 0; d < iMAP->fileNumber; d++) {
				iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			}

			// add to tab button array (for tracking files)
			iMAP->tabButtonArray[iMAP->fileNumber] = file->tabButton;
			file->tabButton->callback((Fl_Callback*)fileTabsCallback);

			for (int g = 0; g < iMAP->fileNumber; g++) {
				iMAP->tabButtonArray[g]->labelfont(iMAP->normalFont);
				file->tabButton->labelcolor(FL_WHITE);
			}

			char topBarTitle [FILENAME_MAX];
			sprintf(topBarTitle,"InferenceMAP - %s",iMAP->fileNameStorage[iMAP->fileNumber]);
			iMAP->mainWindow->label(topBarTitle);

			// adjust sliders
			iMAP->iOffsetSlider->value(file->iOffset);
			iMAP->cMaxSlider->value(file->cMax);
			iMAP->cMinSlider->value(file->cMin);
			iMAP->drawTrajectoriesButton->value(file->drawTrajectories);
			iMAP->animateTrajectoriesButton->value(file->animateTrajectories);
			iMAP->localizationsButton->value(file->drawLocalizations);

			iMAP->startIntervalSlider->activate();
			iMAP->startIntervalSlider->value(file->startInterval);
			iMAP->startIntervalSlider->bounds(file->tMin,file->tMax);
			iMAP->endIntervalSlider->activate();
			iMAP->endIntervalSlider->value(file->endInterval);
			iMAP->endIntervalSlider->bounds(file->tMin,file->tMax);

			iMAP->fpsSlider->activate();
			iMAP->colorbar->activate();
			iMAP->fileInfo->activate();
			iMAP->customSelectionInferenceGui->activate();
			iMAP->iOffsetSlider->activate();
			iMAP->cMinSlider->activate();
			iMAP->cMaxSlider->activate();
			iMAP->drawTrajectoriesButton->activate();
			iMAP->localizationsButton->activate();
			iMAP->animateTrajectoriesButton->activate();

			iMAP->densityGui->fileOpen();

			// adjust zoom level
			file->orthoLimit = 0.6*file->maxRange;

			// update number of files opened
			iMAP->fileNumber++;

			file->tabButton->do_callback();

			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

			iMAP->fileInfo->update();

			Fl::redraw();
		}

	} else { fl_alert("Cannot open more than 10 files."); }
}

void fileOpenExampleCallback(Fl_Widget*, int*v) {

	/* file type legend
	 *  v == 0 -> xyif file
	 *  v == 1 -> multi csv file
	 *  v == 2 -> SLIMfast file
	 *  v == 3 -> SPTrack file
	 *  v == 4 -> tr x y t file
	 *  v == 5 -> tr x y z t file
	 *  v == 6 -> x y t file
	*/

	char exampleString [5];
	sprintf(exampleString,"%i",v);
	const int exampleNumber = atoi(exampleString);

	char exampleFilename[FILENAME_MAX];

	bool fileOpened = false;

	const int fileCount = 1;

	if (iMAP->fileNumber < 10) {

		char **list = NULL;

		switch(exampleNumber) {
			case 0: // Diffusion Example
				// store filename
				strcpy(exampleFilename,iMAP->cwd);

				#ifdef __APPLE__
					strcat(exampleFilename,"diffusion_example.trxyt");
				#elif _WIN32
					strcat(exampleFilename,"/Resources/diffusion_example.trxyt");
				#endif

				list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
				for (int e = 0; e < fileCount; e++) {
					list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
				}

				// copy files to list
				for (int f = 0; f < fileCount; f++) {
					strcpy(list[f],exampleFilename);
					// printf("%s\n",list[f]);
				}

				strcpy(iMAP->fileNameStorage[iMAP->fileNumber],exampleFilename);

				// instantiate plotObject
				iMAP->fileObjectArray[iMAP->fileNumber] = new File(1,list,4);

				// refresh all widgets
				Fl::redraw();
				fileOpened = true;

				//if (iMAP->overlayTiffGui != NULL) {
				//	iMAP->overlayTiffGui->hide();
				//	delete iMAP->overlayTiffGui;
				//	iMAP->overlayTiffGui = NULL;
				//}
				if (iMAP->overlayTiffGui != NULL) {
					iMAP->overlayTiffGui->hide();
					iMAP->overlayTiffGui->window->do_callback();
				}

				break;
			case 1: // Potential Example
				// store filename
				strcpy(exampleFilename,iMAP->cwd);

				#ifdef __APPLE__
					strcat(exampleFilename,"potential_example.trxyt");
				#elif _WIN32
					strcat(exampleFilename,"/Resources/potential_example.trxyt");
				#endif

				list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
				for (int e = 0; e < fileCount; e++) {
					list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
				}

				// copy files to list
				for (int f = 0; f < fileCount; f++) { strcpy(list[f],exampleFilename); }

				strcpy(iMAP->fileNameStorage[iMAP->fileNumber],exampleFilename);

				// instantiate plotObject
				iMAP->fileObjectArray[iMAP->fileNumber] = new File(1,list,4);

				// refresh all widgets
				Fl::redraw();
				fileOpened = true;

				if (iMAP->overlayTiffGui != NULL) {
					iMAP->overlayTiffGui->hide();
					iMAP->overlayTiffGui->window->do_callback();
				}

				break;
			case 2: // Epsilon toxin hopping lipid rafts
				// store filename
				strcpy(exampleFilename,iMAP->cwd);

				#ifdef __APPLE__
					strcat(exampleFilename,"lipid_raft.xyt");
				#elif _WIN32
					strcat(exampleFilename,"/Resources/lipid_raft.xyt");
				#endif

				list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
				for (int e = 0; e < fileCount; e++) {
					list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
				}

				// copy files to list
				for (int f = 0; f < fileCount; f++) {
					strcpy(list[f],exampleFilename);
					// printf("%s\n",list[f]);
				}

				strcpy(iMAP->fileNameStorage[iMAP->fileNumber],exampleFilename);

				// instantiate plotObject
				iMAP->fileObjectArray[iMAP->fileNumber] = new File(1,list,5);

				// refresh all widgets
				Fl::redraw();
				fileOpened = true;

				break;
			case 3: // Glycine receptor in mouse hippocampal neuron
				// store filename
				strcpy(exampleFilename,iMAP->cwd);
				
				#ifdef __APPLE__
					strcat(exampleFilename,"glycine_receptor.trxyt");
				#elif _WIN32
					strcat(exampleFilename,"/Resources/glycine_receptor.trxyt");
				#endif

				list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
				for (int e = 0; e < fileCount; e++) {
					list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
				}

				// copy files to list
				for (int f = 0; f < fileCount; f++) {
					strcpy(list[f],exampleFilename);
					// printf("%s\n",list[f]);
				}

				strcpy(iMAP->fileNameStorage[iMAP->fileNumber],exampleFilename);

				// instantiate plotObject
				iMAP->fileObjectArray[iMAP->fileNumber] = new File(1,list,4);

				// refresh all widgets
				Fl::redraw();
				fileOpened = true;

				if (iMAP->overlayTiffGui != NULL) {
					iMAP->overlayTiffGui->hide();
					iMAP->overlayTiffGui->window->do_callback();
				}

				break;
			case 4: // Transmembrane Domain in COS-7
				// store filename
				strcpy(exampleFilename,iMAP->cwd);
				strcat(exampleFilename,"transmembrane_domain.trxyt");

				list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
				for (int e = 0; e < fileCount; e++) {
					list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
				}

				// copy files to list
				for (int f = 0; f < fileCount; f++) {
					strcpy(list[f],exampleFilename);
					// printf("%s\n",list[f]);
				}

				strcpy(iMAP->fileNameStorage[iMAP->fileNumber],exampleFilename);

				// instantiate plotObject
				iMAP->fileObjectArray[iMAP->fileNumber] = new File(1,list,4);

				// refresh all widgets
				Fl::redraw();
				fileOpened = true;

				break;
		}

		if (fileOpened) {

			// hide any existing gui windows for currentFile
			// hide any existing gui windows for currentFile
			if (iMAP->fileNumber != 0) {
				// hide meshing GUIs if open
				if (file->squareMeshGui != NULL) {
					file->squareMeshGui->close();
					file->squareMeshGui = NULL;
				}
				if (file->voronoiMeshGui != NULL) {
					file->voronoiMeshGui->close();
					file->voronoiMeshGui = NULL;
				}
				if (file->treeMeshGui != NULL) {
					file->treeMeshGui->close();
					file->treeMeshGui = NULL;
				}
				if (iMAP->singleTrajectoryInferenceGui != NULL) {
					iMAP->singleTrajectoryInferenceGui->close();
					iMAP->singleTrajectoryInferenceGui = NULL;
				}
				if (iMAP->overlayTiffGui != NULL) {
					iMAP->overlayTiffGui->hide();
					iMAP->overlayTiffGui->window->do_callback();
				}
			}

			file = iMAP->fileObjectArray[iMAP->fileNumber];
			file->tabButton = new Tabs(iMAP->tabBarWidth*(iMAP->fileNumber)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			file->tabButton->labelfont(iMAP->normalFont);
			file->tabButton->labelcolor(FL_WHITE);
			file->tabButton->copy_label(file->fileName);
			file->tabButton->show();

			iMAP->mainWindow->add(file->tabButton);

			// resize widths of any existing tabs
			for (int d = 0; d < iMAP->fileNumber; d++) {
				iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			}

			// add to tab button array (for tracking files)
			iMAP->tabButtonArray[iMAP->fileNumber] = file->tabButton;
			file->tabButton->callback((Fl_Callback*)fileTabsCallback);

			for (int g = 0; g < iMAP->fileNumber; g++) {
				iMAP->tabButtonArray[g]->labelfont(iMAP->normalFont);
				file->tabButton->labelcolor(FL_WHITE);
			}

			char topBarTitle [FILENAME_MAX];
			sprintf(topBarTitle,"InferenceMAP - %s",iMAP->fileNameStorage[iMAP->fileNumber]);
			iMAP->mainWindow->label(topBarTitle);

			// adjust sliders
			iMAP->iOffsetSlider->value(file->iOffset);
			iMAP->cMaxSlider->value(file->cMax);
			iMAP->cMinSlider->value(file->cMin);
			iMAP->drawTrajectoriesButton->value(file->drawTrajectories);
			iMAP->animateTrajectoriesButton->value(file->animateTrajectories);
			iMAP->localizationsButton->value(file->drawLocalizations);

			iMAP->startIntervalSlider->activate();
			iMAP->startIntervalSlider->value(file->startInterval);
			iMAP->startIntervalSlider->bounds(file->tMin,file->tMax);
			iMAP->endIntervalSlider->activate();
			iMAP->endIntervalSlider->value(file->endInterval);
			iMAP->endIntervalSlider->bounds(file->tMin,file->tMax);

			iMAP->fpsSlider->activate();
			iMAP->colorbar->activate();
			iMAP->fileInfo->activate();
			iMAP->customSelectionInferenceGui->activate();
			iMAP->iOffsetSlider->activate();
			iMAP->cMinSlider->activate();
			iMAP->cMaxSlider->activate();
			iMAP->drawTrajectoriesButton->activate();
			iMAP->localizationsButton->activate();
			iMAP->animateTrajectoriesButton->activate();
			iMAP->densityGui->fileOpen();

			// adjust zoom level
			file->orthoLimit = 0.6*file->maxRange;

			// update number of files opened
			iMAP->fileNumber++;

			file->tabButton->do_callback();

			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

			// overlay TIFF
			switch(exampleNumber) {
				case 0:
				{
					iMAP->drawTrajectoriesButton->value(0);
					iMAP->drawTrajectoriesButton->do_callback();
					iMAP->animateTrajectoriesButton->value(1);
					iMAP->animateTrajectoriesButton->do_callback();
					iMAP->animateDelayedButton->value(1);
					iMAP->animateDelayedButton->do_callback();

					char filename[FILENAME_MAX];
					strcpy(filename,iMAP->cwd);
					
					#ifdef __APPLE__
						strcat(filename,"diffusion_example.tif");
					#elif _WIN32
						strcat(filename,"/Resources/diffusion_example.tif");
					#endif

					iMAP->tiffSequence = new TiffSequence(filename);
					iMAP->overlayTiffGui = new OverlayTiffGui();
					iMAP->tiffOverlay = true;
					iMAP->overlayTiffGui->clampButton->value(1);
					iMAP->overlayTiffGui->normalizeButton->value(1);
					iMAP->overlayTiffGui->normalizeButton->do_callback();
					iMAP->overlayTiffGui->yFlipButton->value(1);
					iMAP->overlayTiffGui->yFlipButton->do_callback();
					iMAP->overlayTiffGui->resolutionSlider->value(20);
				}
				break;
				case 1:
				{
					iMAP->drawTrajectoriesButton->value(0);
					iMAP->drawTrajectoriesButton->do_callback();
					iMAP->animateTrajectoriesButton->value(1);
					iMAP->animateTrajectoriesButton->do_callback();
					iMAP->animateDelayedButton->value(1);
					iMAP->animateDelayedButton->do_callback();

					char filename[FILENAME_MAX];
					strcpy(filename,iMAP->cwd);

					#ifdef __APPLE__
						strcat(filename,"potential_example.tif");
					#elif _WIN32
						strcat(filename,"/Resources/potential_example.tif");
					#endif

					iMAP->tiffSequence = new TiffSequence(filename);
					iMAP->overlayTiffGui = new OverlayTiffGui();
					iMAP->tiffOverlay = true;
					iMAP->overlayTiffGui->clampButton->value(1);
					iMAP->overlayTiffGui->normalizeButton->value(1);
					iMAP->overlayTiffGui->normalizeButton->do_callback();
					iMAP->overlayTiffGui->yFlipButton->value(1);
					iMAP->overlayTiffGui->yFlipButton->do_callback();
					iMAP->overlayTiffGui->resolutionSlider->value(10);
				}
				break;
				case 3:
				{
					iMAP->drawTrajectoriesButton->value(0);
					iMAP->drawTrajectoriesButton->do_callback();
					iMAP->animateTrajectoriesButton->value(1);
					iMAP->animateTrajectoriesButton->do_callback();

					char filename[FILENAME_MAX];
					strcpy(filename,iMAP->cwd);
					
					#ifdef __APPLE__
						strcat(filename,"glycine_receptor_gfp.tif");
					#elif _WIN32
						strcat(filename,"/Resources/glycine_receptor_gfp.tif");
					#endif

					iMAP->tiffSequence = new TiffSequence(filename);
					iMAP->overlayTiffGui = new OverlayTiffGui();
					iMAP->tiffOverlay = true;
					iMAP->overlayTiffGui->xAlignSlider->value(-4000);
					iMAP->overlayTiffGui->yAlignSlider->value(-3800);
					iMAP->overlayTiffGui->clampButton->value(1);
					iMAP->overlayTiffGui->contrastSlider->value(-4.0);
					iMAP->overlayTiffGui->gammaSlider->value(1.3);
					iMAP->overlayTiffGui->normalizeButton->value(1);
					iMAP->overlayTiffGui->normalizeButton->do_callback();
				}
				break;
			}
			iMAP->fileInfo->update();

			Fl::redraw();
		}

	} else { fl_alert("Cannot open more than 10 files."); }
}

void batchVoronoiCallback(Fl_Widget*w, int*v) {

	char **list;

	char fileTypeString [1];
	sprintf(fileTypeString,"%i",v);
	const int fileType = atoi(fileTypeString);

	Fl_Native_File_Chooser native;
	char nativeFilename[FILENAME_MAX];

	bool fileOpened = false;

	native.title("Select tr x y t Files");
	native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
	native.filter("Text Files\t*.trxyt");
	int fileCount;

	// in case save screen is canceled
	if (native.show() == 0) {

		fileCount = native.count();

		strcpy(nativeFilename,native.filename());

		list = new char*[fileCount];
		for (int e = 0; e < fileCount; e++) {
			list[e] = new char[FILENAME_MAX];
		}

		// copy files to list
		for (int f = 0; f < fileCount; f++) {
			strcpy(list[f],native.filename(f));
			// printf("%s\n",list[f]);
		}

		fileOpened = true;
	}

	// loop opening and analysis of files
	if (fileOpened) {
		char **listSingleFile;
		listSingleFile = new char*[1];
		listSingleFile[0] = new char[FILENAME_MAX];

		for (int e = 0; e < fileCount; e++) {

			// store filename
			strcpy(iMAP->fileNameStorage[iMAP->fileNumber], list[e]);
			strcpy(listSingleFile[0],list[e]);

			// instantiate plotObject
			iMAP->fileObjectArray[iMAP->fileNumber] = new File(1, listSingleFile, fileType);

			// hide any existing gui windows for currentFile
			if (iMAP->fileNumber != 0) {
				// hide meshing GUIs if open
				if (file->squareMeshGui != NULL) {
					file->squareMeshGui->close();
				}
				if (file->voronoiMeshGui != NULL) {
					file->voronoiMeshGui->close();
				}
				if (file->treeMeshGui != NULL) {
					file->treeMeshGui->close();
				}
				if (iMAP->singleTrajectoryInferenceGui != NULL) {
					iMAP->singleTrajectoryInferenceGui->close();
				}
			}

			file = iMAP->fileObjectArray[iMAP->fileNumber];
			file->tabButton = new Tabs(iMAP->tabBarWidth*(iMAP->fileNumber)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			file->tabButton->labelfont(iMAP->normalFont);
			file->tabButton->labelcolor(FL_WHITE);
			file->tabButton->copy_label(file->fileName);
			file->tabButton->show();

			iMAP->mainWindow->add(file->tabButton);

			// resize widths of any existing tabs
			for (int d = 0; d < iMAP->fileNumber; d++) {
				iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			}

			// add to tab button array (for tracking files)
			iMAP->tabButtonArray[iMAP->fileNumber] = file->tabButton;
			file->tabButton->callback((Fl_Callback*)fileTabsCallback);

			for (int g = 0; g < iMAP->fileNumber; g++) {
				iMAP->tabButtonArray[g]->labelfont(iMAP->normalFont);
				file->tabButton->labelcolor(FL_WHITE);
			}

			char topBarTitle [FILENAME_MAX];
			sprintf(topBarTitle,"InferenceMAP - %s",iMAP->fileNameStorage[iMAP->fileNumber]);
			iMAP->mainWindow->label(topBarTitle);

			// adjust sliders
			iMAP->iOffsetSlider->value(file->iOffset);
			iMAP->cMaxSlider->value(file->cMax);
			iMAP->cMinSlider->value(file->cMin);
			iMAP->drawTrajectoriesButton->value(file->drawTrajectories);
			iMAP->animateTrajectoriesButton->value(file->animateTrajectories);
			iMAP->localizationsButton->value(file->drawLocalizations);

			iMAP->startIntervalSlider->activate();
			iMAP->startIntervalSlider->value(file->startInterval);
			iMAP->startIntervalSlider->bounds(file->tMin,file->tMax);
			iMAP->endIntervalSlider->activate();
			iMAP->endIntervalSlider->value(file->endInterval);
			iMAP->endIntervalSlider->bounds(file->tMin,file->tMax);

			iMAP->fpsSlider->activate();
			iMAP->colorbar->activate();
			iMAP->fileInfo->activate();
			iMAP->customSelectionInferenceGui->activate();
			iMAP->iOffsetSlider->activate();
			iMAP->cMinSlider->activate();
			iMAP->cMaxSlider->activate();
			iMAP->drawTrajectoriesButton->activate();
			iMAP->localizationsButton->activate();
			iMAP->animateTrajectoriesButton->activate();

			// adjust zoom level
			file->orthoLimit = 0.6*file->maxRange;

			// update number of files opened
			iMAP->fileNumber++;

			file->tabButton->do_callback();

			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

			iMAP->fileInfo->update();

			Fl::redraw();

			/* START VORONOI MESHING */
			file->meshType = 1;
			file->voronoiMeshGui = new VoronoiMeshGui();
			file->voronoiMeshGui->show();

			// set number of zones
//			file->voronoiMeshGui->cellsSlider->value(FMIN(FMAX(15.0,file->localizationCount/50.0),80));
			file->voronoiMeshGui->cellsSlider->value(file->localizationCount/120.0);
			file->voronoiMeshGui->cellsSlider->do_callback();

			// apply mesh
			file->voronoiMeshGui->applyButton->do_callback();

			// adjust minimum points per zone
			file->voronoiMeshGui->minPointsSlider->value(20);
			file->voronoiMeshGui->minPointsSlider->do_callback();

			// adjust localization precision
			file->voronoiMeshGui->localizationPrecisionSlider->value(30);
			file->voronoiMeshGui->localizationPrecisionSlider->do_callback();

			// adjust maximum neighbour distance
			file->voronoiMeshGui->neighbourDistanceSlider->value(2000.0);
			file->voronoiMeshGui->neighbourDistanceSlider->do_callback();

			iMAP->iOffsetSlider->value(0.4);
			file->voronoiMeshGui->overlayAlphaSlider->value(0.65);
			file->voronoiMeshGui->overlayAlphaSlider->do_callback();

			file->voronoiMeshGui->overlayDiffusionButton->value(1);
			file->voronoiMeshGui->overlayDiffusionButton->do_callback();

			// apply prior
//			file->voronoiMeshGui->priorButton->value(1);
//			file->voronoiMeshGui->priorButton->do_callback();

			// select DV inference mode
			file->voronoiMeshGui->inferenceModeChoice->value(0);
			file->voronoiMeshGui->inferenceModeChoice->do_callback();

			// apply inference
			file->voronoiMeshGui->inferButton->do_callback();

			Fl::lock();
			{
				// adjust c axes
//				const double dMax = 0.20;
//				const double dMid = (file->voronoiMesh->getDiffusionMax()-file->voronoiMesh->getDiffusionMin())/2.0;
//
//				iMAP->cMaxSlider->value( dMax/2.0/dMid );
//				iMAP->cMaxSlider->do_callback();
//
//				iMAP->cMinSlider->value( 1.0-dMax/2.0/dMid );
//				iMAP->cMinSlider->do_callback();

//				file->voronoiMeshGui->overlayForceArrowsButton->value(0);
//				file->voronoiMeshGui->overlayForceArrowsButton->do_callback();
//
//				file->voronoiMeshGui->overlayDiffusionButton->value(1);
//				file->voronoiMeshGui->overlayDiffusionButton->do_callback();

//				file->voronoiMeshGui->overlayPotentialButton->value(1);
//				file->voronoiMeshGui->overlayPotentialButton->do_callback();

//				file->voronoiMeshGui->landscapeButton->value(1);
//				file->voronoiMeshGui->landscapeButton->do_callback();
//				file->voronoiMeshGui->landscapeAlphaSlider->value(0.95);
//				file->voronoiMeshGui->landscapeAlphaSlider->do_callback();

//				file->xRotate = -30;
//				file->zRotate = 15;

				file->voronoiMeshGui->displayGridButton->value(0);
				file->voronoiMeshGui->displayGridButton->do_callback();

				file->voronoiMeshGui->overlayButton->value(0);
				file->voronoiMeshGui->overlayButton->do_callback();

//				file->voronoiMeshGui->xLightPositionSlider->value(100);
//				file->voronoiMeshGui->yLightPositionSlider->value(100);

//				file->voronoiMeshGui->scaleLandscapeSlider->value(5.0*file->voronoiMesh->getDiffusionMax());

//				file->voronoiMeshGui->scaleLandscapeSlider->value( (file->voronoiMesh->getPotentialMax()-file->voronoiMesh->getPotentialMin())/vMax );


			}
			Fl::unlock();
			Fl::check();

			iMAP->saveTiffBatch = false;
			void drawVoronoiLandscape() ;

			file->voronoiMesh->exportMeshBatch();

			// save screenshot
			char nativeFilename[FILENAME_MAX];
			char fileNameTemp[FILENAME_MAX];

			const int fileLength = strlen(file->fileName);
			for (int g = 0; g < fileLength-6; g++) { fileNameTemp[g] = file->fileName[g]; }
			fileNameTemp[fileLength-6] = '\0';

			strcpy(nativeFilename,file->filePath);
			strcat(nativeFilename,fileNameTemp);

			captureScreen(nativeFilename);

			writetiffOffscreen(nativeFilename,iMAP->colorbarCanvas,"_colorbar");

			// close voronoi meshing window
			file->voronoiMeshOverlay = false;
			file->voronoiMeshGui->resetButton->do_callback();
			file->voronoiMeshGui->hide();
			file->voronoiMeshGui->saveButton->deactivate();
			file->inferred = false;
			file->landscapePlot = false;
			iMAP->customSelectionInferenceGui->activate();
			if (file->voronoiMeshOverlay) {
				if ( file->voronoiMesh->randomizedOptimizationGui != NULL) {
					file->voronoiMesh->randomizedOptimizationGui->hide();
			//		delete file->voronoiMesh->randomizedOptimizationGui;
				}
			}
			delete file->voronoiMeshGui;
			file->voronoiMeshGui = NULL;
			iMAP->overlayAdjustment = true;

			if (e > 0) {
				closeFile(file->tabButton);
			}

			fprintf(stderr,"File %i / %i\n",e+1,fileCount);
		}
	}

	// refresh all widgets
	Fl::redraw();
}

void batchSingleTrajectoryCallback(Fl_Widget*w, int*v) {

	char **list;

	char fileTypeString [1];
	sprintf(fileTypeString,"%i",v);
	const int fileType = atoi(fileTypeString);

	Fl_Native_File_Chooser native;
	char nativeFilename[FILENAME_MAX];

	bool fileOpened = false;

	native.title("Select tr x y t Files");
	native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
	native.filter("Text Files\t*.trxyt");
	int fileCount;

	// in case save screen is canceled
	if (native.show() == 0) {

		fileCount = native.count();

		strcpy(nativeFilename,native.filename());

		list = new char*[fileCount];
		for (int e = 0; e < fileCount; e++) {
			list[e] = new char[FILENAME_MAX];
		}

		// copy files to list
		for (int f = 0; f < fileCount; f++) {
			strcpy(list[f],native.filename(f));
			// printf("%s\n",list[f]);
		}

		fileOpened = true;
	}

	// loop opening and analysis of files
	if (fileOpened) {
		char **listSingleFile;
		listSingleFile = new char*[1];
		listSingleFile[0] = new char[FILENAME_MAX];

		for (int e = 0; e < fileCount; e++) {

			// store filename
			strcpy(iMAP->fileNameStorage[iMAP->fileNumber], list[e]);
			strcpy(listSingleFile[0],list[e]);

			// instantiate plotObject
			iMAP->fileObjectArray[iMAP->fileNumber] = new File(1, listSingleFile, fileType);

			// hide any existing gui windows for currentFile
			if (iMAP->fileNumber != 0) {
				// hide meshing GUIs if open
				if (file->squareMeshGui != NULL) {
					file->squareMeshGui->close();
				}
				if (file->voronoiMeshGui != NULL) {
					file->voronoiMeshGui->close();
				}
				if (file->treeMeshGui != NULL) {
					file->treeMeshGui->close();
				}
				if (iMAP->singleTrajectoryInferenceGui != NULL) {
					iMAP->singleTrajectoryInferenceGui->close();
				}
			}

			file = iMAP->fileObjectArray[iMAP->fileNumber];
			file->tabButton = new Tabs(iMAP->tabBarWidth*(iMAP->fileNumber)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			file->tabButton->labelfont(iMAP->normalFont);
			file->tabButton->labelcolor(FL_WHITE);
			file->tabButton->copy_label(file->fileName);
			file->tabButton->show();

			iMAP->mainWindow->add(file->tabButton);

			// resize widths of any existing tabs
			for (int d = 0; d < iMAP->fileNumber; d++) {
				iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber+1),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber+1),20);
			}

			// add to tab button array (for tracking files)
			iMAP->tabButtonArray[iMAP->fileNumber] = file->tabButton;
			file->tabButton->callback((Fl_Callback*)fileTabsCallback);

			for (int g = 0; g < iMAP->fileNumber; g++) {
				iMAP->tabButtonArray[g]->labelfont(iMAP->normalFont);
				file->tabButton->labelcolor(FL_WHITE);
			}

			char topBarTitle [FILENAME_MAX];
			sprintf(topBarTitle,"InferenceMAP - %s",iMAP->fileNameStorage[iMAP->fileNumber]);
			iMAP->mainWindow->label(topBarTitle);

			// adjust sliders
			iMAP->iOffsetSlider->value(file->iOffset);
			iMAP->cMaxSlider->value(file->cMax);
			iMAP->cMinSlider->value(file->cMin);
			iMAP->drawTrajectoriesButton->value(file->drawTrajectories);
			iMAP->animateTrajectoriesButton->value(file->animateTrajectories);
			iMAP->localizationsButton->value(file->drawLocalizations);

			iMAP->startIntervalSlider->activate();
			iMAP->startIntervalSlider->value(file->startInterval);
			iMAP->startIntervalSlider->bounds(file->tMin,file->tMax);
			iMAP->endIntervalSlider->activate();
			iMAP->endIntervalSlider->value(file->endInterval);
			iMAP->endIntervalSlider->bounds(file->tMin,file->tMax);

			iMAP->fpsSlider->activate();
			iMAP->colorbar->activate();
			iMAP->fileInfo->activate();
			iMAP->customSelectionInferenceGui->activate();
			iMAP->iOffsetSlider->activate();
			iMAP->cMinSlider->activate();
			iMAP->cMaxSlider->activate();
			iMAP->drawTrajectoriesButton->activate();
			iMAP->localizationsButton->activate();
			iMAP->animateTrajectoriesButton->activate();

			// adjust zoom level
			file->orthoLimit = 0.6*file->maxRange;

			// update number of files opened
			iMAP->fileNumber++;

			file->tabButton->do_callback();

			iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

			iMAP->fileInfo->update();

			Fl::redraw();

			/* START SINGLE TRAJECTORY GUI */
			file->meshType = -1;
			iMAP->singleTrajectoryInferenceGui = new SingleTrajectoryInferenceGui(file->fileName);

			iMAP->singleTrajectoryInferenceGui->inferenceButton->do_callback();

			iMAP->singleTrajectoryInferenceGui->sampleDiffusionButton->do_callback();
			iMAP->singleTrajectoryInferenceGui->saveDPosterioriButton->do_callback();

			iMAP->singleTrajectoryInferenceGui->window->do_callback();

			if (e > 0) {
				closeFile(file->tabButton);
			}

			fprintf(stderr,"File %i / %i\n",e+1,fileCount);
		}
	}

	// refresh all widgets
	Fl::redraw();
}

void stereoVisionCallback(Fl_Widget*w,void*v) {

	if (iMAP->fileNumber > 0) {


		if (file->landscapePlot == true || file->file3D == true) {
			if (iMAP->stereoVision == true) {

				int w = iMAP->screenDimensions[2];//-90;
				int h = iMAP->screenDimensions[3]-90;
				// hide stereo gui
				iMAP->stereoVisionGui->hide();
				iMAP->stereoVision = false;


				if (iMAP->fullScreenButton->value()) {

					const float ar = 730.0f/1110.0f;
					if ((float)w*ar > (float)h) { w = h+380; }
					else { h = w-380; }

					iMAP->mainWindow->fullscreen_off(0,30-iMAP->macOffset,w,h-iMAP->macOffset);
//					iMAP->mainWindow->resize(8,30-iMAP->macOffset,iMAP->screenWidth,iMAP->screenHeight-iMAP->macOffset);
					iMAP->glWindow->resize(0,50-iMAP->macOffset,h-50,h-50);
				} else {
					iMAP->mainWindow->fullscreen_off(0,0,1110,730-iMAP->macOffset);
//					iMAP->mainWindow->resize(8,30-iMAP->macOffset,1150,760-iMAP->macOffset);
					iMAP->glWindow->resize(0,50-iMAP->macOffset,680,680);
				}
			} else {
				// create stereo gui
				if (iMAP->stereoVisionGui == NULL) { iMAP->stereoVisionGui = new StereoVisionGui();	}
//				else { iMAP->stereoVisionGui->show(); }

//				const int w = iMAP->screenDimensions[2]+90;
//				const int h = iMAP->screenDimensions[3]+90;
				int dims[4];
				Fl::screen_xywh(dims[0],dims[1],dims[2],dims[3],0);

				iMAP->mainWindow->resize(0,0,dims[2],dims[3]);
				iMAP->glWindow->resize(0,0,dims[2],dims[3]);
				iMAP->glWindow->fullscreen();
				iMAP->mainWindow->fullscreen();
				iMAP->stereoVision = true;

			}
			Fl::redraw();

		} else { fl_alert("Must be in landscape viewing mode."); }
	} else { fl_alert("No file open."); }
}

void whiteBackgroundCallback(Fl_Widget*w,void*v) {
	if (iMAP->fileNumber > 0) {
		if (iMAP->whiteBackground == false) {
			iMAP->whiteBackground = true;
			iMAP->backgroundRGB[0] = iMAP->backgroundRGB[1] = iMAP->backgroundRGB[2] = 1.0;
		} else {
			iMAP->whiteBackground = false;
			iMAP->backgroundRGB[0] = iMAP->backgroundRGB[1] = iMAP->backgroundRGB[2] = 0.0;
		}
	}
	else {
		fl_alert("No file open.");
	}
}

void movieMakerCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (!iMAP->movieMakerOpened) {
			iMAP->movieMakerOpened = true;
			iMAP->movieMakerGui = new MovieMakerGui();
		}
		else {
			iMAP->movieMakerGui->show();
		}
	}
	else {
		fl_alert("No file open.");
	}
}

void fileTabsCallback(Tabs*w,void*) {

	// find which tab was selected
	int f = 0;
	while (w != iMAP->tabButtonArray[f]) {
		f++;
	}

	if (iMAP->currentFile != f) {
		// hide meshing GUIs if open
		if (file->squareMeshGui != NULL) {	file->squareMeshGui->hide(); }
		if (file->voronoiMeshGui != NULL) {	file->voronoiMeshGui->hide(); }
		if (file->treeMeshGui != NULL) { file->treeMeshGui->hide(); }
		if (iMAP->overlayTiffGui != NULL) {
			iMAP->overlayTiffGui->hide();
			iMAP->overlayTiffGui->window->do_callback();
		}

		// change current file to match tab number
		iMAP->currentFile = f;
		file = iMAP->fileObjectArray[iMAP->currentFile];

		// hide meshing GUIs if open
		if (file->squareMeshGui != NULL) {	file->squareMeshGui->show(); }
		if (file->voronoiMeshGui != NULL) {	file->voronoiMeshGui->show(); }
		if (file->treeMeshGui != NULL) { file->treeMeshGui->show(); }
		
		// turn off make selection button
		iMAP->customSelectionInferenceGui->selectionButton->value(0);
		iMAP->customSelectionInferenceGui->selectionButton->do_callback();

		// unbold unselected tabs
		for (int g = 0; g < iMAP->fileNumber; g++) {
			iMAP->tabButtonArray[g]->labelcolor(FL_DARK3);
			iMAP->tabButtonArray[g]->labelfont(iMAP->normalFont);
		}
		iMAP->tabButtonArray[iMAP->currentFile]->labelcolor(FL_WHITE);
		iMAP->tabButtonArray[iMAP->currentFile]->labelfont(iMAP->boldFont);

		// change title of main window
		char topBarTitle [FILENAME_MAX];
		sprintf(topBarTitle,"InferenceMAP - %s",iMAP->fileNameStorage[iMAP->currentFile]);
		iMAP->mainWindow->label(topBarTitle);

		// adjust visualization widgets
		iMAP->iOffsetSlider->value(file->iOffset);
		iMAP->cMaxSlider->value(file->cMax);
		iMAP->cMinSlider->value(file->cMin);
		iMAP->drawTrajectoriesButton->value(file->drawTrajectories);
		iMAP->animateTrajectoriesButton->value(file->animateTrajectories);
		iMAP->fileInfo->update();

		iMAP->startIntervalSlider->value(file->startInterval);
		iMAP->startIntervalSlider->bounds(file->tMin,file->tMax);
		iMAP->endIntervalSlider->value(file->endInterval);
		iMAP->endIntervalSlider->bounds(file->tMin,file->tMax);

		/* Adjust Initial Visualization Settings */
		// assure ROI described in file is centered and zoomed correctly when loaded
		file->orthoLimit = 0.6*file->maxRange;

		Fl::redraw();
	}
}

void startIntervalCallback(Slider*w,void*v) {
	if (iMAP->fileNumber > 0) {
		if (file->startInterval > file->endInterval) {
			file->startInterval = file->endInterval;
			w->value(file->startInterval);
		} else {
			file->startInterval = w->value();
		}
		file->updateIntervalDisplay();
	}
}

void endIntervalCallback(Slider*w,void*v) {
	if (iMAP->fileNumber > 0) {
		if (file->endInterval < file->startInterval) {
			file->endInterval = file->startInterval;
			w->value(file->endInterval);
		} else {
			file->endInterval = w->value();
		}
		file->updateIntervalDisplay();
	}
}

void plotTypeCallback(Fl_Widget*,void*v) {
	char plotTypeString [5];
	sprintf(plotTypeString,"%i",v);
	int plotType = atoi(plotTypeString);

	if (iMAP->fileNumber > 0) {
		if (plotType == 1 && file->densityCalculated == false) { fl_alert("Density not calculated."); }
		else {
			file->plotType = plotType;
			file->plotChanged = true;
		}
	} else { fl_alert("No file open."); }

}

void exposureTimeButtonCallback(Fl_Widget*,void*) {
	iMAP->exposureTime = atof(iMAP->exposureTimeInput->value());
	iMAP->pixelSize = atof(iMAP->pixelSizeInput->value());
	iMAP->fileParametersWindow->hide();
}

void annotationsGuiCallback(Fl_Widget*w,void*) {
	if (iMAP->fileNumber > 0) {
		if (iMAP->annotationsGui != NULL) {
			iMAP->annotationsGui->show();
		} else {
			iMAP->annotationsGui = new AnnotationsGui();
			iMAP->annotationsGui->show();
		}
	} else {
		fl_alert("No file open.");
	}
}

void iOffsetSliderCallback(Slider*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->iOffset = w->value();
		file->plotChanged = true;
		file->updateIntervalDisplay();
//		file->updateAlpha(w->value());
	}
	Fl::redraw();
}

void cMinSliderCallback(Slider*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->cMin = w->value();
		file->plotChanged = true;
		iMAP->overlayAdjustment = true;
		if (file->squareMeshOverlay) {
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		Fl::redraw();
	}
}

void cMaxSliderCallback(Slider*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->cMax = w->value();
		file->plotChanged = true;
		iMAP->overlayAdjustment = true;
		if (file->squareMeshOverlay) {
			if (file->squareMeshGui->landscapeButton->value()) {
				file->squareMeshGui->landscapeButton->do_callback();
			}
		}
		Fl::redraw();
	}
}

void trajectoriesButtonCallback(Fl_Check_Button*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->drawTrajectories = w->value();
	}
}

void localizationsButtonCallback(Fl_Check_Button*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->drawLocalizations = w->value();
	}
}

void animateButtonCallback(Fl_Check_Button*w,void*v) {
	if (iMAP->fileNumber > 0) {
		file->animateTrajectories = w->value();
		file->animationIndex = 0;

		iMAP->startIntervalSlider->value(file->tMin);
		iMAP->endIntervalSlider->value(file->tMax);

		if (w->value()) {
			iMAP->animateAccumulationButton->activate();
			iMAP->animateDelayedButton->activate();
		} else {
			iMAP->animateAccumulationButton->deactivate();
			iMAP->animateDelayedButton->deactivate();
		}
	}
}

void animateAccumulationCallback(Fl_Check_Button*w,void*v) {
	if (iMAP->fileNumber > 0) {
		if (w->value()) {
			iMAP->animateDelayedButton->value(0);
			iMAP->animateDelayedStepsSlider->deactivate();
		} else {
			iMAP->animateDelayedButton->value(1);
			iMAP->animateDelayedStepsSlider->activate();
		}
	}
}

void animateDelayedCallback(Fl_Check_Button*w,void*v) {
	if (iMAP->fileNumber > 0) {
		if (w->value()) {
			iMAP->animateAccumulationButton->value(0);
			iMAP->animateDelayedStepsSlider->activate();
		} else {
			iMAP->animateAccumulationButton->value(1);
		}
	}
}

void aboutCallback(Fl_Widget*, void*) {
	fl_message("InferenceMAP v1.0\n\nMohamed El Beheiry\nJean-Baptiste Masson\n\nInstitut Curie\nInstitut Pasteur\n\nCopyright (c) 2015\n\nFor academic use only.");
}

void squareMeshingCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (iMAP->selectionButtonPressed && file->selection == NULL) {
			iMAP->customSelectionInferenceGui->selectionButton->value(0);
			iMAP->customSelectionInferenceGui->selectionButton->do_callback();
		}
		if (file->voronoiMeshGui == NULL && file->treeMeshGui == NULL) {
			if (file->squareMeshGui == NULL) {
				file->meshType = 0;
				file->squareMeshGui = new SquareMeshGui();
				file->squareMeshGui->show();
				file->squareMeshOverlay = true;
			} else {
				file->squareMeshGui->hide();
				delete file->squareMeshGui;
				file->meshType = 0;
				file->squareMeshGui = new SquareMeshGui();
				file->squareMeshGui->show();
				file->squareMeshOverlay = true;
			}
			iMAP->pauseCalculation = false;
		} else {
			fl_alert("Cannot have more than one meshing window open.");
		}
	} else {
		fl_alert("No file open.");
	}
}

void voronoiMeshingCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (iMAP->selectionButtonPressed && file->selection == NULL) {
			iMAP->customSelectionInferenceGui->selectionButton->value(0);
			iMAP->customSelectionInferenceGui->selectionButton->do_callback();
		}
		if (file->squareMeshGui == NULL && file->treeMeshGui == NULL) {
			if (file->voronoiMeshGui == NULL) {
				file->meshType = 1;
				file->voronoiMeshGui = new VoronoiMeshGui();
				file->voronoiMeshGui->show();
			} else {
				file->voronoiMeshGui->hide();
				delete file->voronoiMeshGui;
				file->meshType = 1;
				file->voronoiMeshGui = new VoronoiMeshGui();
				file->voronoiMeshGui->show();
			}
			iMAP->pauseCalculation = false;
		} else {
			fl_alert("Cannot have more than one meshing window open.");
		}
	} else {
		fl_alert("No file open.");
	}
}

void quadTreeMeshingCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (iMAP->selectionButtonPressed && file->selection == NULL) {
			iMAP->customSelectionInferenceGui->selectionButton->value(0);
			iMAP->customSelectionInferenceGui->selectionButton->do_callback();
		}
		if (file->voronoiMeshGui == NULL && file->squareMeshGui == NULL) {
			if (file->treeMeshGui == NULL) {
				file->meshType = 2;
				file->treeMeshGui = new TreeMeshGui();
				file->treeMeshGui->show();
			} else {
				file->treeMeshGui->hide();
				delete file->treeMeshGui;
				file->meshType = 2;
				file->treeMeshGui = new TreeMeshGui();
				file->treeMeshGui->show();
			}
			iMAP->pauseCalculation = false;
		} else {
			fl_alert("Cannot have more than one meshing window open.");
		}
	} else {
		fl_alert("No file open.");
	}
}

void singleTrajectoryInferenceCallback(Fl_Widget*,void*) {
	if (iMAP->fileNumber > 0) {
		if (file->numberOfFiles == 1) {
			if (iMAP->singleTrajectoryInferenceGui != NULL) { delete [] iMAP->singleTrajectoryInferenceGui; }
			file->meshType = -1;
			iMAP->singleTrajectoryInferenceGui = new SingleTrajectoryInferenceGui(file->fileName);
			iMAP->drawTrajectoriesButton->value(0);
			iMAP->drawTrajectoriesButton->do_callback();
		} else {
			fl_alert("Only available for single trajectories.");
		}
	} else {
		fl_alert("No files open.");
	}
}

void resizeButtonCallback(Fl_Button*b,void*) {
	Fl::screen_xywh(iMAP->screenDimensions[0],iMAP->screenDimensions[1],iMAP->screenDimensions[2],iMAP->screenDimensions[3],0);

	int w = iMAP->screenDimensions[2]-90;
	int h = iMAP->screenDimensions[3]-90;

	const float ar = 730.0f/1110.0f;

	if (w > 1110 && h > 730) {

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
		iMAP->customSelectionInferenceGui->resize(h-40,285-iMAP->macOffset,410,105);
		iMAP->densityGui->resize(h+120,h-335-iMAP->macOffset,250,45);

		// visualization widgets
		iMAP->intervalBox->resize(h-40,110-iMAP->macOffset,410,55);
		iMAP->startIntervalSlider->resize(h+50,115-iMAP->macOffset,315,20);
		iMAP->endIntervalSlider->resize(h+50,140-iMAP->macOffset,315,20);

		// visualization widgets
		iMAP->visualizationBox->resize(h+167,170-iMAP->macOffset,203,110);
		iMAP->trajectoriesBox->resize(h-40,170-iMAP->macOffset,202,110);
		iMAP->iOffsetSlider->resize(h+172, 185-iMAP->macOffset, 193, 20);
		iMAP->cMinSlider->resize(h+172, 220-iMAP->macOffset, 193, 20);
		iMAP->cMaxSlider->resize(h+172, 255-iMAP->macOffset, 193, 20);
		iMAP->drawTrajectoriesButton->resize(h-35,170-iMAP->macOffset,190,20);
		iMAP->localizationsButton->resize(h-35,190-iMAP->macOffset,190,20);
		iMAP->animateTrajectoriesButton->resize(h-35,210-iMAP->macOffset,190,20);
		iMAP->animateAccumulationButton->resize(h-20,225-iMAP->macOffset,175,20);
		iMAP->animateDelayedButton->resize(h-20,240-iMAP->macOffset,180,20);
		iMAP->animateDelayedStepsSlider->resize(h+45, 242-iMAP->macOffset, 112, 15);
		iMAP->fpsSlider->resize(h+45,260-iMAP->macOffset,112,15);

		iMAP->fullScreenButton->resize(w-35,h-35-iMAP->macOffset,25,25);

		iMAP->colorbar->resize(h-40,h-30-iMAP->macOffset,50,20);
		iMAP->colorbarCanvas->resize(h-45,650-iMAP->macOffset,145,335);

		iMAP->textDisplay->resize(h+121,h-284-iMAP->macOffset,248,243);
		iMAP->textDisplayBox->resize(h+120,h-285-iMAP->macOffset,250,245);

		iMAP->updateDisplayButton->resize(h+120,h-30-iMAP->macOffset,120,20);

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
		iMAP->customSelectionInferenceGui->resize(690,285-iMAP->macOffset,410,105);
		iMAP->densityGui->resize(850,395-iMAP->macOffset,250,45);

		// visualization widgets
		iMAP->intervalBox->resize(690,110-iMAP->macOffset,410,55);
		iMAP->startIntervalSlider->resize(780,115-iMAP->macOffset,315,20);
		iMAP->endIntervalSlider->resize(780,140-iMAP->macOffset,315,20);

		iMAP->visualizationBox->resize(897,170-iMAP->macOffset,203,110);
		iMAP->trajectoriesBox->resize(690,170-iMAP->macOffset,202,110);

		iMAP->iOffsetSlider->resize(902, 185-iMAP->macOffset, 193, 20);
		iMAP->cMinSlider->resize(902, 220-iMAP->macOffset, 193, 20);
		iMAP->cMaxSlider->resize(902, 255-iMAP->macOffset, 193, 20);
		iMAP->drawTrajectoriesButton->resize(693,170-iMAP->macOffset,190,20);
		iMAP->localizationsButton->resize(693,190-iMAP->macOffset,190,20);
		iMAP->animateTrajectoriesButton->resize(693,210-iMAP->macOffset,190,20);
		iMAP->animateAccumulationButton->resize(710,225-iMAP->macOffset,175,20);
		iMAP->animateDelayedButton->resize(710,240-iMAP->macOffset,180,20);
		iMAP->animateDelayedStepsSlider->resize(775, 242-iMAP->macOffset, 112, 15);
		iMAP->fpsSlider->resize(775,260-iMAP->macOffset,112,15);

		iMAP->fullScreenButton->resize(1075,695-iMAP->macOffset,25,25);

		iMAP->colorbar->resize(690,700-iMAP->macOffset,50,20);
		iMAP->colorbarCanvas->resize(685,390-iMAP->macOffset,145,335);

		iMAP->textDisplay->resize(851,446-iMAP->macOffset,248,243);
		iMAP->textDisplayBox->resize(850,445-iMAP->macOffset,250,245);

		iMAP->updateDisplayButton->resize(850,700-iMAP->macOffset,120,20);

		iMAP->screenDimensions[2] = 1110;
		iMAP->screenDimensions[3] = 730;

		iMAP->tabBarWidth = 680;
		// resize widths of any existing tabs
		for (int d = 0; d < iMAP->fileNumber; d++) {
			iMAP->fileObjectArray[d]->tabButton->resize(iMAP->tabBarWidth*(d)/(iMAP->fileNumber),30-iMAP->macOffset,iMAP->tabBarWidth/(iMAP->fileNumber),20);
		}
	}
	
	Fl::redraw();
	iMAP->showLogo = true;
}

void saveScreenCallback(Fl_Widget*,void*) {
    if (iMAP->fileNumber > 0) {
            Fl_Native_File_Chooser native;
            native.title("Save screen as");
            native.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
            native.filter("TIFF Files\t*.tif\n");
            native.options(Fl_Native_File_Chooser::SAVEAS_CONFIRM);
            native.preset_file(file->fileName);

            char capture[FILENAME_MAX];

            if (native.show() == 0) {
				strcpy(capture, native.filename());
				captureScreen((char*)native.filename());

				// save colorbar image in separate tiff file
				writetiffOffscreen(capture,iMAP->colorbarCanvas,"_colorbar");
            }
    }
    else {
		fl_alert("No files open.");
    }
    Fl::redraw();
    Fl::redraw();
}

void overlayTiffCallback(Fl_Widget*,void*) {
    if (iMAP->fileNumber > 0) {

		if (iMAP->overlayTiffGui == NULL) {
			Fl_Native_File_Chooser native;
			native.title("Select TIFF File");
			native.type(Fl_Native_File_Chooser::BROWSE_FILE);
			native.filter("TIFF Files\t*.tif\n");

			char tiffFile[FILENAME_MAX];

			if (native.show() == 0) {
				strcpy(tiffFile, native.filename());
        		iMAP->tiffSequence = new TiffSequence(tiffFile);
        		iMAP->overlayTiffGui = new OverlayTiffGui();

        		iMAP->tiffOverlay = true;
			}
		} else {
			iMAP->overlayTiffGui->show();
		}
    }
    else {
		fl_alert("No files open.");
    }
}

void updateDisplayButtonCallback(Fl_Check_Button*w,void*) {
	if (w->value()) {
		iMAP->updateDisplay = true;
		iMAP->textDisplay->activate();
	}
	else {
		iMAP->updateDisplay = false;
		iMAP->textDisplay->deactivate();
	}
}

int main(int argc, char* argv[]) {
	Fl::lock();
//	Fl::scheme("gtk+");

	// Menu System
	const Fl_Menu_Item menuItems[] = {
		{"&File",0,0,0,FL_SUBMENU},
	    	{"Open",0,0,0,FL_SUBMENU},
//	    		{"x y i t",FL_ALT+'o',(Fl_Callback*)fileOpenCallback,(int*)(0)},
	    		{"tr x y t",0,(Fl_Callback*)fileOpenCallback,(int*)(4)},
//	    		{"Multi CSV",0,(Fl_Callback*)fileOpenCallback,(int*)(1)},
//	    		{"SLIMfast",0,(Fl_Callback*)fileOpenCallback,(int*)(2)},
//	    		{"SPTrack",0,(Fl_Callback*)fileOpenCallback,(int*)(3)},
	    		{"x y t",0,(Fl_Callback*)fileOpenCallback,(int*)(5)},
//	    		{"tr x y t batch",0,(Fl_Callback*)batchVoronoiCallback,(int*)(4)},
//	    		{"tr x y t single batch",0,(Fl_Callback*)batchSingleTrajectoryCallback,(int*)(4)},
//	    		{"tr x y z t",0,(Fl_Callback*)fileOpenCallback,(int*)(6)},
			{0},
			{"Examples",0,0,0,FL_SUBMENU},
				{"Diffusion Example (Simulation)",0,(Fl_Callback*)fileOpenExampleCallback,(int*)(0)},
				{"Potential Example (Simulation)",0,(Fl_Callback*)fileOpenExampleCallback,(int*)(1)},
				{"Toxin Receptor in Lipid Raft",0,(Fl_Callback*)fileOpenExampleCallback,(int*)(2)},
				{"Glycine Receptors on Mouse Hippocampal Neuron",0,(Fl_Callback*)fileOpenExampleCallback,(int*)(3)},
			{0},
			{"Save Screen",0,saveScreenCallback},
			{"Quit", FL_ALT+'q', mainWindowCallback},
	    {0},
		{"&Plot",0,0,0,FL_SUBMENU},
			{"Frame",0,(Fl_Callback*)plotTypeCallback,(int*)(0)},
			{"Density",0,(Fl_Callback*)plotTypeCallback,(int*)(1)},
		{0},
		{"&Inference",0,0,0,FL_SUBMENU},
			{"Single Trajectory",0,(Fl_Callback*)singleTrajectoryInferenceCallback},
//			{"Batch Trajectory",0,(Fl_Callback*)batchTrajectoryInferenceCallback},
			{"Meshing",0,0,0,FL_SUBMENU},
				{"Square",0,(Fl_Callback*)squareMeshingCallback},
				{"Voronoi Tessellation",0,(Fl_Callback*)voronoiMeshingCallback},
				{"Quad-Tree",0,(Fl_Callback*)quadTreeMeshingCallback},
			{0},
		{0},
		{"&View",0,0,0,FL_SUBMENU},
			{"Annotations",0,(Fl_Callback*)annotationsGuiCallback},
			{"Overlay TIFF",0,overlayTiffCallback},
//			{"Stereo Vision",0,stereoVisionCallback},
//			{"Movie Maker",0,movieMakerCallback},
			{"White Background",0,whiteBackgroundCallback},
		{0},
		{"&Help",0,0,0,FL_SUBMENU},
			{"License",0,(Fl_Callback*)licenseCallback},
			{"Manual",0,(Fl_Callback*)userManualCallback},
			{"About",0,(Fl_Callback*)aboutCallback},
		{0},
		{0},
	};

	// setup GL extensions
	setupGLExtensions();

	iMAP = new Globals();

	#ifdef _WIN32

		iMAP->normalFont = FL_HELVETICA;
		iMAP->boldFont = FL_HELVETICA_BOLD;

		// define group to contain all 'widget' windows
		iMAP->mainGroup = new Fl_Group(0,0,iMAP->screenDimensions[2],iMAP->screenDimensions[3],"group");
		iMAP->mainGroup->begin();
		// create main window
		iMAP->mainWindow = new Fl_Window(iMAP->screenDimensions[2],iMAP->screenDimensions[3], "InferenceMAP");
		iMAP->mainWindow->begin();
		iMAP->mainWindow->color(fl_darker(iMAP->bgColor));
		iMAP->mainWindow->callback(mainWindowCallback);
	    // create menubar
		iMAP->menubar = new Fl_Menu_Bar(0,0,iMAP->screenDimensions[2],30,"menu");
		iMAP->menubar->labelsize(12);
		iMAP->menubar->textsize(12);
		iMAP->menubar->menu(menuItems);

		iMAP->macOffset = 0;

		char fontPath[FILENAME_MAX];
		const int pathLen = strlen(argv[0]) - 16; // strlen("InferenceMAP") = 12

		for (int r = 0; r < pathLen; r++) {
			fontPath[r] = argv[0][r];
		}
		fontPath[pathLen] = '\0';

		strcpy(iMAP->cwd,fontPath);
		strcat(fontPath,"/Resources/helvetica48.txf");

		iMAP->fontTex = txfLoadFont(fontPath);

	#elif __APPLE__
		iMAP->macOffset = 30;

		iMAP->normalFont = FL_SCREEN;
		iMAP->boldFont = FL_SCREEN_BOLD;
		
		iMAP->bgColor = fl_darker(iMAP->bgColor);

		// define group to contain all 'widget' windows
		iMAP->mainGroup = new Fl_Group(0,0,iMAP->screenDimensions[2],iMAP->screenDimensions[3]-iMAP->macOffset,"group");
		iMAP->mainGroup->begin();
		// create main window
		iMAP->mainWindow = new Fl_Window(iMAP->screenDimensions[2],iMAP->screenDimensions[3]-iMAP->macOffset, "InferenceMAP");
		iMAP->mainWindow->begin();
		iMAP->mainWindow->color(iMAP->bgColor);
		iMAP->mainWindow->labelfont(iMAP->boldFont);
		iMAP->mainWindow->callback(mainWindowCallback);

		Fl_Sys_Menu_Bar *macMenu = new Fl_Sys_Menu_Bar(0,0,0,0,"menu");
		macMenu->menu(menuItems);
		macMenu->show();
		fl_mac_set_about((Fl_Callback*)aboutCallback,0,0);

		char fontPath[FILENAME_MAX];
		const int pathLen = strlen(argv[0]) - 12; // strlen("InferenceMAP") = 12

		for (int r = 0; r < pathLen; r++) {
			fontPath[r] = argv[0][r];
		}
		fontPath[pathLen] = '\0';

		strcpy(iMAP->cwd,fontPath);
		strcat(fontPath,"helvetica48.txf");

		iMAP->fontTex = txfLoadFont(fontPath);
	#endif

	// create gl window
	iMAP->glWindow = new iMAPGlWindow(0,50-iMAP->macOffset,680,680,"OpenGL Window");

	iMAP->fileInfo = new FileInfo(690,40-iMAP->macOffset,410,65);
	iMAP->fileInfo->deactivate();
	iMAP->fileInfo->show();

	iMAP->customSelectionInferenceGui = new CustomSelectionInferenceGui(690,110-iMAP->macOffset,410,105);
	iMAP->customSelectionInferenceGui->deactivate();
	iMAP->customSelectionInferenceGui->show();

	iMAP->densityGui = new DensityGui(850,395-iMAP->macOffset,250,45);
	iMAP->densityGui->deactivate();
	iMAP->densityGui->show();

	// visualization widgets
	iMAP->visualizationBox = new Fl_Box(895,220-iMAP->macOffset,410,110);
	iMAP->visualizationBox->box(FL_BORDER_FRAME);
	iMAP->visualizationBox->show();

	iMAP->trajectoriesBox = new Fl_Box(690,220-iMAP->macOffset,205,110);
	iMAP->trajectoriesBox->box(FL_BORDER_FRAME);
	iMAP->trajectoriesBox->show();

	iMAP->iOffsetSlider = new Slider(900, 235-iMAP->macOffset, 195, 20, "Intensity Offset");
	iMAP->iOffsetSlider->precision(3);
	iMAP->iOffsetSlider->labelsize(12);
	iMAP->iOffsetSlider->labelcolor(FL_WHITE);
	iMAP->iOffsetSlider->labelfont(iMAP->boldFont);
	iMAP->iOffsetSlider->value(1.0f);
	iMAP->iOffsetSlider->bounds(-1.0f,1.0f);
	iMAP->iOffsetSlider->callback((Fl_Callback*)iOffsetSliderCallback);
	iMAP->iOffsetSlider->deactivate();
	iMAP->iOffsetSlider->show();

	iMAP->cMinSlider = new Slider(900, 270-iMAP->macOffset, 195, 20, "cMin");
	iMAP->cMinSlider->precision(3);
	iMAP->cMinSlider->labelsize(12);
	iMAP->cMinSlider->labelcolor(FL_WHITE);
	iMAP->cMinSlider->labelfont(iMAP->boldFont);
	iMAP->cMinSlider->value(0.0f);
	iMAP->cMinSlider->bounds(0.0f,1.0f);
	iMAP->cMinSlider->callback((Fl_Callback*)cMinSliderCallback);
	iMAP->cMinSlider->deactivate();
	iMAP->cMinSlider->show();

	iMAP->cMaxSlider = new Slider(900, 305-iMAP->macOffset, 195, 20, "cMax");
	iMAP->cMaxSlider->precision(3);
	iMAP->cMaxSlider->labelcolor(FL_WHITE);
	iMAP->cMaxSlider->labelsize(12);
	iMAP->cMaxSlider->labelfont(iMAP->boldFont);
	iMAP->cMaxSlider->value(1.0f);
	iMAP->cMaxSlider->bounds(0.0f,1.0f);
	iMAP->cMaxSlider->callback((Fl_Callback*)cMaxSliderCallback);
	iMAP->cMaxSlider->deactivate();
	iMAP->cMaxSlider->show();

	iMAP->drawTrajectoriesButton = new Fl_Check_Button(695,220-iMAP->macOffset,190,20,"Draw Trajectories");
	iMAP->drawTrajectoriesButton->labelsize(12);
	iMAP->drawTrajectoriesButton->color(FL_BLACK);
	iMAP->drawTrajectoriesButton->labelcolor(FL_WHITE);
	iMAP->drawTrajectoriesButton->labelfont(iMAP->boldFont);
	iMAP->drawTrajectoriesButton->callback((Fl_Callback*)trajectoriesButtonCallback);
	iMAP->drawTrajectoriesButton->value(1);
	iMAP->drawTrajectoriesButton->deactivate();
	iMAP->drawTrajectoriesButton->show();

	iMAP->localizationsButton = new Fl_Check_Button(695,240-iMAP->macOffset,190,20,"Draw Localizations");
	iMAP->localizationsButton->labelsize(12);
	iMAP->localizationsButton->color(FL_BLACK);
	iMAP->localizationsButton->labelcolor(FL_WHITE);
	iMAP->localizationsButton->labelfont(iMAP->boldFont);
	iMAP->localizationsButton->callback((Fl_Callback*)localizationsButtonCallback);
	iMAP->localizationsButton->value(0);
	iMAP->localizationsButton->deactivate();
	iMAP->localizationsButton->show();

	iMAP->animateTrajectoriesButton = new Fl_Check_Button(695,260-iMAP->macOffset,190,20,"Animate Trajectories");
	iMAP->animateTrajectoriesButton->labelsize(12);
	iMAP->animateTrajectoriesButton->color(FL_BLACK);
	iMAP->animateTrajectoriesButton->labelcolor(FL_WHITE);
	iMAP->animateTrajectoriesButton->labelfont(iMAP->boldFont);
	iMAP->animateTrajectoriesButton->callback((Fl_Callback*)animateButtonCallback);
	iMAP->animateTrajectoriesButton->value(0);
	iMAP->animateTrajectoriesButton->deactivate();
	iMAP->animateTrajectoriesButton->show();

	iMAP->animateAccumulationButton = new Fl_Check_Button(710,275-iMAP->macOffset,175,20,"Accumulation");
	iMAP->animateAccumulationButton->labelsize(10);
	iMAP->animateAccumulationButton->color(FL_BLACK);
	iMAP->animateAccumulationButton->labelcolor(FL_WHITE);
	iMAP->animateAccumulationButton->labelfont(iMAP->boldFont);
	iMAP->animateAccumulationButton->callback((Fl_Callback*)animateAccumulationCallback);
	iMAP->animateAccumulationButton->value(1);
	iMAP->animateAccumulationButton->deactivate();
	iMAP->animateAccumulationButton->show();

	iMAP->animateDelayedButton = new Fl_Check_Button(710,290-iMAP->macOffset,180,20,"Interval");
	iMAP->animateDelayedButton->labelsize(10);
	iMAP->animateDelayedButton->color(FL_BLACK);
	iMAP->animateDelayedButton->labelcolor(FL_WHITE);
	iMAP->animateDelayedButton->labelfont(iMAP->boldFont);
	iMAP->animateDelayedButton->callback((Fl_Callback*)animateDelayedCallback);
	iMAP->animateDelayedButton->value(0);
	iMAP->animateDelayedButton->deactivate();
	iMAP->animateDelayedButton->show();

	iMAP->animateDelayedStepsSlider = new Slider(775, 290-iMAP->macOffset, 112, 15);
	iMAP->animateDelayedStepsSlider->precision(0);
	iMAP->animateDelayedStepsSlider->labelcolor(FL_WHITE);
	iMAP->animateDelayedStepsSlider->labelsize(10);
	iMAP->animateDelayedStepsSlider->labelfont(iMAP->boldFont);
	iMAP->animateDelayedStepsSlider->value(20.0f);
	iMAP->animateDelayedStepsSlider->bounds(5.0f,100.0f);
	iMAP->animateDelayedStepsSlider->align(FL_ALIGN_LEFT);
	iMAP->animateDelayedStepsSlider->deactivate();
	iMAP->animateDelayedStepsSlider->show();

	iMAP->fpsSlider = new Slider(775,310-iMAP->macOffset,112,15,"FPS");
	iMAP->fpsSlider->precision(0);
	iMAP->fpsSlider->labelcolor(FL_WHITE);
	iMAP->fpsSlider->labelsize(10);
	iMAP->fpsSlider->labelfont(iMAP->boldFont);
	iMAP->fpsSlider->value(30.0);
	iMAP->fpsSlider->bounds(1,50);
	iMAP->fpsSlider->align(FL_ALIGN_LEFT);
	iMAP->fpsSlider->deactivate();
	iMAP->fpsSlider->show();

	iMAP->intervalBox = new Fl_Box(690,334-iMAP->macOffset,410,55);
	iMAP->intervalBox->box(FL_BORDER_FRAME);
	iMAP->intervalBox->show();

	iMAP->startIntervalSlider = new Slider(780,339-iMAP->macOffset,315,20,"Start Time [s]");
	iMAP->startIntervalSlider->precision(2);
	iMAP->startIntervalSlider->labelcolor(FL_WHITE);
	iMAP->startIntervalSlider->labelsize(10);
	iMAP->startIntervalSlider->labelfont(iMAP->boldFont);
	iMAP->startIntervalSlider->value(0.0f);
	iMAP->startIntervalSlider->bounds(0.0f,100.0f);
	iMAP->startIntervalSlider->align(FL_ALIGN_LEFT);
	iMAP->startIntervalSlider->when(FL_WHEN_RELEASE);
	iMAP->startIntervalSlider->callback((Fl_Callback*)startIntervalCallback);
	iMAP->startIntervalSlider->deactivate();
	iMAP->startIntervalSlider->show();

	iMAP->endIntervalSlider = new Slider(780,364-iMAP->macOffset,315,20,"End Time [s]");
	iMAP->endIntervalSlider->precision(2);
	iMAP->endIntervalSlider->labelcolor(FL_WHITE);
	iMAP->endIntervalSlider->labelsize(10);
	iMAP->endIntervalSlider->labelfont(iMAP->boldFont);
	iMAP->endIntervalSlider->value(0.0f);
	iMAP->endIntervalSlider->bounds(0.0f,100.0f);
	iMAP->endIntervalSlider->align(FL_ALIGN_LEFT);
	iMAP->endIntervalSlider->when(FL_WHEN_RELEASE);
	iMAP->endIntervalSlider->callback((Fl_Callback*)endIntervalCallback);
	iMAP->endIntervalSlider->deactivate();
	iMAP->endIntervalSlider->show();

	iMAP->colorbar = new ColorBar(690,700-iMAP->macOffset,50,20);
	iMAP->colorbar->deactivate();
	iMAP->colorbar->show();

	iMAP->fullScreenButton = new Fl_Button(1075,695-iMAP->macOffset,25,25);
	iMAP->fullScreenButton->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
	iMAP->fullScreenButton->type(FL_TOGGLE_BUTTON);
	iMAP->fullScreenButton->box(FL_BORDER_FRAME);
	iMAP->fullScreenButton->labelsize(12);
	iMAP->fullScreenButton->callback((Fl_Callback*)resizeButtonCallback);
	iMAP->fullScreenButton->show();

	iMAP->colorbarCanvas = new Fl_Box(685,390-iMAP->macOffset,145,335);
	iMAP->colorbarCanvas->show();

	iMAP->textBuffer = new Fl_Text_Buffer();

	iMAP->textDisplay = new Fl_Text_Display(851,316-iMAP->macOffset,248,258);
	iMAP->textDisplay->textfont(iMAP->boldFont);
	iMAP->textDisplay->textsize(12);
	iMAP->textDisplay->textcolor(FL_LIGHT3);
	iMAP->textDisplay->color(iMAP->bgColor);
	iMAP->textDisplay->box(FL_BORDER_FRAME);
	iMAP->textDisplay->buffer(iMAP->textBuffer);
	iMAP->textBuffer->append(iMAP->outputStream);
	iMAP->textDisplay->show_insert_position();
	iMAP->textDisplay->show();

	iMAP->textDisplayBox = new Fl_Box(850,305-iMAP->macOffset,250,300);
	iMAP->textDisplayBox->box(FL_BORDER_FRAME);
	iMAP->textDisplayBox->show();

	iMAP->updateDisplayButton = new Fl_Check_Button(850,700-iMAP->macOffset,120,20,"Update Display");
	iMAP->updateDisplayButton->value(iMAP->updateDisplay);
	iMAP->updateDisplayButton->labelsize(12);
	iMAP->updateDisplayButton->color(iMAP->bgColor);
	iMAP->updateDisplayButton->labelfont(iMAP->boldFont);
	iMAP->updateDisplayButton->labelcolor(FL_WHITE);
	iMAP->updateDisplayButton->box(FL_BORDER_FRAME);
	iMAP->updateDisplayButton->callback((Fl_Callback*)updateDisplayButtonCallback);
	iMAP->updateDisplayButton->show();

	iMAP->mainGroup->end();
	iMAP->mainWindow->end();
	iMAP->mainWindow->show(argc, argv);
	iMAP->glWindow->show();

	// full screen on startup
	iMAP->fullScreenButton->value(1);
	iMAP->fullScreenButton->do_callback();

	Fl::redraw();

	return(Fl::run());
}
