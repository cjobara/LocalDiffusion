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
#include "optimization.h"

extern Globals *iMAP;

File::File(int files, char **list, int fileType) {

	// GUIs
	squareMeshGui = NULL;
	squareMeshOverlay = false;
	voronoiMeshGui = NULL;
	voronoiMeshOverlay = false;
	treeMeshGui = NULL;
	treeMeshOverlay = false;

	this->exposureTime = iMAP->exposureTime/1000.0;
	this->landscapePlot = false;
	this->colormapType = 10; // jet by default
	this->fileLoaded = false;
	this->plotChanged = true;
	this->inferred = false;
	this->averageDx = 0.0;
	this->averageDy = 0.0;
	this->averageDt = 0.0;
	this->tabButton = NULL;
	this->trackNumber = 0;
	this->optimizationMode = 0;
	this->trackLengths = NULL;
	this->trackFrames = NULL;
	this->meshType = 0;
	this->flipColormap = 0;
	this->treeMesh = NULL;
	this->simulation = NULL;
	this->selection = NULL;
	this->xRotate = 0.0;
	this->yRotate = 0.0;
	this->zRotate = 0.0;
	this->randomizedOptimizationTolerance = 0.001;
	this->drawRandomizedOptimizationZones = false;
	this->partitionMatrix = NULL;
	this->xParts = 0;
	this->yParts = 0;
	this->densityCalculated = false;
	this->file3D = false;
	this->xPointerHolder = NULL;
	this->yPointerHolder = NULL;
	this->iPointerHolder = NULL;
	this->tPointerHolder = NULL;
	this->nPointerHolder = NULL;

	this->optimizationFunction = 0;

	/* Optmization Modes
	 * 0 : Global Polynomial Potential
	 * 1 : DF Diffusion and Force
	 * 2 : DV Diffusion and Potential
	 */

	optimizationMode = 0;

	// assign filetype to class variable
	this->fileType = fileType;

	fileNameList = list;
	numberOfFiles = files;

	// truncate fileNameStorage to show only file (without path)
	int fileNameLength = strlen(fileNameList[0]);
	int truncatePoint = fileNameLength;

	#if __APPLE__
	while (fileNameList[0][truncatePoint] != '/') {
		truncatePoint--;
	}
	#elif _WIN32
	while (fileNameList[0][truncatePoint] != '\\') {
		truncatePoint--;
	}
	#endif

	// store path
	for (int p = 0; p < truncatePoint+1; p++) {
		filePath[p] = fileNameList[0][p];
	}
	filePath[truncatePoint+1] = '\0';

	if (files > 1) { sprintf(this->fileName,"%i Trajectories",this->numberOfFiles); }
	else {
		int p = 0;
		for (int t = truncatePoint+1; t < fileNameLength; t++) {
			this->fileName[p] = fileNameList[0][t];
			p++;
		}
		this->fileName[p] = '\0';
	}

	FILE *fileHandle;

	switch(fileType) {
	case 0:
		{
			// count total number of detections
			this->localizationCount = 0;

			// to store number of detections in each file
			localizationCountFile = new int[files];

			for (int v = 0; v < numberOfFiles; v++) {
				// read file
				fileHandle = fopen(fileNameList[v],"r");
				int ch;

				// error if no file read
				if (fileHandle == NULL) {
					printf("Error in opening file\n");
					return;
				}

				localizationCountFile[v] = 0;
				// determine number of lines in file
				do {
					ch = fgetc(fileHandle);
					if (ch=='\n') {
						localizationCountFile[v]++;
						localizationCount++;
					}
				} while( ch != EOF);

				fclose(fileHandle);

				// eliminate last line of file (zero line)
				localizationCountFile[v]--;
				localizationCount--;
			}

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int j = 0;
			for (int v = 0; v < numberOfFiles; v++) {

				// allocate memory for individual tracks
				tracks[v].n = localizationCountFile[v];
				tracks[v].d = new Detection[localizationCountFile[v]];
				tracks[v].xy = new double[2*localizationCountFile[v]];
				tracks[v].averageDt = 0.0;
				tracks[v].averageDx = 0.0;
				tracks[v].averageDy = 0.0;

				fopen(fileNameList[v],"r");

				colormap(tracks[v].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);

				for(int k = 0; k < localizationCountFile[v]; k++) {
					fscanf(fileHandle,"%lf\t%lf\t%lf\t%lf\t",
							&xPointer[j],
							&yPointer[j],
							&iPointer[j],
							&tPointer[j]);

					// tag with channel (file) number
					nPointer[j] = v;
					insideInterval[j] = true;

					// find max/min intensity
					if (iMax < iPointer[j]) { iMax = iPointer[j]; }
					if (iMin > iPointer[j]) { iMin = iPointer[j]; }
					// find max/min frame number
					if (tMax < tPointer[j]) { tMax = tPointer[j]; }
					if (tMin > tPointer[j]) { tMin = tPointer[j]; }
					// find max/min x,y coordinates;
					if (xMax<xPointer[j]) { xMax = xPointer[j]; }
					if (yMax<yPointer[j]) { yMax = yPointer[j]; }
					if (xMin>xPointer[j]) { xMin = xPointer[j]; }
					if (yMin>yPointer[j]) { yMin = yPointer[j]; }

					tracks[v].d[k].x = xPointer[j];
					tracks[v].d[k].y = yPointer[j];
					tracks[v].xy[2*k] = xPointer[j];
					tracks[v].xy[2*k+1] = yPointer[j];
					tracks[v].d[k].t = tPointer[j];
					tracks[v].d[k].i = iPointer[j];

					// total increment
					j++;
				}
				fclose(fileHandle);
			}

			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					tracks[v].averageDt += fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = tracks[h].d[hh].x;
					linesArray[4*l+1] = tracks[h].d[hh].y;
					linesArray[4*l+2] = tracks[h].d[hh+1].x;
					linesArray[4*l+3] = tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;
			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;

			iMAP->exposureTime = averageDt;
			break;
		}
	case 1: // Multi CSV file,
		{
			// count total number of detections
			this->localizationCount = 0;

			// read file
			fileHandle = fopen(fileNameList[0],"r");
			int ch;

			// error if no file read
			if (fileHandle == NULL) {
				printf("Error in opening file\n");
				return;
			}

			// determine number of lines in file
			do {
				ch = fgetc(fileHandle);
				if (ch=='\n') {
					localizationCount++;
				}
			} while( ch != EOF);

			fclose(fileHandle);

			// eliminate last line of file (zero line)
			localizationCount--;

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;
			this->nMin = 10000000;
			this->nMax = -10000000;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

			char comma;

			fopen(fileNameList[0],"r");

			for(int j = 0; j < localizationCount; j++) {
				fscanf(fileHandle,"%lf%c%lf%c%lf%c%i\t",
						&xPointer[j],
						&comma,
						&yPointer[j],
						&comma,
						&tPointer[j],
						&comma,
						&nPointer[j]);

				iPointer[j] = 1.0;

				insideInterval[j] = true;
				xPointer[j] /= 10.0;
				yPointer[j] /= 10.0;

				// find max/min frame number
				if (tMax < tPointer[j]) { tMax = tPointer[j]; }
				if (tMin > tPointer[j]) { tMin = tPointer[j]; }
				if (nMax < nPointer[j]) { nMax = nPointer[j]; }
				if (nMin > nPointer[j]) { nMin = nPointer[j]; }

				// find max/min x,y coordinates;
				if (xMax<xPointer[j]) { xMax = xPointer[j]; }
				if (yMax<yPointer[j]) { yMax = yPointer[j]; }
				if (xMin>xPointer[j]) { xMin = xPointer[j]; }
				if (yMin>yPointer[j]) { yMin = yPointer[j]; }

			}
			fclose(fileHandle);

			// count absolute number of files (trajectories numbers may not be numbered incrementally)
			numberOfFiles = 0;
			for (int bb = 0; bb < localizationCount-1; bb++) {
				if (nPointer[bb] != nPointer[bb+1]) { numberOfFiles++; }
			}
			numberOfFiles++;

			// for fast vertex rendering
			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// to store number of detections in each file
			localizationCountFile = new int[numberOfFiles];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int g = 0;
			int r = 0;
			int id = nPointer[0];

			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
				tracks[h].n = g;
				tracks[h].d = new Detection[g];
				tracks[h].xy = new double[2*g];
				tracks[h].averageDt = 0.0;
				tracks[h].averageDx = 0.0;
				tracks[h].averageDy = 0.0;
				localizationCountFile[h] = g;
				colormap(tracks[h].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);
			}

			r = 0;
			id = nPointer[0];
			// assign points to individual trajectories
			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					tracks[h].d[g].x = xPointer[r];
					tracks[h].d[g].y = yPointer[r];
					tracks[h].xy[2*g] = xPointer[r];
					tracks[h].xy[2*g+1] = yPointer[r];
					tracks[h].d[g].t = tPointer[r];
					tracks[h].d[g].i = iPointer[r];
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
			}

			exposureTime = 1000000.0;
			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					const double dt = fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDt += dt;
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);

					// assume smallest dt is exposure time
					if (exposureTime > dt) { exposureTime = dt; }
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = (float)tracks[h].d[hh].x;
					linesArray[4*l+1] = (float)tracks[h].d[hh].y;
					linesArray[4*l+2] = (float)tracks[h].d[hh+1].x;
					linesArray[4*l+3] = (float)tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;

			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;

			break;
		}
	case 2: // SLIMfast file,
		{
			// count total number of detections
			this->localizationCount = 0;

			// read file
			fileHandle = fopen(fileNameList[0],"r");
			int ch;

			// error if no file read
			if (fileHandle == NULL) {
				printf("Error in opening file\n");
				return;
			}

			// determine number of lines in file
			do {
				ch = fgetc(fileHandle);
				if (ch=='\n') {
					localizationCount++;
				}
			} while( ch != EOF);

			fclose(fileHandle);

			// eliminate last line of file (zero line)
//			localizationCount--;

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;
			this->nMin = 10000000;
			this->nMax = -10000000;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

			double *dummy = new double[1];
			int previousChannel = -1;
			double *currentChannel = new double[1];
			int actualChannel = 0;
			fopen(fileNameList[0],"r");

			iMAP->exposureTime /= 1000.0;

			for(int j = 0; j < localizationCount; j++) {
				fscanf(fileHandle,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t",
						&xPointer[j],
						&yPointer[j],
						&tPointer[j],
						&currentChannel[0],
						&iPointer[j],
						&dummy[0],
						&dummy[0],
						&dummy[0]);

				if (j == 0) {
					previousChannel = (int)currentChannel[0];
					actualChannel = 0;
				}

//				fprintf(stderr,"%i\n",previousChannel);

				if ((int)currentChannel[0] != previousChannel) {
//					fprintf(stderr,"%i\t%i\n",previousChannel,currentChannel[0]);
					actualChannel++;
					previousChannel = (int)currentChannel[0];
				}
				nPointer[j] = actualChannel;
				insideInterval[j] = true;

				tPointer[j] *= iMAP->exposureTime;
				xPointer[j] *= iMAP->pixelSize/1000.0;
				yPointer[j] *= iMAP->pixelSize/1000.0;

				// find max/min intensity
				if (iMax < iPointer[j]) { iMax = iPointer[j]; }
				if (iMin > iPointer[j]) { iMin = iPointer[j]; }

				// find max/min frame number
				if (tMax < tPointer[j]) { tMax = tPointer[j]; }
				if (tMin > tPointer[j]) { tMin = tPointer[j]; }
				if (nMax < nPointer[j]) { nMax = nPointer[j]; }
				if (nMin > nPointer[j]) { nMin = nPointer[j]; }

				// find max/min x,y coordinates;
				if (xMax<xPointer[j]) { xMax = xPointer[j]; }
				if (yMax<yPointer[j]) { yMax = yPointer[j]; }
				if (xMin>xPointer[j]) { xMin = xPointer[j]; }
				if (yMin>yPointer[j]) { yMin = yPointer[j]; }

			}
			fclose(fileHandle);

			// count absolute number of files (trajectories numbers may not be numbered incrementally)
			numberOfFiles = 0;
			for (int bb = 0; bb < localizationCount-1; bb++) {
				if (nPointer[bb] != nPointer[bb+1]) { numberOfFiles++; }
			}
			numberOfFiles++;

			// for fast vertex rendering
			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// to store number of detections in each file
			localizationCountFile = new int[numberOfFiles];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int g = 0;
			int r = 0;
			int id = nPointer[0];

			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
				tracks[h].n = g;
				tracks[h].d = new Detection[g];
				tracks[h].xy = new double[2*g];
				tracks[h].averageDt = 0.0;
				tracks[h].averageDx = 0.0;
				tracks[h].averageDy = 0.0;
				localizationCountFile[h] = g;
				colormap(tracks[h].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);
			}

			r = 0;
			id = nPointer[0];
			// assign points to individual trajectories
			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					tracks[h].d[g].x = xPointer[r];
					tracks[h].d[g].y = yPointer[r];
					tracks[h].xy[2*g] = xPointer[r];
					tracks[h].xy[2*g+1] = yPointer[r];
					tracks[h].d[g].t = tPointer[r];
					tracks[h].d[g].i = iPointer[r];
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
			}


			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					tracks[v].averageDt += fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = (float)tracks[h].d[hh].x;
					linesArray[4*l+1] = (float)tracks[h].d[hh].y;
					linesArray[4*l+2] = (float)tracks[h].d[hh+1].x;
					linesArray[4*l+3] = (float)tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;

			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;

			break;
		}
	case 3: // SPTrack file,
		{
			// count total number of detections
			this->localizationCount = 0;

			int spaceCount = 0;
			bool connected = false;

			// read file
			fileHandle = fopen(fileNameList[0],"r");
			int ch;

			// error if no file read
			if (fileHandle == NULL) {
				printf("Error in opening file\n");
				return;
			}

			// count number of spaces (15: connected trajectories, 18: unconnected trajectories)

			do {
				ch = fgetc(fileHandle);
				if (ch==' ') { spaceCount++; }
			} while (ch != '\n');

			if (spaceCount == 15) { connected = true; }

			rewind(fileHandle);

			// determine number of lines in file
			do {
				ch = fgetc(fileHandle);
				if (ch=='\n') {
					localizationCount++;
				}
			} while ( ch != EOF);

			fclose(fileHandle);

			// eliminate last line of file (zero line)
			localizationCount--;

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;
			this->nMin = 10000000;
			this->nMax = -10000000;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

			double *dummy = new double[1];
			int previousChannel = -1;
			double *currentChannel = new double[1];
			int actualChannel = 0;
			fopen(fileNameList[0],"r");

			for(int j = 0; j < localizationCount; j++) {

				if (connected) {
					fscanf(fileHandle,"%lf\t%lf\t%lf\t%lf\t%lf\t",
							&currentChannel[0],
							&tPointer[j],
							&xPointer[j],
							&yPointer[j],
							&dummy[0]);
				} else {
					fscanf(fileHandle,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
							&currentChannel[0],
							&tPointer[j],
							&xPointer[j],
							&yPointer[j],
							&dummy[0],
							&dummy[0]);
				}

				if (j == 0) {
					previousChannel = (int)currentChannel[0];
					actualChannel = 0;
				}

				if ((int)currentChannel[0] != previousChannel) {
//					fprintf(stderr,"%i\t%i\n",previousChannel,currentChannel[0]);
					actualChannel++;
					previousChannel = (int)currentChannel[0];
				}

				iPointer[j] = 1.0;
				insideInterval[j] = true;

				nPointer[j] = actualChannel;

				tPointer[j] *= iMAP->exposureTime/1000.0;
				xPointer[j] *= iMAP->pixelSize/1000.0;
				yPointer[j] *= iMAP->pixelSize/1000.0;

				// find max/min intensity
				if (iMax < iPointer[j]) { iMax = iPointer[j]; }
				if (iMin > iPointer[j]) { iMin = iPointer[j]; }

				// find max/min frame number
				if (tMax < tPointer[j]) { tMax = tPointer[j]; }
				if (tMin > tPointer[j]) { tMin = tPointer[j]; }
				if (nMax < nPointer[j]) { nMax = nPointer[j]; }
				if (nMin > nPointer[j]) { nMin = nPointer[j]; }

				// find max/min x,y coordinates;
				if (xMax<xPointer[j]) { xMax = xPointer[j]; }
				if (yMax<yPointer[j]) { yMax = yPointer[j]; }
				if (xMin>xPointer[j]) { xMin = xPointer[j]; }
				if (yMin>yPointer[j]) { yMin = yPointer[j]; }

			}
			fclose(fileHandle);

			// count absolute number of files (trajectories numbers may not be numbered incrementally)
			numberOfFiles = 0;
			for (int bb = 0; bb < localizationCount-1; bb++) {
				if (nPointer[bb] != nPointer[bb+1]) { numberOfFiles++; }
			}
			numberOfFiles++;

			// for fast vertex rendering
			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// to store number of detections in each file
			localizationCountFile = new int[numberOfFiles];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int g = 0;
			int r = 0;
			int id = nPointer[0];

			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
				tracks[h].n = g;
				tracks[h].d = new Detection[g];
				tracks[h].xy = new double[2*g];
				tracks[h].averageDt = 0.0;
				tracks[h].averageDx = 0.0;
				tracks[h].averageDy = 0.0;
				localizationCountFile[h] = g;
				colormap(tracks[h].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);
			}

			r = 0;
			id = nPointer[0];
			// assign points to individual trajectories
			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					tracks[h].d[g].x = xPointer[r];
					tracks[h].d[g].y = yPointer[r];
					tracks[h].xy[2*g] = xPointer[r];
					tracks[h].xy[2*g+1] = yPointer[r];
					tracks[h].d[g].t = tPointer[r];
					tracks[h].d[g].i = iPointer[r];
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
			}


			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					tracks[v].averageDt += fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = (float)tracks[h].d[hh].x;
					linesArray[4*l+1] = (float)tracks[h].d[hh].y;
					linesArray[4*l+2] = (float)tracks[h].d[hh+1].x;
					linesArray[4*l+3] = (float)tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;

			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;
			iMAP->exposureTime /= 1000.0;

			break;
		}

	case 4: // tr x y t file,
		{
			// count total number of detections
			this->localizationCount = 0;

			// read file
			fileHandle = fopen(fileNameList[0],"r");
			int ch;

			// error if no file read
			if (fileHandle == NULL) {
				printf("Error in opening file\n");
				return;
			}

			// determine number of lines in file
			do {
				ch = fgetc(fileHandle);
				if (ch=='\n') {
					localizationCount++;
				}
			} while( ch != EOF);

			fclose(fileHandle);

			// eliminate last line of file (zero line)
//			localizationCount--;

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;
			this->nMin = 10000000;
			this->nMax = -10000000;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

//			this->linesColorArray = new float[];
			fopen(fileNameList[0],"r");

			for(int j = 0; j < localizationCount; j++) {
				fscanf(fileHandle,"%i\t%lf\t%lf\t%lf\t",
						&nPointer[j],
						&xPointer[j],
						&yPointer[j],
						&tPointer[j]);

//				tFloat[j] = (float)tPointer[j];
				iPointer[j] = 1.0;
				insideInterval[j] = true;

				// find max/min frame number
				if (tMax < tPointer[j]) { tMax = tPointer[j]; }
				if (tMin > tPointer[j]) { tMin = tPointer[j]; }
				if (nMax < nPointer[j]) { nMax = nPointer[j]; }
				if (nMin > nPointer[j]) { nMin = nPointer[j]; }

				// find max/min x,y coordinates;
				if (xMax < xPointer[j]) { xMax = xPointer[j]; }
				if (yMax < yPointer[j]) { yMax = yPointer[j]; }
				if (xMin > xPointer[j]) { xMin = xPointer[j]; }
				if (yMin > yPointer[j]) { yMin = yPointer[j]; }
			}
			fclose(fileHandle);

			// count absolute number of files (trajectories numbers may not be numbered incrementally)
			numberOfFiles = 0;
			for (int bb = 0; bb < localizationCount-1; bb++) {
				if (nPointer[bb] != nPointer[bb+1]) { numberOfFiles++; }
			}
			numberOfFiles++;

			// for fast vertex rendering
			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// to store number of detections in each file
			localizationCountFile = new int[numberOfFiles];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int g = 0;
			int r = 0;
			int id = nPointer[0];

			int maxLength = -10000;
			int minLength = 10000;

			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					g++;
					r++;
				}
//				do {
//					g++;
//					r++;
//				} while (id == nPointer[r] && r < localizationCount);

				if ( r >= localizationCount) { break; }
				id = nPointer[r];

				if (g == 1) { g++; } // to correct for issue where two-point trajectory is at very end of file

				tracks[h].n = g;
				tracks[h].d = new Detection[g];
				tracks[h].xy = new double[2*g];
				tracks[h].averageDt = 0.0;
				tracks[h].averageDx = 0.0;
				tracks[h].averageDy = 0.0;
				localizationCountFile[h] = g;
				colormap(tracks[h].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);

				if (g > maxLength) { maxLength = g; }
				if (g < minLength) { minLength = g; }
			}

			r = 0;
			id = nPointer[0];
			// assign points to individual trajectories
			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					tracks[h].d[g].x = xPointer[r];
					tracks[h].d[g].y = yPointer[r];
					tracks[h].xy[2*g] = xPointer[r];
					tracks[h].xy[2*g+1] = yPointer[r];
					tracks[h].d[g].t = tPointer[r];
					tracks[h].d[g].i = iPointer[r];
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
			}

//			fprintf(stderr,"[%i,%i]\n",minLength,maxLength);

			exposureTime = 1000000.0;
			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					const double dt = fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDt += dt;
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);

					// assume smallest dt is exposure time
					if (exposureTime > dt) { exposureTime = dt; }
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

//			fprintf(stderr,"%i:\t%f\t%f\t%f\n",numberOfFiles,averageDx,averageDy,averageDt);
			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;

//				xFloat[w] = (float)xPointer[w];
//				yFloat[w] = (float)yPointer[w];
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = (float)tracks[h].d[hh].x;
					linesArray[4*l+1] = (float)tracks[h].d[hh].y;
					linesArray[4*l+2] = (float)tracks[h].d[hh+1].x;
					linesArray[4*l+3] = (float)tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;

			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;
			iMAP->exposureTime = exposureTime;

			break;
		}
	case 5: // x y t
		{
			// count total number of detections
			this->localizationCount = 0;

			// to store number of detections in each file
			localizationCountFile = new int[files];

			for (int v = 0; v < numberOfFiles; v++) {
				// read file
				fileHandle = fopen(fileNameList[v],"r");
				int ch;

				// error if no file read
				if (fileHandle == NULL) {
					printf("Error in opening file\n");
					return;
				}

				localizationCountFile[v] = 0;
				// determine number of lines in file
				do {
					ch = fgetc(fileHandle);
					if (ch=='\n') {
						localizationCountFile[v]++;
						localizationCount++;
					}
				} while( ch != EOF);

				fclose(fileHandle);

				// eliminate last line of file (zero line)
				localizationCountFile[v]--;
				localizationCount--;
			}

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

			this->linesArray = new float[2*2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];
			
			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int j = 0;
			for (int v = 0; v < numberOfFiles; v++) {

				// allocate memory for individual tracks
				tracks[v].n = localizationCountFile[v];
				tracks[v].d = new Detection[localizationCountFile[v]];
				tracks[v].xy = new double[2*localizationCountFile[v]];
				tracks[v].averageDt = 0.0;
				tracks[v].averageDx = 0.0;
				tracks[v].averageDy = 0.0;

				fopen(fileNameList[v],"r");

				colormap(tracks[v].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);

				for(int k = 0; k < localizationCountFile[v]; k++) {
					fscanf(fileHandle,"%lf\t%lf\t%lf\t",
							&xPointer[j],
							&yPointer[j],
							&tPointer[j]);

					// tag with channel (file) number
					nPointer[j] = v;
					insideInterval[j] = true;

					// dummy value for intensity
					iPointer[j] = 1.0;

					// find max/min intensity
					if (iMax < iPointer[j]) { iMax = iPointer[j]; }
					if (iMin > iPointer[j]) { iMin = iPointer[j]; }
					// find max/min frame number
					if (tMax < tPointer[j]) { tMax = tPointer[j]; }
					if (tMin > tPointer[j]) { tMin = tPointer[j]; }
					// find max/min x,y coordinates;
					if (xMax<xPointer[j]) { xMax = xPointer[j]; }
					if (yMax<yPointer[j]) { yMax = yPointer[j]; }
					if (xMin>xPointer[j]) { xMin = xPointer[j]; }
					if (yMin>yPointer[j]) { yMin = yPointer[j]; }

					tracks[v].d[k].x = xPointer[j];
					tracks[v].d[k].y = yPointer[j];
					tracks[v].xy[2*k] = xPointer[j];
					tracks[v].xy[2*k+1] = yPointer[j];
					tracks[v].d[k].t = tPointer[j];
					tracks[v].d[k].i = iPointer[j];

					// total increment
					j++;
				}
				fclose(fileHandle);
			}

			exposureTime = 1000000.0;
			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					const double dt = fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDt += dt;
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);

					// assume smallest dt is exposure time
					if (exposureTime > dt) { exposureTime = dt; }
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].xy[2*f] -= xFileOffset;
					tracks[e].xy[2*f+1] -= yFileOffset;
				}
			}

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[4*l] = tracks[h].d[hh].x;
					linesArray[4*l+1] = tracks[h].d[hh].y;
					linesArray[4*l+2] = tracks[h].d[hh+1].x;
					linesArray[4*l+3] = tracks[h].d[hh+1].y;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					l++;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = 0.0; // to change for 3d files
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;

			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;

			break;
		}
	case 6: // tr x y z t
		{
			// 3D file format
			file3D = true;

			// count total number of detections
			this->localizationCount = 0;

			// read file
			fileHandle = fopen(fileNameList[0],"r");
			int ch;

			// error if no file read
			if (fileHandle == NULL) {
				printf("Error in opening file\n");
				return;
			}

			// determine number of lines in file
			do {
				ch = fgetc(fileHandle);
				if (ch=='\n') {
					localizationCount++;
				}
			} while( ch != EOF);

			fclose(fileHandle);

			// eliminate last line of file (zero line)
	//			localizationCount--;

			// initialize extrema
			this->xMax = -1000000.0;
			this->xMin = 1000000.0;
			this->yMax = -1000000.0;
			this->yMin = 1000000.0;
			this->zMax = -1000000.0;
			this->zMin = 1000000.0;
			this->iMax = -1000000.0;
			this->iMin = 1000000.0;
			this->tMax = -1000000.0;
			this->tMin = 1000000.0;
			this->nMin = 10000000;
			this->nMax = -10000000;

			// initialize storage arrays
			this->xPointer = new double[localizationCount];
			this->yPointer = new double[localizationCount];
			this->zPointer = new double[localizationCount];
			this->iPointer = new double[localizationCount];
			this->tPointer = new double[localizationCount];
			this->nPointer = new int[localizationCount];
			this->insideInterval = new bool[localizationCount];

	//			this->linesColorArray = new float[];
			fopen(fileNameList[0],"r");

			for(int j = 0; j < localizationCount; j++) {
				fscanf(fileHandle,"%i\t%lf\t%lf\t%lf\t%lf\t",
						&nPointer[j],
						&xPointer[j],
						&yPointer[j],
						&zPointer[j],
						&tPointer[j]);

	//				tFloat[j] = (float)tPointer[j];
				iPointer[j] = 1.0;
				insideInterval[j] = true;

				// find max/min frame number
				if (tMax < tPointer[j]) { tMax = tPointer[j]; }
				if (tMin > tPointer[j]) { tMin = tPointer[j]; }
				if (nMax < nPointer[j]) { nMax = nPointer[j]; }
				if (nMin > nPointer[j]) { nMin = nPointer[j]; }

				// find max/min x,y coordinates;
				if (xMax < xPointer[j]) { xMax = xPointer[j]; }
				if (yMax < yPointer[j]) { yMax = yPointer[j]; }
				if (xMin > xPointer[j]) { xMin = xPointer[j]; }
				if (yMin > yPointer[j]) { yMin = yPointer[j]; }
				if (zMax < zPointer[j]) { zMax = zPointer[j]; }
				if (zMin > zPointer[j]) { zMin = zPointer[j]; }
			}
			fclose(fileHandle);

			// count absolute number of files (trajectories numbers may not be numbered incrementally)
			numberOfFiles = 0;
			for (int bb = 0; bb < localizationCount-1; bb++) {
				if (nPointer[bb] != nPointer[bb+1]) { numberOfFiles++; }
			}
			numberOfFiles++;

			// for fast vertex rendering
			this->linesArray = new float[3*2*(localizationCount-numberOfFiles)];
			this->linesIndexArray = new int[2*(localizationCount-numberOfFiles)];
			this->linesColorArray = new float[4*2*(localizationCount-numberOfFiles)];

			// to store number of detections in each file
			localizationCountFile = new int[numberOfFiles];

			// initialize trajectory struct array
			tracks = new Trajectory[numberOfFiles];//(Trajectory*)calloc(numberOfFiles,sizeof(Trajectory));

			int g = 0;
			int r = 0;
			int id = nPointer[0];

			int maxLength = -10000;
			int minLength = 10000;

			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					g++;
					r++;
				}
	//				do {
	//					g++;
	//					r++;
	//				} while (id == nPointer[r] && r < localizationCount);

				if ( r >= localizationCount) { break; }
				id = nPointer[r];

				if (g == 1) { g++; } // to correct for issue where two-point trajectory is at very end of file

				tracks[h].n = g;
				tracks[h].d = new Detection[g];
//				tracks[h].xy = new double[2*g];
				tracks[h].xyz = new double[3*g];
				tracks[h].averageDt = 0.0;
				tracks[h].averageDx = 0.0;
				tracks[h].averageDy = 0.0;
				tracks[h].averageDz = 0.0;
				localizationCountFile[h] = g;
				colormap(tracks[h].rgb,colormapType-1,((float)ran1(&iMAP->seed)),0.0,1.0,false);

				if (g > maxLength) { maxLength = g; }
				if (g < minLength) { minLength = g; }
			}

			r = 0;
			id = nPointer[0];
			// assign points to individual trajectories
			for (int h = 0; h < numberOfFiles; h++) {
				g = 0;
				while (id == nPointer[r] && r < localizationCount-1) {
					tracks[h].d[g].x = xPointer[r];
					tracks[h].d[g].y = yPointer[r];
					tracks[h].d[g].z = zPointer[r];
					tracks[h].xyz[3*g] = xPointer[r];
					tracks[h].xyz[3*g+1] = yPointer[r];
					tracks[h].xyz[3*g+2] = zPointer[r];
					tracks[h].d[g].t = tPointer[r];
					tracks[h].d[g].i = iPointer[r];
					g++;
					r++;
				}
				if ( r >= localizationCount) { break; }
				id = nPointer[r];
			}

	//			fprintf(stderr,"[%i,%i]\n",minLength,maxLength);

			exposureTime = 1000000.0;
			for (int v = 0; v < numberOfFiles; v++) {
				for (int k = 0; k < tracks[v].n-1; k++) {
					const double dt = fabs(tracks[v].d[k+1].t-tracks[v].d[k].t);
					tracks[v].averageDt += dt;
					tracks[v].averageDx += fabs(tracks[v].d[k+1].x-tracks[v].d[k].x);
					tracks[v].averageDy += fabs(tracks[v].d[k+1].y-tracks[v].d[k].y);
					tracks[v].averageDz += fabs(tracks[v].d[k+1].z-tracks[v].d[k].z);

					// assume smallest dt is exposure time
					if (exposureTime > dt) { exposureTime = dt; }
				}
				tracks[v].averageDt /= tracks[v].n-1;
				tracks[v].averageDx /= tracks[v].n-1;
				tracks[v].averageDy /= tracks[v].n-1;
				tracks[v].averageDz /= tracks[v].n-1;

				averageDx += tracks[v].averageDx;
				averageDy += tracks[v].averageDy;
				averageDz += tracks[v].averageDz;
				averageDt += tracks[v].averageDt;
			}

			averageDx /= (double)numberOfFiles;
			averageDy /= (double)numberOfFiles;
			averageDz /= (double)numberOfFiles;
			averageDt /= (double)numberOfFiles;

	//			fprintf(stderr,"%i:\t%f\t%f\t%f\n",numberOfFiles,averageDx,averageDy,averageDt);
			// offset x,y,z coordinates to centre of screen
			xFileOffset = (xMin+xMax)/2.0;
			yFileOffset = (yMin+yMax)/2.0;
			zFileOffset = (zMin+zMax)/2.0;

			for (int w = 0; w < localizationCount; w++) {
				xPointer[w] -= xFileOffset;
				yPointer[w] -= yFileOffset;
				zPointer[w] -= zFileOffset;
			}

			for (int e = 0; e < numberOfFiles; e++) {
				for (int f = 0; f < tracks[e].n; f++) {
					tracks[e].d[f].x -= xFileOffset;
					tracks[e].d[f].y -= yFileOffset;
					tracks[e].d[f].z -= zFileOffset;
					tracks[e].xyz[3*f] -= xFileOffset;
					tracks[e].xyz[3*f+1] -= yFileOffset;
					tracks[e].xyz[3*f+2] -= zFileOffset;
				}
			}

			int l = 0;
			for (int h = 0; h < numberOfFiles; h++) {
				for (int hh = 0; hh < tracks[h].n-1; hh++) {
					linesArray[6*l] = (float)tracks[h].d[hh].x;
					linesArray[6*l+1] = (float)tracks[h].d[hh].y;
					linesArray[6*l+2] = (float)tracks[h].d[hh].z;
					linesArray[6*l+3] = (float)tracks[h].d[hh+1].x;
					linesArray[6*l+4] = (float)tracks[h].d[hh+1].y;
					linesArray[6*l+5] = (float)tracks[h].d[hh+1].z;
					linesColorArray[8*l] = linesColorArray[8*l+4] = tracks[h].rgb[0];
					linesColorArray[8*l+1] = linesColorArray[8*l+5] = tracks[h].rgb[1];
					linesColorArray[8*l+2] = linesColorArray[8*l+6] = tracks[h].rgb[2];
					linesColorArray[8*l+3] = linesColorArray[8*l+7] = iMAP->iOffsetSlider->value();
					linesIndexArray[l] = l;
					l++;
				}
			}

			xMax -= xFileOffset;
			xMin -= xFileOffset;
			yMax -= yFileOffset;
			yMin -= yFileOffset;
			zMax -= zFileOffset;
			zMin -= zFileOffset;

			// determine ranges, and max range
			xRange = xMax - xMin;
			yRange = yMax - yMin;
			zRange = zMax - zMin;
			maxRange = FMAX(xRange,yRange);

			// construct element array for rendering
			pointsArray = new float[3*localizationCount];
			indexArray = new int[localizationCount];
			colorArray = new float[4*localizationCount];

			// detection thresholding arrays
			tMinVisibility = new int[localizationCount];
			tMaxVisibility = new int[localizationCount];

			int q = 0;
			for (int p = 0; p < localizationCount; p++) {
				pointsArray[q] = (float) xPointer[p];
				q++;
				pointsArray[q] = (float) yPointer[p];
				q++;
				pointsArray[q] = (float) zPointer[p];
				q++;

				indexArray[p] = p;
				tMinVisibility[p] = true;
				tMaxVisibility[p] = true;
			}

			// fixed slider values
			tMinInitial = tMin;
			tMaxInitial = tMax;

			fRange = tMaxInitial-tMinInitial;
			iMAP->exposureTime = exposureTime;

			break;
		}

	}

	// display parameters
	plotType = 0; // frame plot by default
	colormapFlip = false;
	cMin = 0.0;
	cMax = 1.0;
	iOffset = 1.0;
	detectionSize = 10.0;
	orthoLimit = 1.0;
	xTranslate = 0.0;
	yTranslate = 0.0;
	zTranslate = 0.0;
	drawTrajectories = true;
	drawLocalizations = false;
	animationIndex = 0;
	animationProgress = tMin;
	animateTrajectories = false;
	boundingBoxEnable = true;
	dimensionsEnable = true;
	gridEnable = false;
	ticksEnable = true;
	gridSpacing = maxRange/25.0;

	// purely for diplay purposes
	localizationCountDisplay = localizationCount;
	localizationCountVoronoi = localizationCount;

	startInterval = tMin;
	endInterval = tMax;

	// meshing classes
	squareMesh = NULL;
	voronoiMesh = NULL;
	squareMeshApply = false;
	voronoiMeshApply = false;
	treeMeshApply = false;

	// activate OpenGL display
	this->fileLoaded = true;
}

void File::updateIntervalDisplay() {
	int l = 0;
	for (int h = 0; h < numberOfFiles; h++) {
		for (int hh = 0; hh < tracks[h].n-1; hh++) {
			if (tracks[h].d[hh].t <= startInterval || tracks[h].d[hh].t >= endInterval) {
				// adjust alpha
				linesColorArray[8*l+3] = linesColorArray[8*l+7] = 0.0;
			} else {
				linesColorArray[8*l+3] = linesColorArray[8*l+7] = iOffset;
			}
			l++;
		}
	}

	int u = 0;
	for (int l = 0; l < localizationCountDisplay; l++) {
		if (tPointer[l] <= startInterval || tPointer[l] >= endInterval) {
			file->colorArray[4*l+3] = 0.0;
			file->insideInterval[l] = false;
		} else {
			file->colorArray[4*l+3] = iOffset;
			file->insideInterval[l] = true;
			u++;
		}
	}

	localizationCountVoronoi = u;
	// hack: update cell numbers if voronoi tessellation gui open
	if (file->voronoiMeshGui != NULL && file->voronoiMeshApply == false) {
		file->voronoiMeshGui->cellsSlider->value(localizationCountVoronoi/80.0);
	}

	iMAP->fileInfo->updatePoints();
}

File::~File() {
	delete [] pointsArray;
	delete [] indexArray;
	delete [] colorArray;
	delete [] xPointer;
	delete [] yPointer;
	delete [] iPointer;
	delete [] tPointer;
	delete [] tMinVisibility;
	delete [] tMaxVisibility;
	delete [] tracks;
	delete [] nPointer;
}

void File::updateAlpha(float alpha) {
	int l = 0;
	for (int h = 0; h < numberOfFiles; h++) {
		for (int hh = 0; hh < tracks[h].n-1; hh++) {
			linesColorArray[8*l+3] = linesColorArray[8*l+7] = alpha;
			l++;
		}
	}
}

int File::getVisibility(int index) {
	return tMinVisibility[index] * tMaxVisibility[index];
}

void File::setMinFrame(float minimum) {
	for (int c = 0; c < localizationCount; c++) {
		if (tPointer[c] < minimum) {
			// mark detections to be blacked out
			tMinVisibility[c] = 0;
		}
		else {
			tMinVisibility[c] = 1;
		}
	}
	tMin = minimum;
	plotChanged = true;
}

void File::setMaxFrame(float maximum) {
	for (int c = 0; c < localizationCount; c++) {
		if (tPointer[c] > maximum) {
			// mark detections to be blacked out
			tMaxVisibility[c] = 0;
		}
		else {
			tMaxVisibility[c] = 1;
		}
	}
	tMax = maximum;
	plotChanged = true;
}

int File::updateIntervalData() {
	// count localizations inside interval
	int localizationsInsideInterval = 0;
	for (int l = 0; l < localizationCount; l++) {
		if (tPointer[l] >= startInterval && tPointer[l] <= endInterval) {
			localizationsInsideInterval++;
		}
	}

	if (localizationsInsideInterval <= 0) { return 0; }

	xPointerHolder = xPointer;
	yPointerHolder = yPointer;
	iPointerHolder = iPointer;
	tPointerHolder = tPointer;
	nPointerHolder = nPointer;

	// create new arrays
	xPointer = new double[localizationsInsideInterval];
	yPointer = new double[localizationsInsideInterval];
	iPointer = new double[localizationsInsideInterval];
	tPointer = new double[localizationsInsideInterval];
	nPointer = new int[localizationsInsideInterval];

	// store date inside interval
	int u = 0;
	for (int l = 0; l < localizationCount; l++) {
		if (tPointerHolder[l] >= startInterval && tPointerHolder[l] <= endInterval) {
			xPointer[u] = xPointerHolder[l];
			yPointer[u] = yPointerHolder[l];
			iPointer[u] = iPointerHolder[l];
			tPointer[u] = tPointerHolder[l];
			nPointer[u] = nPointerHolder[l];
			u++;
		}
	}
	localizationCount = localizationsInsideInterval;

	return 1;
}

void File::reloadOriginalData() {
	if (xPointerHolder != NULL && yPointerHolder != NULL && iPointerHolder != NULL && tPointerHolder != NULL && nPointerHolder != NULL) {
		delete [] xPointer;
		delete [] yPointer;
		delete [] iPointer;
		delete [] tPointer;
		delete [] nPointer;

		xPointer = xPointerHolder;
		yPointer = yPointerHolder;
		iPointer = iPointerHolder;
		tPointer = tPointerHolder;
		nPointer = nPointerHolder;

		xPointerHolder = NULL;
		yPointerHolder = NULL;
		iPointerHolder = NULL;
		tPointerHolder = NULL;
		nPointerHolder = NULL;

		localizationCount = this->localizationCountDisplay;
	}
}
