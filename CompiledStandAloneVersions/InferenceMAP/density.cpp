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

#include "draw.h"
#include "graphics.h"
#include "globals.h"
#include "tessellation.h"
#include "density.h"

#define PI 3.1415926535897932384626433832795028841971693993751

extern Globals *iMAP;
extern File *file;

void createPartitions() {

	// set initial density values to zero
	for (int c = 0; c < file->localizationCount; c++) {
		file->iPointer[c] = 0;
	}

	// delete any previously stored partitionMatrix to prevent memory leak
	if (file->partitionMatrix != NULL) {
		file->xParts = (int)ceil((file->xRange)/file->unitLength);
		file->yParts = (int)ceil((file->yRange)/file->unitLength);
		for (int i = 0; i < file->xParts; ++i) {
		  delete [] file->partitionMatrix[i];
		}
		delete [] file->partitionMatrix;
	}

	file->unitLength = iMAP->densityGui->getRadius();

	// calculate number of partitions per dimension
	file->xParts = (int)ceil((file->xRange)/file->unitLength);
	file->yParts = (int)ceil((file->yRange)/file->unitLength);

	// allocate memory for partitionMatrix
	file->partitionMatrix = new Partition*[file->xParts];
	for (int a = 0; a < file->xParts; a++) {
		file->partitionMatrix[a] = new Partition[file->yParts];
	}

	// count how many coordinates in each partition
	int x_pos,y_pos,z_pos;
	for (int u = 0; u < file->localizationCount; u++) {
		x_pos = (int)( (file->pointsArray[u*3]   + fabsf(file->xMin))/ceil(file->xRange)*(float)file->xParts );
		y_pos = (int)( (file->pointsArray[u*3+1] + fabsf(file->yMin))/ceil(file->yRange)*(float)file->yParts );

		if (x_pos >= file->xParts) { x_pos = file->xParts-1; }
		if (x_pos < 0) { x_pos = 0; }
		if (y_pos >= file->yParts) { y_pos = file->yParts-1; }
		if (y_pos < 0) { y_pos = 0; }

		file->partitionMatrix[x_pos][y_pos].count++;
	}

	// allocate memory to hold coordinates, in each partition
	for (int a = 0; a < file->xParts; a++) {
		for (int b = 0; b < file->yParts; b++) {
			file->partitionMatrix[a][b].coords = new Coordinate[file->partitionMatrix[a][b].count];
		}
	}

	// populate partitions with the coordinates
	for (int u = 0; u < file->localizationCount; u++) {
		x_pos = (int)( (file->pointsArray[u*3]   + fabsf(file->xMin))/ceil(file->xRange)*(float)file->xParts );
		y_pos = (int)( (file->pointsArray[u*3+1] + fabsf(file->yMin))/ceil(file->yRange)*(float)file->yParts );
		if (x_pos >= file->xParts) { x_pos = file->xParts-1; }
		if (x_pos < 0) { x_pos = 0; }
		if (y_pos >= file->yParts) { y_pos = file->yParts-1; }
		if (y_pos < 0) { y_pos = 0; }

		file->partitionMatrix[x_pos][y_pos].coords[file->partitionMatrix[x_pos][y_pos].progress].x = file->pointsArray[u*3];
		file->partitionMatrix[x_pos][y_pos].coords[file->partitionMatrix[x_pos][y_pos].progress].y = file->pointsArray[u*3+1];
		file->partitionMatrix[x_pos][y_pos].coords[file->partitionMatrix[x_pos][y_pos].progress].index = u;

		file->partitionMatrix[x_pos][y_pos].progress++;

		// reset progress when all coordinates filled
		if (file->partitionMatrix[x_pos][y_pos].count == file->partitionMatrix[x_pos][y_pos].progress) {
			file->partitionMatrix[x_pos][y_pos].progress = 0;
		}
	}

}

void calculateDensity() {
	// make sure computation is only done once
	if (file->densityCalculated == false) {

		createPartitions();

		// calculate densities for individual coordinates
		char progress[100];
		float distance = 0;
		for (int a = 0; a < file->xParts; a++) {
			for (int b = 0; b < file->yParts; b++) {
				// search points in surrounding partitions (27 total -- including own)
				for (int aa = -1; aa <= 1; aa++) {
					for (int bb = -1; bb <= 1; bb++) {
						for (int cc = -1; cc <= 1; cc++) {
							if ( a+aa >= 0 && b+bb >= 0 && a+aa < file->xParts && b+bb < file->yParts ) {
								for (int d = 0; d < file->partitionMatrix[a][b].count; d++) {
									for (int e = 0; e < file->partitionMatrix[a+aa][b+bb].count; e++) {
										distance = sqrt( (file->partitionMatrix[a][b].coords[d].x-file->partitionMatrix[a+aa][b+bb].coords[e].x)*(file->partitionMatrix[a][b].coords[d].x-file->partitionMatrix[a+aa][b+bb].coords[e].x) +
														 (file->partitionMatrix[a][b].coords[d].y-file->partitionMatrix[a+aa][b+bb].coords[e].y)*(file->partitionMatrix[a][b].coords[d].y-file->partitionMatrix[a+aa][b+bb].coords[e].y) );
										if ( distance < file->unitLength && distance != 0.0 ) {
											file->iPointer[file->partitionMatrix[a][b].coords[d].index]++;
										}
									}
								}
							}
						}
					}
				}
			}
			if (iMAP->updateDisplay) {
				iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
				char progress[70];
				sprintf(progress,"Density Calculation: %2.0f%%\n\n",(float)a/(float)file->xParts*100);
				textDisplayUpdate(progress);
			}
		}

		iMAP->textBuffer->remove(0,iMAP->textBuffer->length());
		textDisplayUpdate("Density Calculation 100%\n");
		textDisplayUpdate("Density Calculation Complete");

		// determine min/max density values
		file->iMin = 10000000;
		file->iMax = -10000000;

		for (int c = 0; c < file->localizationCount; c++) {
			if (file->iMax < file->iPointer[c]) {
				file->iMax = file->iPointer[c];
			}
			if (file->iMin > file->iPointer[c]) {
				file->iMin = file->iPointer[c];
			}
		}

		file->densityCalculated = true;

        Fl::redraw();
	}
}



