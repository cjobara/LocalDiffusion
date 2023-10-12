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

#define PI 3.1415926535897932384626433832795028841971693993751

extern Globals *iMAP;
extern File *file;


const char* uniformSpriteVS = {

"#version 120\n"

"uniform float angle2pixels;\n"
"uniform float detectionSize;\n"

"void main(void) {\n"

"    gl_Position = ftransform();\n"
"    gl_FrontColor = gl_Color;\n"
"    gl_FogFragCoord = length(gl_Position);\n"

"    float d = length(vec3(gl_Vertex)-vec3(0.0,0.0,15.0)); // hypotenuse to point\n"
"	float angle = asin(detectionSize/700.0/d); // angular size of point, in radians\n"

"	gl_PointSize = angle * angle2pixels; // ==fraction of screen * screen size\n"

"}"

};

const char* uniformSpriteFS = {

"#version 120\n"

"void main(void){\n"
"    float dist = distance( gl_PointCoord, vec2(0.5) );\n"
"    float alpha = 1.0 - smoothstep(0.0, 0.5, dist);\n"
"    gl_FragColor = gl_Color * vec4(1.0,1.0,1.0,alpha);\n"
"}"

};

void updateView() {
	iMAP->glWindow->redraw();

	if (file->landscapePlot == true || file->file3D == true) {
		glRotatef(file->xRotate,1,0,0);
		glRotatef(file->yRotate,0,1,0);
		glRotatef(file->zRotate,0,0,1);
		if (fabsf(file->xRotate) >= 360) { file->xRotate -= 360; }
		if (fabsf(file->yRotate) >= 360) { file->yRotate -= 360; }
		if (fabsf(file->zRotate) >= 360) { file->zRotate -= 360; }
		// for panning
		glTranslatef(iMAP->mv[0]*file->xTranslate+iMAP->mv[1]*file->yTranslate,
					 iMAP->mv[4]*file->xTranslate+iMAP->mv[5]*file->yTranslate,
					 iMAP->mv[8]*file->xTranslate+iMAP->mv[9]*file->yTranslate);

	} else {
		glTranslatef(file->xTranslate,file->yTranslate,file->zTranslate);
	}

}

void updateViewAxis() {
	glRotatef(file->xRotate,1.0,0,0);
	glRotatef(file->yRotate,0,1.0,0);
	glRotatef(file->zRotate,0,0,1.0);
}

void drawLocalizations() {
//	glShadeModel(GL_SMOOTH);
	glEnable(GL_POINT_SPRITE);
//	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glEnable(GL_POINT_SMOOTH);
//	glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);	// Really Nice Perspective Calculations
//	glDisable(GL_DEPTH_TEST);							// Disable Depth Testing

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

	if (file->landscapePlot) { glPointSize(35.0/iMAP->fovy); }
	else { glPointSize(35.0f/file->orthoLimit); }

	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	// Make sure updateView is within push & pop statements
	glPushMatrix();
		updateView();
		if (file->file3D) { zSortLocalizations(); }
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
		glColorPointer(4,GL_FLOAT,0,file->colorArray);
		glVertexPointer(3,GL_FLOAT,0,file->pointsArray);
		glDrawElements(GL_POINTS, file->localizationCountDisplay, GL_UNSIGNED_INT, file->indexArray);
	glPopMatrix();

	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glBindTexture(GL_TEXTURE_2D,0);

//	glEnable(GL_DEPTH_TEST);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_POINT_SPRITE);
	glDisable(GL_BLEND);
//	glDisable(GL_POINT_SMOOTH);
}

void drawFrame() {

	float relIntensity;
	float relFrame;

	if (file->plotChanged) {
		for (int r=0; r < file->localizationCountDisplay; r++) {
			relFrame = (float)(file->tPointer[r] - file->tMin) / (float)(file->tMax-file->tMin);
			relIntensity = (float)(file->iPointer[r]/file->iMax);
			colormap(file,r,relFrame);
			if (file->insideInterval[r]) {
				file->colorArray[r*4+3] = (relIntensity+file->iOffset)*file->getVisibility(r);
			}
		}
		file->plotChanged = false;
	}

	drawLocalizations();
}

void drawIntensity() {

	float relIntensity;

	if (file->plotChanged) {
		for (int r=0; r < file->localizationCountDisplay; r++) {
			relIntensity = (float)(file->iPointer[r]/file->iMax);
			colormap(file,r,relIntensity);
			if (file->insideInterval[r]) {
				file->colorArray[r*4+3] = (relIntensity+file->iOffset)*file->getVisibility(r);
			}
		}
		file->plotChanged = false;
	}

	drawLocalizations();
}

void drawTrajectories() {
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLineWidth(2.0);

	if (file->squareMeshOverlay && file->landscapePlot) {
		if (file->squareMeshGui->fogButton->value()) {
			glEnable(GL_FOG);
			const float fogCol[3]={0.0,0.0,0.0};
			glFogfv(GL_FOG_COLOR,fogCol);
			glFogi(GL_FOG_MODE, GL_LINEAR);
			glFogf(GL_FOG_DENSITY, 0.5f);
			glHint(GL_FOG_HINT,GL_NICEST);
			glFogf(GL_FOG_START, (float)file->squareMeshGui->fogStartSlider->value());
			glFogf(GL_FOG_END, (float)file->squareMeshGui->fogEndSlider->value());
		}
	}
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glPushMatrix();
		updateView();
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
		glVertexPointer(2,GL_FLOAT,0,file->linesArray);
		glColorPointer(4,GL_FLOAT,0,file->linesColorArray);
		glDrawArrays(GL_LINES,0,2*(file->localizationCountDisplay-file->numberOfFiles-1));
	glPopMatrix();

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glDisable(GL_FOG);

	glDisable(GL_BLEND);
}

void animateTrajectories() {
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(3.0);

	glPushMatrix();
		updateView();
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
		if (file->numberOfFiles == 1) {
			glBegin(GL_LINES);
			float *rgb = new float[3];
			for (int p = 0; p < file->animationIndex-1; p++) {
				rgb = colormap(9,(double)(p+1)/(double)(file->animationIndex),file->cMin,file->cMax,false);
				if (iMAP->animateAccumulationButton->value()) {
					glColor4f(rgb[0],rgb[1],rgb[2],(double)(p+1)/(double)(file->animationIndex));
					glVertex2f(file->xPointer[p],file->yPointer[p]);
					glVertex2f(file->xPointer[p+1],file->yPointer[p+1]);
				} else {
					if (file->animationIndex-p < iMAP->animateDelayedStepsSlider->value()) {
						glColor4f(rgb[0],rgb[1],rgb[2],1.0);
						glVertex2f(file->xPointer[p],file->yPointer[p]);
						glVertex2f(file->xPointer[p+1],file->yPointer[p+1]);
					}
				}
			}
			glEnd();
		} else {
			int m = 0;
			for (int v = 0; v < file->numberOfFiles-1; v++) {
				glBegin(GL_LINES);
					for (int k = 0; k < file->tracks[v].n-1; k++) {
						if (file->tracks[v].d[k].t <= (double)file->animationIndex*file->exposureTime) {
							if (file->drawTrajectories) {
								glColor4f(1.0f,1.0f,1.0f,1.0f);
								glVertex2f(file->tracks[v].d[k].x,file->tracks[v].d[k].y);
								glVertex2f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y);
							}
							else {
								if (iMAP->animateAccumulationButton->value()) {
									glColor4f(file->tracks[v].rgb[0],file->tracks[v].rgb[1],file->tracks[v].rgb[2],pow(2000.0*file->tracks[v].d[0].t/(float)file->animationIndex*iMAP->exposureTime,2.0));
//									glColor4f(file->tracks[v].rgb[0],file->tracks[v].rgb[1],file->tracks[v].rgb[2],1500.0*file->tracks[v].d[0].t/(float)file->animationIndex*iMAP->exposureTime);
									glVertex2f(file->tracks[v].d[k].x,file->tracks[v].d[k].y);
									glVertex2f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y);
								} else {
									if (file->animationIndex-(int)(file->tracks[v].d[k].t/file->exposureTime) < iMAP->animateDelayedStepsSlider->value()) {
										glColor4f(file->tracks[v].rgb[0],file->tracks[v].rgb[1],file->tracks[v].rgb[2],1.0);
										glVertex2f(file->tracks[v].d[k].x,file->tracks[v].d[k].y);
										glVertex2f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y);
									}
								}
							}
						}
						m++;
					}
				glEnd();
			}
		}
	glPopMatrix();

	glDisable(GL_BLEND);
}

void drawSquareMesh() {
	// mesh width and height
	const float dx = file->squareMeshGui->getDx();

	// number of lines
	int xLines,yLines;
	if (file->squareMeshApplied()) {
		if (file->squareMesh->selectionMode()) {
			xLines = (int)(file->squareMesh->selection.xRange/dx)+1;
			yLines = (int)(file->squareMesh->selection.yRange/dx)+1;
		} else {
			xLines = (int)(file->xRange/dx)+1;
			yLines = (int)(file->yRange/dx)+1;
		}
	} else {
		xLines = (int)(file->xRange/dx)+1;
		yLines = (int)(file->yRange/dx)+1;
	}

	// minimum line spacing
	float dMin = dx;

	// label parameters
	const float gridSpacing = (dMin)/25.0f;
	const float scale = gridSpacing/4.0f;
	txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);

	if (!file->squareMeshApplied()) {
		
		int **pointArray = new int*[xLines];
		for (int a = 0; a < xLines; a++) {
			pointArray[a] = new int[yLines];
			for (int b = 0; b < yLines; b++) {
				pointArray[a][b] = 0;
			}
		}

		// count how many coordinates in each partition
		for (int u = 0; u < file->localizationCount; u++) {
			const int x_pos = (int)( (file->pointsArray[u*3]   + fabsf(file->xMin))/dx );
			const int y_pos = (int)( (file->pointsArray[u*3+1] + fabsf(file->yMin))/dx );
			pointArray[x_pos][y_pos]++;
		}

		if (iMAP->overlayAdjustment) {

			if (iMAP->linesOverlay != NULL) {
				delete [] iMAP->linesOverlay;
				iMAP->linesOverlay = NULL;
			}

			// grid lines
			iMAP->linesOverlay = new float[4*xLines+4*yLines];

			for (int a = 0; a < xLines; a++) {
				iMAP->linesOverlay[4*a] = file->xMin+(float)a*dx;
				iMAP->linesOverlay[4*a+1] = file->yMin;
				iMAP->linesOverlay[4*a+2] = file->xMin+(float)a*dx;
				iMAP->linesOverlay[4*a+3] = file->yMax;
			}
			for (int b = 0; b < yLines; b++) {
				iMAP->linesOverlay[4*xLines+4*b] = file->xMin;
				iMAP->linesOverlay[4*xLines+4*b+1] = file->yMin+(float)b*dx;
				iMAP->linesOverlay[4*xLines+4*b+2] = file->xMax;
				iMAP->linesOverlay[4*xLines+4*b+3] = file->yMin+(float)b*dx;
			}

			iMAP->overlayAdjustment = false;
		}

		// draw lines + number of points annotation
		glPushMatrix();
			updateView();
			glLineWidth(1.0);

			glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[0],1.0f);

			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(2,GL_FLOAT,0,iMAP->linesOverlay);
			glDrawArrays(GL_LINES,0,2*xLines+2*yLines);
			glDisableClientState(GL_VERTEX_ARRAY);

			if (xLines < 50) {
				// draw string for each box of mesh
				char label[20];
				glEnable(GL_TEXTURE_2D);
				glColor4f(1.0f,1.0f,1.0f,1.0f);
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER,0.5f);

				// offset-less boxes
				if (file->squareMeshGui->localizationNumberLabelButton->value()) {
					for (int a = 0; a < xLines; a++) {
						for (int b = 0; b < yLines; b++) {
							if (pointArray[a][b] > 0) {
								glPushMatrix();
									sprintf(label, "%i", pointArray[a][b]);
//									sprintf(label,"%i,%i",a,b);
									glTranslatef(file->xMin+a*dx,file->yMin+b*dx,0.0f);
									glScalef(scale, scale, scale);
									txfRenderString(iMAP->fontTex, label, strlen(label));
								glPopMatrix();
							}
						}
					}
				}

				glDisable(GL_ALPHA_TEST);
				glBindTexture(GL_TEXTURE_2D,0);
				glDisable(GL_TEXTURE_2D);
			}
		glPopMatrix();

		for (int rr = 0; rr < xLines; rr++) {
			delete [] pointArray[rr];
		}
		delete [] pointArray;
	}
	// mesh has been applied
	else {
		if (iMAP->overlayAdjustment) {

			if (iMAP->linesOverlay != NULL) {
				delete [] iMAP->linesOverlay;
				iMAP->linesOverlay = NULL;
			}

			iMAP->activeZones = 0;
			for (int a = 0; a < xLines; a++) {
				for (int b = 0; b < yLines; b++) {
					if (file->squareMesh->active(a,b)) {
						iMAP->activeZones++;
					}
				}
			}
			// grid lines
			iMAP->linesOverlay = new float[16*iMAP->activeZones];

			float xMin,xMax,yMin,yMax;

			if (file->squareMesh->selectionMode()) {
				xMin = (float)file->squareMesh->selection.xMin;//-file->xTranslate;
				xMax = (float)file->squareMesh->selection.xMax;//-file->xTranslate;
				yMin = (float)file->squareMesh->selection.yMin;//-file->yTranslate;
				yMax = (float)file->squareMesh->selection.yMax;//-file->yTranslate;
			} else {
				xMin = (float)file->xMin;//-file->xTranslate;
				xMax = (float)file->xMax;//-file->xTranslate;
				yMin = (float)file->yMin;//-file->yTranslate;
				yMax = (float)file->yMax;//-file->yTranslate;
			}

			int g = 0;
			for (int a = 0; a < xLines; a++) {
				for (int b = 0; b < yLines; b++) {
					if (file->squareMesh->active(a,b)) {
						// bottom segment
						iMAP->linesOverlay[16*g] = (float)a*dx+xMin;
						iMAP->linesOverlay[16*g+1] = (float)b*dx+yMin;
						iMAP->linesOverlay[16*g+2] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->linesOverlay[16*g+3] = (float)b*dx+yMin;
						// right segment
						iMAP->linesOverlay[16*g+4] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->linesOverlay[16*g+5] = (float)b*dx+yMin;
						iMAP->linesOverlay[16*g+6] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->linesOverlay[16*g+7] = (float)FMIN(((float)b+1)*dx+yMin,yMax);
						// top segment
						iMAP->linesOverlay[16*g+8] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->linesOverlay[16*g+9] = (float)FMIN(((float)b+1)*dx+yMin,yMax);
						iMAP->linesOverlay[16*g+10] = (float)a*dx+xMin;
						iMAP->linesOverlay[16*g+11] = (float)FMIN(((float)b+1)*dx+yMin,yMax);
						// left segment
						iMAP->linesOverlay[16*g+12] = (float)a*dx+xMin;
						iMAP->linesOverlay[16*g+13] = (float)FMIN(((float)b+1)*dx+yMin,yMax);
						iMAP->linesOverlay[16*g+14] = (float)a*dx+xMin;
						iMAP->linesOverlay[16*g+15] = (float)b*dx+yMin;
						g++;
					}
				}
			}

			const float alpha = file->squareMeshGui->overlayAlphaSlider->value();

			if (iMAP->quadsOverlay != NULL) {
				delete [] iMAP->quadsOverlay;
				iMAP->quadsOverlay = NULL;
				delete [] iMAP->quadsColor;
				iMAP->quadsColor = NULL;
			}

			iMAP->quadsOverlay = new float[8*iMAP->activeZones];
			iMAP->quadsColor = new float[4*4*iMAP->activeZones];

			int k = 0;
			for (int a = 0; a < xLines; a++) {
				for (int b = 0; b < yLines; b++) {
					if (file->squareMesh->active(a,b)) {
						iMAP->quadsOverlay[8*k] = (float)a*dx+xMin;
						iMAP->quadsOverlay[8*k+1] = (float)b*dx+yMin;
						iMAP->quadsOverlay[8*k+2] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->quadsOverlay[8*k+3] = (float)b*dx+yMin;
						iMAP->quadsOverlay[8*k+4] = (float)FMIN(((float)a+1)*dx+xMin,xMax);
						iMAP->quadsOverlay[8*k+5] = (float)FMIN(((float)b+1)*dx+yMin,yMax);
						iMAP->quadsOverlay[8*k+6] = (float)a*dx+xMin;
						iMAP->quadsOverlay[8*k+7] = (float)FMIN(((float)b+1)*dx+yMin,yMax);

						// diffusion overlay
						if (file->squareMeshGui->overlayDiffusionButton->value() && !file->squareMeshGui->spotVisualizationButton->value()) {
							const float dValue = (float) ( file->squareMesh->getCell(a,b)->getDiffusion() - file->squareMesh->getDiffusionMin() ) /
												 (float) ( file->squareMesh->getDiffusionMax() - file->squareMesh->getDiffusionMin() );
							if (dValue >= 0.0f && dValue <= 1.0f) {
								colormap(&iMAP->quadsColor[16*k],file->squareMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
								iMAP->quadsColor[16*k+3] = alpha;
							} else {
								iMAP->quadsColor[16*k] = iMAP->quadsColor[16*k+1] = iMAP->quadsColor[16*k+2] = iMAP->quadsColor[16*k+3] = 0.0;
							}
						}
						// potential overlay
						else if (file->squareMeshGui->overlayPotentialButton->value() && !file->squareMeshGui->spotVisualizationButton->value()) {
							const float pValue = (float) ( file->squareMesh->getCell(a,b)->getPotential() - file->squareMesh->getPotentialMin() ) /
												 (float) ( file->squareMesh->getPotentialMax() - file->squareMesh->getPotentialMin() );
							if (pValue >= 0.0f && pValue <= 1.0f) {
								colormap(&iMAP->quadsColor[16*k],file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
								iMAP->quadsColor[16*k+3] = alpha;
							} else {
								iMAP->quadsColor[16*k] = iMAP->quadsColor[16*k+1] = iMAP->quadsColor[16*k+2] = iMAP->quadsColor[16*k+3] = 0.0;
							}
						}
						// force magnitude overlay
						else if (file->squareMeshGui->overlayForceMagnitudeButton->value() && !file->squareMeshGui->spotVisualizationButton->value()) {
							const float pValue = (float) ( file->squareMesh->getCell(a,b)->getForceMagnitude() - file->squareMesh->getForceMin() ) /
												 (float) ( file->squareMesh->getForceMax() - file->squareMesh->getForceMin() );
							if (pValue >= 0.0f && pValue <= 1.0f) {
								colormap(&iMAP->quadsColor[16*k],file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
								iMAP->quadsColor[16*k+3] = alpha;
							} else {
								iMAP->quadsColor[16*k] = iMAP->quadsColor[16*k+1] = iMAP->quadsColor[16*k+2] = iMAP->quadsColor[16*k+3] = 0.0;
							}
						}
						// number of points overlay
						else if (file->squareMeshGui->overlayPointNumberButton->value() && !file->squareMeshGui->spotVisualizationButton->value()) {
							const float pValue = ( (float)file->squareMesh->getCell(a,b)->getCount() - (float)file->squareMesh->getMinCell()->getCount() ) /
												 ( (float)file->squareMesh->getMaxCell()->getCount() - (float)file->squareMesh->getMinCell()->getCount() );
							if (pValue >= 0.0f && pValue <= 1.0f) {
								colormap(&iMAP->quadsColor[16*k],file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
								iMAP->quadsColor[16*k+3] = alpha;
							} else {
								iMAP->quadsColor[16*k] = iMAP->quadsColor[16*k+1] = iMAP->quadsColor[16*k+2] = iMAP->quadsColor[16*k+3] = 0.0;
							}
						}
						// green (default)
						else if (!file->squareMeshGui->spotVisualizationButton->value()){
							iMAP->quadsColor[16*k] = 0.2;
							iMAP->quadsColor[16*k+1] = 1.0;
							iMAP->quadsColor[16*k+2] = 0.2;
							iMAP->quadsColor[16*k+3] = alpha;
						}
						iMAP->quadsColor[16*k+4] = iMAP->quadsColor[16*k+8] = iMAP->quadsColor[16*k+12] = iMAP->quadsColor[16*k];
						iMAP->quadsColor[16*k+5] = iMAP->quadsColor[16*k+9] = iMAP->quadsColor[16*k+13] = iMAP->quadsColor[16*k+1];
						iMAP->quadsColor[16*k+6] = iMAP->quadsColor[16*k+10] = iMAP->quadsColor[16*k+14] = iMAP->quadsColor[16*k+2];
						iMAP->quadsColor[16*k+7] = iMAP->quadsColor[16*k+11] = iMAP->quadsColor[16*k+15] = iMAP->quadsColor[16*k+3];
						k++;
					}
				}
			}
			iMAP->overlayAdjustment = false;
		}

		if (file->squareMeshGui->displayGridButton->value()) {
			// draw function
			glPushMatrix();
				updateView();
				glLineWidth(1.0);

				glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[0],file->squareMeshGui->gridAlphaSlider->value());

				glEnableClientState(GL_VERTEX_ARRAY);
				glVertexPointer(2,GL_FLOAT,0,iMAP->linesOverlay);
				glDrawArrays(GL_LINES,0,8*iMAP->activeZones);
				glDisableClientState(GL_VERTEX_ARRAY);

				if (xLines < 50) {
					// draw string for each box of mesh
					char label[20];
					glEnable(GL_TEXTURE_2D);
					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glEnable(GL_ALPHA_TEST);
					glAlphaFunc(GL_GREATER,0.5f);

					// offset-less boxes
					if (file->squareMeshGui->localizationNumberLabelButton->value()) {
						for (int a = 0; a < xLines; a++) {
							for (int b = 0; b < yLines; b++) {
								if (file->squareMesh->active(a,b)) {
									glPushMatrix();
										sprintf(label,"%i",file->squareMesh->getCount(a,b));
		//								sprintf(label,"%i,%i",a,b);
										glTranslatef(file->xMin+a*dx,file->yMin+b*dx,0.0f);
										glScalef(scale, scale, scale);
										txfRenderString(iMAP->fontTex, label, strlen(label));
									glPopMatrix();
								}
							}
						}
					}

					glDisable(GL_ALPHA_TEST);
					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
				}
			glPopMatrix();

		}

		// color boxes when mesh is applied
		if (file->squareMeshGui->overlay()) {
			glPushMatrix();
				updateView();
				glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
				glEnableClientState(GL_VERTEX_ARRAY);
				glEnableClientState(GL_COLOR_ARRAY);
				glVertexPointer(2,GL_FLOAT,0,iMAP->quadsOverlay);
				glColorPointer(4,GL_FLOAT,0,iMAP->quadsColor);
				glDrawArrays(GL_QUADS,0,4*iMAP->activeZones);
				glDisableClientState(GL_COLOR_ARRAY);
				glDisableClientState(GL_VERTEX_ARRAY);
			glPopMatrix();
		}
		glDisable(GL_DEPTH_TEST);

		if (file->squareMeshGui->neighbourDistanceViewButton->value()) {
			glColor4f(1.0,1.0,0.0,0.7);
			glLineWidth(2.0);
			glPushMatrix();
			updateView();
			float am1dist,ap1dist,bm1dist,bp1dist;
			const float maxDist = file->squareMesh->getMaximumNeighbourDistance();
			for (int a = 0; a < xLines; a++) {
				for (int b = 0; b < yLines; b++) {
					if (file->squareMesh->active(a,b)) {
						if (a > 0) {
							am1dist = sqrt( (file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a-1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
											 (file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a-1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
						} else { am1dist = 1.0e5; }
						if (a < xLines-1) {
							ap1dist = sqrt( (file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a+1,b)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
											 (file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a+1,b)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
						} else { ap1dist = 1.0e5; }
						if (b > 0) {
							bm1dist = sqrt( (file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b-1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
											 (file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b-1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
						} else { bm1dist = 1.0e5; }
						if (b < yLines-1) {
							bp1dist = sqrt( (file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid())*(file->squareMesh->getCell(a,b+1)->getXCentroid()-file->squareMesh->getCell(a,b)->getXCentroid()) +
											 (file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid())*(file->squareMesh->getCell(a,b+1)->getYCentroid()-file->squareMesh->getCell(a,b)->getYCentroid()) );
						} else { bp1dist = 1.0e5; }

						if (a == xLines-1) {
							if (file->squareMesh->active(a-1,b) && am1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a-1,b)->getXCentroid(),file->squareMesh->getCell(a-1,b)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}
						else if (a == 0) {
							if (file->squareMesh->active(a+1,b) && ap1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a+1,b)->getXCentroid(),file->squareMesh->getCell(a+1,b)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}
						else {
							if (file->squareMesh->active(a-1,b) && am1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a-1,b)->getXCentroid(),file->squareMesh->getCell(a-1,b)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
							if (file->squareMesh->active(a+1,b) && ap1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a+1,b)->getXCentroid(),file->squareMesh->getCell(a+1,b)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}

						if (b == yLines-1) {
							if (file->squareMesh->active(a,b-1) && bm1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a,b-1)->getXCentroid(),file->squareMesh->getCell(a,b-1)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}
						else if (b == 0) {
							if (file->squareMesh->active(a,b+1) && bp1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a,b+1)->getXCentroid(),file->squareMesh->getCell(a,b+1)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}
						else {
							if (file->squareMesh->active(a,b-1) && bm1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a,b-1)->getXCentroid(),file->squareMesh->getCell(a,b-1)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
							if (file->squareMesh->active(a,b+1) && bp1dist < maxDist) {
								glBegin(GL_LINES);
									glVertex2f(file->squareMesh->getCell(a,b+1)->getXCentroid(),file->squareMesh->getCell(a,b+1)->getYCentroid());
									glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
								glEnd();
							}
						}


					}
				}
			}

			glPopMatrix();
		}

	}

	glDisable(GL_BLEND);

	if (file->squareMeshGui->spotVisualizationButton->value()) {
		glShadeModel(GL_SMOOTH);
		glEnable(GL_POINT_SPRITE);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_POINT_SMOOTH);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
		glDisable(GL_DEPTH_TEST);							// Disable Depth Testing

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

		for (int a = 0; a < xLines; a++) {
			for (int b = 0; b < yLines; b++) {
				const float sizeFactor = file->squareMeshGui->spotVisualizationScaleSlider->value();
				glPointSize(sizeFactor*50.0f);
				glPushMatrix();
					updateView();
					if (file->squareMesh->active(a,b)) {
						float *rgb = new float[3];
						// diffusion overlay
						if (file->squareMeshGui->overlayDiffusionButton->value()) {
							const float dValue = (float) ( file->squareMesh->getCell(a,b)->getDiffusion() - file->squareMesh->getDiffusionMin() ) /
												 (float) ( file->squareMesh->getDiffusionMax() - file->squareMesh->getDiffusionMin() );
							colormap(rgb,file->squareMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
							glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
							glBegin(GL_POINTS);
								glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
							glEnd();
						}
						// potential overlay
						else if (file->squareMeshGui->overlayPotentialButton->value()) {
							const float pValue = (float) ( file->squareMesh->getCell(a,b)->getPotential() - file->squareMesh->getPotentialMin() ) /
												 (float) ( file->squareMesh->getPotentialMax() - file->squareMesh->getPotentialMin() );
							colormap(rgb,file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
							glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
							glBegin(GL_POINTS);
								glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
							glEnd();
						}
						// force magnitude overlay
						else if (file->squareMeshGui->overlayForceMagnitudeButton->value()) {
							const float pValue = (float) ( file->squareMesh->getCell(a,b)->getForceMagnitude() - file->squareMesh->getForceMin() ) /
												 (float) ( file->squareMesh->getForceMax() - file->squareMesh->getForceMin() );
							colormap(rgb,file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
							glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
							glBegin(GL_POINTS);
								glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
							glEnd();
						}
						// point number overlay
						else if (file->squareMeshGui->overlayPointNumberButton->value()) {
							const float pValue = (float) ( file->squareMesh->getCell(a,b)->getCount() - file->squareMesh->getMinCell()->getCount() ) /
												 (float) ( file->squareMesh->getMaxCell()->getCount() - file->squareMesh->getMinCell()->getCount() );
							colormap(rgb,file->squareMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
							glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
							glBegin(GL_POINTS);
								glVertex2f(file->squareMesh->getCell(a,b)->getXCentroid(),file->squareMesh->getCell(a,b)->getYCentroid());
							glEnd();
						}
					}
				glPopMatrix();
			}
		}
		glBindTexture(GL_TEXTURE_2D,0);

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_POINT_SPRITE);
		glDisable(GL_BLEND);
		glDisable(GL_POINT_SMOOTH);
	}

	if (file->squareMeshGui->overlayForceArrowsButton->value()) {
		clipPlanes();
		drawSquareForces();
		unclipPlanes();
	}
	
}

void drawSquareRandomizedOptimizationZones() {
	const float dx = file->squareMeshGui->getDx();
//	const int sgdZones = iMAP->regularMeshGui->roZoneNumberSlider->value();
	const int roZones = file->squareMesh->roZones;

	float xMin,xMax,yMin,yMax;
	if (file->squareMesh->selectionMode()) {
		xMin = (float)file->squareMesh->selection.xMin-file->xTranslate;
		xMax = (float)file->squareMesh->selection.xMax-file->xTranslate;
		yMin = (float)file->squareMesh->selection.yMin-file->yTranslate;
		yMax = (float)file->squareMesh->selection.yMax-file->yTranslate;
	} else {
		xMin = (float)file->xMin-file->xTranslate;
		xMax = (float)file->xMax-file->xTranslate;
		yMin = (float)file->yMin-file->yTranslate;
		yMax = (float)file->yMax-file->yTranslate;
	}

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	for (int k = 0; k < roZones; k++) {

		const int a = file->squareMesh->roXCoords[k];
		const int b = file->squareMesh->roYCoords[k];

		glPushMatrix();
			updateView();

			glColor4f(1.0f,1.0f,1.0f,0.8f);

			glBegin(GL_QUADS);

			glVertex2f((float)a*dx+xMin,(float)b*dx+yMin);
			glVertex2f((float)FMIN((a+1)*dx+xMin,xMax),(float)b*dx+yMin);
			glVertex2f((float)FMIN((a+1)*dx+xMin,xMax),(float)FMIN((b+1)*dx+yMin,yMax));
			glVertex2f((float)a*dx+xMin,(float)FMIN((b+1)*dx+yMin,yMax));

			glEnd();

		glPopMatrix();
	}
	glDisable(GL_BLEND);

}

void drawSquareForces() {
	const float dx = file->squareMeshGui->getDx();

	float xMin,xMax,yMin,yMax;
	if (file->squareMesh->selectionMode()) {
		xMin = (float)file->squareMesh->selection.xMin;
		xMax = (float)file->squareMesh->selection.xMax;
		yMin = (float)file->squareMesh->selection.yMin;
		yMax = (float)file->squareMesh->selection.yMax;
	} else {
		xMin = (float)file->xMin;
		xMax = (float)file->xMax;
		yMin = (float)file->yMin;
		yMax = (float)file->yMax;
	}

	// minimum line spacing
	float dMin = dx;
	dMin *= 50.0f;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
	// draw arrows signifying force magnitude and direction
	glPushMatrix();
		updateView();

		glLineWidth(1.0);

		float arrowColor[] = {1.0,1.0,0.0,0.7};
		if (file->optimizationMode == 2) { arrowColor[0] = 0.0; arrowColor[1] = 1.0; arrowColor[2] = 0.0; }

		for (int a = 0; a < file->squareMesh->getXCells(); a++) {
			for (int b = 0; b < file->squareMesh->getYCells(); b++) {
				if (file->squareMesh->getCell(a,b)->getForceMagnitude() != 0.0f) {
//					const double l = dMin*file->regularMesh->getCell(a,b)->getForceMagnitude()/file->regularMesh->getForceMax();
					const double l = dMin*file->squareMesh->getCell(a,b)->getForceMagnitude()/file->squareMesh->getForceMax();

					const double fx = file->squareMesh->getCell(a,b)->getForceX();
					const double fy = file->squareMesh->getCell(a,b)->getForceY();

					const double xc = xMin+((double)a+0.5f)*dx;
					const double yc = yMin+((double)b+0.5f)*dx;

					glPushMatrix();
						glTranslatef(xc,yc,0.0f);
						glRotatef(180.0f/PI*atan2(fy,fx),0.0f,0.0f,1.0f);
						glTranslatef(-l/140.0f,0.0f,0.0f);
						glRotatef(90.0f,0.0f,0.0f,1.0f);
						glBegin(GL_TRIANGLES);
							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
							glVertex2f(l*0.005f,0.0f);

							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
							glVertex2f(-l*0.005f,0.0f);

							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],1.0);
							glVertex2f(0.0f,-l*0.015f);
						glEnd();
						glBegin(GL_LINE_LOOP);
							glColor4f(1.0f,1.0f,1.0f,0.0f);
							glVertex2f(l*0.005f,0.0f);

							glColor4f(1.0f,1.0f,1.0f,0.0f);
							glVertex2f(-l*0.005f,0.0f);

							glColor4f(1.0f,1.0f,1.0f,0.0f);
							glVertex2f(0.0f,-l*0.015f);
						glEnd();
					glPopMatrix();
				}
			}
		}

	glPopMatrix();
	glDisable(GL_BLEND);

}

void drawSquarePosteriorHighlight() {
	if (file->squareMeshGui->posterioriGroup->visible() && file->squareMeshGui->posterioriGroup->active()) {

		const float dx = file->squareMeshGui->getDx();

		const int a = file->squareMesh->getCurrentZoneX();
		const int b = file->squareMesh->getCurrentZoneY();
		const int a2 = file->squareMesh->getCurrentZone2X();
		const int b2 = file->squareMesh->getCurrentZone2Y();

		float xMin,xMax,yMin,yMax;

		if (file->squareMesh->selectionMode()) {
			xMin = (float)file->squareMesh->selection.xMin-file->xTranslate;
			xMax = (float)file->squareMesh->selection.xMax-file->xTranslate;
			yMin = (float)file->squareMesh->selection.yMin-file->yTranslate;
			yMax = (float)file->squareMesh->selection.yMax-file->yTranslate;
		} else {
			xMin = (float)file->xMin-file->xTranslate;
			xMax = (float)file->xMax-file->xTranslate;
			yMin = (float)file->yMin-file->yTranslate;
			yMax = (float)file->yMax-file->yTranslate;
		}

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		clipPlanes();

		// draw arrows signifying force magnitude and direction
		glPushMatrix();
			updateView();

			glLineWidth(3.0);

			glColor4f(1.0f,1.0f,1.0f,1.0f);

			glBegin(GL_LINE_LOOP);
				glVertex2f(a*dx+xMin,b*dx+yMin);
				glVertex2f((a+1)*dx+xMin,b*dx+yMin);
				glVertex2f((a+1)*dx+xMin,(b+1)*dx+yMin);
				glVertex2f(a*dx+xMin,(b+1)*dx+yMin);
			glEnd();

			glColor4f(1.0f,1.0f,1.0f,0.6f);
			glBegin(GL_QUADS);
				glVertex2f(a*dx+xMin,b*dx+yMin);
				glVertex2f((a+1)*dx+xMin,b*dx+yMin);
				glVertex2f((a+1)*dx+xMin,(b+1)*dx+yMin);
				glVertex2f(a*dx+xMin,(b+1)*dx+yMin);
			glEnd();

			if (file->squareMeshGui->potentialReferenceButton->value()) {
				glBegin(GL_LINES);
					glLineWidth(4.0);
					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glVertex2f(file->squareMesh->getCentroidX(a,b),file->squareMesh->getCentroidY(a,b));
					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glVertex2f(file->squareMesh->getCentroidX(a2,b2),file->squareMesh->getCentroidY(a2,b2));
				glEnd();

				glBegin(GL_LINE_LOOP);
					glVertex2f(a2*dx+xMin,b2*dx+yMin);
					glVertex2f((a2+1)*dx+xMin,b2*dx+yMin);
					glVertex2f((a2+1)*dx+xMin,(b2+1)*dx+yMin);
					glVertex2f(a2*dx+xMin,(b2+1)*dx+yMin);
				glEnd();

				glColor4f(1.0f,1.0f,1.0f,0.6f);
				glBegin(GL_QUADS);
					glVertex2f(a2*dx+xMin,b2*dx+yMin);
					glVertex2f((a2+1)*dx+xMin,b2*dx+yMin);
					glVertex2f((a2+1)*dx+xMin,(b2+1)*dx+yMin);
					glVertex2f(a2*dx+xMin,(b2+1)*dx+yMin);
				glEnd();

			}
		glPopMatrix();

		unclipPlanes();

		glDisable(GL_BLEND);
	}
}

void drawVoronoiMesh() {

	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glBlendEquation(GL_MAX);

	// draw voronoi diagram
	glEnableClientState(GL_VERTEX_ARRAY);

	glLineWidth(2.0f);
//	iMAP->activeZones = 0;
	for (int h = 0; h < file->voronoiMesh->getNumberOfClusters(); h++) {
		glPushMatrix();
			updateView();
			glVertexPointer(2,GL_FLOAT,0,file->voronoiMesh->getCell(h)->vertices);
			if (file->voronoiMesh->active(h) && file->voronoiMeshGui->overlay()) {
//				iMAP->activeZones++;
				float *rgb = new float[3];
				// diffusion overlay
				if (file->voronoiMeshGui->overlayDiffusionButton->value()) {
					const float dValue = (float) ( file->voronoiMesh->getCell(h)->getDiffusion() - file->voronoiMesh->getDiffusionMin() ) /
										 (float) ( file->voronoiMesh->getDiffusionMax() - file->voronoiMesh->getDiffusionMin() );
					colormap(rgb,file->voronoiMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
					glColor4f(rgb[0],rgb[1],rgb[2],file->voronoiMeshGui->overlayAlphaSlider->value());
					glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(h)->nVertices);
				}
				// potential overlay
				else if (file->voronoiMeshGui->overlayPotentialButton->value()) {
					const float pValue = (float) ( file->voronoiMesh->getCell(h)->getPotential() - file->voronoiMesh->getPotentialMin() ) /
										 (float) ( file->voronoiMesh->getPotentialMax() - file->voronoiMesh->getPotentialMin() );
					colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
					glColor4f(rgb[0],rgb[1],rgb[2],file->voronoiMeshGui->overlayAlphaSlider->value());
					glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(h)->nVertices);
				}
				// force magnitude overlay
				else if (file->voronoiMeshGui->overlayForceMagnitudeButton->value()) {
					const float pValue = (float) ( file->voronoiMesh->getCell(h)->getForceMagnitude() - file->voronoiMesh->getForceMin() ) /
										 (float) ( file->voronoiMesh->getForceMax() - file->voronoiMesh->getForceMin() );
					colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
					glColor4f(rgb[0],rgb[1],rgb[2],file->voronoiMeshGui->overlayAlphaSlider->value());
					glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(h)->nVertices);
				}
				// point number overlay
				else if (file->voronoiMeshGui->overlayPointNumberButton->value()) {
					const float pValue = (float) ( file->voronoiMesh->getCell(h)->getCount() - file->voronoiMesh->getMinCell()->getCount() ) /
										 (float) ( file->voronoiMesh->getMaxCell()->getCount() - file->voronoiMesh->getMinCell()->getCount() );
					colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
					glColor4f(rgb[0],rgb[1],rgb[2],file->voronoiMeshGui->overlayAlphaSlider->value());
					glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(h)->nVertices);
				}
			}
			else {
				glColor4f(0.0f,0.0f,0.0f,0.0f);
				glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(h)->nVertices);
			}
			if (file->voronoiMeshGui->displayGridButton->value()) {
				glColor4f(file->voronoiMeshGui->gridRGB[0],file->voronoiMeshGui->gridRGB[1],file->voronoiMeshGui->gridRGB[2],file->voronoiMeshGui->gridAlphaSlider->value());
				glDrawArrays(GL_LINE_LOOP,0,file->voronoiMesh->getCell(h)->nVertices);
			}
		glPopMatrix();
	}
	glDisableClientState(GL_VERTEX_ARRAY);

	// draw neighbour connections
	if (file->voronoiMeshGui->neighbourDistanceViewButton->value()) {
		glColor4f(1.0,1.0,0.0,0.7);
		glLineWidth(2);
		glPushMatrix();
		updateView();

		for (int h = 0; h < file->voronoiMesh->getNumberOfClusters(); h++) {
			if (file->voronoiMesh->getCell(h)->active()) {
				for (int l = 0; l < file->voronoiMesh->getCell(h)->getNLeftNeighbours(); l++) {
					if (file->voronoiMesh->getCell(h)->getLeftNeighbours(l)->active()) {
						glBegin(GL_LINES);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
							glVertex2f(file->voronoiMesh->getCell(h)->getLeftNeighbours(l)->getXCentroid(),
									   file->voronoiMesh->getCell(h)->getLeftNeighbours(l)->getYCentroid());
						glEnd();
					}
				}
				for (int r = 0; r < file->voronoiMesh->getCell(h)->getNRightNeighbours(); r++) {
					if (file->voronoiMesh->getCell(h)->getRightNeighbours(r)->active()) {
						glBegin(GL_LINES);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
							glVertex2f(file->voronoiMesh->getCell(h)->getRightNeighbours(r)->getXCentroid(),
									   file->voronoiMesh->getCell(h)->getRightNeighbours(r)->getYCentroid());
						glEnd();
					}
				}
				for (int t = 0; t < file->voronoiMesh->getCell(h)->getNTopNeighbours(); t++) {
					if (file->voronoiMesh->getCell(h)->getTopNeighbours(t)->active()) {
						glBegin(GL_LINES);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
							glVertex2f(file->voronoiMesh->getCell(h)->getTopNeighbours(t)->getXCentroid(),
									   file->voronoiMesh->getCell(h)->getTopNeighbours(t)->getYCentroid());
						glEnd();
					}
				}
				for (int b = 0; b < file->voronoiMesh->getCell(h)->getNBottomNeighbours(); b++) {
					if (file->voronoiMesh->getCell(h)->getBottomNeighbours(b)->active()) {
						glBegin(GL_LINES);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
							glVertex2f(file->voronoiMesh->getCell(h)->getBottomNeighbours(b)->getXCentroid(),
									   file->voronoiMesh->getCell(h)->getBottomNeighbours(b)->getYCentroid());
						glEnd();
					}
				}
			}
		}
		glPopMatrix();
	}

	if (file->voronoiMeshGui->spotVisualizationButton->value()) {
		glShadeModel(GL_SMOOTH);
		glEnable(GL_POINT_SPRITE);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
//		glEnable(GL_BLEND);
//		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_POINT_SMOOTH);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
		glDisable(GL_DEPTH_TEST);							// Disable Depth Testing

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, iMAP->gaussianTex);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

		for (int h = 0; h < file->voronoiMesh->getNumberOfClusters(); h++) {
//			const float sizeFactor = file->voronoiMeshGui->spotVisualizationScaleSlider->value()*file->voronoiMesh->getCell(h)->getVariance()*file->voronoiMesh->getCell(h)->getVariance();
			const float sizeFactor = file->voronoiMeshGui->spotVisualizationScaleSlider->value()*file->voronoiMesh->getCell(h)->getArea()/file->voronoiMesh->getMaxArea();
			glPointSize(sizeFactor*50.0f);
			glPushMatrix();
				updateView();
				if (file->voronoiMesh->active(h)) {
					float *rgb = new float[3];
					// diffusion overlay
					if (file->voronoiMeshGui->overlayDiffusionButton->value()) {
						const float dValue = (float) ( file->voronoiMesh->getCell(h)->getDiffusion() - file->voronoiMesh->getDiffusionMin() ) /
											 (float) ( file->voronoiMesh->getDiffusionMax() - file->voronoiMesh->getDiffusionMin() );
						colormap(rgb,file->voronoiMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
						glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
						glBegin(GL_POINTS);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
						glEnd();
					}
					// potential overlay
					else if (file->voronoiMeshGui->overlayPotentialButton->value()) {
						const float pValue = (float) ( file->voronoiMesh->getCell(h)->getPotential() - file->voronoiMesh->getPotentialMin() ) /
											 (float) ( file->voronoiMesh->getPotentialMax() - file->voronoiMesh->getPotentialMin() );
						colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
						glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
						glBegin(GL_POINTS);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
						glEnd();
					}
					// force magnitude overlay
					else if (file->voronoiMeshGui->overlayForceMagnitudeButton->value()) {
						const float pValue = (float) ( file->voronoiMesh->getCell(h)->getForceMagnitude() - file->voronoiMesh->getForceMin() ) /
											 (float) ( file->voronoiMesh->getForceMax() - file->voronoiMesh->getForceMin() );
						colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
						glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
						glBegin(GL_POINTS);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
						glEnd();
					}
					// point number overlay
					else if (file->voronoiMeshGui->overlayPointNumberButton->value()) {
						const float pValue = (float) ( file->voronoiMesh->getCell(h)->getCount() - file->voronoiMesh->getMinCell()->getCount() ) /
											 (float) ( file->voronoiMesh->getMaxCell()->getCount() - file->voronoiMesh->getMinCell()->getCount() );
						colormap(rgb,file->voronoiMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
						glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
						glBegin(GL_POINTS);
							glVertex2f(file->voronoiMesh->getCell(h)->getXCentroid(),file->voronoiMesh->getCell(h)->getYCentroid());
						glEnd();
					}
				}
			glPopMatrix();
		}
		glBindTexture(GL_TEXTURE_2D,0);

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_POINT_SPRITE);
//					glDisable(GL_BLEND);
		glDisable(GL_POINT_SMOOTH);
	}

	if (file->voronoiMeshGui->localizationNumberLabelButton->value()) {
		glEnable(GL_TEXTURE_2D);
		txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

		char label[20];
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER,0.5f);
		glColor4f(1.0f,1.0f,1.0f,1.0f);

		const float scale = file->voronoiMeshGui->fontScaleSlider->value()*file->voronoiMesh->getCharDistance()/250.0f;
		for (int h = 0; h < file->voronoiMesh->getNumberOfClusters(); h++) {
			glPushMatrix();
				updateView();
				sprintf(label,"%i",h);
//				sprintf(label,"%i",file->voronoiMesh->getCell(h)->getCount());
				glTranslatef((float)file->voronoiMesh->getCell(h)->getXMean()-20.0f*scale,(float)file->voronoiMesh->getCell(h)->getYMean()-20.0f*scale,0.0f);
				glScalef(scale,scale,scale);
				if (file->voronoiMesh->active(h)) { glColor4f(1.0f,1.0f,1.0f,1.0f); }
				else { glColor4f(1.0f,0.0f,0.0f,1.0f); }
				txfRenderString(iMAP->fontTex, label, strlen(label));
			glPopMatrix();
		}

		glDisable(GL_ALPHA_TEST);
		glBindTexture(GL_TEXTURE_2D,0);
		glDisable(GL_TEXTURE_2D);
	}
	glDisable(GL_BLEND);

	if (file->voronoiMeshGui->overlayForceArrowsButton->value()) {
		clipPlanes();
		drawVoronoiForces();
		unclipPlanes();
	}
}

void drawVoronoiForces() {

//	const float scale = 60.0*file->voronoiMesh->getCharDistance();

	// draw arrows signifying force magnitude and direction
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);

	glPushMatrix();
		updateView();

		float arrowColor[] = {1.0,1.0,0.0,0.7};
		if (file->optimizationMode == 2) { arrowColor[0] = 0.0; arrowColor[1] = 1.0; arrowColor[2] = 0.0; }

		glLineWidth(1.0);
		for (int a = 0; a < file->voronoiMesh->getNumberOfClusters(); a++) {
			if (file->voronoiMesh->getCell(a)->getForceMagnitude() != 0.0f) {
				const double l = 8.0*file->voronoiMesh->getCell(a)->getPerimeter();
				const double fx = file->voronoiMesh->getCell(a)->getForceX();
				const double fy = file->voronoiMesh->getCell(a)->getForceY();

				const double xc = file->voronoiMesh->getCell(a)->getXMean();
				const double yc = file->voronoiMesh->getCell(a)->getYMean();

				glPushMatrix();
					glTranslatef(xc,yc,0.0f);
					glRotatef(180.0f/PI*atan2(fy,fx),0.0f,0.0f,1.0f);
					glTranslatef(-l/140.0f,0.0f,0.0f);
					glRotatef(90.0f,0.0f,0.0f,1.0f);
					glBegin(GL_TRIANGLES);
						glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
						glVertex2f(l*0.005f,0.0f);

						glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
						glVertex2f(-l*0.005f,0.0f);

						glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],1.0);
						glVertex2f(0.0f,-l*0.015f);
					glEnd();
					glBegin(GL_LINE_LOOP);
						glColor4f(1.0f,1.0f,1.0f,0.0f);
						glVertex2f(l*0.005f,0.0f);

						glColor4f(1.0f,1.0f,1.0f,0.0f);
						glVertex2f(-l*0.005f,0.0f);

						glColor4f(1.0f,1.0f,1.0f,0.0f);
						glVertex2f(0.0f,-l*0.015f);
					glEnd();
				glPopMatrix();
			}
		}

	glPopMatrix();
	glDisable(GL_BLEND);
}

void drawVoronoiPosteriorHighlight() {
	if (file->voronoiMeshGui->posterioriGroup->visible() && file->voronoiMeshGui->posterioriGroup->active()) {

		const int a = file->voronoiMesh->getCurrentZone();
		const int a2 = file->voronoiMesh->getCurrentZone2();

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2,GL_FLOAT,0,file->voronoiMesh->getCell(a)->vertices);

		clipPlanes();

		// draw arrows signifying force magnitude and direction
		glPushMatrix();
			updateView();

			glLineWidth(3.0);

			glColor4f(1.0f,1.0f,1.0f,1.0f);
			glDrawArrays(GL_LINE_LOOP,0,file->voronoiMesh->getCell(a)->nVertices);

			glColor4f(1.0f,1.0f,1.0f,0.6f);
			glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(a)->nVertices);

			if (file->voronoiMeshGui->potentialReferenceButton->value()) {
				glBegin(GL_LINES);
					glLineWidth(4.0);
					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glVertex2f(file->voronoiMesh->getCentroidX(a),file->voronoiMesh->getCentroidY(a));
					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glVertex2f(file->voronoiMesh->getCentroidX(a2),file->voronoiMesh->getCentroidY(a2));
				glEnd();

				glVertexPointer(2,GL_FLOAT,0,file->voronoiMesh->getCell(a2)->vertices);

				glColor4f(1.0f,1.0f,1.0f,1.0f);
				glDrawArrays(GL_LINE_LOOP,0,file->voronoiMesh->getCell(a2)->nVertices);

				glColor4f(1.0f,1.0f,1.0f,0.5f);
				glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(a2)->nVertices);
			}
		glPopMatrix();
		glDisableClientState(GL_VERTEX_ARRAY);

		unclipPlanes();

		glDisable(GL_BLEND);

	}
}

void drawVoronoiLandscape() {
	const GLfloat position [] = { file->voronoiMeshGui->xLightPositionSlider->value(),
								  file->voronoiMeshGui->yLightPositionSlider->value(),
								  file->voronoiMeshGui->zLightPositionSlider->value(), 1.0 };
	const GLfloat ambient [] = { 0.10f,0.10f,0.10f,0.10f };
	const GLfloat diffuse [] = { 1.0f,1.0f,1.0f,1.0f };
	const GLfloat specular [] = { 0.6f,0.6f,0.6f,1.0f };

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
//	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE|GL_SPECULAR);

	glMaterialfv(GL_FRONT, GL_AMBIENT,ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR,specular);
	glMaterialf(GL_FRONT,  GL_SHININESS, 128);

	glEnable(GL_DEPTH_TEST);

//	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (file->voronoiMeshGui->fogButton->value()) {
		glEnable(GL_FOG);
		const float fogCol[3]={0.0,0.0,0.0};
		glFogfv(GL_FOG_COLOR,fogCol);
		glFogi(GL_FOG_MODE, GL_LINEAR);
		glFogf(GL_FOG_DENSITY, 0.5f);
		glHint(GL_FOG_HINT,GL_NICEST);
		glFogf(GL_FOG_START, (float)file->voronoiMeshGui->fogStartSlider->value());
		glFogf(GL_FOG_END, (float)file->voronoiMeshGui->fogEndSlider->value());
    }

	glPushMatrix();
		updateView();

		clipPlanes();

		glScalef(1.0,1.0,file->voronoiMeshGui->scaleLandscapeSlider->value());

		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);

		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glNormalPointer(GL_FLOAT,0,file->voronoiMesh->landscapeNormals);
		glVertexPointer(3,GL_FLOAT,0,file->voronoiMesh->landscapeTriangles);
		glColorPointer(4,GL_FLOAT,0,file->voronoiMesh->landscapeColors);

		glDrawArrays(GL_TRIANGLES,0,file->voronoiMesh->landscapeVertices);

		glDisableClientState(GL_COLOR_ARRAY);

		glDisable(GL_DEPTH_TEST);

		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		unclipPlanes();
	glPopMatrix();

	glDisable(GL_FOG);

	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

//	if (iMAP->saveTiffBatch == false) {
//
//		// save screenshot
//		char nativeFilename[FILENAME_MAX];
//		char fileNameTemp[FILENAME_MAX];
//
//		const int fileLength = strlen(file->fileName);
//		for (int g = 0; g < fileLength-6; g++) { fileNameTemp[g] = file->fileName[g]; }
//		fileNameTemp[fileLength-6] = '\0';
//
//		strcpy(nativeFilename,file->filePath);
//		strcat(nativeFilename,fileNameTemp);
//
//		captureScreen(nativeFilename);
//
//		iMAP->saveTiffBatch = true;
//	}
}

void drawVoronoiRandomizedOptimizationZones() {
	const int roZones = file->voronoiMesh->roZones;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// draw voronoi diagram
	glEnableClientState(GL_VERTEX_ARRAY);

	glLineWidth(1.0f);
	for (int h = 0; h < roZones; h++) {
		glPushMatrix();
			updateView();
			glColor4f(1.0f,1.0f,1.0f,0.8f);
			glVertexPointer(2,GL_FLOAT,0,file->voronoiMesh->getCell(file->voronoiMesh->roClusters[h])->vertices);
			glDrawArrays(GL_POLYGON,0,file->voronoiMesh->getCell(file->voronoiMesh->roClusters[h])->nVertices);
		glPopMatrix();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
}

void drawQuadTreeMesh() {

	if (iMAP->overlayAdjustment) {
		// allocate for variables (active zones in quad tree)
		if (file->treeMesh->linesOverlay != NULL) {	delete [] file->treeMesh->linesOverlay;	}
		if (file->treeMesh->quadsOverlay != NULL) {	delete [] file->treeMesh->quadsOverlay;	}
		if (file->treeMesh->quadsColor != NULL) { delete [] file->treeMesh->quadsColor; }

		file->treeMesh->linesOverlay = new float[4*2*2*file->treeMesh->totalVariables];
		file->treeMesh->quadsOverlay = new float[4*2*file->treeMesh->totalVariables];
		file->treeMesh->quadsColor = new float[4*4*2*file->treeMesh->totalVariables];

		// accumulate vertices to respective arrays
		int lineProgress = 0;
		int quadProgress = 0;
		int quadColorProgress = 0;
		createQuadTreeOverlay(file->treeMesh->quadTree,&lineProgress,&quadProgress,&quadColorProgress);

		iMAP->overlayAdjustment = false;
	}

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glPushMatrix();
    	updateView();
		if (file->treeMeshGui->overlayButton->value()) {
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);
			glVertexPointer(2,GL_FLOAT,0,file->treeMesh->quadsOverlay);
			glColorPointer(4,GL_FLOAT,0,file->treeMesh->quadsColor);
			glDrawArrays(GL_QUADS,0,4*file->treeMesh->totalVariables);
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_VERTEX_ARRAY);
		}

		if (file->treeMeshGui->displayGridButton->value()) {
			glLineWidth(1.0);
			glColor4f(file->treeMeshGui->gridRGB[0],file->treeMeshGui->gridRGB[1],file->treeMeshGui->gridRGB[2],file->treeMeshGui->gridAlphaSlider->value());
//			glColor4f(1.0-iMAP->backgroundRGB[0],1.0-iMAP->backgroundRGB[1],1.0-iMAP->backgroundRGB[2],file->treeMeshGui->gridAlphaSlider->value());
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(2,GL_FLOAT,0,file->treeMesh->linesOverlay);
			glDrawArrays(GL_LINES,0,4*2*file->treeMesh->totalVariables);
			glDisableClientState(GL_VERTEX_ARRAY);
		}

		if (file->treeMeshGui->localizationNumberLabelButton->value()) {
			if (file->treeMesh->totalVariables < 200) {
				glEnable(GL_TEXTURE_2D);
				glColor4f(1.0f,1.0f,1.0f,1.0f);
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER,0.5f);
				txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);

				drawQuadTreeLabels(file->treeMesh->quadTree);

				glDisable(GL_ALPHA_TEST);
				glDisable(GL_TEXTURE_2D);
			}
		}

		// draw neighbour connections
		if (file->treeMeshGui->neighbourDistanceViewButton->value()) {
			glColor4f(1.0,1.0,0.0,0.7);
			glLineWidth(2);
			drawQuadTreeNeighbourConnections(file->treeMesh->quadTree);
		}

		if (file->treeMeshGui->overlayForceArrowsButton->value()) {
			clipPlanes();
			drawQuadTreeForces(file->treeMesh->quadTree);
			unclipPlanes();
		}

	glPopMatrix();

    glDisable(GL_BLEND);
}

void drawQuadTreeNeighbourConnections(QuadTree *tree) {
    if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
    	if (tree->active()) {
			// draw connections for each quad of mesh
			for (int l = 0; l < tree->getNumberOfLeftNeighbours(); l++) {
				glBegin(GL_LINES);
					glVertex2f(tree->getXCentroid(),tree->getYCentroid());
					glVertex2f(tree->getLeftNeighbour(l)->getXCentroid(),tree->getLeftNeighbour(l)->getYCentroid());
				glEnd();
			}
			for (int r = 0; r < tree->getNumberOfRightNeighbours(); r++) {
				glBegin(GL_LINES);
					glVertex2f(tree->getXCentroid(),tree->getYCentroid());
					glVertex2f(tree->getRightNeighbour(r)->getXCentroid(),tree->getRightNeighbour(r)->getYCentroid());
				glEnd();
			}
			for (int t = 0; t < tree->getNumberOfTopNeighbours(); t++) {
				glBegin(GL_LINES);
					glVertex2f(tree->getXCentroid(),tree->getYCentroid());
					glVertex2f(tree->getTopNeighbour(t)->getXCentroid(),tree->getTopNeighbour(t)->getYCentroid());
				glEnd();
			}
			for (int b = 0; b < tree->getNumberOfBottomNeighbours(); b++) {
				glBegin(GL_LINES);
					glVertex2f(tree->getXCentroid(),tree->getYCentroid());
					glVertex2f(tree->getBottomNeighbour(b)->getXCentroid(),tree->getBottomNeighbour(b)->getYCentroid());
				glEnd();
			}
    	}
		return;
    }
    drawQuadTreeNeighbourConnections(tree->nw);
    drawQuadTreeNeighbourConnections(tree->ne);
    drawQuadTreeNeighbourConnections(tree->sw);
    drawQuadTreeNeighbourConnections(tree->se);
    return;

}

void createQuadTreeOverlay(QuadTree* tree, int *lineProgress, int *quadProgress, int *quadColorProgress) {
    if (tree == 0) { return; }

    if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
    	if (tree->active()) {
			/*** LINES ***/
			// bottom segment
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.x ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.y ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.y ;
			// right segment
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.y ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;
			// top segment
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.x ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;
			// left segment
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.x ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.x ;
			file->treeMesh->linesOverlay[(*lineProgress)++] = tree->bounds.y ;

			/*** QUADS ***/
			file->treeMesh->quadsOverlay[(*quadProgress)++] = tree->bounds.x ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = tree->bounds.y ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = tree->bounds.y ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = FMIN(tree->bounds.x+tree->bounds.w,file->xMax) ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = tree->bounds.x ;
			file->treeMesh->quadsOverlay[(*quadProgress)++] = FMIN(tree->bounds.y+tree->bounds.w,file->yMax) ;

			/*** QUAD COLOR ***/
			// diffusion overlay
			if (file->treeMeshGui->overlayDiffusionButton->value()) {
				const float dValue = (float) ( tree->getD() - file->treeMesh->getDiffusionMin() ) /
									 (float) ( file->treeMesh->getDiffusionMax() - file->treeMesh->getDiffusionMin() );
				colormap(&file->treeMesh->quadsColor[*quadColorProgress],file->treeMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
				*quadColorProgress += 3;
			}
			// potential overlay
			else if (file->treeMeshGui->overlayPotentialButton->value()) {
				const float pValue = (float) ( tree->getEp() - file->treeMesh->getPotentialMin() ) /
									 (float) ( file->treeMesh->getPotentialMax() - file->treeMesh->getPotentialMin() );
				colormap(&file->treeMesh->quadsColor[*quadColorProgress],file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
				*quadColorProgress += 3;
			}
			// force magnitude overlay
			else if (file->treeMeshGui->overlayForceMagnitudeButton->value()) {
				const float pValue = (float) ( tree->getFmag() - file->treeMesh->getForceMin() ) /
									 (float) ( file->treeMesh->getForceMax() - file->treeMesh->getForceMin() );
				colormap(&file->treeMesh->quadsColor[*quadColorProgress],file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
				*quadColorProgress += 3;
			}
			// # points overlay
			else if (file->treeMeshGui->overlayPointNumberButton->value()) {
				const float pValue = (float) ( tree->getCount() - file->treeMesh->getMinQuadTree()->getCount() ) /
									 (float) ( file->treeMesh->getMaxQuadTree()->getCount() - file->treeMesh->getMinQuadTree()->getCount() );
				colormap(&file->treeMesh->quadsColor[*quadColorProgress],file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
				*quadColorProgress += 3;
			} else {
				file->treeMesh->quadsColor[(*quadColorProgress)++] = 0.0;
				file->treeMesh->quadsColor[(*quadColorProgress)++] = 1.0;
				file->treeMesh->quadsColor[(*quadColorProgress)++] = 0.0;
			}
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMeshGui->overlayAlphaSlider->value(); // 3
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 4
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 5
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 6
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 7
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 8
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 9
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 10
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 11
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 12
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 13
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 14
			file->treeMesh->quadsColor[(*quadColorProgress)++] = file->treeMesh->quadsColor[(*quadColorProgress)-4]; // 15
    	}
		return;
    }

    createQuadTreeOverlay(tree->nw,lineProgress,quadProgress,quadColorProgress);
    createQuadTreeOverlay(tree->ne,lineProgress,quadProgress,quadColorProgress);
    createQuadTreeOverlay(tree->sw,lineProgress,quadProgress,quadColorProgress);
    createQuadTreeOverlay(tree->se,lineProgress,quadProgress,quadColorProgress);

    return;
}

void drawQuadTreeLabels(QuadTree *tree) {
    if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
    	if (tree->active()) {
			// label parameters
			const float scale = file->treeMeshGui->fontScaleSlider->value()*tree->bounds.w/90.0f;
			char label[20];

			// draw string for each quad of mesh
			glPushMatrix();
				updateView();
				sprintf(label,"%i",tree->count);
				glTranslated(tree->bounds.x-file->xTranslate,tree->bounds.y-file->yTranslate,0.0);
				glScalef(scale, scale, scale);
				txfRenderString(iMAP->fontTex, label, strlen(label));
			glPopMatrix();
    	}
		return;
    }
    drawQuadTreeLabels(tree->nw);
    drawQuadTreeLabels(tree->ne);
    drawQuadTreeLabels(tree->sw);
    drawQuadTreeLabels(tree->se);
    return;
}

void drawQuadTreeForces(QuadTree* tree) {
    if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
    	if (tree->active()) {
			// minimum line spacing
			const float dMin = tree->bounds.w * 50.0;

			float arrowColor[] = {1.0,1.0,0.0,0.7};
			if (file->optimizationMode == 2) { arrowColor[0] = 0.0; arrowColor[1] = 1.0; arrowColor[2] = 0.0; }

			// draw arrows signifying force magnitude and direction
			glPushMatrix();
				glLineWidth(1.0);

				if (tree->getFmag() != 0.0f) {
//					const double l = dMin*tree->getFmag()/file->treeMesh->getForceMax();

					const double l = dMin;

					const double fx = tree->getFx();
					const double fy = tree->getFy();

					const double xc = tree->bounds.x+tree->bounds.w/2.0;
					const double yc = tree->bounds.y+tree->bounds.w/2.0;

					glPushMatrix();
						glTranslatef(xc,yc,0.0f);
						glRotatef(180.0f/PI*atan2(fy,fx),0.0f,0.0f,1.0f);
						glTranslatef(-l/140.0f,0.0f,0.0f);
						glRotatef(90.0f,0.0f,0.0f,1.0f);

						glBegin(GL_TRIANGLES);
							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
							glVertex2f(l*0.005f,0.0f);

							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],arrowColor[3]);
							glVertex2f(-l*0.005f,0.0f);

							glColor4f(arrowColor[0],arrowColor[1],arrowColor[2],1.0);
							glVertex2f(0.0f,-l*0.015f);
						glEnd();

						glBegin(GL_LINE_LOOP);
							glColor4f(0.0f,0.0f,0.0f,1.0f);
							glVertex2f(l*0.005f,0.0f);

							glColor4f(0.0f,0.0f,0.0f,1.0f);
							glVertex2f(-l*0.005f,0.0f);

							glColor4f(0.0f,0.0f,0.0f,1.0f);
							glVertex2f(0.0f,-l*0.015f);
						glEnd();

					glPopMatrix();
				}
			glPopMatrix();
    	}
    	return;
    }
    drawQuadTreeForces(tree->nw);
    drawQuadTreeForces(tree->ne);
    drawQuadTreeForces(tree->sw);
    drawQuadTreeForces(tree->se);
    return;
}

void drawQuadTreeSpotVisualization(QuadTree *tree) {
	if (tree == 0) { return; }

	const float sizeFactor = file->treeMeshGui->spotVisualizationScaleSlider->value()/tree->power/tree->power*file->treeMesh->maxQuadTreePower*file->treeMesh->maxQuadTreePower;
	glPointSize(sizeFactor*50.0f);
	if (tree->active()) {
		float *rgb = new float[3];
		// diffusion overlay
		if (file->treeMeshGui->overlayDiffusionButton->value()) {
			const float dValue = (float) ( tree->diffusion - file->treeMesh->getDiffusionMin() ) /
								 (float) ( file->treeMesh->getDiffusionMax() - file->treeMesh->getDiffusionMin() );
			colormap(rgb,file->treeMeshGui->getColormap(),dValue,file->cMin,file->cMax,file->flipColormap);
			glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
			glBegin(GL_POINTS);
				glVertex2f(tree->xAverage,tree->yAverage);
			glEnd();
		}
		// potential overlay
		else if (file->treeMeshGui->overlayPotentialButton->value()) {
			const float pValue = (float) ( tree->potential - file->treeMesh->getPotentialMin() ) /
								 (float) ( file->treeMesh->getPotentialMax() - file->treeMesh->getPotentialMin() );
			colormap(rgb,file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
			glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
			glBegin(GL_POINTS);
				glVertex2f(tree->xAverage,tree->yAverage);
			glEnd();
		}
		// force magnitude overlay
		else if (file->treeMeshGui->overlayForceMagnitudeButton->value()) {
			const float pValue = (float) ( tree->forceMagnitude - file->treeMesh->getForceMin() ) /
								 (float) ( file->treeMesh->getForceMax() - file->treeMesh->getForceMin() );
			colormap(rgb,file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
			glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
			glBegin(GL_POINTS);
				glVertex2f(tree->xAverage,tree->yAverage);
			glEnd();
		}
		// point number overlay
		else if (file->treeMeshGui->overlayPointNumberButton->value()) {
			const float pValue = (float) ( tree->getCount() - file->treeMesh->minCount ) /
								 (float) ( file->treeMesh->maxCount - file->treeMesh->minCount );
			colormap(rgb,file->treeMeshGui->getColormap(),pValue,file->cMin,file->cMax,file->flipColormap);
			glColor4f(rgb[0],rgb[1],rgb[2],1.0f);
			glBegin(GL_POINTS);
				glVertex2f(tree->xAverage,tree->yAverage);
			glEnd();
		}
	}

	drawQuadTreeSpotVisualization(tree->nw);
	drawQuadTreeSpotVisualization(tree->ne);
	drawQuadTreeSpotVisualization(tree->sw);
	drawQuadTreeSpotVisualization(tree->se);
}
void drawQuadtreePosteriorHighlight() {
	if (file->treeMeshGui->posterioriGroup->visible() && file->treeMeshGui->posterioriGroup->active()) {

		if (file->treeMesh->leafSelected) {

			if (file->treeMesh->selectedQuadLeaf->active()) {
				const float x = (float)file->treeMesh->selectedQuadLeaf->bounds.x;
				const float y = (float)file->treeMesh->selectedQuadLeaf->bounds.y;
				const float w = (float)file->treeMesh->selectedQuadLeaf->bounds.w;
				const float h = (float)file->treeMesh->selectedQuadLeaf->bounds.h;

				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

				clipPlanes();

				// draw arrows signifying force magnitude and direction
				glPushMatrix();
					updateView();

					glLineWidth(3.0);

					glColor4f(1.0f,1.0f,1.0f,1.0f);
					glBegin(GL_LINE_LOOP);
						glVertex2f(x,y);
						glVertex2f(x+w,y);
						glVertex2f(x+w,y+h);
						glVertex2f(x,y+h);
					glEnd();

					glColor4f(1.0f,1.0f,1.0f,0.8f);
					glBegin(GL_QUADS);
						glVertex2f(x,y);
						glVertex2f(x+w,y);
						glVertex2f(x+w,y+h);
						glVertex2f(x,y+h);
					glEnd();

					if (file->treeMeshGui->potentialReferenceButton->value() && file->treeMesh->leaf2Selected) {

						const float x2 = (float)file->treeMesh->selectedQuadLeaf2->bounds.x;
						const float y2 = (float)file->treeMesh->selectedQuadLeaf2->bounds.y;
						const float w2 = (float)file->treeMesh->selectedQuadLeaf2->bounds.w;
						const float h2 = (float)file->treeMesh->selectedQuadLeaf2->bounds.h;

						glBegin(GL_LINES);
							glLineWidth(4.0);
							glColor4f(1.0f,1.0f,1.0f,1.0f);
							glVertex2f(file->treeMesh->selectedQuadLeaf->getXCentroid(),file->treeMesh->selectedQuadLeaf->getYCentroid());
							glColor4f(1.0f,1.0f,1.0f,1.0f);
							glVertex2f(file->treeMesh->selectedQuadLeaf2->getXCentroid(),file->treeMesh->selectedQuadLeaf2->getYCentroid());
						glEnd();

						glColor4f(1.0f,1.0f,1.0f,1.0f);
						glBegin(GL_LINE_LOOP);
							glVertex2f(x2,y2);
							glVertex2f(x2+w2,y2);
							glVertex2f(x2+w2,y2+h2);
							glVertex2f(x2,y2+h2);
						glEnd();

						glColor4f(1.0f,1.0f,1.0f,0.20f);
						glBegin(GL_QUADS);
							glVertex2f(x2,y2);
							glVertex2f(x2+w2,y2);
							glVertex2f(x2+w2,y2+h2);
							glVertex2f(x2,y2+h2);
						glEnd();

					}
				glPopMatrix();

				unclipPlanes();

				glDisable(GL_BLEND);
			}
		}
	}
}

void drawQuadTreeRandomizedOptimizationZones(QuadTree *tree) {
    if (tree->nw == 0 && tree->ne == 0 && tree->sw == 0 && tree->se == 0) {
    	if (tree->roActive()) {
			// draw arrows signifying force magnitude and direction
    		glPushMatrix();
    			updateView();
    			glColor4f(1.0f,1.0f,1.0f,0.8f);
    			glBegin(GL_QUADS);
					glVertex2f((float)tree->bounds.x,(float)tree->bounds.y);
					glVertex2f((float)FMIN(tree->bounds.x+tree->bounds.w,file->treeMesh->getXMax()),tree->bounds.y);
					glVertex2f((float)FMIN(tree->bounds.x+tree->bounds.w,file->treeMesh->getXMax()),(float)FMIN(tree->bounds.y+tree->bounds.w,file->treeMesh->getYMax()));
					glVertex2f((float)tree->bounds.x,(float)FMIN(tree->bounds.y+tree->bounds.w,file->treeMesh->getYMax()));
    			glEnd();
    		glPopMatrix();
    	}
    	return;
    }
    drawQuadTreeRandomizedOptimizationZones(tree->nw);
    drawQuadTreeRandomizedOptimizationZones(tree->ne);
    drawQuadTreeRandomizedOptimizationZones(tree->sw);
    drawQuadTreeRandomizedOptimizationZones(tree->se);
    return;
}

void drawQuadTreeRandomizedOptimizationZones() {
	QuadTree *tree;
	for (int h = 0; h < file->treeMesh->getRoZones(); h++) {
		tree = file->treeMesh->getROQuadTree(h);
		glPushMatrix();
			updateView();
			glColor4f(1.0f,1.0f,1.0f,0.8f);
			glBegin(GL_QUADS);
				glVertex2f((float)tree->bounds.x,(float)tree->bounds.y);
				glVertex2f((float)FMIN(tree->bounds.x+tree->bounds.w,file->treeMesh->getXMax()),tree->bounds.y);
				glVertex2f((float)FMIN(tree->bounds.x+tree->bounds.w,file->treeMesh->getXMax()),(float)FMIN(tree->bounds.y+tree->bounds.w,file->treeMesh->getYMax()));
				glVertex2f((float)tree->bounds.x,(float)FMIN(tree->bounds.y+tree->bounds.w,file->treeMesh->getYMax()));
			glEnd();
		glPopMatrix();
	}
}

void billboardBegin() {

	glDisable(GL_LIGHTING);

	float modelview[16];
	int i,j;

	// save the current modelview matrix
	glPushMatrix();

	// get the current modelview matrix
	glGetFloatv(GL_MODELVIEW_MATRIX , modelview);

	// undo all rotations
	// beware all scaling is lost as well
	for( i=0; i<3; i++ )
	    for( j=0; j<3; j++ ) {
		if ( i==j )
		    modelview[i*4+j] = 1.0;
		else
		    modelview[i*4+j] = 0.0;
	    }

	// set the modelview with no rotations
	glLoadMatrixf(modelview);
}

void billboardEnd() {
	// restore the previously
	// stored modelview matrix
	glPopMatrix();
	glEnable(GL_LIGHTING);
}

void drawVariables() {
	glViewport(0,iMAP->tabBarWidth/1.2,iMAP->tabBarWidth/4,iMAP->tabBarWidth/4);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_TEXTURE_2D);
	txfEstablishTexture(iMAP->fontTex, 0, GL_TRUE);
	char label[100];
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER,0.5f);

	const float scale = 0.1f*file->orthoLimit;
	printf("scale = %f\n",scale);
	glPushMatrix();
//		billboardBegin();
		glColor4f(1.0f,1.0f,1.0f,1.0f);
		sprintf(label,"%i Variables",iMAP->activeZones);
		glTranslatef((float)-iMAP->tabBarWidth/220,0,0.0f);
		glScalef(scale,scale,scale);
		txfRenderString(iMAP->fontTex, label, strlen(label));
//		billboardEnd();
	glPopMatrix();

	glDisable(GL_ALPHA_TEST);
	glBindTexture(GL_TEXTURE_2D,0);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glViewport(0,0,iMAP->tabBarWidth,iMAP->tabBarWidth);

}

void drawCircle(float cx, float cy, float r, int num_segments) {
	glBegin(GL_LINE_LOOP);
	for(int ii = 0; ii < num_segments; ii++)
	{
		const float theta = 2.0f * PI * float(ii) / float(num_segments);//get the current angle

		const float x = r * cosf(theta);//calculate the x component
		const float y = r * sinf(theta);//calculate the y component

		glVertex2f(x + cx, y + cy);//output vertex

	}
	glEnd();
}

void drawSquareLandscape() {
	const GLfloat position [] = { file->squareMeshGui->xLightPositionSlider->value(),
								  file->squareMeshGui->yLightPositionSlider->value(),
								  file->squareMeshGui->zLightPositionSlider->value(), 1.0 };
	const GLfloat ambient [] = { 0.10f,0.10f,0.10f,0.10f };
	const GLfloat diffuse [] = { 1.0f,1.0f,1.0f,1.0f };
	const GLfloat specular [] = { 0.6f,0.6f,0.6f,1.0f };

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
//	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE|GL_SPECULAR);

	glMaterialfv(GL_FRONT, GL_AMBIENT,ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR,specular);
	glMaterialf(GL_FRONT,  GL_SHININESS, 128);

	glEnable(GL_DEPTH_TEST);

//	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (file->squareMeshGui->fogButton->value()) {
		glEnable(GL_FOG);
		const float fogCol[3]={0.0,0.0,0.0};
		glFogfv(GL_FOG_COLOR,fogCol);
		glFogi(GL_FOG_MODE, GL_LINEAR);
		glFogf(GL_FOG_DENSITY, 0.5f);
		glHint(GL_FOG_HINT,GL_NICEST);
		glFogf(GL_FOG_START, (float)file->squareMeshGui->fogStartSlider->value());
		glFogf(GL_FOG_END, (float)file->squareMeshGui->fogEndSlider->value());
    }

	glPushMatrix();
		updateView();

		clipPlanes();

		glScalef(1.0,1.0,file->squareMeshGui->scaleLandscapeSlider->value());
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);

		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glNormalPointer(GL_FLOAT,0,file->squareMesh->landscapeNormals);
		glVertexPointer(3,GL_FLOAT,0,file->squareMesh->landscapeTriangles);
		glColorPointer(4,GL_FLOAT,0,file->squareMesh->landscapeColors);

		glDrawArrays(GL_TRIANGLES,0,4*3*file->squareMesh->getXCells()*file->squareMesh->getYCells());

		glDisableClientState(GL_COLOR_ARRAY);

		glDisable(GL_DEPTH_TEST);

		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		unclipPlanes();
	glPopMatrix();

	glDisable(GL_FOG);

//	glDisable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
//	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

}

void drawQuadTreeLandscape() {
	const GLfloat position [] = { file->treeMeshGui->xLightPositionSlider->value(),
								  file->treeMeshGui->yLightPositionSlider->value(),
								  file->treeMeshGui->zLightPositionSlider->value(), 1.0 };
	const GLfloat ambient [] = { 0.10f,0.10f,0.10f,0.10f };
	const GLfloat diffuse [] = { 1.0f,1.0f,1.0f,1.0f };
	const GLfloat specular [] = { 0.6f,0.6f,0.6f,1.0f };

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
//	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
//	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE|GL_SPECULAR);

	glMaterialfv(GL_FRONT, GL_AMBIENT,ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,diffuse);
//	glMaterialfv(GL_FRONT, GL_SPECULAR,specular);
//	glMaterialf(GL_FRONT,  GL_SHININESS, 128);

	glEnable(GL_DEPTH_TEST);

//	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (file->treeMeshGui->fogButton->value()) {
		glEnable(GL_FOG);
		const float fogCol[3]={0.0,0.0,0.0};
		glFogfv(GL_FOG_COLOR,fogCol);
		glFogi(GL_FOG_MODE, GL_LINEAR);
		glFogf(GL_FOG_DENSITY, 0.5f);
		glHint(GL_FOG_HINT,GL_NICEST);
		glFogf(GL_FOG_START, (float)file->treeMeshGui->fogStartSlider->value());
		glFogf(GL_FOG_END, (float)file->treeMeshGui->fogEndSlider->value());
    }

	glPushMatrix();
		updateView();

		clipPlanes();

		glScalef(1.0,1.0,file->treeMeshGui->scaleLandscapeSlider->value());

		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		glNormalPointer(GL_FLOAT,0,file->treeMesh->landscapeNormals);
		glVertexPointer(3,GL_FLOAT,0,file->treeMesh->landscapeTriangles);
		glColorPointer(4,GL_FLOAT,0,file->treeMesh->landscapeColors);
		glDrawArrays(GL_TRIANGLES,0,file->treeMesh->landscapeVertices);

		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_DEPTH_TEST);

		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		unclipPlanes();

	glPopMatrix();

	glDisable(GL_FOG);

//	glDisable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_NORMALIZE);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

}

void clipPlanes() {
	// clip the 6 faces of the cube
	GLdouble plane[4] ;

	plane[0] = +1. ;  plane[1] =  0. ;  plane[2] =  0. ;  plane[3] = -file->xMin ;
	glEnable( GL_CLIP_PLANE0 ) ;
	glClipPlane( GL_CLIP_PLANE0, plane ) ;

	plane[0] = -1. ;  plane[1] =  0. ;  plane[2] =  0. ;  plane[3] = file->xMax ;
	glEnable( GL_CLIP_PLANE1 ) ;
	glClipPlane( GL_CLIP_PLANE1, plane ) ;

	plane[0] =  0. ;  plane[1] = +1. ;  plane[2] =  0. ;  plane[3] = -file->yMin ;
	glEnable( GL_CLIP_PLANE2 ) ;
	glClipPlane( GL_CLIP_PLANE2, plane ) ;

	plane[0] =  0. ;  plane[1] = -1. ;  plane[2] =  0. ;  plane[3] = file->yMax ;
	glEnable( GL_CLIP_PLANE3 ) ;
	glClipPlane( GL_CLIP_PLANE3, plane ) ;
}

void unclipPlanes() {
	  // disable cube clip plane
	  glDisable( GL_CLIP_PLANE0 ) ;
	  glDisable( GL_CLIP_PLANE1 ) ;
	  glDisable( GL_CLIP_PLANE2 ) ;
	  glDisable( GL_CLIP_PLANE3 ) ;
}

void drawTrajectories3D() {
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLineWidth(2.0);
	glEnable(GL_DEPTH_TEST);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	if (file->squareMeshOverlay) {
		if (file->squareMeshGui->fogButton->value()) {
			glEnable(GL_FOG);
			const float fogCol[3]={0.0,0.0,0.0};
			glFogfv(GL_FOG_COLOR,fogCol);
			glFogi(GL_FOG_MODE, GL_LINEAR);
			glFogf(GL_FOG_DENSITY, 0.5f);
			glHint(GL_FOG_HINT,GL_NICEST);
			glFogf(GL_FOG_START, (float)file->squareMeshGui->fogStartSlider->value());
			glFogf(GL_FOG_END, (float)file->squareMeshGui->fogEndSlider->value());
		}
	} else if (file->voronoiMeshOverlay) {
		if (file->voronoiMeshGui->fogButton->value()) {
			glEnable(GL_FOG);
			const float fogCol[3]={0.0,0.0,0.0};
			glFogfv(GL_FOG_COLOR,fogCol);
			glFogi(GL_FOG_MODE, GL_LINEAR);
			glFogf(GL_FOG_DENSITY, 0.5f);
			glHint(GL_FOG_HINT,GL_NICEST);
			glFogf(GL_FOG_START, (float)file->voronoiMeshGui->fogStartSlider->value());
			glFogf(GL_FOG_END, (float)file->voronoiMeshGui->fogEndSlider->value());
		}
	} else if (file->treeMeshOverlay) {
		if (file->treeMeshGui->fogButton->value()) {
			glEnable(GL_FOG);
			const float fogCol[3]={0.0,0.0,0.0};
			glFogfv(GL_FOG_COLOR,fogCol);
			glFogi(GL_FOG_MODE, GL_LINEAR);
			glFogf(GL_FOG_DENSITY, 0.5f);
			glHint(GL_FOG_HINT,GL_NICEST);
			glFogf(GL_FOG_START, (float)file->treeMeshGui->fogStartSlider->value());
			glFogf(GL_FOG_END, (float)file->treeMeshGui->fogEndSlider->value());
		}
	}
	glPushMatrix();
		updateView();
		zSortLines();
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
		glVertexPointer(3,GL_FLOAT,0,file->linesArray);
		glColorPointer(4,GL_FLOAT,0,file->linesColorArray);
		glDrawArrays(GL_LINES,0,2*(file->localizationCountDisplay-file->numberOfFiles-1));

	glPopMatrix();

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisable(GL_FOG);

	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
}

void animateTrajectories3D() {
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(3);

	glPushMatrix();
		updateView();
		glGetDoublev(GL_MODELVIEW_MATRIX,iMAP->mv);
		if (file->numberOfFiles == 1) {
			glBegin(GL_LINES);
			float *rgb = new float[3];
			for (int p = 0; p < file->animationIndex-1; p++) {
				rgb = colormap(9,(double)(p+1)/(double)(file->animationIndex),file->cMin,file->cMax,false);
				if (iMAP->animateAccumulationButton->value()) {
					glColor4f(rgb[0],rgb[1],rgb[2],(double)(p+1)/(double)(file->animationIndex));
					glVertex3f(file->xPointer[p],file->yPointer[p],file->zPointer[p]);
					glVertex3f(file->xPointer[p+1],file->yPointer[p+1],file->zPointer[p+1]);
				} else {
					if (file->animationIndex-p < iMAP->animateDelayedStepsSlider->value()) {
						glColor4f(rgb[0],rgb[1],rgb[2],1.0);
						glVertex3f(file->xPointer[p],file->yPointer[p],file->zPointer[p]);
						glVertex3f(file->xPointer[p+1],file->yPointer[p+1],file->zPointer[p+1]);
					}
				}
			}
			glEnd();
		} else {
			int m = 0;
			for (int v = 0; v < file->numberOfFiles-1; v++) {
				glBegin(GL_LINES);
					for (int k = 0; k < file->tracks[v].n-1; k++) {
						if (file->tracks[v].d[k].t <= (double)file->animationIndex*file->exposureTime) {
							if (file->drawTrajectories) {
								glColor4f(1.0f,1.0f,1.0f,1.0f);
								glVertex3f(file->tracks[v].d[k].x,file->tracks[v].d[k].y,file->tracks[v].d[k].z);
								glVertex3f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y,file->tracks[v].d[k+1].z);
							}
							else {
								if (iMAP->animateAccumulationButton->value()) {
									glColor4f(file->tracks[v].rgb[0],file->tracks[v].rgb[1],file->tracks[v].rgb[2],pow(2000.0*file->tracks[v].d[0].t/(float)file->animationIndex*iMAP->exposureTime,2.0));
									glVertex3f(file->tracks[v].d[k].x,file->tracks[v].d[k].y,file->tracks[v].d[k].z);
									glVertex3f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y,file->tracks[v].d[k+1].z);
								} else {
									if (file->animationIndex-(int)(file->tracks[v].d[k].t/file->exposureTime) < iMAP->animateDelayedStepsSlider->value()) {
										glColor4f(file->tracks[v].rgb[0],file->tracks[v].rgb[1],file->tracks[v].rgb[2],1.0);
										glVertex3f(file->tracks[v].d[k].x,file->tracks[v].d[k].y,file->tracks[v].d[k].z);
										glVertex3f(file->tracks[v].d[k+1].x,file->tracks[v].d[k+1].y,file->tracks[v].d[k+1].z);
									}
								}
							}
						}
						m++;
					}
				glEnd();
			}
		}
	glPopMatrix();

	glDisable(GL_BLEND);
}

void zSortLocalizations() {
	// adjust zcoordinates according to modelview matrix
	file->coordsTemp = new float[file->localizationCountDisplay];
	for (int l = 0; l < file->localizationCountDisplay; l++) {
		file->coordsTemp[l] = iMAP->mv[2]*(float)file->xPointer[l] + iMAP->mv[6]*file->yPointer[l] + iMAP->mv[10]*(float)file->zPointer[l] + iMAP->mv[14];
		file->indexArray[l] = l;
	}
	qsort(file->indexArray,file->localizationCountDisplay,sizeof(int),compareZ);
	delete [] file->coordsTemp;
}

void zSortLines() {
	// adjust zcoordinates according to modelview matrix
	file->coordsTemp = new float[(file->localizationCountDisplay-file->numberOfFiles)];
	for (int l = 0; l < (file->localizationCountDisplay-file->numberOfFiles); l++) {
		file->coordsTemp[l] = iMAP->mv[2]*(float)file->linesArray[6*l] + iMAP->mv[6]*file->linesArray[6*l+1] + iMAP->mv[10]*(float)file->linesArray[6*l+2] + iMAP->mv[14];
		file->linesIndexArray[l] = l;
	}
	qsort(file->linesIndexArray,(file->localizationCountDisplay-file->numberOfFiles),sizeof(int),compareZ);
	delete [] file->coordsTemp;
}

int compareZ(const void* a, const void* b) {
	float fa = file->coordsTemp[*(int*)a];
	float fb = file->coordsTemp[*(int*)b];
	return (fa > fb) - (fa < fb);
}
