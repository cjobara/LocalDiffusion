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

#include "selection.h"

#include "file.h"
#include "globals.h"
#include "draw.h"
#include "text.h"
#include "inference.h"

#define GTOL 1.e-16
#define K_R 4.e-7
#define PI 3.1415926535897932384626433832795028841971693993751

extern File *file;
extern Globals *iMAP;

Selection::Selection(double x, double y) {
	head = new Node();
	head->xCoord = x-(double)file->xTranslate;
	head->yCoord = y-(double)file->yTranslate;
	head->next = NULL;
	currentNode = NULL;
	tail = head;
	vertices = NULL;
	n = 0;
	optimizationArray[0] = 0.0;
	optimizationArray[1] = 0.0;
	optimizationArray[2] = 0.0;
	//cell = NULL;
	tobj = NULL;
	xVertexCentre = 0.0;
	yVertexCentre = 0.0;
}

void Selection::add(double x, double y) {
	Node *newNode = new Node();
	newNode->xCoord = x-(double)file->xTranslate;
	newNode->yCoord = y-(double)file->yTranslate;
	newNode->next = NULL;

	currentNode = head;
	while (currentNode) {
		if (currentNode->next == NULL) {
			currentNode->next = newNode;
			tail = currentNode->next;
			n++;
			return;
		}
		currentNode = currentNode->next;
	}
}

void Selection::createArray() {
	if (n > 0) {
		n++;
		n++;

		vertices = NULL;
		vertices = new double[2*n];
		
		Node *temp = head;
		int a = 0;
		while (temp != NULL) {
			vertices[2*a] = temp->xCoord;
			vertices[2*a+1] = temp->yCoord;
			temp = temp->next;
			a++;
		}
		vertices[2*(n-1)] = vertices[0];
		vertices[2*(n-1)+1] = vertices[1];
		//vertices[2*a+2] = 0.0;
	}
}

bool Selection::farEnough(double x, double y) {
	const double minDistance = 0.01;
	if ( sqrt((tail->xCoord-x)*(tail->xCoord-x)+(tail->yCoord-y)*(tail->yCoord-y)) - minDistance > 0.000001) { return true; }
	else { return false; }
}

void Selection::assignPoints() {
	createArray();
	
	// count number of points
	int points = 0;
	int i,j,c;
	for (int b = 0; b < file->localizationCount; b++) {
		c = 0;
		const double x = file->xPointer[b];
		const double y = file->yPointer[b];
		for (i = 0, j = n-1; i < n; j = i++) {
			if ( ((vertices[2*i+1] >= y) != (vertices[2*j+1] >= y)) &&
				(x <= (vertices[2*j]-vertices[2*i]) * (y-vertices[2*i+1]) / (vertices[2*j+1]-vertices[2*i+1]) + vertices[2*i]) ) {
				c = !c;
			}
		}
		points += c;
	}
	
	// assign values	
	//cell = new SelectionCell();
	//SelectionCell cell;
	cell.count = points;
	
	//Sleep(50);
	cell.xPointer = new double[points];
	cell.yPointer = new double[points];
	cell.indexArray = new int[points];
	cell.iPointer = new double[points];
	cell.tPointer = new double[points];
	
	int f = 0;
	for (int b = 0; b < file->localizationCount; b++) {
		c = 0;
		const double x = file->xPointer[b];
		const double y = file->yPointer[b];
		for (i = 0, j = n-1; i < n; j = i++) {
			if ( ((vertices[2*i+1] >= y) != (vertices[2*j+1] >= y)) &&
				(x <= (vertices[2*j]-vertices[2*i]) * (y-vertices[2*i+1]) / (vertices[2*j+1]-vertices[2*i+1]) + vertices[2*i]) ) {
				c = !c;
			}
		}
		if (c) {
			cell.xPointer[f] = x;
			cell.yPointer[f] = y;
			cell.iPointer[f] = file->iPointer[b];
			cell.indexArray[f] = b;
			cell.tPointer[f] = file->tPointer[b];
			f++;
		}
	}
	
	// determine initial values
	cell.averageDt = 0.0;
	cell.averageDx = 0.0;
	cell.averageDy = 0.0;
	for (int g = 0; g < points-1; g++) {
		cell.averageDt += fabs(cell.tPointer[g+1]-cell.tPointer[g]);
		cell.averageDx += fabs(cell.xPointer[g+1]-cell.xPointer[g]);
		cell.averageDy += fabs(cell.yPointer[g+1]-cell.yPointer[g]);
	}
	cell.averageDt /= points-1;
	cell.averageDx /= points-1;
	cell.averageDy /= points-1;
	
	// calculate area and perimter
	// calculate perimeter and area of each cell
	cell.perimeter = 0.0;
	for (int w = 0; w < n; w++) {
		if (w < n-1) {
			cell.perimeter += sqrt( (vertices[2*w]-vertices[2*(w+1)])*(vertices[2*w]-vertices[2*(w+1)]) +
							   (vertices[2*w+1]-vertices[2*(w+1)+1])*(vertices[2*w+1]-vertices[2*(w+1)+1]) );
		} else {
			cell.perimeter += sqrt( (vertices[2*w]-vertices[0])*(vertices[2*w]-vertices[0]) +
							   (vertices[2*w+1]-vertices[1])*(vertices[2*w+1]-vertices[1]) );
		}


		xVertexCentre += vertices[2*w];
		yVertexCentre += vertices[2*w+1];
	}

	cell.area = 0.0;

	int m = n-1;
	for (int i=0; i<n ; i++){
		cell.area += (vertices[2*m] + vertices[2*i])*(vertices[2*m+1] - vertices[2*i+1]);
		m = i;
	}
	cell.area *= 0.5;
	cell.area = fabs(cell.area);
	
	xVertexCentre /= (double)n;
	yVertexCentre /= (double)n;

	// calculate bounds based on selection region
	for (int j = 0; j < cell.count; j++) {
		if (cell.xPointer[j] > cell.xMax) { cell.xMax = cell.xPointer[j]; }
		if (cell.xPointer[j] < cell.xMin) { cell.xMin = cell.xPointer[j]; }
		if (cell.yPointer[j] > cell.yMax) { cell.yMax = cell.yPointer[j]; }
		if (cell.yPointer[j] < cell.yMin) { cell.yMin = cell.yPointer[j]; }
	}

	cell.xRange = cell.xMax-cell.xMin;
	cell.yRange = cell.yMax-cell.yMin;
	
	tessellate();

	char label[100];
	sprintf(label,"%i",cell.count);
	iMAP->customSelectionInferenceGui->pointsBox->labelfont(0);
	iMAP->customSelectionInferenceGui->pointsBox->copy_label(label);

}

void Selection::inferDF() {
	iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

	iMAP->updateDisplayButton->value(1);
	iMAP->updateDisplayButton->do_callback();

	const double D_eff_x = cell.averageDx*cell.averageDx/cell.averageDt;
	const double D_eff_y = cell.averageDy*cell.averageDy/cell.averageDt; // valeur initial de diffusion

	optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);
	optimizationArray[1] = optimizationArray[2] = 0.0;

	dfpmin(optimizationArray, 3, GTOL, iterations, fret, customSelectionDFPosterior, (dfunc));

	cell.diffusion = optimizationArray[0];
	cell.xForce = optimizationArray[1];
	cell.yForce = optimizationArray[2];
	cell.forceMagnitude = sqrt(optimizationArray[1]*optimizationArray[1] + optimizationArray[2]*optimizationArray[2]);

	char label[100];

	sprintf(label,"%.4f",optimizationArray[0]);
	iMAP->customSelectionInferenceGui->diffusionBox->labelfont(0);
	iMAP->customSelectionInferenceGui->diffusionBox->copy_label(label);

	sprintf(label,"{%.2f,%.2f}",optimizationArray[1],optimizationArray[2]);
	iMAP->customSelectionInferenceGui->forceBox->labelfont(0);
	iMAP->customSelectionInferenceGui->forceBox->copy_label(label);

	sprintf(label,"%.3f",sqrt(optimizationArray[1]*optimizationArray[1]+optimizationArray[2]*optimizationArray[2]));
	iMAP->customSelectionInferenceGui->forceMagnitudeBox->labelfont(0);
	iMAP->customSelectionInferenceGui->forceMagnitudeBox->copy_label(label);

	sprintf(label,"%.3f",cell.area);
	iMAP->customSelectionInferenceGui->areaBox->labelfont(0);
	iMAP->customSelectionInferenceGui->areaBox->copy_label(label);

	iMAP->customSelectionInferenceGui->forceLabelBox->copy_label("Force [pN]");
	iMAP->customSelectionInferenceGui->forceMagnitudeLabelBox->copy_label("Force Magnitude [pN]");
}

void Selection::inferDDr() {
	iMAP->textBuffer->remove(0,iMAP->textBuffer->length());

	iMAP->updateDisplayButton->value(1);
	iMAP->updateDisplayButton->do_callback();

	const double D_eff_x = cell.averageDx*cell.averageDx/cell.averageDt;
	const double D_eff_y = cell.averageDy*cell.averageDy/cell.averageDt; // valeur initial de diffusion

	optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);
	optimizationArray[1] = optimizationArray[2] = 0.0;

	dfpmin(optimizationArray, 3, GTOL, iterations, fret, customSelectionDDrPosterior, (dfunc));

	cell.diffusion = optimizationArray[0];
	cell.xForce = optimizationArray[1];
	cell.yForce = optimizationArray[2];
	cell.forceMagnitude = sqrt(optimizationArray[1]*optimizationArray[1] + optimizationArray[2]*optimizationArray[2]);

	char label[100];

	sprintf(label,"%.4f",optimizationArray[0]);
	iMAP->customSelectionInferenceGui->diffusionBox->labelfont(0);
	iMAP->customSelectionInferenceGui->diffusionBox->copy_label(label);

	sprintf(label,"{%.2f,%.2f}",optimizationArray[1],optimizationArray[2]);
	iMAP->customSelectionInferenceGui->forceBox->labelfont(0);
	iMAP->customSelectionInferenceGui->forceBox->copy_label(label);

	sprintf(label,"%.3f",sqrt(optimizationArray[1]*optimizationArray[1]+optimizationArray[2]*optimizationArray[2]));
	iMAP->customSelectionInferenceGui->forceMagnitudeBox->labelfont(0);
	iMAP->customSelectionInferenceGui->forceMagnitudeBox->copy_label(label);

	sprintf(label,"%.3f",cell.area);
	iMAP->customSelectionInferenceGui->areaBox->labelfont(0);
	iMAP->customSelectionInferenceGui->areaBox->copy_label(label);

	iMAP->customSelectionInferenceGui->forceLabelBox->copy_label("Drift [um/s]");
	iMAP->customSelectionInferenceGui->forceMagnitudeLabelBox->copy_label("Drift Magnitude [um/s]");
}

void Selection::tessellate() {

	// Create a new tessellation object
	tobj = gluNewTess();

	// Set callback functions
	gluTessCallback(tobj, GLU_TESS_BEGIN, (void (CALLBACK *)()) &beginCallback);
	gluTessCallback(tobj, GLU_TESS_END, (void (CALLBACK *)()) &endCallback);
	gluTessCallback(tobj, GLU_TESS_VERTEX, (void (CALLBACK *)()) &vertexCallback);
	gluTessCallback(tobj, GLU_TESS_COMBINE, (void (CALLBACK *)()) &combineCallback);
	gluTessCallback(tobj, GLU_TESS_ERROR, (void (CALLBACK *)()) &errorCallback);

	// Set winding rule
	gluTessProperty(tobj, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_ODD);
	gluTessNormal(tobj, 0.0, 0.0, 1.0);
}

void Selection::draw() {
	// Begin polygon
	gluTessBeginPolygon(tobj, NULL);

		// Begin contour
		gluTessBeginContour(tobj);

			// Render contour
			for (int i = 0; i < n-1; i++)	{
				gluTessVertex(tobj,&vertices[2*i],(void*)&vertices[2*i]);
			}

		// End contour
		gluTessEndContour(tobj);

	// End polygon
	gluTessEndPolygon(tobj);
}

void Selection::drawForce() {

	float xc = (float)xVertexCentre;
	float yc = (float)yVertexCentre;
	float zc = -0.5;

	float fx = (float)this->cell.xForce;
	float fy = (float)this->cell.yForce;
	float fz = 0.0;
//	float mag = (float)this->cell.forceMagnitude*(file->maxRange*file->maxRange/750.0);

	const GLfloat mat_ambient[]    = { 0.20f,0.20f,0.20f,1.0f };  // RGBA
	const GLfloat mat_diffuse[]    = { 1.0f,1.0f,1.0,1.0f };  // RGBA
	const GLfloat mat_specular[]   = { 1.0f,1.0f,1.0,1.0f };  // RGBA
	const GLfloat light_position0[] = { 0.0, 0.0, 50.0, 0.0 };  // XYZ

	glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);

	glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
	glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
	glMaterialf(GL_FRONT,  GL_SHININESS, 100.0);

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, mat_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, mat_ambient);
	glEnable(GL_LIGHTING);

	if (iMAP->customSelectionInferenceGui->DDrinferred) { glColor4f(0.0,1.0,0.0,0.8f); }
	else if (iMAP->customSelectionInferenceGui->DFinferred) { glColor4f(1.0,1.0,0.0,0.8f); }


	const double RADPERDEG = 0.0174533;

	GLUquadricObj *quadObj;

	glPushMatrix ();

	const float mag = sqrt(cell.xRange*cell.xRange+cell.yRange*cell.yRange)/4.0;

	glTranslated(xc,yc,zc);

		if ((fx!=0.)||(fy!=0.)) {
			glRotated(atan2(fy,fx)/RADPERDEG,0.,0.,1.);
			glRotated(atan2(sqrt(fx*fx+fy*fy),fz)/RADPERDEG,0.,1.,0.);
		}
		else if (fz<0){ glRotated(180,1.,0.,0.); }

		quadObj = gluNewQuadric ();
		gluQuadricDrawStyle (quadObj, GLU_FILL);
		gluQuadricNormals (quadObj, GLU_SMOOTH);
		gluCylinder(quadObj, 0.5*mag, 0.0, mag, 20, 20);
		gluDeleteQuadric(quadObj);

		quadObj = gluNewQuadric ();
		gluQuadricDrawStyle (quadObj, GLU_FILL);
		gluQuadricNormals (quadObj, GLU_SMOOTH);
		gluDisk(quadObj, 0.0, 0.5*mag, 32, 1);
		gluDeleteQuadric(quadObj);

		glTranslatef(0.0,0.0,-mag);
		glPushMatrix ();
			quadObj = gluNewQuadric ();
			gluQuadricDrawStyle (quadObj, GLU_FILL);
			gluQuadricNormals (quadObj, GLU_SMOOTH);
			gluCylinder(quadObj, mag/4.0, mag/4.0, mag, 20, 20);
			gluDeleteQuadric(quadObj);
		glPopMatrix ();

	glPopMatrix ();

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_DEPTH_TEST);

}

void CALLBACK combineCallback(GLdouble coords[3], GLdouble *vertex_data[4], GLdouble weight[4], GLdouble **dataOut) {
	GLdouble *vertex;

	vertex = (GLdouble *) malloc(6 * sizeof(GLdouble));

	vertex[0] = coords[0];
	vertex[1] = coords[1];
	vertex[2] = 0.0;

	vertex[3] = 0.0;
	vertex[4] = 0.0;
	vertex[5] = 0.0;

	*dataOut = vertex;
}

void CALLBACK vertexCallback(GLvoid *vertex) {
	const GLdouble *ptr = (const GLdouble*)vertex;
//	glColor3dv(ptr+3);
	glVertex2dv(ptr);
}

void CALLBACK errorCallback(GLenum errorCode) {
	   const GLubyte *estring;
	   estring = gluErrorString(errorCode);
	   fprintf(stderr, "%s\n", estring);
}

void CALLBACK beginCallback(GLenum which) {
    glBegin(which);
}

void CALLBACK endCallback() {
    glEnd();
}

Selection::~Selection() {
	//if (head->next != NULL) {
	//	Node *temp = head->next;
	//	while (temp != NULL) {
	//		head->next = temp->next;
	//		temp->next = NULL;
	//		delete temp;
	//		temp = head->next;
	//	}
	//	if (vertices != NULL) {
	//		delete [] vertices;
	//		vertices = NULL;
	//	}
	//}
}
