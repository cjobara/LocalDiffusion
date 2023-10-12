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

#ifndef SELECTION_H_
#define SELECTION_H_

#include <stdio.h>
#include <iostream>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>

// OpenGL Libraries
#ifdef __APPLE__
	#include <GL/glew.h>
	#include <OpenGL/OpenGL.h>
//	#include <OpenGL/gl.h>
//	#include <OpenGL/glu.h>
	#include <GLUT/glut.H>
	#include <FL/gl.h>
#elif _WIN32
	#include <GL/glew.h>
	#include <GL/glext.h>
	#include <glut.h>
#endif

class Selection;

void CALLBACK errorCallback(GLenum errorCode);
void CALLBACK combineCallback(GLdouble coords[3], GLdouble *vertex_data[4], GLdouble weight[4], GLdouble **dataOut);
void CALLBACK vertexCallback(GLvoid *vertex);
void CALLBACK beginCallback(GLenum which);
void CALLBACK endCallback();

struct Node {
	double xCoord;
	double yCoord;
	Node *next;
	Node() {
		xCoord = 0.0;
		yCoord = 0.0;
		next = NULL;
	}
};

struct SelectionCell {
	int count;
	double *xPointer;
	double *yPointer;
	int *indexArray;
	double *iPointer;
	double *tPointer;
	float area;
	float perimeter;
	double xForce;
	double yForce;
	double forceMagnitude;
	double diffusion;
	double averageDx;
	double averageDy;
	double averageDt;
	double xMin,xMax;
	double yMin,yMax;
	double xRange,yRange;
	SelectionCell() {
		count = 0;
		xPointer = NULL;
		yPointer = NULL;
		indexArray = NULL;
		iPointer = NULL;
		tPointer = NULL;
		area = 0.0;
		perimeter = 0.0;
		xForce = 0.0;
		yForce = 0.0;
		forceMagnitude = 0.0;
		diffusion = 0.0;
		averageDx = 0.0;
		averageDy = 0.0;
		averageDt = 0.0;
		xMin = yMin = 10000000.0;
		xMax = yMax = -10000000.0;
		xRange = yRange = 0.0;
	}
};

class Selection {
	public:
		Selection(double x, double y);

		// linked list treatment
		Node *head;
		Node *tail;
		Node *currentNode;
		void add(double x, double y);
		void createArray();
		bool farEnough(double x, double y);
		void assignPoints();

		// inference
		double xVertexCentre,yVertexCentre;
		double optimizationArray[3];
		void inferDF();
		void inferDDr();
		int getIndex(int a) { return cell.indexArray[a]; }
		double fret[1];
		int iterations[1];

		int n;
		double *vertices;

		// tessellation
		GLUtesselator *tobj;
		void tessellate();
		void draw();
		void drawForce();

		SelectionCell cell;

		virtual ~Selection();
};

#endif /* SELECTION_H_ */
